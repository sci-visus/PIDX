#include "../PIDX_io.h"

static PIDX_return_code rst_init(PIDX_hybrid_idx_io file, int gi, int svi, int evi);
static void set_default_patch_size(PIDX_hybrid_idx_io file, unsigned long long* process_bounds, int nprocs);
static int getPowerOftwo(int x);


PIDX_return_code restructure_init(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
  int ret = 0;
  int start_index = 0, end_index = 0;

  ret = rst_init(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (start_index = svi; start_index < evi; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (start_index + file->idx_d->var_pipe_length);

    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, svi, start_index, end_index);

    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_rst_aggregate_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;
  }

  return PIDX_success;
}



PIDX_return_code restructure(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
  int ret = 0;
  int start_index = 0, end_index = 0;

  for (start_index = svi; start_index < evi; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (start_index + file->idx_d->var_pipe_length);

    /* Perform data restructuring */
    if (file->idx_dbg->debug_do_rst == 1)
    {
      ret = PIDX_rst_write(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_WRITE);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure_cleanup(PIDX_hybrid_idx_io file, int gi)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  /* Destroy buffers allocated during restructuring phase */
  ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  ret = PIDX_rst_meta_data_destroy(file->rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  /* Deleting the restructuring ID */
  PIDX_rst_finalize(file->rst_id);


  free(var_grp->rank_r_offset);
  var_grp->rank_r_offset = 0;

  free(var_grp->rank_r_count);
  var_grp->rank_r_count = 0;

  return PIDX_success;
}


static PIDX_return_code rst_init(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
   int nprocs = 1;
   PIDX_variable_group var_grp = file->idx->variable_grp[gi];
   PIDX_variable var = var_grp->variable[svi];

#if PIDX_HAVE_MPI
  file->comm = file->global_comm;

  if (file->idx_d->parallel_mode == 1)
    MPI_Comm_size(file->comm,  &nprocs);
#endif

  var_grp->rank_r_offset = malloc(sizeof (unsigned long long) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(var_grp->rank_r_offset, 0, (sizeof (unsigned long long) * nprocs * PIDX_MAX_DIMENSIONS));

  var_grp->rank_r_count =  malloc(sizeof (unsigned long long) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(var_grp->rank_r_count, 0, (sizeof (unsigned long long) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(var->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, var_grp->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->comm);

    MPI_Allgather(var->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, var_grp->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(var_grp->rank_r_offset, var->sim_patch[0]->offset, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
    memcpy(var_grp->rank_r_count, var->sim_patch[0]->size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  }
#endif

  set_default_patch_size(file, var_grp->rank_r_count, nprocs);

  return PIDX_success;
}


static void set_default_patch_size(PIDX_hybrid_idx_io file, unsigned long long* process_bounds, int nprocs)
{
  int i = 0, j = 0;
  unsigned long long average_count = 0;
  int check_bit = 0;
  unsigned long long max_dim_length[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int equal_partiton = 0;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * 0 + i];
    for (j = 0; j < nprocs; j++)
    {
      if (max_dim_length[i] <= process_bounds[PIDX_MAX_DIMENSIONS * j + i])
        max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * j + i];
    }
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    average_count = average_count + max_dim_length[i];
  }
  average_count = average_count / PIDX_MAX_DIMENSIONS;
  average_count = getPowerOftwo(average_count);

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    check_bit = check_bit || ((double) file->idx->bounds[i] / average_count > (double) file->idx->bounds[i] / max_dim_length[i]);

  while (check_bit)
  {
    average_count = average_count * 2;
    check_bit = 0;
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      check_bit = check_bit || ((double) file->idx->bounds[i] / average_count > (double) file->idx->bounds[i] / max_dim_length[i]);
  }
  unsigned long long reg_patch_size[PIDX_MAX_DIMENSIONS];
  //reg_patch_size =  average_count;
  if (equal_partiton == 1)
  {
    reg_patch_size[0] = average_count / 1;
    reg_patch_size[1] = average_count / 1;
    reg_patch_size[2] = average_count / 1;
  }
  else
  {
    reg_patch_size[0] = getPowerOftwo(max_dim_length[0]) * 1;
    reg_patch_size[1] = getPowerOftwo(max_dim_length[1]) * 1;
    reg_patch_size[2] = getPowerOftwo(max_dim_length[2]) * 1;
  }

  memcpy(file->idx->reg_patch_size, reg_patch_size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  //reg_patch_size = reg_patch_size * 4;
}


static int getPowerOftwo(int x)
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}
