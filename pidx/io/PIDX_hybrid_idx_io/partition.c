#include "../PIDX_io.h"

///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
static PIDX_return_code restructure(PIDX_hybrid_idx_io file, int gi, int svi, int evi);


///
/// \brief rst_init
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
static PIDX_return_code rst_init(PIDX_hybrid_idx_io file, int gi, int svi, int evi);


///
/// \brief partition_setup
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
static PIDX_return_code partition_setup(PIDX_hybrid_idx_io file, int gi, int svi, int evi);


///
/// \brief getPowerOftwo
/// \param x
/// \return
///
static int getPowerOftwo(int x);


///
/// \brief set_default_patch_size
/// \param file
/// \param process_bounds
/// \param nprocs
///
static void set_default_patch_size(PIDX_hybrid_idx_io file, unsigned long long* process_bounds, int nprocs);


///
/// \brief intersectNDChunk
/// \param A
/// \param B
/// \return
///
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);



PIDX_return_code partition(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{

  int ret = 0;
  // Restructuring into power two domain
  ret = restructure(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // calculate number of partitions
  int d = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    file->idx_d->idx_count[d] = file->idx->bounds[d] / file->idx_d->idx_size[d];
    if (file->idx->bounds[d] % file->idx_d->idx_size[d] != 0)
      file->idx_d->idx_count[d]++;

    file->idx_d->idx_count[d] = pow(2, (int)ceil(log2(file->idx_d->idx_count[d])));
  }

  // perform partitioning, make sub-communicators
  ret = partition_setup(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int grank = 0, gnprocs = 1;
  MPI_Comm_rank(file->global_comm, &grank);
  MPI_Comm_size(file->global_comm, &gnprocs);
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  var_grp->rank_buffer = malloc(gnprocs * sizeof(*var_grp->rank_buffer));
  memset(var_grp->rank_buffer, 0, gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&grank, 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->global_comm);

  return PIDX_success;
}


static PIDX_return_code restructure(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
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

  //set_default_patch_size(file, var_grp->rank_r_count, nprocs);

  return PIDX_success;
}


static PIDX_return_code partition_setup(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
  int *colors;
  int index_i = 0, index_j = 0, index_k = 0;
  int i = 0, j = 0, k = 0, d = 0;
  //PIDX_return_code ret;

  int rank = 0;
  MPI_Comm_rank(file->global_comm, &rank);

  colors = malloc(sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  memset(colors, 0, sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  file->idx_d->color = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]) + 1;

  for (k = 0; k < file->idx_d->idx_count[2]; k++)
    for (j = 0; j < file->idx_d->idx_count[1]; j++)
      for (i = 0; i < file->idx_d->idx_count[0]; i++)
      {
        colors[(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i] = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i;
      }

  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[svi];

  if (var->patch_group_count == 1)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = var->rst_patch_group[0]->reg_patch->offset[d];
      local_proc_patch->size[d] = var->rst_patch_group[0]->reg_patch->size[d];
    }
    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    PointND bounds_point;
    int maxH = 0;
    bounds_point.x = (int) file->idx_d->idx_count[0];
    bounds_point.y = (int) file->idx_d->idx_count[1];
    bounds_point.z = (int) file->idx_d->idx_count[2];
    bounds_point.u = (int) 1;
    bounds_point.v = (int) 1;
    char bitSequence[512];
    char bitPattern[512];
    GuessBitmaskPattern(bitSequence, bounds_point);
    maxH = strlen(bitSequence);

    for (i = 0; i <= maxH; i++)
      bitPattern[i] = RegExBitmaskBit(bitSequence, i);

    int z_order = 0;
    int number_levels = maxH - 1;

    for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx_d->idx_size[0], index_i++)
    {
      for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx_d->idx_size[1], index_j++)
      {
        for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx_d->idx_size[2], index_k++)
        {
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = file->idx_d->idx_size[0];
          reg_patch->size[1] = file->idx_d->idx_size[1];
          reg_patch->size[2] = file->idx_d->idx_size[2];

         if (intersectNDChunk(reg_patch, local_proc_patch))
          {
            PointND xyzuv_Index;
            xyzuv_Index.x = index_i;
            xyzuv_Index.y = index_j;
            xyzuv_Index.z = index_k;
            xyzuv_Index.u = 0;
            xyzuv_Index.v = 0;

            z_order = 0;
            PointND zero;
            zero.x = 0;
            zero.y = 0;
            zero.z = 0;
            zero.u = 0;
            zero.v = 0;
            memset(&zero, 0, sizeof (PointND));

            int cnt = 0;
            for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
            {
              int bit = bitPattern[number_levels];
              z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            file->idx_d->color = colors[z_order];
            //printf("[%d] ---> %d\n", rank, file->idx_d->color);

            assert(var->sim_patch_count == 1);
            //var->sim_patch_count = 1;
            var->sim_patch[0]->offset[0] = var->rst_patch_group[0]->reg_patch->offset[0];
            var->sim_patch[0]->offset[1] = var->rst_patch_group[0]->reg_patch->offset[1];
            var->sim_patch[0]->offset[2] = var->rst_patch_group[0]->reg_patch->offset[2];
            var->sim_patch[0]->offset[3] = 0;
            var->sim_patch[0]->offset[4] = 0;

            var->sim_patch[0]->size[0] = var->rst_patch_group[0]->reg_patch->size[0];
            var->sim_patch[0]->size[1] = var->rst_patch_group[0]->reg_patch->size[1];
            var->sim_patch[0]->size[2] = var->rst_patch_group[0]->reg_patch->size[2];
            var->sim_patch[0]->size[3] = 1;
            var->sim_patch[0]->size[4] = 1;

            break;
          }
        }
      }
    }
    free(reg_patch);
  }
  else
  {
    file->idx->bounds[0] = 0;//reg_patch->size[0];
    file->idx->bounds[1] = 0;//reg_patch->size[1];
    file->idx->bounds[2] = 0;//reg_patch->size[2];

    var->sim_patch_count = 0;
    var->sim_patch[0]->offset[0] = 0;
    var->sim_patch[0]->offset[1] = 0;
    var->sim_patch[0]->offset[2] = 0;
    var->sim_patch[0]->offset[3] = 0;
    var->sim_patch[0]->offset[4] = 0;

    var->sim_patch[0]->size[0] = 0;//-1;
    var->sim_patch[0]->size[1] = 0;//-1;
    var->sim_patch[0]->size[2] = 0;//-1;
    var->sim_patch[0]->size[3] = 0;//-1;
    var->sim_patch[0]->size[4] = 0;//-1;
  }


  free(colors);

  //
  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';

  if (file->idx_d->idx_count[0] != 1 || file->idx_d->idx_count[1] != 1 || file->idx_d->idx_count[2] != 1)
  {
    sprintf(file->idx->filename_partition, "%s_%d.idx", file_name_skeleton, file->idx_d->color);
    sprintf(file->idx->filename_file_zero, "%s_%s.idx", file_name_skeleton, "file_zero");
  }
  else
  {
    strcpy(file->idx->filename_partition, file->idx->filename);
    strcpy(file->idx->filename_file_zero, file->idx->filename);
  }

  strcpy(file->idx->filename_global, file->idx->filename);

  free(local_proc_patch);
  local_proc_patch = 0;

  return PIDX_success;
}


PIDX_return_code create_local_comm(PIDX_hybrid_idx_io file)
{
  int rank = 0;
  int nprocs = 1;
  int ret;
  MPI_Comm_rank(file->global_comm, &rank);
  MPI_Comm_size(file->global_comm, &nprocs);

  ret = MPI_Comm_split(file->global_comm, file->idx_d->color, rank, &(file->comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  return PIDX_success;
}


PIDX_return_code destroy_local_comm(PIDX_hybrid_idx_io file)
{
  int ret;
  ret = MPI_Comm_free(&(file->comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  return PIDX_success;
}


PIDX_return_code partition_cleanup(PIDX_hybrid_idx_io file, int gi)
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


static void set_default_patch_size(PIDX_hybrid_idx_io file, unsigned long long* process_bounds, int nprocs)
{
  int i = 0, j = 0;
  unsigned long long average_count = 0;
  int check_bit = 0;
  unsigned long long max_dim_length[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
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
    reg_patch_size[3] = 1;
    reg_patch_size[4] = 1;
  }
  else
  {
    reg_patch_size[0] = getPowerOftwo(max_dim_length[0]) * 1;
    reg_patch_size[1] = getPowerOftwo(max_dim_length[1]) * 1;
    reg_patch_size[2] = getPowerOftwo(max_dim_length[2]) * 1;
    //printf("Box size = %d %d %d\n", reg_patch_size[0], reg_patch_size[1], reg_patch_size[2]);
    //reg_patch_size[0] = getPowerOftwo(process_bounds[0]) * 1;
    //reg_patch_size[1] = getPowerOftwo(process_bounds[1]) * 1;
    //reg_patch_size[2] = getPowerOftwo(process_bounds[2]) * 1;
    reg_patch_size[3] = 1;//getPowerOftwo(process_bounds[3]) * 1;
    reg_patch_size[4] = 1;//getPowerOftwo(process_bounds[4]) * 1;
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

static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}
