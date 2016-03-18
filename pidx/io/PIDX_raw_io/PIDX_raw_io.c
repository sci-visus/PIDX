#include "../PIDX_io.h"

static int maximum_neighbor_count = 256;

static int intersectNDChunk(Ndim_patch A, Ndim_patch B);


struct PIDX_raw_io_descriptor
{

#if PIDX_HAVE_MPI
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
#endif

  PIDX_header_io_id header_io_id;              ///< IDX metadata id
  PIDX_rst_id rst_id;                          ///< Restructuring phase id

  //int flush_used;
  //int write_on_close;                          ///< HPC Writes
  int one_time_initializations;                ///<


  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;

  //PIDX_time time;
};


PIDX_raw_io PIDX_raw_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_raw_io raw_io_id;
  raw_io_id = malloc(sizeof (*raw_io_id));
  memset(raw_io_id, 0, sizeof (*raw_io_id));

  raw_io_id->idx = idx_meta_data;
  raw_io_id->idx_d = idx_derived_ptr;
  raw_io_id->idx_dbg = idx_dbg;

  return (raw_io_id);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_raw_io_set_communicator(PIDX_raw_io id, MPI_Comm comm)
{
  if (id == NULL)
    return PIDX_err_id;

  id->comm = comm;

  return PIDX_success;
}
#endif

PIDX_return_code PIDX_raw_write(PIDX_raw_io file, int start_var_index, int end_var_index)
{
  file->idx_d->var_pipe_length = file->idx->variable_count - 1;
  if (file->idx_d->var_pipe_length == 0)
    file->idx_d->var_pipe_length = 1;

  PIDX_return_code ret;
  int nprocs = 1;

  PIDX_time time = file->idx_d->time;
  time->populate_idx_start_time = PIDX_get_time();

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
    MPI_Comm_size(file->comm,  &nprocs);
#endif

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char *directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", directory_path);
  sprintf(size_path, "%s_SIZE", directory_path);
  free(directory_path);


  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);

    file->idx->enable_rst = 0;
  }
#endif


  if (file->idx_d->parallel_mode == 1)
  {
    PIDX_rst_id temp_id;
    temp_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_var_index, start_var_index);

#if PIDX_HAVE_MPI
    ret = PIDX_rst_set_communicator(temp_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_rst;
#endif

    ret = PIDX_rst_meta_data_create(temp_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    PIDX_variable var0 = file->idx->variable[start_var_index];
    int max_patch_count;
    int patch_count =var0->patch_group_count;
    MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->comm);

    int64_t *local_patch_offset = malloc(sizeof(int64_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));
    memset(local_patch_offset, 0, sizeof(int64_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));

    int64_t *local_patch_size = malloc(sizeof(int64_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));
    memset(local_patch_size, 0, sizeof(int64_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));

    int pcounter = 0;
    int i = 0, d = 0;
    local_patch_offset[0] = patch_count;
    local_patch_size[0] = patch_count;
    for (i = 0; i < patch_count; i++)
    {
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        local_patch_offset[i * PIDX_MAX_DIMENSIONS + d + 1] = var0->rst_patch_group[i]->reg_patch->offset[d];
        local_patch_size[i * PIDX_MAX_DIMENSIONS + d + 1] = var0->rst_patch_group[i]->reg_patch->size[d];
      }
      pcounter++;
    }

    file->idx_d->global_patch_offset = malloc((nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t));
    memset(file->idx_d->global_patch_offset, 0,(nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t));

    file->idx_d->global_patch_size = malloc((nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t));
    memset(file->idx_d->global_patch_size, 0, (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t));

#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode == 1)
    {
      MPI_Allgather(local_patch_offset, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_LONG_LONG, file->idx_d->global_patch_offset + 2, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_LONG_LONG, file->comm);

      MPI_Allgather(local_patch_size, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_LONG_LONG, file->idx_d->global_patch_size + 2, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_LONG_LONG, file->comm);
    }
    else
    {
      memcpy(file->idx_d->global_patch_offset, local_patch_offset, sizeof(uint64_t) * (PIDX_MAX_DIMENSIONS * max_patch_count + 1));
      memcpy(file->idx_d->global_patch_size, local_patch_size, sizeof(uint64_t) * (PIDX_MAX_DIMENSIONS * max_patch_count + 1));

      file->idx->enable_rst = 0;
    }
    file->idx_d->global_patch_size[0] = nprocs;
    file->idx_d->global_patch_offset[0] = nprocs;
    file->idx_d->global_patch_size[1] = max_patch_count;
    file->idx_d->global_patch_offset[1] = max_patch_count;

    int rank = 0;
    MPI_Comm_rank(file->comm, &rank);
    if (rank == 0)
    {
      int fp = open(offset_path, O_CREAT | O_WRONLY, 0664);
      ssize_t write_count = pwrite(fp, file->idx_d->global_patch_offset, (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t), 0);
      if (write_count != (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      close(fp);

      fp = open(size_path, O_CREAT | O_WRONLY, 0664);
      write_count = pwrite(fp, file->idx_d->global_patch_size, (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t), 0);
      if (write_count != (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(int64_t))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      close(fp);

      /*
      uint64_t *buffer = malloc((nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 1) * sizeof(int64_t));
      memset(buffer, 0, (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 1) * sizeof(int64_t));

      fp = open("Size_File", O_RDONLY);
      write_count = pread(fp, buffer, (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 1) * sizeof(int64_t), 0);
      if (write_count != (nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 1) * sizeof(int64_t))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      close(fp);
      printf("np %lld mp %lld:  %lld %lld %lld %lld %lld\n", buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6]);
      */
    }

    free(local_patch_offset);
    free(local_patch_size);

    ret = PIDX_rst_meta_data_destroy(temp_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    PIDX_rst_finalize(temp_id);
  }
#else
  file->idx->enable_rst = 0;
#endif

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, 0);

    /*
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;
    */

    file->one_time_initializations = 1;
  }

  time->populate_idx_end_time = PIDX_get_time();

  time->write_init_start[time->header_counter] = PIDX_get_time();

#if !SIMULATE_IO
  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif

    ret = PIDX_header_io_enable_raw_dump (file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
#endif

  time->write_init_end[time->header_counter] = PIDX_get_time();


  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);


    /*--------------------------------------Create RST IDs [start]------------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*----------------------------------------Create RST IDs [end]------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    /* Attaching the communicator to the restructuring phase */
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/



    /*--------------------------------------------RST [start]------------------------------------------------*/
    time->rst_start[start_index] = PIDX_get_time();

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
    }
    time->rst_end[start_index] = PIDX_get_time();

    time->rst_io_start[start_index] = PIDX_get_time();
    if (file->idx_dbg->debug_do_io == 1)
    {
      ret = PIDX_rst_buf_aggregate_write(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->idx_dbg->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
    time->rst_io_end[start_index] = PIDX_get_time();

    /*--------------------------------------------RST [end]---------------------------------------------------*/



    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[start_index] = PIDX_get_time();

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    time->finalize_end[start_index] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/
  }

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  free(file->idx_d->global_patch_offset);
  file->idx_d->global_patch_offset = 0;

  free(file->idx_d->global_patch_size);
  file->idx_d->global_patch_size = 0;

  return PIDX_success;
}


  static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
  {
    int d = 0, check_bit = 0;
    for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
      check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

    return !(check_bit);
  }


PIDX_return_code PIDX_forced_raw_read(PIDX_raw_io file, int start_var_index, int end_var_index)
{
  file->idx_d->var_pipe_length = file->idx->variable_count - 1;
  if (file->idx_d->var_pipe_length == 0)
    file->idx_d->var_pipe_length = 1;

  int nprocs = 1, rank = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  char *idx_directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  idx_directory_path = malloc(sizeof(*idx_directory_path) * PATH_MAX);
  memset(idx_directory_path, 0, sizeof(*idx_directory_path) * PATH_MAX);
  strncpy(idx_directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", idx_directory_path);
  sprintf(size_path, "%s_SIZE", idx_directory_path);
  free(idx_directory_path);

  int64_t number_cores = 0;
  int fp = open(size_path, O_RDONLY);
  ssize_t write_count = pread(fp, &number_cores, sizeof(int64_t), 0);
  if (write_count != sizeof(int64_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  int64_t max_patch_count = 0;
  write_count = pread(fp, &max_patch_count, sizeof(int64_t), sizeof(int64_t));
  if (write_count != sizeof(int64_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  uint64_t *size_buffer = malloc((number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t));
  memset(size_buffer, 0, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t));

  write_count = pread(fp, size_buffer, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t), 2 * sizeof(int64_t));
  if (write_count != (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  close(fp);


  uint64_t *offset_buffer = malloc((number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t));
  memset(offset_buffer, 0, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t));

  int fp1 = open(offset_path, O_RDONLY);
  write_count = pread(fp1, offset_buffer, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t), 2 * sizeof(int64_t));
  if (write_count != (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(int64_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }
  close(fp1);


  PIDX_time time = file->idx_d->time;
  time->populate_idx_start_time = PIDX_get_time();

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, 0);
    file->one_time_initializations = 1;
  }

  time->populate_idx_end_time = PIDX_get_time();

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, file->idx->current_time_step);


  int n = 0, m = 0, d = 0;
  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_proc_patch->offset[d] = file->idx->variable[start_var_index]->sim_patch[0]->offset[d];
    local_proc_patch->size[d] = file->idx->variable[start_var_index]->sim_patch[0]->size[d];
  }

  // PC - - - - - - - - - -     PC - - - - - - - - - -            PC - - - - -
  // 0  1 2 3 4 5 6 7 8 9 10   11 12 13 14 15 16 17 18 19 20 21  22
  Ndim_patch n_proc_patch = (Ndim_patch)malloc(sizeof (*n_proc_patch));
  memset(n_proc_patch, 0, sizeof (*n_proc_patch));
  int p_counter = 1;
  int pc = 0, pc_index = 0;

  Ndim_patch_group patch_grp;
  patch_grp = malloc(sizeof(*(patch_grp)));
  memset(patch_grp, 0, sizeof(*(patch_grp)));

  int *source_patch_id = malloc(sizeof(int) * maximum_neighbor_count);
  patch_grp->source_patch_rank = (int*)malloc(sizeof(int) * maximum_neighbor_count);
  patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
  patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
  memset(source_patch_id, 0, sizeof(int) * maximum_neighbor_count);
  memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
  memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
  memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

  int patch_count = 0;
  for (n = 0; n < number_cores; n++)
  {
    pc = offset_buffer[n * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)];
    pc_index = n * (max_patch_count * PIDX_MAX_DIMENSIONS + 1);
    for (m = 0; m < pc; m++)
    {
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        n_proc_patch->offset[d] = offset_buffer[pc_index + m * PIDX_MAX_DIMENSIONS + d + 1];
        n_proc_patch->size[d] = size_buffer[pc_index + m * PIDX_MAX_DIMENSIONS + d + 1];
      }

      if (intersectNDChunk(local_proc_patch, n_proc_patch))
      {
        //sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, n, m);
        patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
        memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));
        for (d = 0; d < 3; d++)
        {
          //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
          if (n_proc_patch->offset[d] <= local_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
          {
            patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
            patch_grp->patch[patch_count]->size[d] = (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) - local_proc_patch->offset[d] + 1;
          }
          else if (local_proc_patch->offset[d] <= n_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) >= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
          {
            patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
            patch_grp->patch[patch_count]->size[d] = (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) - n_proc_patch->offset[d] + 1;
          }
          else if (( local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) <= (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) && local_proc_patch->offset[d] >= n_proc_patch->offset[d])
          {
            patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
            patch_grp->patch[patch_count]->size[d] = local_proc_patch->size[d];
          }
          else if (( n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) && n_proc_patch->offset[d] >= local_proc_patch->offset[d])
          {
            patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
            patch_grp->patch[patch_count]->size[d] = n_proc_patch->size[d];
          }
        }
        patch_grp->source_patch_rank[patch_count] = n;
        source_patch_id[patch_count] = m;

        patch_count++;
      }
      p_counter++;
    }
  }

  free(local_proc_patch);
  local_proc_patch = 0;
  free(n_proc_patch);
  n_proc_patch = 0;

  uint64_t k1 = 0, j1 = 0, i1 = 0, i = 0, j = 0;
  int count1 = 0, send_o = 0, send_c = 0, index = 0;
  int64_t sim_patch_offsetx[5];
  int64_t sim_patch_countx[5];

  unsigned char ***temp_patch_buffer;
  temp_patch_buffer = malloc(sizeof(*temp_patch_buffer) * (end_var_index - start_var_index));
  memset(temp_patch_buffer, 0, sizeof(*temp_patch_buffer) * (end_var_index - start_var_index));
  //printf("start_var_index = %d end_var_index %d\n", start_var_index, end_var_index);
  for (i = 0; i < (end_var_index - start_var_index); i++)
  {
    PIDX_variable var = file->idx->variable[i + start_var_index];
    temp_patch_buffer[i] = malloc(sizeof(*(temp_patch_buffer[i])) * patch_count);
    memset(temp_patch_buffer[i], 0, sizeof(*(temp_patch_buffer[i])) * patch_count);

    for (j = 0; j < patch_count; j++)
    {
      temp_patch_buffer[i][j] = malloc(sizeof(*(temp_patch_buffer[i][j])) * patch_grp->patch[j]->size[0] * patch_grp->patch[j]->size[1] * patch_grp->patch[j]->size[2] * var->bits_per_value/8 * var->values_per_sample);
      memset(temp_patch_buffer[i][j], 0, sizeof(*(temp_patch_buffer[i][j])) * patch_grp->patch[j]->size[0] * patch_grp->patch[j]->size[1] * patch_grp->patch[j]->size[2] * var->bits_per_value/8 * var->values_per_sample);
    }
  }

  for (i = 0; i < patch_count; i++)
  {
    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, patch_grp->source_patch_rank[i], source_patch_id[i]);
    int fpx = open(file_name, O_RDONLY);
    //patch_grp->patch[i]->buffer = malloc(patch_grp->patch[i]->size[0] * patch_grp->patch[i]->size[1] * patch_grp->patch[i]->size[2] * sizeof(double));

    pc_index = patch_grp->source_patch_rank[i] * (max_patch_count * PIDX_MAX_DIMENSIONS + 1);
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      sim_patch_offsetx[d] = offset_buffer[pc_index + source_patch_id[i] * PIDX_MAX_DIMENSIONS + d + 1];
      sim_patch_countx[d] = size_buffer[pc_index + source_patch_id[i] * PIDX_MAX_DIMENSIONS + d + 1];
    }

    count1 = 0;

    for (k1 = patch_grp->patch[i]->offset[2]; k1 < patch_grp->patch[i]->offset[2] + patch_grp->patch[i]->size[2]; k1++)
    {
      for (j1 = patch_grp->patch[i]->offset[1]; j1 < patch_grp->patch[i]->offset[1] + patch_grp->patch[i]->size[1]; j1++)
      {
        for (i1 = patch_grp->patch[i]->offset[0]; i1 < patch_grp->patch[i]->offset[0] + patch_grp->patch[i]->size[0]; i1 = i1 + patch_grp->patch[i]->size[0])
        {
          int64_t *sim_patch_offset = sim_patch_offsetx;// file->idx->variable[start_var_index]->sim_patch[0]->offset;
          int64_t *sim_patch_count = sim_patch_countx;// file->idx->variable[start_var_index]->sim_patch[0]->size;

          index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                  (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                  (i1 - sim_patch_offset[0]);

          int start_index = 0, other_offset = 0, v1 = 0;
          for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + 1)
          {
            other_offset = 0;
            for (v1 = 0; v1 < start_index; v1++)
            {
              PIDX_variable var1 = file->idx->variable[v1];
              other_offset = other_offset + ((var1->bits_per_value/8) * var1->values_per_sample * sim_patch_countx[0] * sim_patch_countx[1] * sim_patch_countx[2]);
              //printf("[%d] ----> %d [%d %d %d]\n", start_index, other_offset, sim_patch_countx[0], sim_patch_countx[1], sim_patch_countx[2]);
            }

            PIDX_variable var = file->idx->variable[start_index];
            send_o = index * var->values_per_sample;
            send_c = patch_grp->patch[i]->size[0] * var->values_per_sample;

            //size_t preadc = pread(fp, patch_grp->patch[i]->buffer + (count1 * send_c * var->bits_per_value/8), send_c * var->bits_per_value/8, send_o * var->bits_per_value/8);

            size_t preadc = pread(fp, temp_patch_buffer[start_index - start_var_index][i] + (count1 * send_c * var->bits_per_value/8), send_c * var->bits_per_value/8, other_offset + (send_o * var->bits_per_value/8));

            if (preadc != send_c * var->bits_per_value/8)
            {
              fprintf(stderr, "[%s] [%d] Error in pread [%d %d]\n", __FILE__, __LINE__, (int)preadc, (int)send_c * var->bits_per_value/8);
              return PIDX_err_rst;
            }
          }
          count1++;
        }
      }
    }
    close(fpx);
  }

  int r, recv_o;
  for (r = 0; r < patch_count; r++)
  {
    //Ndim_patch out_patch = file->idx->variable[start_var_index]->sim_patch[0];// var->rst_patch_group[g]->reg_patch;
    //int nx = out_patch->size[0];
    //int ny = out_patch->size[1];

    for (k1 = patch_grp->patch[r]->offset[2]; k1 < patch_grp->patch[r]->offset[2] + patch_grp->patch[r]->size[2]; k1++)
    {
      for (j1 = patch_grp->patch[r]->offset[1]; j1 < patch_grp->patch[r]->offset[1] + patch_grp->patch[r]->size[1]; j1++)
      {
        for (i1 = patch_grp->patch[r]->offset[0]; i1 < patch_grp->patch[r]->offset[0] + patch_grp->patch[r]->size[0]; i1 = i1 + patch_grp->patch[r]->size[0])
        {
          index = ((patch_grp->patch[r]->size[0])* (patch_grp->patch[r]->size[1]) * (k1 - patch_grp->patch[r]->offset[2])) + ((patch_grp->patch[r]->size[0]) * (j1 - patch_grp->patch[r]->offset[1])) + (i1 - patch_grp->patch[r]->offset[0]);

          int start_index;
          for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + 1)
          {
            PIDX_variable var = file->idx->variable[start_index];

            send_o = index * var->values_per_sample * (var->bits_per_value/8);
            send_c = (patch_grp->patch[r]->size[0]);
            recv_o = (file->idx->variable[start_index]->sim_patch[0]->size[0] * file->idx->variable[start_index]->sim_patch[0]->size[1] * (k1 - file->idx->variable[start_index]->sim_patch[0]->offset[2])) + (file->idx->variable[start_index]->sim_patch[0]->size[0] * (j1 - file->idx->variable[start_index]->sim_patch[0]->offset[1])) + (i1 - file->idx->variable[start_index]->sim_patch[0]->offset[0]);

#if !SIMULATE_IO
            //memcpy(file->idx->variable[start_index]->sim_patch[0]->buffer + (recv_o * var->values_per_sample * (var->bits_per_value/8)), patch_grp->patch[r]->buffer + send_o, send_c * var->values_per_sample * (var->bits_per_value/8));
            memcpy(file->idx->variable[start_index]->sim_patch[0]->buffer + (recv_o * var->values_per_sample * (var->bits_per_value/8)), temp_patch_buffer[start_index - start_var_index][r] + send_o, send_c * var->values_per_sample * (var->bits_per_value/8));
          }
#endif
        }
      }
    }
    free(patch_grp->patch[r]);
  }

  for (i = 0; i < (end_var_index - start_var_index); i++)
  {
    for (j = 0; j < patch_count; j++)
      free(temp_patch_buffer[i][j]);

    free(temp_patch_buffer[i]);
  }
  free(temp_patch_buffer);

  free(patch_grp->source_patch_rank);
  free(patch_grp->patch);
  free(patch_grp->reg_patch);
  free(patch_grp);
  free(source_patch_id);

  free(file_name);
  free(data_set_path);
  free(directory_path);

  free(offset_buffer);
  free(size_buffer);

  return PIDX_success;
}



PIDX_return_code PIDX_raw_read(PIDX_raw_io file, int start_var_index, int end_var_index)
{
  file->idx_d->var_pipe_length = file->idx->variable_count - 1;
  if (file->idx_d->var_pipe_length == 0)
    file->idx_d->var_pipe_length = 1;


  PIDX_return_code ret;
  int nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
    MPI_Comm_size(file->comm,  &nprocs);
#endif

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  }
#endif

  PIDX_time time = file->idx_d->time;
  time->populate_idx_start_time = PIDX_get_time();

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, 0);

    /*
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;
    */

    file->one_time_initializations = 1;
  }

  time->populate_idx_end_time = PIDX_get_time();

  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);


    /*--------------------------------------Create RST IDs [start]------------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*----------------------------------------Create RST IDs [end]------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    /* Attaching the communicator to the restructuring phase */
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/



    /*--------------------------------------------RST [start]------------------------------------------------*/
    time->rst_start[start_index] = PIDX_get_time();

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    if (file->idx_dbg->debug_do_io == 1)
    {
      ret = PIDX_rst_buf_aggregate_read(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->idx_dbg->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    /* Perform data restructuring */
    if (file->idx_dbg->debug_do_rst == 1)
    {
      ret = PIDX_rst_read(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    time->rst_end[start_index] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/



    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[start_index] = PIDX_get_time();

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    time->finalize_end[start_index] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

  }

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_raw_io_finalize(PIDX_raw_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}
