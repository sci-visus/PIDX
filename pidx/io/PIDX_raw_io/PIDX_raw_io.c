#include "../PIDX_io.h"
//#include "PIDX_raw_io.h"

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

    file->idx->enable_rst = 0;
  }
#else
  file->idx->enable_rst = 0;
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

    if (file->idx_dbg->debug_do_io == 1)
    {
      PIDX_rst_buf_aggregate_write(file->rst_id);
      //PIDX_rst_buf_aggregate_read(file->rst_id);
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->idx_dbg->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
    time->rst_end[start_index] = PIDX_get_time();
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
