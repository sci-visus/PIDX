#include "../PIDX_io.h"


PIDX_return_code PIDX_global_async_io(PIDX_hybrid_idx_io file, PIDX_file_io_id **io_id, Agg_buffer **agg_buffer, PIDX_block_layout* block_layout_by_level,  MPI_File *fp, MPI_Request *request,  int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level, int file_zero)
{
  int j;
  int rank = 0, nprocs = 1;
  int ret = 0;
  int j_1 = 0;
  int async_status = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;

  /*------------------------------------Create ALL the IDs [start]---------------------------------------*/

  for(j = layout_start ; j < layout_end; j++)
  {
    j_1 = j - layout_start;
    io_id[var_index][j] = PIDX_file_io_init(file->idx, file->idx_d, init_index, var_index, var_index);

    ret = PIDX_file_io_set_communicator(io_id[var_index][j], file->comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  //MPI_Request *req = *request;
  for(j = layout_start ; j < agg_io_level; j++)
  {
    j_1 = j - layout_start;
    Agg_buffer temp_agg = agg_buffer[var_index][j];
    PIDX_block_layout temp_layout = block_layout_by_level[j_1];// file->idx->variable[init_index]->block_layout_by_level_files[j_1];

    if (file->idx_dbg->debug_do_io == 1)
    {
      time->io_start[var_index][j] = PIDX_get_time();
      if (file_zero == 1)
      {
        ret = PIDX_async_aggregated_io(io_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE, /*&(req[j_1])*/&(request[j_1]), &(fp[j_1]), file->idx->filename_template_file_zero, async_status);
      }
      else
      {
        if (index == 0)
          ret = PIDX_async_aggregated_io(io_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE, /*&(req[j_1])*/&(request[j_1]), &(fp[j_1]), file->idx->filename_template_partition, async_status);
        else
          ret = PIDX_async_aggregated_io(io_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE, /*&(req[j_1])*/&(request[j_1]), &(fp[j_1]), file->idx->filename_template_global, async_status);
      }

      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      time->io_end[var_index][j] = PIDX_get_time();
    }

  }

  if (file->idx_dbg->debug_do_io == 1)
  {
    for(j = agg_io_level ; j < layout_end; j++)
    {
      time->io_per_process_start[var_index][j] = PIDX_get_time();
      //ret = PIDX_file_io_per_process(io_id[var_index][j - agg_io_level], file->idx->variable[init_index]->block_layout_by_level_files[j - agg_io_level], PIDX_WRITE);
      ret = PIDX_file_io_per_process(io_id[var_index][j - agg_io_level], block_layout_by_level[j - agg_io_level], PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      time->io_per_process_end[var_index][j] = PIDX_get_time();
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_global_aggregate(PIDX_hybrid_idx_io file, PIDX_agg_id** agg_id, Agg_buffer** agg_buffer, PIDX_block_layout* block_layout_by_level, PIDX_block_layout global_block_layout_files, MPI_Comm comm, int init_index, int var_index, int layout_start, int agg_io_level, int agg_factor, int file_status)
{
  int j;
  int rank = 0, nprocs = 1;
  int ret = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(comm,  &nprocs);
    MPI_Comm_rank(comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;

  /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
  /* Create the aggregation ID */
  /* Create the I/O ID */
  int j_1 = 0;

  for(j = layout_start ; j < agg_io_level; j++)
  {
    j_1 = j - layout_start;

    agg_id[var_index][j] = PIDX_agg_init(file->idx, file->idx_d, init_index, var_index, var_index);

    agg_buffer[var_index][j] = malloc(sizeof(*(agg_buffer[var_index][j])) );
    memset(agg_buffer[var_index][j], 0, sizeof(*(agg_buffer[var_index][j])) );

    agg_buffer[var_index][j]->file_number = -1;
    agg_buffer[var_index][j]->var_number = -1;
    agg_buffer[var_index][j]->sample_number = -1;

    agg_buffer[var_index][j]->no_of_aggregators = 0;
    agg_buffer[var_index][j]->aggregator_interval = 0;
    agg_buffer[var_index][j]->aggregation_factor = 1;
    //agg_buffer[var_index][j]->aggregation_factor = agg_factor;//(int)pow(2, (agg_io_level - j));

    ret = PIDX_agg_set_global_communicator(agg_id[var_index][j], file->global_comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    ret = PIDX_agg_set_communicator(agg_id[var_index][j], comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    PIDX_agg_id temp_id = agg_id[var_index][j];
    Agg_buffer temp_agg = agg_buffer[var_index][j];
    PIDX_block_layout temp_layout = block_layout_by_level[j_1];
    PIDX_block_layout temp_global_layout = global_block_layout_files;

    time->agg_meta_start[var_index][j] = PIDX_get_time();
    ret = PIDX_agg_meta_data_create(temp_id, temp_agg, temp_global_layout, temp_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    time->agg_meta_end[var_index][j] = PIDX_get_time();
    time->agg_buf_start[var_index][j] = PIDX_get_time();
    ret = PIDX_agg_buf_create_multiple_level(temp_id, temp_agg, temp_layout, temp_global_layout, var_index, j, file_status);
    //ret = PIDX_agg_buf_create(temp_id, temp_agg, temp_layout, temp_global_layout, var_index, j);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    time->agg_buf_end[var_index][j] = PIDX_get_time();


    if (file->idx_dbg->debug_do_agg == 1)
    {
      time->agg_start[var_index][j] = PIDX_get_time();
      ret = PIDX_agg_global_and_local(temp_id, temp_agg, j, temp_layout, temp_global_layout, PIDX_WRITE, var_index, j_1);

      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_end[var_index][j] = PIDX_get_time();
    }

    ret = PIDX_agg_meta_data_destroy(temp_id, temp_layout, temp_global_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

  }

  return PIDX_success;
}
