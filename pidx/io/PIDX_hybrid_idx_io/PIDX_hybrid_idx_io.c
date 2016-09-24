#include "../PIDX_io.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"


PIDX_hybrid_idx_io PIDX_hybrid_idx_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_hybrid_idx_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_d = idx_derived_ptr;
  idx_io_id->idx_dbg = idx_dbg;

  return (idx_io_id);
}



PIDX_return_code PIDX_hybrid_idx_io_set_communicator(PIDX_hybrid_idx_io id, MPI_Comm comm)
{
  if (id == NULL)
    return PIDX_err_id;

  id->global_comm = comm;

  return PIDX_success;
}



PIDX_return_code PIDX_hybrid_idx_write(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index)
{
  PIDX_time time = file->idx_d->time;
  time->SX = PIDX_get_time();

  PIDX_return_code ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  int i = 0;
  int nprocs = 1, rank = 0;
  int start_index = 0;

  time->partition_start_time = PIDX_get_time();

  ret = partition_domain(file, group_index, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->partition_end_time = PIDX_get_time();

  int grank = 0, gnprocs = 1;
  MPI_Comm_rank(file->global_comm, &grank);
  MPI_Comm_size(file->global_comm, &gnprocs);
  var_grp->rank_buffer = malloc(gnprocs * sizeof(*var_grp->rank_buffer));
  memset(var_grp->rank_buffer, 0, gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&grank, 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->global_comm);

  ret = populate_idx_file_structure(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = one_time_initialize(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  time->hz_s_time = PIDX_get_time();
  ret = create_hz_buffers(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->hz_e_time = PIDX_get_time();


  file->idx_d->shared_block_level = (int)log2(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]) + file->idx->bits_per_block + 1;
  if (file->idx_d->shared_block_level >= file->idx_d->maxh)
    file->idx_d->shared_block_level = file->idx_d->maxh;

  int partion_level = (int) log2(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  int hz_from_file_zero = 0, hz_from_shared = 0, hz_from_non_shared = 0;
  int hz_to_file_zero = 0, hz_to_shared = 0, hz_to_non_shared = 0;

  ret = init_agg_io_buffer(file, group_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int modes;
  switch(modes)
  {
    case 0:
      hz_from_file_zero = 0;
      hz_to_file_zero =  file->idx_d->shared_block_level;

      hz_from_shared = file->idx_d->shared_block_level;
      hz_to_shared =  file->idx_d->total_partiton_level;

      hz_from_non_shared = file->idx_d->total_partiton_level;
      hz_to_non_shared =  file->idx_d->maxh;

      if (hz_from_file_zero == hz_to_file_zero)
      {
        var_grp->f0_start_layout_index = 0;
        var_grp->f0_end_layout_index = 0;
        var_grp->f0_layout_count = 0;
      }
      else
      {
        ret = populate_idx_dataset(file,
                                   var_grp->f0_block_layout,
                                   var_grp->f0_block_layout_by_level,
                                   &(var_grp->f0_start_layout_index),
                                   &(var_grp->f0_end_layout_index),
                                   &(var_grp->f0_layout_count),
                                   group_index,
                                   start_var_index, end_var_index,
                                   hz_from_file_zero, hz_to_file_zero);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
        }
      }

      break;

    case 1:

      hz_from_file_zero = 0;
      hz_to_file_zero =  0;

      hz_from_shared = 0;
      hz_to_shared =  file->idx_d->total_partiton_level;

      hz_from_non_shared = file->idx_d->total_partiton_level;
      hz_to_non_shared =  file->idx_d->maxh;

      break;

    case 2:

      hz_from_file_zero = 0;
      hz_to_file_zero =  0;

      hz_from_shared = 0;
      hz_to_shared =  0;

      hz_from_non_shared = 0;
      hz_to_non_shared =  file->idx_d->maxh;

      break;
  }

  if (hz_from_file_zero == hz_to_file_zero)
  {
    var_grp->f0_start_layout_index = 0;
    var_grp->f0_end_layout_index = 0;
  }

  if (hz_from_shared == hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
  }

  if (hz_from_non_shared == hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
  }

  ret = partition_communicator(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }

  time->populate_idx_start_time_s = PIDX_get_time();

  memset(var_grp->rank_buffer, 0, gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&rank, 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->global_comm);


  ret = populate_idx_dataset(file,
                             var_grp->shared_block_layout, var_grp->shared_block_layout_by_level,
                             &(var_grp->shared_start_layout_index), &(var_grp->shared_end_layout_index),
                             &(var_grp->shared_layout_count),
                             group_index,  start_var_index, end_var_index,
                             hz_from_shared, hz_to_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->populate_idx_end_time_s = PIDX_get_time();


  time->populate_idx_start_time_ns = PIDX_get_time();


  ret = populate_idx_dataset(file,
                             var_grp->nshared_block_layout, var_grp->nshared_block_layout_by_level,
                             &(var_grp->nshared_start_layout_index), &(var_grp->nshared_end_layout_index),
                             &(var_grp->nshared_layout_count),
                             group_index,  start_var_index, end_var_index,
                             hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->populate_idx_end_time_ns = PIDX_get_time();

  ret = write_headers(file, group_index, start_var_index, end_var_index, 0);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  int agg_io_level_non_shared = 0, no_of_aggregators = 0, agg_io_level_shared = 0, agg_io_level_file_zero = 0;
  if (file->idx->enable_agg == 0)
    agg_io_level_non_shared = var_grp->nshared_start_layout_index;
  else
  {
    for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->nshared_block_layout_by_level[i - var_grp->nshared_start_layout_index]->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level_non_shared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    agg_io_level_shared = var_grp->shared_start_layout_index;
  else
  {
    for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->shared_block_layout_by_level[i - var_grp->shared_start_layout_index]->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level_shared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    agg_io_level_file_zero = var_grp->f0_start_layout_index;
  else
  {
    for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->f0_block_layout_by_level[i - var_grp->f0_start_layout_index]->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level_file_zero = i + 1;
    }
  }

  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + 1)
  {
    ret = PIDX_global_aggregate(file,
                                file->nshared_agg_id,
                                file->idx_d->nshared_agg_buffer,
                                var_grp->nshared_block_layout_by_level, var_grp->nshared_block_layout,
                                file->comm,
                                start_var_index, start_index, 1,
                                var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index,
                                var_grp->nshared_layout_count,
                                agg_io_level_non_shared,
                                file->idx_d->aggregator_multiplier, 1);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = PIDX_global_aggregate(file,
                                file->shared_agg_id,
                                file->idx_d->shared_agg_buffer,
                                var_grp->shared_block_layout_by_level, var_grp->shared_block_layout,
                                file->comm,
                                start_var_index, start_index, 0,
                                var_grp->shared_start_layout_index, var_grp->shared_end_layout_index,
                                var_grp->shared_layout_count,
                                agg_io_level_shared,
                                file->idx_d->aggregator_multiplier, 0);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

#if 1
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (/*file->idx_d->var_pipe_length + */1))
  {

    create_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);
    create_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);


    ret = PIDX_global_async_io(file, file->nshared_io_id,
                               file->idx_d->nshared_agg_buffer,
                               var_grp->nshared_block_layout_by_level,
                               file->idx_d->fp_non_shared,
                               file->idx_d->request_non_shared,
                               start_var_index, start_index, 1,
                               var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index,
                               var_grp->nshared_layout_count,
                               agg_io_level_non_shared, 0, file->idx_d->async_io);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = PIDX_global_async_io(file, file->shared_io_id,
                               file->idx_d->shared_agg_buffer,
                               var_grp->shared_block_layout_by_level,
                               file->idx_d->fp_shared,
                               file->idx_d->request_shared,
                               start_var_index, start_index, 0,
                               var_grp->shared_start_layout_index, var_grp->shared_end_layout_index,
                               var_grp->shared_layout_count,
                               agg_io_level_shared, 0, file->idx_d->async_io);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (file->idx_d->file_zero == 1)
    {
      ret = PIDX_global_aggregate(file, file->f0_agg_id,
                                  file->idx_d->f0_agg_buffer,
                                  var_grp->f0_block_layout_by_level,
                                  var_grp->f0_block_layout,
                                  file->global_comm,
                                  start_var_index, start_index, 0,
                                  var_grp->f0_start_layout_index, var_grp->f0_end_layout_index,
                                  var_grp->f0_layout_count,
                                  agg_io_level_file_zero, 1, 2);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      create_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);

      ret = PIDX_global_async_io(file, file->f0_io_id,
                                 file->idx_d->f0_agg_buffer,
                                 var_grp->f0_block_layout_by_level,
                                 file->idx_d->fp_file_zero, file->idx_d->request_file_zero,
                                 start_var_index, start_index, 0,
                                 var_grp->f0_start_layout_index, var_grp->f0_end_layout_index,
                                 var_grp->f0_layout_count,
                                 agg_io_level_file_zero, 1, file->idx_d->async_io);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);

      destroy_file_zero_ids_and_buffers(file, start_index, var_grp->f0_start_layout_index, var_grp->f0_end_layout_index, agg_io_level_file_zero);
    }

    if (file->idx_d->async_io == 1)
    {
      wait_and_destroy_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);
      wait_and_destroy_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);
    }

    free(file->idx_d->status_shared);
    free(file->idx_d->request_shared);
    free(file->idx_d->fp_shared);
    free(file->idx_d->status_non_shared);
    free(file->idx_d->request_non_shared);
    free(file->idx_d->fp_non_shared);

    destroy_non_shared_ids_and_buffers(file, start_index, var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index, agg_io_level_non_shared);

    destroy_shared_ids_and_buffers(file, start_index, var_grp->shared_start_layout_index, var_grp->shared_end_layout_index, agg_io_level_shared);
  }
#endif

  time->buffer_cleanup_start = PIDX_get_time();


  free(var_grp->rank_buffer);
  ret = destroy_hz_buffers(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = destroy_agg_io_buffer(file, group_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = partition_destroy(file, group_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = partition_communicator_destroy(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->buffer_cleanup_end = PIDX_get_time();


  time->EX = PIDX_get_time();
  return PIDX_success;
}



PIDX_return_code PIDX_hybrid_idx_io_finalize(PIDX_hybrid_idx_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}
