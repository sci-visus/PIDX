/*
 * 1. One global dataset.
 * 2. Partitions in local index space.
 * 3. Partitions in global index space.
 * 4. Partition divided into non-shared, shared, file 0 (global index)
 * 5. Partition divided into non-shared, shared, file 0 (local index)
 */


#include "../PIDX_io.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"

static int hz_from_file_zero = 0, hz_from_shared = 0, hz_from_non_shared = 0;
static int hz_to_file_zero = 0, hz_to_shared = 0, hz_to_non_shared = 0;
static int agg_io_level_non_shared = 0, agg_io_level_shared = 0, agg_io_level_file_zero = 0;
static PIDX_return_code init(PIDX_hybrid_idx_io file, int gi);
static PIDX_return_code find_partition_count(PIDX_hybrid_idx_io file, int gi);
static PIDX_return_code populate_block_layouts(PIDX_hybrid_idx_io file, int gi, int svi, int gvi);
static PIDX_return_code data_aggregate(PIDX_hybrid_idx_io file, int gi, int svi, int start_index);
static PIDX_return_code data_io(PIDX_hybrid_idx_io file, int gi, int svi, int start_index);
static PIDX_return_code create_async_buffers(PIDX_hybrid_idx_io file, int gi);
static PIDX_return_code wait_and_destroy_async_buffers(PIDX_hybrid_idx_io file, int gi);
static PIDX_return_code finalize_aggregation(PIDX_hybrid_idx_io file, int start_index, int gi);

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



PIDX_return_code PIDX_hybrid_idx_write(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  time->SX = PIDX_get_time();

  PIDX_return_code ret;
  int start_index = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  ret = init(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Restructuring the grid into power two blocks
  // After this step every process has got a power two block
  // 15 x 31 x 10 ---> 16 x 32 x 16
  ret = restructure(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  time->partition_start_time = PIDX_get_time();

  // Calculate the number of partititons
  ret = find_partition_count(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Setting up for partition
  // group processes (assign same color) to
  // processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_end_time = PIDX_get_time();


  ret = populate_bit_string(file);
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
  ret = create_hz_buffers(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->hz_e_time = PIDX_get_time();

  ret = create_local_comm(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  time->populate_idx_start_time_s = PIDX_get_time();
  ret = populate_block_layouts(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->populate_idx_end_time_s = PIDX_get_time();



  ret = write_headers(file, gi, svi, evi, 0);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = create_agg_io_buffer(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, svi, start_index);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

#if 1
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    create_async_buffers(file, gi);

    ret = data_io(file, gi, svi, start_index);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    wait_and_destroy_async_buffers(file, gi);

    finalize_aggregation(file, gi, start_index);
  }
#endif

  time->buffer_cleanup_start = PIDX_get_time();

  ret = destroy_agg_io_buffer(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = delete_block_layout(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  free(var_grp->rank_buffer);

  ret = destroy_local_comm(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = destroy_hz_buffers(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->buffer_cleanup_end = PIDX_get_time();


  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->EX = PIDX_get_time();
  return PIDX_success;
}



PIDX_return_code PIDX_hybrid_idx_io_finalize(PIDX_hybrid_idx_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}


static PIDX_return_code data_aggregate(PIDX_hybrid_idx_io file, int gi, int svi, int start_index)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  ret = PIDX_global_aggregate(file, file->f0_agg_id,
                              idx->f0_agg_buffer,
                              var_grp->f0_block_layout_by_level,
                              var_grp->f0_block_layout,
                              file->global_comm,
                              svi, start_index,
                              var_grp->f0_start_layout_index,
                              agg_io_level_file_zero, 1, 2);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = PIDX_global_aggregate(file,
                              file->nshared_agg_id,
                              idx->nshared_agg_buffer,
                              var_grp->nshared_block_layout_by_level, var_grp->nshared_block_layout,
                              file->comm,
                              svi, start_index,
                              var_grp->nshared_start_layout_index,
                              agg_io_level_non_shared,
                              idx->aggregator_multiplier, 1);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = PIDX_global_aggregate(file,
                              file->shared_agg_id,
                              idx->shared_agg_buffer,
                              var_grp->shared_block_layout_by_level, var_grp->shared_block_layout,
                              file->comm,
                              svi, start_index,
                              var_grp->shared_start_layout_index,
                              agg_io_level_shared,
                              idx->aggregator_multiplier, 0);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}

static PIDX_return_code find_partition_count(PIDX_hybrid_idx_io file, int gi)
{

  int d = 0;
  idx_dataset_derived_metadata idx = file->idx_d;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  // calculate number of partitions
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
      file->idx_d->partition_count[d] = file->idx->bounds[d] / file->idx_d->partition_size[d];
      if (file->idx->bounds[d] % file->idx_d->partition_size[d] != 0)
          file->idx_d->partition_count[d]++;

      file->idx_d->partition_count[d] = pow(2, (int)ceil(log2(file->idx_d->partition_count[d])));
  }

  idx->shared_block_level = (int)log2(idx->partition_count[0] * idx->partition_count[1] * idx->partition_count[2]) + file->idx->bits_per_block + 1;
  if (idx->shared_block_level >= idx->maxh)
    idx->shared_block_level = idx->maxh;

  int partion_level = (int) log2(idx->partition_count[0] * idx->partition_count[1] * idx->partition_count[2]);
  idx->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (idx->total_partiton_level >= idx->maxh)
    idx->total_partiton_level = idx->maxh;

  if (idx->io_mode == 0)
  {
    hz_from_file_zero = 0;
    hz_to_file_zero =  idx->shared_block_level;

    hz_from_shared = idx->shared_block_level;
    hz_to_shared =  idx->total_partiton_level;

    hz_from_non_shared = idx->total_partiton_level;
    hz_to_non_shared =  idx->maxh;
  }
  else if (idx->io_mode == 1)
  {
    hz_from_file_zero = 0;
    hz_to_file_zero =  0;

    hz_from_shared = 0;
    hz_to_shared =  idx->total_partiton_level;

    hz_from_non_shared = idx->total_partiton_level;
    hz_to_non_shared =  idx->maxh;
  }
  else if (idx->io_mode == 2)
  {
    hz_from_file_zero = 0;
    hz_to_file_zero =  0;

    hz_from_shared = 0;
    hz_to_shared =  0;

    hz_from_non_shared = 0;
    hz_to_non_shared =  idx->maxh;
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

  return PIDX_success;
}

static PIDX_return_code init(PIDX_hybrid_idx_io file, int gi)
{
  int d = 0;
  int grank = 0, gnprocs = 1;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    if (file->idx->bounds[d] % file->idx->chunk_size[d] == 0)
      file->idx->chunked_bounds[d] = (int) file->idx->bounds[d] / file->idx->chunk_size[d];
    else
      file->idx->chunked_bounds[d] = (int) (file->idx->bounds[d] / file->idx->chunk_size[d]) + 1;
  }

  MPI_Comm_rank(file->global_comm, &grank);
  MPI_Comm_size(file->global_comm, &gnprocs);

  var_grp->rank_buffer = malloc(gnprocs * sizeof(*var_grp->rank_buffer));
  memset(var_grp->rank_buffer, 0, gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&grank, 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->global_comm);

  return PIDX_success;
}


static PIDX_return_code populate_block_layouts(PIDX_hybrid_idx_io file, int gi, int svi, int evi)
{
  int ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  create_file_zero_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero);
  ret = populate_idx_block_layout(file,
                             var_grp->f0_block_layout,
                             var_grp->f0_block_layout_by_level,
                             var_grp->f0_start_layout_index,
                             var_grp->f0_end_layout_index,
                             var_grp->f0_layout_count,
                             gi,
                             svi, evi,
                             hz_from_file_zero, hz_to_file_zero);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  create_shared_block_layout(file, gi, hz_from_shared, hz_to_shared);
  ret = populate_idx_block_layout(file,
                             var_grp->shared_block_layout, var_grp->shared_block_layout_by_level,
                             var_grp->shared_start_layout_index, var_grp->shared_end_layout_index,
                             var_grp->shared_layout_count,
                             gi,  svi, evi,
                             hz_from_shared, hz_to_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  create_shared_block_layout(file, gi, hz_from_non_shared, hz_to_non_shared);
  ret = populate_idx_block_layout(file,
                             var_grp->nshared_block_layout, var_grp->nshared_block_layout_by_level,
                             var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index,
                             var_grp->nshared_layout_count,
                             gi,  svi, evi,
                             hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int i = 0;
  int no_of_aggregators = 0;
  int nprocs = 1;

  MPI_Comm_size(file->comm, &nprocs);
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


  return PIDX_success;
}


static PIDX_return_code create_async_buffers(PIDX_hybrid_idx_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  create_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);
  create_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);
  create_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);

  return PIDX_success;
}


static PIDX_return_code data_io(PIDX_hybrid_idx_io file, int gi, int svi, int start_index)
{
  int ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  ret = PIDX_global_async_io(file, file->f0_io_id,
                             idx->f0_agg_buffer,
                             var_grp->f0_block_layout_by_level,
                             idx->fp_file_zero, idx->request_file_zero,
                             svi, start_index, 0,
                             var_grp->f0_start_layout_index, var_grp->f0_end_layout_index,
                             var_grp->f0_layout_count,
                             agg_io_level_file_zero, 1);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = PIDX_global_async_io(file, file->nshared_io_id,
                             idx->nshared_agg_buffer,
                             var_grp->nshared_block_layout_by_level,
                             idx->fp_non_shared,
                             idx->request_non_shared,
                             svi, start_index, 1,
                             var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index,
                             var_grp->nshared_layout_count,
                             agg_io_level_non_shared, 0);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = PIDX_global_async_io(file, file->shared_io_id,
                             idx->shared_agg_buffer,
                             var_grp->shared_block_layout_by_level,
                             idx->fp_shared,
                             idx->request_shared,
                             svi, start_index, 0,
                             var_grp->shared_start_layout_index, var_grp->shared_end_layout_index,
                             var_grp->shared_layout_count,
                             agg_io_level_shared, 0);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}

static PIDX_return_code wait_and_destroy_async_buffers(PIDX_hybrid_idx_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  wait_and_destroy_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);
  wait_and_destroy_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);
  wait_and_destroy_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);

  return PIDX_success;
}


static PIDX_return_code finalize_aggregation(PIDX_hybrid_idx_io file, int start_index, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  destroy_file_zero_ids_and_buffers(file, start_index, var_grp->f0_start_layout_index, var_grp->f0_end_layout_index, agg_io_level_file_zero);

  destroy_non_shared_ids_and_buffers(file, start_index, var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index, agg_io_level_non_shared);

  destroy_shared_ids_and_buffers(file, start_index, var_grp->shared_start_layout_index, var_grp->shared_end_layout_index, agg_io_level_shared);

  return PIDX_success;
}
