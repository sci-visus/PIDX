#include "../PIDX_inc.h"
#include "restructure.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"
#include "timming.h"

static int hz_from_file_zero = 0, hz_from_shared = 0, hz_from_non_shared = 0;
static int hz_to_file_zero = 0, hz_to_shared = 0, hz_to_non_shared = 0;
static int agg_l_nshared = 0, agg_l_shared = 0, agg_l_f0 = 0;
static PIDX_return_code init(PIDX_io file, int gi);
static PIDX_return_code find_agg_level(PIDX_io file, int gi);
static PIDX_return_code select_io_mode(PIDX_io file, int gi);

PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;

  int start_index = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  // Calculate bounds with compression and
  // populate rank buffer
  //time->idx_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_init_end = MPI_Wtime();


  // Restructuring the grid into power two blocks
  // After this step every process has got a power two block
  // 15 x 31 x 10 ---> 16 x 32 x 16
  //time->idx_rst_start = MPI_Wtime();
  ret = restructure_setup(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_rst_end = MPI_Wtime();


  // Calculates the number of partititons
  //time->idx_partition_start = PIDX_get_time();
  ret = find_partition_count(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  // Setting up for partition
  // group processes (assign same color) to
  // processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_partition_end = PIDX_get_time();


  // calculates maxh and bitstring
  //time->idx_bit_string_start = PIDX_get_time();
  ret = populate_bit_string(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  // selects levels based on mode and maxh
  select_io_mode(file, gi);
  //time->idx_bit_string_end = PIDX_get_time();


  // Reorders data using HZ index scheme
  // Chunks data and compresses using ZFP
  //time->idx_hz_start = PIDX_get_time();
  // TODO
  //ret = create_hz_buffers(file, svi, evi);
  //if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_hz_end = PIDX_get_time();


  // Splits the global communicator into local communicators
  // Essential for MPI scaling
  //time->idx_comm_create_start = PIDX_get_time();
  ret = create_local_comm(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_comm_create_end = PIDX_get_time();


  // Populates the idx block layout
  // individually for file zero, shared and non-sharef file
  //time->idx_layout_start = PIDX_get_time();
  ret = populate_block_layouts(file, gi, svi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_LOCAL_PARTITION_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_layout_end = PIDX_get_time();


  // Creates the file heirarchy
  // Also writes the header info for all binary files
  //time->header_write_start = PIDX_get_time();
  ret = write_headers(file, gi, svi, evi, 0);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->header_write_end = PIDX_get_time();


  // Creates the agg and io ids
  //time->agg_buffer_start = PIDX_get_time();
  ret = create_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->agg_buffer_end = PIDX_get_time();


  // Performs data aggregation
  //time->idx_agg_start = PIDX_get_time();
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 0, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  //time->idx_agg_end = PIDX_get_time();

  // Performs data io
  //time->idx_io_start = PIDX_get_time();
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    create_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);

    ret = data_io(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    wait_and_destroy_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);
    finalize_aggregation(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared);
  }
  //time->idx_io_end = PIDX_get_time();


  // Cleanup all buffers nd ids
  //time->buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);

  ret = destroy_local_comm(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = hz_encode_cleanup(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->buffer_cleanup_end = PIDX_get_time();

  return PIDX_success;
}


PIDX_return_code PIDX_local_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_return_code ret;
  int start_index = 0;

  PIDX_time time = file->idx_d->time;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  // Calculate bounds with compression and
  // populate rank buffer
  //time->idx_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_init_end = MPI_Wtime();


  // Restructuring the grid into power two blocks
  // After this step every process has got a power two block
  // 15 x 31 x 10 ---> 16 x 32 x 16
  //time->idx_rst_start = MPI_Wtime();
  ret = restructure_setup(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  //ret = restructure(file, gi, svi, evi);
  //if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_rst_end = MPI_Wtime();



  // calculates maxh and bitstring
  //time->idx_bit_string_start = PIDX_get_time();
  ret = populate_bit_string(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  // selects levels based on mode and maxh
  select_io_mode(file, gi);
  //time->idx_bit_string_end = PIDX_get_time();



  // Calculates the number of partititons
  //time->idx_partition_start = PIDX_get_time();
  ret = find_partition_count(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /*
  // Setting up for partition
  // group processes (assign same color) to
  // processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_partition_end = PIDX_get_time();
  */



  // Reorders data using HZ index scheme
  // Chunks data and compresses using ZFP
  //time->idx_hz_start = PIDX_get_time();
  ret = hz_encode_setup(file, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_hz_end = PIDX_get_time();


  /*
  // Splits the global communicator into local communicators
  // Essential for MPI scaling
  //time->idx_comm_create_start = PIDX_get_time();
  ret = create_local_comm(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_comm_create_end = PIDX_get_time();
  */


  // Populates the idx block layout
  // individually for file zero, shared and non-sharef file
  //time->idx_layout_start = PIDX_get_time();
  ret = populate_block_layouts(file, gi, svi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_LOCAL_PARTITION_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_layout_end = PIDX_get_time();


  /*
   // Creates the file heirarchy
  // Also writes the header info for all binary files
  //time->header_write_start = PIDX_get_time();
  ret = write_headers(file, gi, svi, evi, 0);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->header_write_end = PIDX_get_time();
  */


  // Creates the agg and io ids
  //time->agg_buffer_start = PIDX_get_time();
  ret = create_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->agg_buffer_end = PIDX_get_time();

  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 1, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }


  // Performs data io
  //time->idx_io_start = PIDX_get_time();
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    //create_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);

    ret = data_io(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    //wait_and_destroy_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);
  }
  //time->idx_io_end = PIDX_get_time();

  // Performs data aggregation
  //time->idx_agg_start = PIDX_get_time();
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 2, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  finalize_aggregation(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared);
  //time->idx_agg_end = PIDX_get_time();

  //hz_encode(file, svi, evi, PIDX_READ);
  hz_encode(file, PIDX_READ);
  restructure(file, PIDX_READ);

  // Cleanup all buffers nd ids
  //time->buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);

  ret = destroy_local_comm(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = hz_encode_cleanup(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->buffer_cleanup_end = PIDX_get_time();

  //time->EX = PIDX_get_time();

  return PIDX_success;

}


static PIDX_return_code init(PIDX_io file, int gi)
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

  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, file->idx_d->perm_layout_count);
    file->one_time_initializations = 1;
  }

  return PIDX_success;
}


static PIDX_return_code find_agg_level(PIDX_io file, int gi)
{
  int i = 0;
  int no_of_aggregators = 0;
  int nprocs = 1;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  MPI_Comm_size(file->comm, &nprocs);
  if (file->idx->enable_agg == 0)
    agg_l_nshared = var_grp->nshared_start_layout_index;
  else
  {
    for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->nshared_block_layout_by_level[i - var_grp->nshared_start_layout_index]->efc;
      if (no_of_aggregators <= nprocs)
        agg_l_nshared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    agg_l_shared = var_grp->shared_start_layout_index;
  else
  {
    for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->shared_block_layout_by_level[i - var_grp->shared_start_layout_index]->efc;
      if (no_of_aggregators <= nprocs)
        agg_l_shared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    agg_l_f0 = var_grp->f0_start_layout_index;
  else
  {
    for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->f0_block_layout_by_level[i - var_grp->f0_start_layout_index]->efc;
      if (no_of_aggregators <= nprocs)
        agg_l_f0 = i + 1;
    }
  }

  return PIDX_success;
}


static PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  //file->idx->io_type = 2;
  if (file->idx->io_type == PIDX_LOCAL_PARTITION_IDX_IO)
  {
    hz_from_file_zero = 0;
    hz_to_file_zero =  idx->shared_block_level;

    hz_from_shared = idx->shared_block_level;
    hz_to_shared =  idx->total_partiton_level;

    hz_from_non_shared = idx->total_partiton_level;
    hz_to_non_shared =  idx->maxh;
  }
  else if (file->idx->io_type == PIDX_GLOBAL_PARTITION_IDX_IO)
  {
    hz_from_file_zero = 0;
    hz_to_file_zero =  0;

    hz_from_shared = 0;
    hz_to_shared =  idx->total_partiton_level;

    hz_from_non_shared = idx->total_partiton_level;
    hz_to_non_shared =  idx->maxh;
  }
  else if (file->idx->io_type == PIDX_IDX_IO)
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
