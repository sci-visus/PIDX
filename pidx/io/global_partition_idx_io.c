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



/// Global Partitioned IDX Write Steps
/****************************************************************
*  Step 1                                                       *
*  Restructure initialization  ---->  restructure_setup()        *
*                                                               *
*  Step 2                                                       *
*  Restrucure                  ---->  restructure()             *
*                                                               *
*  Step 3                                                       *
*  Partition setup           ---->    partition_setup()         *
*                                                               *
*  Step 4                                                       *
*  create partition (comm)    ---->  create_local_comm()        *
*                                                               *
*  Step 5                                                       *
*  create bit string           ---->  populate_bit_string()     *
*                                                               *
*  Step 6                                                       *
*  HZ encode initialization    ---->  hz_encode_setup()        *
*                                                               *
*  Step 7                                                       *
*  HZ encoding                 ---->  hz_encode()     *
*                                                               *
*  Step 8                                                       *
*  create block layout         ---->  populate_idx_layout()     *
*                                                               *
*  Step 9                                                       *
*  write headers               ---->  write_headers()           *
*                                                               *
*  Step 10                                                      *
*  Aggregation Initialization  ---->  data_aggregate(init mode) *
*                                                               *
*  Step 11                                                      *
*  Aggregation                 ---->  data_aggregate(perform)   *
*                                                               *
*  Step 12                                                      *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 13                                                      *
*  Cleanup                                                      *
*****************************************************************/

PIDX_return_code PIDX_global_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;

  int start_index = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  /// Step 0 [Start]
  /// Calculate bounds with compression and
  /// populate rank buffer
  time->global_idx_write_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_init_end = MPI_Wtime();
  /// Step 0 [End]


  /// Step 1 [Start]
  time->global_idx_write_rst_init_start = MPI_Wtime();
  /// initialization of restructuring phase
  /// calculate the imposed box size and relevant metadata
  ret = restructure_setup(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_rst_init_end = MPI_Wtime();
  /// Step 1 [End]


  /// Step 2 [Start]
  time->global_idx_write_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_rst_end = MPI_Wtime();
  /// Step 2 [End]


  /// Step 3 [Start]
  time->global_idx_write_partition_setup_start = MPI_Wtime();
  /// Calculates the number of partititons
  ret = find_partition_count(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Setting up for partition
  /// group processes (assign same color) to
  /// processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_partition_setup_end = PIDX_get_time();
  /// Step 3 [End]


  /// Step 4 [Start]
  time->global_idx_write_comm_create_start = PIDX_get_time();
  /// Splits the global communicator into local communicators
  /// Essential for MPI scaling
  ret = create_local_comm(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  if (var_grp->variable[svi]->patch_group_count == 0)
    goto cleanup;
  time->global_idx_write_comm_create_end = PIDX_get_time();
  /// Step 4 [End]


  /// Step 5 [Start]
  time->global_idx_write_bitstring_start = PIDX_get_time();
  /// calculates maxh and bitstring
  /// bitstring depends on restrcuring box size
  ret = populate_bit_string(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// selects levels based on mode and maxh
  select_io_mode(file, gi);
  time->global_idx_write_bitstring_end = PIDX_get_time();
  /// Step 5 [End]


  /// Step 6 [Start]
  time->global_idx_write_hz_init_start = PIDX_get_time();
  /// Encoding initialization and meta data formation
  ret = hz_encode_setup(file, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_hz_init_end = PIDX_get_time();
  /// Step 6 [End]


  /// Step 7 [Start]
  time->global_idx_write_hz_start = PIDX_get_time();
  /// Reorders data using HZ index scheme
  /// Chunks data and compresses using ZFP
  //ret = hz_encode(file, svi, evi, PIDX_WRITE);
  ret = hz_encode(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_hz_end = PIDX_get_time();
  /// Step 7 [End]


  /// Step 8 [Start]
  time->global_idx_write_layout_start = PIDX_get_time();
  /// Populates the idx block layout
  /// individually for file zero, shared and non-sharef file
  ret = populate_block_layouts(file, gi, svi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_GLOBAL_PARTITION_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_layout_end = PIDX_get_time();
  /// Step 8 [End]


  /// Step 9 [Start]
  time->global_idx_write_header_start = PIDX_get_time();
  /// Creates the file heirarchy
  /// Also writes the header info for all binary files
  ret = write_headers(file, gi, svi, evi, 0);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_write_header_end = PIDX_get_time();
  /// Step 9 [End]


  /// Step 10 [Start]
  /// Aggregation initialization
  time->global_idx_write_agg_init_start = PIDX_get_time();
  /// Creates the agg and io ids
  ret = create_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Perform aggregation initialization creation of relevant meta data
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 1, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->global_idx_write_agg_init_end = PIDX_get_time();
  /// Step 10 [End]


  /// Step 11 [Start]
  time->global_idx_write_agg_start = PIDX_get_time();
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 2, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->global_idx_write_agg_end = PIDX_get_time();
  /// Step 11 [End]


  /// Step 12 [Start]
  time->global_idx_write_io_start = PIDX_get_time();
  /// performs actual file io
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    create_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);

    ret = data_io(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    wait_and_destroy_async_buffers(file, gi, agg_l_f0, agg_l_shared, agg_l_nshared);
    finalize_aggregation(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared);
  }
  time->global_idx_write_io_end = PIDX_get_time();
  /// Step 12 [End]


  /// Step 13 [Start]
  /// Cleanup all buffers and ids
  time->global_idx_write_buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = hz_encode_cleanup(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

cleanup:
  ret = destroy_local_comm(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);
  time->global_idx_write_buffer_cleanup_end = PIDX_get_time();
  /// Step 13 [End]

  return PIDX_success;
}




/// Global Partitioned IDX Write Steps
/****************************************************************
*  Step 1                                                       *
*  Restructure initialization  ---->  restructure_setup()        *
*                                                               *
*  Step 2                                                       *
*  Partition setup           ---->    partition_setup()         *
*                                                               *
*  Step 3                                                       *
*  create partition (comm)    ---->  create_local_comm()        *
*                                                               *
*  Step 4                                                       *
*  create bit string           ---->  populate_bit_string()     *
*                                                               *
*  Step 5                                                       *
*  HZ encode initialization    ---->  hz_encode_setup()        *
*                                                               *
*  Step 6                                                       *
*  create block layout         ---->  populate_idx_layout()     *
*                                                               *
*  Step 7                                                       *
*  Aggregation Initialization  ---->  data_aggregate(init mode) *
*                                                               *
*  Step 8                                                       *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 9                                                       *
*  Aggregation                 ---->  data_aggregate(perform)   *
*                                                               *
*  Step 10                                                      *
*  HZ encoding                 ---->  hz_encode()     *
*                                                               *
*  Step 11                                                      *
*  Cleanup                                                      *
*                                                               *
*  Step 12                                                      *
*  Restrucure                  ---->  restructure()             *
*****************************************************************/

PIDX_return_code PIDX_global_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;

  int start_index = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  /// Step 0 [Start]
  /// Calculate bounds with compression and
  /// populate rank buffer
  time->global_idx_read_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_init_end = MPI_Wtime();
  /// Step 0 [End]


  /// Step 1 [Start]
  time->global_idx_read_rst_init_start = MPI_Wtime();
  /// initialization of restructuring phase
  /// calculate the imposed box size and relevant metadata
  ret = restructure_setup(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_rst_init_end = MPI_Wtime();
  /// Step 1 [End]


  /// Step 2 [Start]
  time->global_idx_read_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_rst_end = MPI_Wtime();
  /// Step 2 [End]


  /// Step 3 [Start]
  time->global_idx_read_comm_create_start = PIDX_get_time();
  /// Splits the global communicator into local communicators
  /// Essential for MPI scaling
  ret = create_local_comm(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  if (var_grp->variable[svi]->patch_group_count == 0)
    goto cleanup;
  time->global_idx_read_comm_create_end = MPI_Wtime();
  /// Step 3 [End]



  /// Step 4 [Start]
  time->global_idx_read_bitstring_start = PIDX_get_time();
  /// calculates maxh and bitstring
  /// bitstring depends on restrcuring box size
  ret = populate_bit_string(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// selects levels based on mode and maxh
  select_io_mode(file, gi);
  time->global_idx_read_bitstring_end = PIDX_get_time();
  /// Step 4 [End]


  /// Step 5 [Start]
  time->global_idx_read_hz_init_start = PIDX_get_time();
  /// Encoding initialization and meta data formation
  ret = hz_encode_setup(file, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_hz_init_end = PIDX_get_time();
  /// Step 5 [End]


  /// Step 6 [Start]
  time->global_idx_read_layout_start = PIDX_get_time();
  /// Populates the idx block layout
  /// individually for file zero, shared and non-sharef file
  ret = populate_block_layouts(file, gi, svi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_GLOBAL_PARTITION_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  /// Step 6 [End]


  /// Step 7 [Start]
  time->global_idx_read_agg_init_start = PIDX_get_time();
  /// Aggregation initialization
  ret = create_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 1, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->global_idx_read_agg_init_end = PIDX_get_time();
  /// Step 7 [End]


  /// Step 8 [Start]
  time->global_idx_read_io_start = PIDX_get_time();
  /// performs actual file io
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    ret = data_io(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->global_idx_read_io_end = PIDX_get_time();
  /// Step 8 [End]


  /// Step 9 [Start]
  time->global_idx_read_agg_start = PIDX_get_time();
  /// Performs data aggregation
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared, 2, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    finalize_aggregation(file, gi, start_index, agg_l_f0, agg_l_shared, agg_l_nshared);
  }
  time->global_idx_read_agg_end = PIDX_get_time();
  /// Step 9 [End]


  /// Step 10 [Start]
  time->global_idx_read_hz_start = PIDX_get_time();
  /// Reorders data using HZ index scheme
  /// Chunks data and compresses using ZFP
  //ret = hz_encode(file, svi, evi, PIDX_READ);
  ret = hz_encode(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_hz_end = PIDX_get_time();
  /// Step 10 [End]


  /// Step 11 [Start]
  /// Cleanup all buffers and ids
  time->global_idx_read_buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = hz_encode_cleanup(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = destroy_local_comm(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->global_idx_read_buffer_cleanup_end = PIDX_get_time();
  /// Step 11 [End]


  /// Step 12 [Start]
  time->global_idx_read_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  cleanup:
  ret = restructure(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);
  time->global_idx_read_rst_end = MPI_Wtime();
  /// Step 12 [End]

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

  hz_from_file_zero = 0;
  hz_to_file_zero =  0;

  hz_from_shared = 0;
  hz_to_shared =  idx->total_partiton_level;

  hz_from_non_shared = idx->total_partiton_level;
  hz_to_non_shared =  idx->maxh;

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
