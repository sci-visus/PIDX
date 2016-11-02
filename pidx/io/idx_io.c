#include "../PIDX_inc.h"
#include "restructure.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"

static int hz_from_non_shared = 0;
static int hz_to_non_shared = 0;
static int agg_l_nshared = 0;
static PIDX_return_code init(PIDX_io file, int gi);
static PIDX_return_code find_agg_level(PIDX_io file, int gi);
static PIDX_return_code select_io_mode(PIDX_io file, int gi);

/// IDX Write Steps
/****************************************************************
*  Step 1                                                       *
*  Restructure initialization  ---->  restructure_init()        *
*                                                               *
*  Step 2                                                       *
*  Restrucure                  ---->  restructure()             *
*                                                               *
*  Step 3                                                       *
*  create bit string           ---->  populate_bit_string()     *
*  create block layout         ---->  populate_idx_layout()     *
*                                                               *
*  Step 4                                                       *
*  write headers               ---->  write_headers()           *
*                                                               *
*  Step 5                                                       *
*  HZ encode initialization    ---->  setup_hz_buffers()        *
*                                                               *
*  Step 6                                                       *
*  HZ encoding                 ---->  populate_hz_buffers()     *
*                                                               *
*  Step 7                                                       *
*  Aggregation Initialization  ---->  data_aggregate(init mode) *
*                                                               *
*  Step 8                                                       *
*  Aggregation                 ---->  data_aggregate(perform)   *
*                                                               *
*  Step 9                                                       *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 10                                                      *
*  Cleanup                                                      *
*****************************************************************/


PIDX_return_code PIDX_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int start_index = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  /// Calculate bounds with compression and
  /// populate rank buffer
  time->idx_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_init_end = MPI_Wtime();


  /// Step 1 [Start]
  time->idx_write_rst_init_start = MPI_Wtime();
  /// initialization of restructuring phase
  /// calculate the imposed box size and relevant metadata
  ret = restructure_init(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_rst_init_end = MPI_Wtime();
  /// Step 1 [End]


  /// Step 2 [Start]
  time->idx_write_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_rst_end = MPI_Wtime();
  /// Step 2 [End]


  /// Step 3 [Start]
  time->idx_write_layout_start = PIDX_get_time();
  /// calculates maxh and bitstring
  /// bitstring depends on restrcuring box size
  ret = populate_bit_string(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// selects layout levels based on maxh
  select_io_mode(file, gi);

  /// calculates the block layout, given this is pure IDX only non-share block layout is populated
  ret = populate_block_layouts(file, gi, svi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared, PIDX_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_layout_end = PIDX_get_time();
  /// Step 3 [End]


  /// Step 4 [Start]
  time->idx_write_header_start = PIDX_get_time();
  /// finds the filename temlates
  ret = one_time_initialize(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Creates the file heirarchy
  /// Also writes the header info for all binary files
  ret = write_headers(file, gi, svi, evi, PIDX_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_header_end = PIDX_get_time();
  /// Step 4 [End]


  /// Step 5 [Start]
  time->idx_write_hz_init_start = PIDX_get_time();
  /// Encoding initialization and meta data formation
  ret = setup_hz_buffers(file, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_hz_init_end = PIDX_get_time();
  /// Step 5 [End]


  /// Step 6 [Start]
  time->idx_write_hz_start = PIDX_get_time();
  /// Reorders data using HZ index scheme
  /// Chunks data and compresses using ZFP
  ret = populate_hz_buffers(file, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_hz_end = PIDX_get_time();
  /// Step 6 [End]


  /// Step 7 [Start]
  time->idx_write_agg_init_start = PIDX_get_time();
  /// Creates the agg and io ids
  ret = create_agg_io_buffer(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}


  /// Perform aggregation initialization creation of relevant meta data
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, svi, start_index, 0, 0, agg_l_nshared, 1, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->idx_write_agg_init_end = PIDX_get_time();
  /// Step 7 [End]


  /// Step 8 [Start]
  time->idx_write_agg_start = PIDX_get_time();
  /// Performs data aggregation
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, svi, start_index, 0, 0, agg_l_nshared, 2, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->idx_write_agg_end = PIDX_get_time();
  /// Step 8 [End]


  /// Step 9 [Start]
  time->idx_write_io_start = PIDX_get_time();
  /// performs actual file io
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    create_async_buffers(file, gi, 0, 0, agg_l_nshared);

    ret = data_io(file, gi, svi, start_index, 0, 0, agg_l_nshared, PIDX_WRITE);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    wait_and_destroy_async_buffers(file, gi, 0, 0, agg_l_nshared);
    finalize_aggregation(file, gi, start_index, 0, 0, agg_l_nshared);
  }
  time->idx_write_io_end = PIDX_get_time();
  /// Step 9 [End]


  /// Step 10 [Start]
  /// Cleanup all buffers and ids
  time->idx_write_buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);

  ret = destroy_hz_buffers(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_write_buffer_cleanup_end = PIDX_get_time();
  /// Step 10 [End]
#if 0
#endif
  return PIDX_success;
}



/// IDX Read Steps
/****************************************************************
*  Step 1                                                       *
*  Restructure initialization  ---->  restructure_init()        *
*                                                               *
*  Step 2                                                       *
*  create bit string           ---->  populate_bit_string()     *
*  create block layout         ---->  populate_idx_layout()     *
*                                                               *
*  Step 3                                                       *
*  read headers                ---->  write_headers()           *
*                                                               *
*  Step 4                                                       *
*  HZ encode initialization    ---->  setup_hz_buffers()        *
*                                                               *
*  Step 5                                                       *
*  Aggregation Initialization  ---->  data_aggregate(init mode) *
*                                                               *
*  Step 6                                                       *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 7                                                       *
*  Aggregation                 ---->  data_aggregate(perform)   *
*                                                               *
*  Step 8                                                       *
*  HZ encoding                 ---->  populate_hz_buffers()     *
*                                                               *
*  Step 9                                                       *
*  Restrucure                  ---->  restructure()             *
*                                                               *
*  Step 10                                                      *
*  Cleanup                                                      *
*****************************************************************/
PIDX_return_code PIDX_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  int start_index = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  /// Calculate bounds with compression and
  /// populate rank buffer
  time->idx_init_start = MPI_Wtime();
  ret = init(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_init_end = MPI_Wtime();


  /// Step 1 [Start]
  time->idx_read_rst_init_start = MPI_Wtime();
  /// initialization of restructuring phase
  /// calculate the imposed box size and relevant metadata
  ret = restructure_init(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_read_rst_init_end = MPI_Wtime();
  /// Step 1 [End]

  /// Step 2 [Start]
  time->idx_read_layout_start = PIDX_get_time();
  /// calculates maxh and bitstring
  /// bitstring depends on restrcuring box size
  ret = populate_bit_string(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// selects layout levels based on maxh
  select_io_mode(file, gi);

  /// calculates the block layout, given this is pure IDX only non-share block layout is populated
  ret = populate_block_layouts(file, gi, svi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared, PIDX_IDX_IO);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  time->idx_read_layout_end = PIDX_get_time();
  /// Step 2 [End]


  /// Step 3 [Start]
  time->idx_read_header_start = PIDX_get_time();
  /// finds the filename temlates
  ret = one_time_initialize(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  /// Step 3 [End]


  /// Step 4 [Start]
  time->idx_read_hz_init_start = PIDX_get_time();
  /// Encoding initialization and meta data formation
  ret = setup_hz_buffers(file, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_read_hz_init_end = PIDX_get_time();
  /// Step 4 [End]


  /// Step 5 [Start]
  time->idx_read_agg_init_start = PIDX_get_time();
  /// Creates the agg and io ids
  ret = create_agg_io_buffer(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  /// Perform aggregation initialization creation of relevant meta data
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, svi, start_index, 0, 0, agg_l_nshared, 1, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->idx_read_agg_init_end = PIDX_get_time();
  /// Step 5 [End]

#if 1
  /// Step 6 [Start]
  time->idx_read_io_start = PIDX_get_time();
  /// performs actual file io
  for (start_index = svi; start_index < evi; start_index = start_index + (/*idx->var_pipe_length + */1))
  {
    ret = data_io(file, gi, svi, start_index, 0, 0, agg_l_nshared, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->idx_read_io_end = PIDX_get_time();
  /// Step 6 [End]


  /// Step 7 [Start]
  time->idx_read_agg_start = PIDX_get_time();
  /// Performs data aggregation
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = data_aggregate(file, gi, svi, start_index, 0, 0, agg_l_nshared, 2, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  for (start_index = svi; start_index < evi; start_index = start_index + 1)
  {
    ret = finalize_aggregation(file, gi, start_index, 0, 0, agg_l_nshared);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  time->idx_read_agg_end = PIDX_get_time();
  /// Step 7 [End]


  /// Step 8 [Start]
  time->idx_read_hz_start = PIDX_get_time();
  /// Reorders data using HZ index scheme
  /// Chunks data and compresses using ZFP
  ret = populate_hz_buffers(file, svi, evi, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_read_hz_end = PIDX_get_time();
  /// Step 8 [End]


  /// Step 9 [Start]
  time->idx_read_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  ret = restructure(file, PIDX_READ);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_read_rst_end = MPI_Wtime();
  /// Step 9 [End]


  /// Step 10 [Start]
  /// Cleanup all buffers and ids
  time->idx_read_buffer_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = delete_block_layout(file, gi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  free(var_grp->rank_buffer);

  ret = destroy_hz_buffers(file);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->idx_read_buffer_cleanup_end = PIDX_get_time();
  /// Step 10 [End]
#endif
  time->EX = PIDX_get_time();
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

  return PIDX_success;
}


static PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  var_grp->f0_start_layout_index = 0;
  var_grp->f0_end_layout_index = 0;
  var_grp->shared_start_layout_index = 0;
  var_grp->shared_end_layout_index = 0;

  hz_from_non_shared = 0;
  hz_to_non_shared =  idx->maxh;

  return PIDX_success;
}
