#include "../PIDX_inc.h"
#include "restructure.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"


/// Raw Write Steps
/****************************************************************
*  Step 1                                                       *
*  write headers               ---->  init_raw_headers_layout() *
*                                                               *
*  Step 2                                                       *
*  Restructure initialization  ---->  restructure_init()        *
*                                                               *
*  Step 3                                                       *
*  Restrucure                  ---->  restructure()             *
*                                                               *
*  Step 4                                                       *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 5                                                       *
*  Cleanup                                                      *
*****************************************************************/

PIDX_return_code PIDX_raw_write(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;

  /// Step 1 [Start]
  time->raw_write_header_start = MPI_Wtime();
  /// Creates the directory heirarchy
  ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->raw_write_header_end = MPI_Wtime();
  /// Step 1 [End]


  /// Step 2 [Start]
  time->raw_write_rst_init_start = MPI_Wtime();
  /// initialization of restructuring phase
  /// calculate the imposed box size and relevant metadata
  ret = restructure_init(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->raw_write_rst_init_end = MPI_Wtime();
  /// Step 2 [End]


  /// Step 3 [Start]
  time->raw_write_rst_start = MPI_Wtime();
  /// Restructuring the grid into power two blocks
  /// After this step every process has got a power two block
  /// 15 x 31 x 10 ---> 16 x 32 x 16
  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->raw_write_rst_end = MPI_Wtime();
  /// Step 3 [End]


  /// Step 4 [Start]
  time->raw_write_io_start = MPI_Wtime();
  /// Writes data out
  ret = restructure_io(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->raw_write_io_end = MPI_Wtime();
  /// Step 4 [End]


  /// Step 5 [Start]
  time->raw_write_buffer_cleanup_start = PIDX_get_time();
  /// Cleanup all buffers and ids
  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = write_and_close_raw_headers(file, file->idx->filename);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  time->raw_write_buffer_cleanup_end = PIDX_get_time();
  /// Step 5 [End]

  return PIDX_success;
}



/// Raw Read Steps
/****************************************************************
*  Step 1                                                       *
*  Restructure initialization  ---->  restructure_init()        *
*                                                               *
*  Step 2                                                       *
*  File I/O                    ---->  data_io()                 *
*                                                               *
*  Step 3                                                       *
*  Restrucure                  ---->  restructure()             *
*                                                               *
*  Step 4                                                       *
*  Cleanup                                                      *
*****************************************************************/

PIDX_return_code PIDX_raw_read(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;
  int nprocs = 1;
  int rst_case_type = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  int patch_count = var0->sim_patch_count;
  int max_patch_count = 0;

  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->comm);
  if (max_patch_count > 1)
    rst_case_type = 1;
  MPI_Comm_size(file->comm, &nprocs);

  ///  Find if this is amulti patch per process read or not
  ///  If multi patch per process then use raw forced io
  ///  else use the regular read routine
  if (file->idx_d->data_core_count == nprocs && rst_case_type == 0)
  {
    time->raw_read_rst_init_start = MPI_Wtime();
    /// initialization of restructuring phase
    /// calculate the imposed box size and relevant metadata
    ret = restructure_init(file, gi, svi, evi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    time->raw_read_rst_init_end = MPI_Wtime();
    /// Step 1 [End]

    /// Step 4 [Start]
    time->raw_read_io_start = MPI_Wtime();
    /// Writes data out
    ret = restructure_io(file, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    time->raw_read_io_end = MPI_Wtime();
    /// Step 4 [End]

    /// Step 2 [Start]
    time->raw_read_rst_start = MPI_Wtime();
    /// Restructuring the grid into power two blocks
    /// After this step every process has got a power two block
    /// 15 x 31 x 10 ---> 16 x 32 x 16
    ret = restructure(file, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    time->raw_read_rst_end = MPI_Wtime();
    /// Step 2 [End]

    /// Step 5 [Start]
    time->raw_read_buffer_cleanup_start = PIDX_get_time();
    /// Cleanup all buffers and ids
    ret = restructure_cleanup(file, gi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    time->raw_read_buffer_cleanup_end = PIDX_get_time();
    /// Step 5 [End]
  }
  else
  {
    /// Step 1 [Start]
    time->raw_forced_read_start = PIDX_get_time();
    ret = restructure_forced_read(file, svi, evi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    time->raw_forced_read_end = PIDX_get_time();
    /// Step 1 [End]
  }

  return PIDX_success;
}
