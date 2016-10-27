#include "../PIDX_inc.h"
#include "restructure.h"
#include "partition.h"
#include "local_buffer.h"
#include "headers.h"
#include "blocks.h"
#include "hz_buffer.h"
#include "agg_io.h"

PIDX_return_code PIDX_raw_write(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;
  PIDX_return_code ret;

  // Creates the file heirarchy
  // Also writes the header info for all binary files
  //time->header_write_start = PIDX_get_time();
  ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->header_write_end = PIDX_get_time();

  // Restructuring the grid into power two blocks
  // After this step every process has got a power two block
  // 15 x 31 x 10 ---> 16 x 32 x 16
  //time->idx_rst_start = MPI_Wtime();
  ret = restructure_init(file, gi, svi, evi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  ret = restructure_io(file, PIDX_WRITE);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->idx_rst_end = MPI_Wtime();

  // Cleanup all buffers and ids
  //time->buffer_cleanup_start = PIDX_get_time();
  ret = restructure_cleanup(file, gi);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  //time->buffer_cleanup_end = PIDX_get_time();

  ret = write_and_close_raw_headers(file, file->idx->filename);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  //time->EX = PIDX_get_time();

  return PIDX_success;
}



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

  // Restructuring the grid into power two blocks
  // After this step every process has got a power two block
  // 15 x 31 x 10 ---> 16 x 32 x 16
  //time->idx_rst_start = MPI_Wtime();

  MPI_Comm_size(file->comm, &nprocs);
  if (file->idx_d->data_core_count == nprocs && rst_case_type == 0)
  {
    ret = restructure_init(file, gi, svi, evi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    ret = restructure_io(file, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    ret = restructure(file, PIDX_READ);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

    // Cleanup all buffers nd ids
    //time->buffer_cleanup_start = PIDX_get_time();
    ret = restructure_cleanup(file, gi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
    //time->buffer_cleanup_end = PIDX_get_time();
  }
  else
  {
    ret = restructure_forced_read(file, svi, evi);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}
  }
  //time->idx_rst_end = MPI_Wtime();

  return PIDX_success;
}
