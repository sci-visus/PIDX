#include "PIDX_io.h"

/// Returns elapsed time
double PIDX_get_time()
{
#if PIDX_HAVE_MPI
  return MPI_Wtime();
#else
  struct timeval temp;
  gettimeofday(&temp, NULL);
  return (double)(temp.tv_sec) + (double)(temp.tv_usec)/1000000.0;
#endif
}

void PIDX_init_timming_buffers1(PIDX_time time, int variable_count)
{
  time->init_start = malloc (sizeof(double) * variable_count);            memset(time->init_start, 0, sizeof(double) * variable_count);
  time->init_end = malloc (sizeof(double) * variable_count);              memset(time->init_end, 0, sizeof(double) * variable_count);
  time->write_init_start = malloc (sizeof(double) * variable_count);      memset(time->write_init_start, 0, sizeof(double) * variable_count);
  time->write_init_end = malloc (sizeof(double) * variable_count);        memset(time->write_init_end, 0, sizeof(double) * variable_count);
}



void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count)
{
  time->header_counter = 0;
  time->startup_start = malloc (sizeof(double) * variable_count);                 memset(time->startup_start, 0, sizeof(double) * variable_count);
  time->startup_end = malloc (sizeof(double) * variable_count);                   memset(time->startup_end, 0, sizeof(double) * variable_count);

  time->rst_start = malloc (sizeof(double) * variable_count);                     memset(time->rst_start, 0, sizeof(double) * variable_count);
  time->rst_end = malloc (sizeof(double) * variable_count);                       memset(time->rst_end, 0, sizeof(double) * variable_count);
  time->hz_start = malloc (sizeof(double) * variable_count);                      memset(time->hz_start, 0, sizeof(double) * variable_count);
  time->hz_end = malloc (sizeof(double) * variable_count);                        memset(time->hz_end, 0, sizeof(double) * variable_count);

  time->cleanup_start = malloc (sizeof(double) * variable_count);                 memset(time->cleanup_start, 0, sizeof(double) * variable_count);
  time->cleanup_end = malloc (sizeof(double) * variable_count);                   memset(time->cleanup_end, 0, sizeof(double) * variable_count);
  time->finalize_start = malloc (sizeof(double) * variable_count);                memset(time->finalize_start, 0, sizeof(double) * variable_count);
  time->finalize_end = malloc (sizeof(double) * variable_count);                  memset(time->finalize_end, 0, sizeof(double) * variable_count);

  time->buffer_start = malloc (sizeof(double) * variable_count);                  memset(time->buffer_start, 0, sizeof(double) * variable_count);
  time->buffer_end = malloc (sizeof(double) * variable_count);                    memset(time->buffer_end, 0, sizeof(double) * variable_count);

  time->chunk_start =  malloc (sizeof(double) * variable_count);                  memset(time->chunk_start, 0, sizeof(double) * variable_count);
  time->chunk_end =  malloc (sizeof(double) * variable_count);                    memset(time->chunk_end, 0, sizeof(double) * variable_count);
  time->compression_start =  malloc (sizeof(double) * variable_count);            memset(time->compression_start, 0, sizeof(double) * variable_count);
  time->compression_end =  malloc (sizeof(double) * variable_count);              memset(time->compression_end, 0, sizeof(double) * variable_count);

  time->agg_start = malloc (sizeof(double*) * variable_count);       memset(time->agg_start, 0, sizeof(double*) * variable_count);
  time->agg_end = malloc (sizeof(double*) * variable_count);         memset(time->agg_end, 0, sizeof(double*) * variable_count);
  time->agg_buf_start = malloc (sizeof(double*) * variable_count);   memset(time->agg_buf_start, 0, sizeof(double*) * variable_count);
  time->agg_buf_end = malloc (sizeof(double*) * variable_count);     memset(time->agg_buf_end, 0, sizeof(double*) * variable_count);
  time->io_start = malloc (sizeof(double*) * variable_count);        memset(time->io_start, 0, sizeof(double*) * variable_count);
  time->io_end = malloc (sizeof(double*) * variable_count);          memset(time->io_end, 0, sizeof(double*) * variable_count);
  time->io_per_process_start = malloc (sizeof(double*) * variable_count);          memset(time->io_per_process_start, 0, sizeof(double*) * variable_count);
  time->io_per_process_end = malloc (sizeof(double*) * variable_count);          memset(time->io_per_process_end, 0, sizeof(double*) * variable_count);

  int i = 0;
  for (i = 0; i < variable_count; i++)
  {
    time->agg_start[i] = malloc (sizeof(double) * layout_count);       memset(time->agg_start[i], 0, sizeof(double) * layout_count);
    time->agg_end[i] = malloc (sizeof(double) * layout_count);         memset(time->agg_end[i], 0, sizeof(double) * layout_count);
    time->agg_buf_start[i] = malloc (sizeof(double) * layout_count);   memset(time->agg_buf_start[i], 0, sizeof(double) * layout_count);
    time->agg_buf_end[i] = malloc (sizeof(double) * layout_count);     memset(time->agg_buf_end[i], 0, sizeof(double) * layout_count);
    time->io_start[i] = malloc (sizeof(double) * layout_count);        memset(time->io_start[i], 0, sizeof(double) * layout_count);
    time->io_end[i] = malloc (sizeof(double) * layout_count);          memset(time->io_end[i], 0, sizeof(double) * layout_count);

    time->io_per_process_start[i] = malloc (sizeof(double) * layout_count);        memset(time->io_per_process_start[i], 0, sizeof(double) * layout_count);
    time->io_per_process_end[i] = malloc (sizeof(double) * layout_count);          memset(time->io_per_process_end[i], 0, sizeof(double) * layout_count);
  }
}


void PIDX_delete_timming_buffers1(PIDX_time time)
{
  free(time->init_start);
  free(time->init_end);
  free(time->write_init_start);
  free(time->write_init_end);
}



void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count)
{

  free(time->startup_start);
  free(time->startup_end);

  free(time->rst_start);
  free(time->rst_end);
  free(time->hz_start);
  free(time->hz_end);

  free(time->cleanup_start);
  free(time->cleanup_end);
  free(time->finalize_start);
  free(time->finalize_end);

  free(time->buffer_start);
  free(time->buffer_end);

  free(time->chunk_start);
  free(time->chunk_end);
  free(time->compression_start);
  free(time->compression_end);

  int i = 0;
  for (i = 0; i < variable_count; i++)
  {
    free(time->agg_start[i]);
    free(time->agg_end[i]);
    free(time->agg_buf_start[i]);
    free(time->agg_buf_end[i]);
    free(time->io_start[i]);
    free(time->io_end[i]);

    free(time->io_per_process_start[i]);
    free(time->io_per_process_end[i]);
  }
  free(time->agg_start);
  free(time->agg_end);
  free(time->agg_buf_start);
  free(time->agg_buf_end);
  free(time->io_start);
  free(time->io_end);
  free(time->io_per_process_start);
  free(time->io_per_process_end);

}
