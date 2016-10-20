#include "../PIDX_inc.h"

void PIDX_init_timming_buffers1(PIDX_time time, int variable_count)
{
  variable_count = variable_count * 2;
  time->init_start = malloc (sizeof(double) * variable_count);            memset(time->init_start, 0, sizeof(double) * variable_count);
  time->init_end = malloc (sizeof(double) * variable_count);              memset(time->init_end, 0, sizeof(double) * variable_count);
  time->write_init_start = malloc (sizeof(double) * variable_count);      memset(time->write_init_start, 0, sizeof(double) * variable_count);
  time->write_init_end = malloc (sizeof(double) * variable_count);        memset(time->write_init_end, 0, sizeof(double) * variable_count);
}



void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count)
{
  variable_count = variable_count * 2;
  time->header_counter = 0;
  time->variable_counter = 0;
  time->startup_start = malloc (sizeof(double) * variable_count);                 memset(time->startup_start, 0, sizeof(double) * variable_count);
  time->startup_end = malloc (sizeof(double) * variable_count);                   memset(time->startup_end, 0, sizeof(double) * variable_count);

  time->rst_meta_data_start_io = malloc (sizeof(double) * variable_count);        memset(time->rst_meta_data_start_io, 0, sizeof(double) * variable_count);
  time->rst_meta_data_end_io = malloc (sizeof(double) * variable_count);          memset(time->rst_meta_data_end_io, 0, sizeof(double) * variable_count);

  time->rst_start = malloc (sizeof(double) * variable_count);                     memset(time->rst_start, 0, sizeof(double) * variable_count);
  time->rst_end = malloc (sizeof(double) * variable_count);                       memset(time->rst_end, 0, sizeof(double) * variable_count);

  time->rst_io_start = malloc (sizeof(double) * variable_count);                  memset(time->rst_io_start, 0, sizeof(double) * variable_count);
  time->rst_io_end = malloc (sizeof(double) * variable_count);                    memset(time->rst_io_end, 0, sizeof(double) * variable_count);

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

  time->agg_meta_start = malloc (sizeof(double*) * variable_count);   memset(time->agg_meta_start, 0, sizeof(double*) * variable_count);
  time->agg_meta_end = malloc (sizeof(double*) * variable_count);     memset(time->agg_meta_end, 0, sizeof(double*) * variable_count);

  time->windows_start = malloc (sizeof(double*) * variable_count);       memset(time->windows_start, 0, sizeof(double*) * variable_count);
  time->windows_end = malloc (sizeof(double*) * variable_count);     memset(time->windows_end, 0, sizeof(double*) * variable_count);

  time->first_fence_end = malloc (sizeof(double*) * variable_count);       memset(time->first_fence_end, 0, sizeof(double*) * variable_count);
  time->first_fence_start = malloc (sizeof(double*) * variable_count);     memset(time->first_fence_start, 0, sizeof(double*) * variable_count);

  time->second_fence_end = malloc (sizeof(double*) * variable_count);       memset(time->second_fence_end, 0, sizeof(double*) * variable_count);
  time->second_fence_start = malloc (sizeof(double*) * variable_count);     memset(time->second_fence_start, 0, sizeof(double*) * variable_count);

  time->fence_free = malloc (sizeof(double*) * variable_count);       memset(time->fence_free, 0, sizeof(double*) * variable_count);

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

    time->agg_meta_start[i] = malloc (sizeof(double) * layout_count);   memset(time->agg_meta_start[i], 0, sizeof(double) * layout_count);
    time->agg_meta_end[i] = malloc (sizeof(double) * layout_count);     memset(time->agg_meta_end[i], 0, sizeof(double) * layout_count);

    time->first_fence_end[i] = malloc (sizeof(double) * layout_count);   memset(time->first_fence_end[i], 0, sizeof(double) * layout_count);
    time->first_fence_start[i] = malloc (sizeof(double) * layout_count); memset(time->first_fence_start[i], 0, sizeof(double) * layout_count);

    time->windows_end[i] = malloc (sizeof(double) * layout_count);   memset(time->windows_end[i], 0, sizeof(double) * layout_count);
    time->windows_start[i] = malloc (sizeof(double) * layout_count); memset(time->windows_start[i], 0, sizeof(double) * layout_count);

    time->fence_free[i] = malloc (sizeof(double) * layout_count);   memset(time->fence_free[i], 0, sizeof(double) * layout_count);

    time->second_fence_end[i] = malloc (sizeof(double) * layout_count);   memset(time->second_fence_end[i], 0, sizeof(double) * layout_count);
    time->second_fence_start[i] = malloc (sizeof(double) * layout_count); memset(time->second_fence_start[i], 0, sizeof(double) * layout_count);

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

  free(time->rst_io_end);
  free(time->rst_io_start);

  free(time->hz_start);
  free(time->hz_end);

  free(time->rst_meta_data_end_io);
  free(time->rst_meta_data_start_io);

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
  variable_count = variable_count * 2;
  for (i = 0; i < variable_count; i++)
  {
    free(time->agg_start[i]);
    free(time->agg_end[i]);
    free(time->agg_buf_start[i]);
    free(time->agg_buf_end[i]);

    free(time->agg_meta_start[i]);
    free(time->agg_meta_end[i]);

    free(time->first_fence_start[i]);
    free(time->first_fence_end[i]);


    free(time->windows_start[i]);
    free(time->windows_end[i]);

    free(time->second_fence_start[i]);
    free(time->second_fence_end[i]);

    free(time->io_start[i]);
    free(time->io_end[i]);

    free(time->io_per_process_start[i]);
    free(time->io_per_process_end[i]);

    free(time->fence_free[i]);
  }
  free(time->agg_start);
  free(time->agg_end);
  free(time->agg_buf_start);
  free(time->agg_buf_end);

  free(time->agg_meta_start);
  free(time->agg_meta_end);

  free(time->first_fence_start);
  free(time->first_fence_end);
  free(time->fence_free);

  free(time->windows_end);
  free(time->windows_start);

  free(time->second_fence_start);
  free(time->second_fence_end);

  free(time->io_start);
  free(time->io_end);
  free(time->io_per_process_start);
  free(time->io_per_process_end);

}
