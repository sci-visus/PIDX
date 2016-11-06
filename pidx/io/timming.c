#include "../PIDX_inc.h"

void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count)
{
  // Restructuring phase timings
  time->rst_init_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_init_start, 0, sizeof(double) * variable_count);
  time->rst_init_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_init_end, 0, sizeof(double) * variable_count);

  time->rst_meta_data_create_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_meta_data_create_start, 0, sizeof(double) * variable_count);
  time->rst_meta_data_create_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_meta_data_create_end, 0, sizeof(double) * variable_count);

  time->rst_meta_data_io_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_meta_data_io_start, 0, sizeof(double) * variable_count);
  time->rst_meta_data_io_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_meta_data_io_end, 0, sizeof(double) * variable_count);

  time->rst_buffer_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_buffer_start, 0, sizeof(double) * variable_count);
  time->rst_buffer_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_buffer_end, 0, sizeof(double) * variable_count);

  time->rst_write_read_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_write_read_start, 0, sizeof(double) * variable_count);
  time->rst_write_read_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_write_read_end, 0, sizeof(double) * variable_count);

  time->rst_buff_agg_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_start, 0, sizeof(double) * variable_count);
  time->rst_buff_agg_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_end, 0, sizeof(double) * variable_count);

  time->rst_buff_agg_free_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_free_start, 0, sizeof(double) * variable_count);
  time->rst_buff_agg_free_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_free_end, 0, sizeof(double) * variable_count);

  time->rst_buff_agg_io_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_io_start, 0, sizeof(double) * variable_count);
  time->rst_buff_agg_io_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_buff_agg_io_end, 0, sizeof(double) * variable_count);

  time->rst_cleanup_start = malloc (sizeof(double) * variable_count);
  memset(time->rst_cleanup_start, 0, sizeof(double) * variable_count);
  time->rst_cleanup_end = malloc (sizeof(double) * variable_count);
  memset(time->rst_cleanup_end, 0, sizeof(double) * variable_count);


  // HZ encoding phase timings
  time->hz_init_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_init_start, 0, sizeof(double) * variable_count);
  time->hz_init_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_init_end, 0, sizeof(double) * variable_count);

  time->hz_meta_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_meta_start, 0, sizeof(double) * variable_count);
  time->hz_meta_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_meta_end, 0, sizeof(double) * variable_count);

  time->hz_buffer_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_start, 0, sizeof(double) * variable_count);
  time->hz_buffer_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_end, 0, sizeof(double) * variable_count);

  time->hz_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_start, 0, sizeof(double) * variable_count);
  time->hz_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_end, 0, sizeof(double) * variable_count);

  time->hz_buffer_free_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_free_start, 0, sizeof(double) * variable_count);
  time->hz_buffer_free_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_free_end, 0, sizeof(double) * variable_count);

  time->hz_cleanup_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_cleanup_start, 0, sizeof(double) * variable_count);
  time->hz_cleanup_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_cleanup_end, 0, sizeof(double) * variable_count);


  // Chunking timings
  time->chunk_init_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_init_start, 0, sizeof(double) * variable_count);
  time->chunk_init_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_init_end, 0, sizeof(double) * variable_count);

  time->chunk_meta_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_meta_start, 0, sizeof(double) * variable_count);
  time->chunk_meta_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_meta_end, 0, sizeof(double) * variable_count);

  time->chunk_buffer_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_buffer_start, 0, sizeof(double) * variable_count);
  time->chunk_buffer_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_buffer_end, 0, sizeof(double) * variable_count);

  time->chunk_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_start, 0, sizeof(double) * variable_count);
  time->chunk_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_end, 0, sizeof(double) * variable_count);

  time->chunk_buffer_free_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_buffer_free_start, 0, sizeof(double) * variable_count);
  time->chunk_buffer_free_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_buffer_free_end, 0, sizeof(double) * variable_count);

  time->chunk_cleanup_start =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_cleanup_start, 0, sizeof(double) * variable_count);
  time->chunk_cleanup_end =  malloc (sizeof(double) * variable_count);
  memset(time->chunk_cleanup_end, 0, sizeof(double) * variable_count);


  // Compression timings
  time->compression_init_start =  malloc (sizeof(double) * variable_count);
  memset(time->compression_init_start, 0, sizeof(double) * variable_count);
  time->compression_init_end =  malloc (sizeof(double) * variable_count);
  memset(time->compression_init_end, 0, sizeof(double) * variable_count);

  time->compression_start =  malloc (sizeof(double) * variable_count);
  memset(time->compression_start, 0, sizeof(double) * variable_count);
  time->compression_end =  malloc (sizeof(double) * variable_count);
  memset(time->compression_end, 0, sizeof(double) * variable_count);


  // Aggregation phase timings
  time->agg_init_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_init_start, 0, sizeof(double*) * variable_count);
  time->agg_init_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_init_end, 0, sizeof(double*) * variable_count);

  time->agg_meta_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_start, 0, sizeof(double*) * variable_count);
  time->agg_meta_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_end, 0, sizeof(double*) * variable_count);

  time->agg_buf_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_buf_start, 0, sizeof(double*) * variable_count);
  time->agg_buf_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_buf_end, 0, sizeof(double*) * variable_count);

  time->agg_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_start, 0, sizeof(double*) * variable_count);
  time->agg_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_end, 0, sizeof(double*) * variable_count);

  time->agg_meta_cleanup_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_cleanup_start, 0, sizeof(double*) * variable_count);
  time->agg_meta_cleanup_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_cleanup_end, 0, sizeof(double*) * variable_count);

  int i = 0;
  for (i = 0; i < variable_count; i++)
  {
    time->agg_init_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_init_start[i], 0, sizeof(double) * layout_count);
    time->agg_init_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_init_end[i], 0, sizeof(double) * layout_count);

    time->agg_meta_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_start[i], 0, sizeof(double) * layout_count);
    time->agg_meta_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_end[i], 0, sizeof(double) * layout_count);

    time->agg_buf_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_buf_start[i], 0, sizeof(double) * layout_count);
    time->agg_buf_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_buf_end[i], 0, sizeof(double) * layout_count);

    time->agg_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_start[i], 0, sizeof(double) * layout_count);
    time->agg_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_end[i], 0, sizeof(double) * layout_count);

    time->agg_meta_cleanup_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_cleanup_start[i], 0, sizeof(double) * layout_count);
    time->agg_meta_cleanup_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_cleanup_end[i], 0, sizeof(double) * layout_count);
  }

  // File io phase timings
  time->io_start = malloc (sizeof(double) * variable_count);
  memset(time->io_start, 0, sizeof(double) * variable_count);
  time->io_end = malloc (sizeof(double) * variable_count);
  memset(time->io_end, 0, sizeof(double) * variable_count);
}


void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count)
{

  free(time->rst_init_start);
  free(time->rst_init_end);
  free(time->rst_meta_data_create_start);
  free(time->rst_meta_data_create_end);
  free(time->rst_meta_data_io_start);
  free(time->rst_meta_data_io_end);
  free(time->rst_buffer_start);
  free(time->rst_buffer_end);
  free(time->rst_write_read_start);
  free(time->rst_write_read_end);
  free(time->rst_buff_agg_start);
  free(time->rst_buff_agg_end);
  free(time->rst_buff_agg_free_start);
  free(time->rst_buff_agg_free_end);
  free(time->rst_buff_agg_io_start);
  free(time->rst_buff_agg_io_end);
  free(time->rst_cleanup_start);
  free(time->rst_cleanup_end);

  free(time->hz_init_start);
  free(time->hz_init_end);
  free(time->hz_meta_start);
  free(time->hz_meta_end);
  free(time->hz_buffer_start);
  free(time->hz_buffer_end);
  free(time->hz_start);
  free(time->hz_end);
  free(time->hz_buffer_free_start);
  free(time->hz_buffer_free_end);
  free(time->hz_cleanup_start);
  free(time->hz_cleanup_end);

  free(time->chunk_init_start);
  free(time->chunk_init_end);
  free(time->chunk_meta_start);
  free(time->chunk_meta_end);
  free(time->chunk_buffer_start);
  free(time->chunk_buffer_end);
  free(time->chunk_start);
  free(time->chunk_end);
  free(time->chunk_buffer_free_start);
  free(time->chunk_buffer_free_end);
  free(time->chunk_cleanup_start);
  free(time->chunk_cleanup_end);

  free(time->compression_init_start);
  free(time->compression_init_end);
  free(time->compression_start);
  free(time->compression_end);

  int i = 0;
  for (i = 0; i < variable_count; i++)
  {
    free(time->agg_init_start[i]);
    free(time->agg_init_end[i]);
    free(time->agg_meta_start[i]);
    free(time->agg_meta_end[i]);
    free(time->agg_buf_start[i]);
    free(time->agg_buf_end[i]);
    free(time->agg_start[i]);
    free(time->agg_end[i]);
    free(time->agg_meta_cleanup_start[i]);
    free(time->agg_meta_cleanup_end[i]);
  }

  free(time->agg_init_start);
  free(time->agg_init_end);
  free(time->agg_meta_start);
  free(time->agg_meta_end);
  free(time->agg_buf_start);
  free(time->agg_buf_end);
  free(time->agg_start);
  free(time->agg_end);
  free(time->agg_meta_cleanup_start);
  free(time->agg_meta_cleanup_end);

  free(time->io_start);
  free(time->io_end);
}
