#include "../../PIDX_inc.h"


void PIDX_init_timming_buffers1(PIDX_time time, int group_count, int variable_count, int layout_count, int wavelet_level_count)
{
  int i = 0;
  int g = 0;

  // Restructuring phase timings
  time->rst_init_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_init_start, 0, sizeof(double*) * group_count);
  time->rst_init_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_init_end, 0, sizeof(double*) * group_count);

  time->rst_meta_data_create_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_meta_data_create_start, 0, sizeof(double*) * group_count);
  time->rst_meta_data_create_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_meta_data_create_end, 0, sizeof(double*) * group_count);

  time->rst_meta_data_io_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_meta_data_io_start, 0, sizeof(double*) * group_count);
  time->rst_meta_data_io_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_meta_data_io_end, 0, sizeof(double*) * group_count);

  time->rst_buffer_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_buffer_start, 0, sizeof(double*) * group_count);
  time->rst_buffer_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_buffer_end, 0, sizeof(double*) * group_count);

  time->rst_write_read_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_write_read_start, 0, sizeof(double*) * group_count);
  time->rst_write_read_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_write_read_end, 0, sizeof(double*) * group_count);

  time->rst_buff_agg_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_start, 0, sizeof(double*) * group_count);
  time->rst_buff_agg_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_end, 0, sizeof(double*) * group_count);

  time->rst_buff_agg_free_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_free_start, 0, sizeof(double*) * group_count);
  time->rst_buff_agg_free_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_free_end, 0, sizeof(double*) * group_count);

  time->rst_buff_agg_io_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_io_start, 0, sizeof(double*) * group_count);
  time->rst_buff_agg_io_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_buff_agg_io_end, 0, sizeof(double*) * group_count);

  time->rst_cleanup_start = malloc (sizeof(double*) * group_count);
  memset(time->rst_cleanup_start, 0, sizeof(double*) * group_count);
  time->rst_cleanup_end = malloc (sizeof(double*) * group_count);
  memset(time->rst_cleanup_end, 0, sizeof(double*) * group_count);

  // wavelet phase timings
  time->w_stencil_comm_x_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_x_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_x_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_x_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comm_y_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_y_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_y_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_y_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comm_z_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_z_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_z_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_z_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comm_x_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_x_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_x_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_x_even_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comm_y_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_y_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_y_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_y_even_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comm_z_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_z_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comm_z_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comm_z_even_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_x_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_x_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_x_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_x_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_y_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_y_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_y_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_y_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_z_odd_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_z_odd_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_z_odd_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_z_odd_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_x_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_x_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_x_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_x_even_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_y_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_y_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_y_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_y_even_end, 0, sizeof(double**) * group_count);

  time->w_stencil_comp_z_even_start = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_z_even_start, 0, sizeof(double**) * group_count);
  time->w_stencil_comp_z_even_end = malloc (sizeof(double**) * group_count);
  memset(time->w_stencil_comp_z_even_end, 0, sizeof(double**) * group_count);

  time->w_rst_comp_x_start = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_x_start, 0, sizeof(double**) * group_count);
  time->w_rst_comp_x_end = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_x_end, 0, sizeof(double**) * group_count);

  time->w_rst_comp_y_start = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_y_start, 0, sizeof(double**) * group_count);
  time->w_rst_comp_y_end = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_y_end, 0, sizeof(double**) * group_count);

  time->w_rst_comp_z_start = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_z_start, 0, sizeof(double**) * group_count);
  time->w_rst_comp_z_end = malloc (sizeof(double**) * group_count);
  memset(time->w_rst_comp_z_end, 0, sizeof(double**) * group_count);

  // HZ encoding phase timings
  time->hz_init_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_init_start, 0, sizeof(double*) * group_count);
  time->hz_init_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_init_end, 0, sizeof(double*) * group_count);

  time->hz_meta_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_meta_start, 0, sizeof(double*) * group_count);
  time->hz_meta_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_meta_end, 0, sizeof(double*) * group_count);

  time->hz_buffer_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_buffer_start, 0, sizeof(double*) * group_count);
  time->hz_buffer_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_buffer_end, 0, sizeof(double*) * group_count);

  time->hz_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_start, 0, sizeof(double*) * group_count);
  time->hz_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_end, 0, sizeof(double*) * group_count);

  time->hz_compress_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_compress_start, 0, sizeof(double*) * group_count);
  time->hz_compress_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_compress_end, 0, sizeof(double*) * group_count);

  time->hz_buffer_free_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_buffer_free_start, 0, sizeof(double*) * group_count);
  time->hz_buffer_free_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_buffer_free_end, 0, sizeof(double*) * group_count);

  time->hz_io_start = malloc (sizeof(double**) * group_count);
  memset(time->hz_io_start, 0, sizeof(double**) * group_count);
  time->hz_io_end = malloc (sizeof(double**) * group_count);
  memset(time->hz_io_end, 0, sizeof(double**) * group_count);

  time->hz_cleanup_start = malloc (sizeof(double*) * group_count);
  memset(time->hz_cleanup_start, 0, sizeof(double*) * group_count);
  time->hz_cleanup_end = malloc (sizeof(double*) * group_count);
  memset(time->hz_cleanup_end, 0, sizeof(double*) * group_count);


  // Chunking timings
  time->chunk_init_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_init_start, 0, sizeof(double*) * group_count);
  time->chunk_init_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_init_end, 0, sizeof(double*) * group_count);

  time->chunk_meta_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_meta_start, 0, sizeof(double*) * group_count);
  time->chunk_meta_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_meta_end, 0, sizeof(double*) * group_count);

  time->chunk_buffer_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_buffer_start, 0, sizeof(double*) * group_count);
  time->chunk_buffer_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_buffer_end, 0, sizeof(double*) * group_count);

  time->chunk_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_start, 0, sizeof(double*) * group_count);
  time->chunk_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_end, 0, sizeof(double*) * group_count);

  time->chunk_buffer_free_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_buffer_free_start, 0, sizeof(double*) * group_count);
  time->chunk_buffer_free_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_buffer_free_end, 0, sizeof(double*) * group_count);

  time->chunk_cleanup_start =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_cleanup_start, 0, sizeof(double*) * group_count);
  time->chunk_cleanup_end =  malloc (sizeof(double*) * group_count);
  memset(time->chunk_cleanup_end, 0, sizeof(double*) * group_count);


  // Compression timings
  time->compression_init_start =  malloc (sizeof(double*) * group_count);
  memset(time->compression_init_start, 0, sizeof(double*) * group_count);
  time->compression_init_end =  malloc (sizeof(double*) * group_count);
  memset(time->compression_init_end, 0, sizeof(double*) * group_count);

  time->compression_start =  malloc (sizeof(double*) * group_count);
  memset(time->compression_start, 0, sizeof(double*) * group_count);
  time->compression_end =  malloc (sizeof(double*) * group_count);
  memset(time->compression_end, 0, sizeof(double*) * group_count);


  // File io phase timings
  time->io_start = malloc (sizeof(double*) * group_count);
  memset(time->io_start, 0, sizeof(double*) * group_count);
  time->io_end = malloc (sizeof(double*) * group_count);
  memset(time->io_end, 0, sizeof(double*) * group_count);


  // Aggregation phase timings
  time->agg_init_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_init_start, 0, sizeof(double**) * group_count);
  time->agg_init_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_init_end, 0, sizeof(double**) * group_count);

  time->agg_meta_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_meta_start, 0, sizeof(double**) * group_count);
  time->agg_meta_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_meta_end, 0, sizeof(double**) * group_count);

  time->agg_buf_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_buf_start, 0, sizeof(double**) * group_count);
  time->agg_buf_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_buf_end, 0, sizeof(double**) * group_count);

  time->agg_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_start, 0, sizeof(double**) * group_count);
  time->agg_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_end, 0, sizeof(double**) * group_count);

  time->agg_compress_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_compress_start, 0, sizeof(double**) * group_count);
  time->agg_compress_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_compress_end, 0, sizeof(double**) * group_count);

  time->agg_meta_cleanup_start = malloc (sizeof(double**) * group_count);
  memset(time->agg_meta_cleanup_start, 0, sizeof(double**) * group_count);
  time->agg_meta_cleanup_end = malloc (sizeof(double**) * group_count);
  memset(time->agg_meta_cleanup_end, 0, sizeof(double**) * group_count);


  // Restructuring phase timings
  for (i = 0; i < group_count; i++)
  {
    time->rst_init_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_init_start[i], 0, sizeof(double) * variable_count);
    time->rst_init_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_init_end[i], 0, sizeof(double) * variable_count);

    time->rst_meta_data_create_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_meta_data_create_start[i], 0, sizeof(double) * variable_count);
    time->rst_meta_data_create_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_meta_data_create_end[i], 0, sizeof(double) * variable_count);

    time->rst_meta_data_io_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_meta_data_io_start[i], 0, sizeof(double) * variable_count);
    time->rst_meta_data_io_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_meta_data_io_end[i], 0, sizeof(double) * variable_count);

    time->rst_buffer_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buffer_start[i], 0, sizeof(double) * variable_count);
    time->rst_buffer_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buffer_end[i], 0, sizeof(double) * variable_count);

    time->rst_write_read_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_write_read_start[i], 0, sizeof(double) * variable_count);
    time->rst_write_read_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_write_read_end[i], 0, sizeof(double) * variable_count);

    time->rst_buff_agg_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_start[i], 0, sizeof(double) * variable_count);
    time->rst_buff_agg_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_end[i], 0, sizeof(double) * variable_count);

    time->rst_buff_agg_free_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_free_start[i], 0, sizeof(double) * variable_count);
    time->rst_buff_agg_free_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_free_end[i], 0, sizeof(double) * variable_count);

    time->rst_buff_agg_io_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_io_start[i], 0, sizeof(double) * variable_count);
    time->rst_buff_agg_io_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_buff_agg_io_end[i], 0, sizeof(double) * variable_count);

    time->rst_cleanup_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_cleanup_start[i], 0, sizeof(double) * variable_count);
    time->rst_cleanup_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->rst_cleanup_end[i], 0, sizeof(double) * variable_count);
  }

  // Wavelet phase timings
  for (g = 0; g < group_count; g++)
  {
    time->w_stencil_comm_x_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_x_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_x_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_x_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comm_y_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_y_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_y_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_y_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comm_z_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_z_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_z_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_z_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comm_x_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_x_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_x_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_x_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comm_y_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_y_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_y_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_y_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comm_z_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_z_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comm_z_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comm_z_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_x_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_x_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_x_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_x_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_y_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_y_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_y_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_y_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_z_odd_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_z_odd_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_z_odd_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_z_odd_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_x_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_x_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_x_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_x_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_y_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_y_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_y_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_y_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_stencil_comp_z_even_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_z_even_start[g], 0, sizeof(double*) * variable_count);
    time->w_stencil_comp_z_even_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_stencil_comp_z_even_end[g], 0, sizeof(double*) * variable_count);

    time->w_rst_comp_x_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_x_start[g], 0, sizeof(double*) * variable_count);
    time->w_rst_comp_x_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_x_end[g], 0, sizeof(double*) * variable_count);

    time->w_rst_comp_y_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_y_start[g], 0, sizeof(double*) * variable_count);
    time->w_rst_comp_y_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_y_end[g], 0, sizeof(double*) * variable_count);

    time->w_rst_comp_z_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_z_start[g], 0, sizeof(double*) * variable_count);
    time->w_rst_comp_z_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->w_rst_comp_z_end[g], 0, sizeof(double*) * variable_count);
  }


  // HZ encoding phase timings
  for (i = 0; i < group_count; i++)
  {
    time->hz_init_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_init_start[i], 0, sizeof(double) * variable_count);
    time->hz_init_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_init_end[i], 0, sizeof(double) * variable_count);

    time->hz_meta_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_meta_start[i], 0, sizeof(double) * variable_count);
    time->hz_meta_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_meta_end[i], 0, sizeof(double) * variable_count);

    time->hz_buffer_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_buffer_start[i], 0, sizeof(double) * variable_count);
    time->hz_buffer_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_buffer_end[i], 0, sizeof(double) * variable_count);

    time->hz_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_start[i], 0, sizeof(double) * variable_count);
    time->hz_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_end[i], 0, sizeof(double) * variable_count);

    time->hz_compress_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_compress_start[i], 0, sizeof(double) * variable_count);
    time->hz_compress_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_compress_end[i], 0, sizeof(double) * variable_count);

    time->hz_buffer_free_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_buffer_free_start[i], 0, sizeof(double) * variable_count);
    time->hz_buffer_free_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_buffer_free_end[i], 0, sizeof(double) * variable_count);

    time->hz_cleanup_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_cleanup_start[i], 0, sizeof(double) * variable_count);
    time->hz_cleanup_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->hz_cleanup_end[i], 0, sizeof(double) * variable_count);

    time->hz_io_start[i] = malloc (sizeof(double*) * variable_count);
    memset(time->hz_io_start[i], 0, sizeof(double*) * variable_count);
    time->hz_io_end[i] = malloc (sizeof(double*) * variable_count);
    memset(time->hz_io_end[i], 0, sizeof(double*) * variable_count);
  }


  // Chunking timings
  for (i = 0; i < group_count; i++)
  {
    time->chunk_init_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_init_start[i], 0, sizeof(double) * variable_count);
    time->chunk_init_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_init_end[i], 0, sizeof(double) * variable_count);

    time->chunk_meta_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_meta_start[i], 0, sizeof(double) * variable_count);
    time->chunk_meta_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_meta_end[i], 0, sizeof(double) * variable_count);

    time->chunk_buffer_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_buffer_start[i], 0, sizeof(double) * variable_count);
    time->chunk_buffer_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_buffer_end[i], 0, sizeof(double) * variable_count);

    time->chunk_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_start[i], 0, sizeof(double) * variable_count);
    time->chunk_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_end[i], 0, sizeof(double) * variable_count);

    time->chunk_buffer_free_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_buffer_free_start[i], 0, sizeof(double) * variable_count);
    time->chunk_buffer_free_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_buffer_free_end[i], 0, sizeof(double) * variable_count);

    time->chunk_cleanup_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_cleanup_start[i], 0, sizeof(double) * variable_count);
    time->chunk_cleanup_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->chunk_cleanup_end[i], 0, sizeof(double) * variable_count);
  }

  // Compression timings
  for (i = 0; i < group_count; i++)
  {
    time->compression_init_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->compression_init_start[i], 0, sizeof(double) * variable_count);
    time->compression_init_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->compression_init_end[i], 0, sizeof(double) * variable_count);

    time->compression_start[i] =  malloc (sizeof(double) * variable_count);
    memset(time->compression_start[i], 0, sizeof(double) * variable_count);
    time->compression_end[i] =  malloc (sizeof(double) * variable_count);
    memset(time->compression_end[i], 0, sizeof(double) * variable_count);
  }


  // File io phase timings
  for (i = 0; i < group_count; i++)
  {
    time->io_start[i] = malloc (sizeof(double) * variable_count);
    memset(time->io_start[i], 0, sizeof(double) * variable_count);
    time->io_end[i] = malloc (sizeof(double) * variable_count);
    memset(time->io_end[i], 0, sizeof(double) * variable_count);
  }

  //fprintf(stderr, "gc vc lc %d %d %d\n", group_count, variable_count, layout_count);
  // Aggregation phase timings
  for (g = 0; g < group_count; g++)
  {
    time->agg_init_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_init_start[g], 0, sizeof(double*) * variable_count);
    time->agg_init_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_init_end[g], 0, sizeof(double*) * variable_count);

    time->agg_meta_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_meta_start[g], 0, sizeof(double*) * variable_count);
    time->agg_meta_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_meta_end[g], 0, sizeof(double*) * variable_count);

    time->agg_buf_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_buf_start[g], 0, sizeof(double*) * variable_count);
    time->agg_buf_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_buf_end[g], 0, sizeof(double*) * variable_count);

    time->agg_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_start[g], 0, sizeof(double*) * variable_count);
    time->agg_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_end[g], 0, sizeof(double*) * variable_count);

    time->agg_compress_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_compress_start[g], 0, sizeof(double*) * variable_count);
    time->agg_compress_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_compress_end[g], 0, sizeof(double*) * variable_count);

    time->agg_meta_cleanup_start[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_meta_cleanup_start[g], 0, sizeof(double*) * variable_count);
    time->agg_meta_cleanup_end[g] = malloc (sizeof(double*) * variable_count);
    memset(time->agg_meta_cleanup_end[g], 0, sizeof(double*) * variable_count);
  }

  // Wavelet phase timings
  for (g = 0; g < group_count; g++)
  {
    for (i = 0; i < variable_count; i++)
    {
      time->w_stencil_comm_x_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_x_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_x_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_x_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comm_y_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_y_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_y_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_y_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comm_z_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_z_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_z_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_z_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comm_x_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_x_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_x_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_x_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comm_y_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_y_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_y_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_y_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comm_z_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_z_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comm_z_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comm_z_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_x_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_x_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_x_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_x_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_y_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_y_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_y_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_y_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_z_odd_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_z_odd_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_z_odd_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_z_odd_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_x_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_x_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_x_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_x_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_y_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_y_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_y_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_y_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_stencil_comp_z_even_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_z_even_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_stencil_comp_z_even_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_stencil_comp_z_even_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_rst_comp_x_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_x_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_rst_comp_x_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_x_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_rst_comp_y_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_y_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_rst_comp_y_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_y_end[g][i], 0, sizeof(double) * wavelet_level_count);

      time->w_rst_comp_z_start[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_z_start[g][i], 0, sizeof(double) * wavelet_level_count);
      time->w_rst_comp_z_end[g][i] = malloc (sizeof(double) * wavelet_level_count);
      memset(time->w_rst_comp_z_end[g][i], 0, sizeof(double) * wavelet_level_count);
    }
  }



  // Aggregation phase timings
  for (g = 0; g < group_count; g++)
  {
    for (i = 0; i < variable_count; i++)
    {
      time->agg_init_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_init_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_init_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_init_end[g][i], 0, sizeof(double) * layout_count);

      time->agg_meta_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_meta_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_meta_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_meta_end[g][i], 0, sizeof(double) * layout_count);

      time->agg_buf_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_buf_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_buf_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_buf_end[g][i], 0, sizeof(double) * layout_count);

      time->agg_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_end[g][i], 0, sizeof(double) * layout_count);

      time->agg_compress_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_compress_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_compress_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_compress_end[g][i], 0, sizeof(double) * layout_count);

      time->agg_meta_cleanup_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_meta_cleanup_start[g][i], 0, sizeof(double) * layout_count);
      time->agg_meta_cleanup_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->agg_meta_cleanup_end[g][i], 0, sizeof(double) * layout_count);

      time->hz_io_start[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->hz_io_start[g][i], 0, sizeof(double) * layout_count);
      time->hz_io_end[g][i] = malloc (sizeof(double) * layout_count);
      memset(time->hz_io_end[g][i], 0, sizeof(double) * layout_count);
    }
  }
}



void PIDX_delete_timming_buffers1(PIDX_time time, int group_count, int variable_count)
{
  int i = 0;
  int g = 0;

  for (g = 0; g < group_count; g++)
  {
    for (i = 0; i < variable_count; i++)
    {
      free(time->w_stencil_comm_x_even_start[g][i]);
      free(time->w_stencil_comm_x_even_end[g][i]);
      free(time->w_stencil_comm_y_even_start[g][i]);
      free(time->w_stencil_comm_y_even_end[g][i]);
      free(time->w_stencil_comm_z_even_start[g][i]);
      free(time->w_stencil_comm_z_even_end[g][i]);
      free(time->w_stencil_comm_x_odd_start[g][i]);
      free(time->w_stencil_comm_x_odd_end[g][i]);
      free(time->w_stencil_comm_y_odd_start[g][i]);
      free(time->w_stencil_comm_y_odd_end[g][i]);
      free(time->w_stencil_comm_z_odd_start[g][i]);
      free(time->w_stencil_comm_z_odd_end[g][i]);


      free(time->w_stencil_comp_x_even_start[g][i]);
      free(time->w_stencil_comp_x_even_end[g][i]);
      free(time->w_stencil_comp_y_even_start[g][i]);
      free(time->w_stencil_comp_y_even_end[g][i]);
      free(time->w_stencil_comp_z_even_start[g][i]);
      free(time->w_stencil_comp_z_even_end[g][i]);
      free(time->w_stencil_comp_x_odd_start[g][i]);
      free(time->w_stencil_comp_x_odd_end[g][i]);
      free(time->w_stencil_comp_y_odd_start[g][i]);
      free(time->w_stencil_comp_y_odd_end[g][i]);
      free(time->w_stencil_comp_z_odd_start[g][i]);
      free(time->w_stencil_comp_z_odd_end[g][i]);

      free(time->w_rst_comp_x_start[g][i]);
      free(time->w_rst_comp_x_end[g][i]);
      free(time->w_rst_comp_y_start[g][i]);
      free(time->w_rst_comp_y_end[g][i]);
      free(time->w_rst_comp_z_start[g][i]);
      free(time->w_rst_comp_z_end[g][i]);

      free(time->agg_init_start[g][i]);
      free(time->agg_init_end[g][i]);
      free(time->agg_meta_start[g][i]);
      free(time->agg_meta_end[g][i]);
      free(time->agg_buf_start[g][i]);
      free(time->agg_buf_end[g][i]);
      free(time->agg_start[g][i]);
      free(time->agg_end[g][i]);
      free(time->agg_compress_start[g][i]);
      free(time->agg_compress_end[g][i]);
      free(time->agg_meta_cleanup_start[g][i]);
      free(time->agg_meta_cleanup_end[g][i]);
      free(time->hz_io_start[g][i]);
      free(time->hz_io_end[g][i]);
    }

    free(time->rst_init_start[g]);
    free(time->rst_init_end[g]);
    free(time->rst_meta_data_create_start[g]);
    free(time->rst_meta_data_create_end[g]);
    free(time->rst_meta_data_io_start[g]);
    free(time->rst_meta_data_io_end[g]);
    free(time->rst_buffer_start[g]);
    free(time->rst_buffer_end[g]);
    free(time->rst_write_read_start[g]);
    free(time->rst_write_read_end[g]);
    free(time->rst_buff_agg_start[g]);
    free(time->rst_buff_agg_end[g]);
    free(time->rst_buff_agg_free_start[g]);
    free(time->rst_buff_agg_free_end[g]);
    free(time->rst_buff_agg_io_start[g]);
    free(time->rst_buff_agg_io_end[g]);
    free(time->rst_cleanup_start[g]);
    free(time->rst_cleanup_end[g]);

    free(time->w_stencil_comm_x_even_start[g]);
    free(time->w_stencil_comm_x_even_end[g]);
    free(time->w_stencil_comm_y_even_start[g]);
    free(time->w_stencil_comm_y_even_end[g]);
    free(time->w_stencil_comm_z_even_start[g]);
    free(time->w_stencil_comm_z_even_end[g]);
    free(time->w_stencil_comm_x_odd_start[g]);
    free(time->w_stencil_comm_x_odd_end[g]);
    free(time->w_stencil_comm_y_odd_start[g]);
    free(time->w_stencil_comm_y_odd_end[g]);
    free(time->w_stencil_comm_z_odd_start[g]);
    free(time->w_stencil_comm_z_odd_end[g]);

    free(time->w_stencil_comp_x_even_start[g]);
    free(time->w_stencil_comp_x_even_end[g]);
    free(time->w_stencil_comp_y_even_start[g]);
    free(time->w_stencil_comp_y_even_end[g]);
    free(time->w_stencil_comp_z_even_start[g]);
    free(time->w_stencil_comp_z_even_end[g]);
    free(time->w_stencil_comp_x_odd_start[g]);
    free(time->w_stencil_comp_x_odd_end[g]);
    free(time->w_stencil_comp_y_odd_start[g]);
    free(time->w_stencil_comp_y_odd_end[g]);
    free(time->w_stencil_comp_z_odd_start[g]);
    free(time->w_stencil_comp_z_odd_end[g]);

    free(time->w_rst_comp_x_start[g]);
    free(time->w_rst_comp_x_end[g]);
    free(time->w_rst_comp_y_start[g]);
    free(time->w_rst_comp_y_end[g]);
    free(time->w_rst_comp_z_start[g]);
    free(time->w_rst_comp_z_end[g]);

    free(time->hz_init_start[g]);
    free(time->hz_init_end[g]);
    free(time->hz_meta_start[g]);
    free(time->hz_meta_end[g]);
    free(time->hz_buffer_start[g]);
    free(time->hz_buffer_end[g]);
    free(time->hz_start[g]);
    free(time->hz_end[g]);
    free(time->hz_compress_start[g]);
    free(time->hz_compress_end[g]);
    free(time->hz_io_start[g]);
    free(time->hz_io_end[g]);
    free(time->hz_buffer_free_start[g]);
    free(time->hz_buffer_free_end[g]);
    free(time->hz_cleanup_start[g]);
    free(time->hz_cleanup_end[g]);

    free(time->chunk_init_start[g]);
    free(time->chunk_init_end[g]);
    free(time->chunk_meta_start[g]);
    free(time->chunk_meta_end[g]);
    free(time->chunk_buffer_start[g]);
    free(time->chunk_buffer_end[g]);
    free(time->chunk_start[g]);
    free(time->chunk_end[g]);
    free(time->chunk_buffer_free_start[g]);
    free(time->chunk_buffer_free_end[g]);
    free(time->chunk_cleanup_start[g]);
    free(time->chunk_cleanup_end[g]);

    free(time->compression_init_start[g]);
    free(time->compression_init_end[g]);
    free(time->compression_start[g]);
    free(time->compression_end[g]);

    free(time->io_start[g]);
    free(time->io_end[g]);

    free(time->agg_init_start[g]);
    free(time->agg_init_end[g]);
    free(time->agg_meta_start[g]);
    free(time->agg_meta_end[g]);
    free(time->agg_buf_start[g]);
    free(time->agg_buf_end[g]);
    free(time->agg_start[g]);
    free(time->agg_end[g]);
    free(time->agg_compress_start[g]);
    free(time->agg_compress_end[g]);
    free(time->agg_meta_cleanup_start[g]);
    free(time->agg_meta_cleanup_end[g]);
  }

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

  free(time->w_stencil_comm_x_even_start);
  free(time->w_stencil_comm_x_even_end);
  free(time->w_stencil_comm_y_even_start);
  free(time->w_stencil_comm_y_even_end);
  free(time->w_stencil_comm_z_even_start);
  free(time->w_stencil_comm_z_even_end);
  free(time->w_stencil_comm_x_odd_start);
  free(time->w_stencil_comm_x_odd_end);
  free(time->w_stencil_comm_y_odd_start);
  free(time->w_stencil_comm_y_odd_end);
  free(time->w_stencil_comm_z_odd_start);
  free(time->w_stencil_comm_z_odd_end);

  free(time->w_stencil_comp_x_even_start);
  free(time->w_stencil_comp_x_even_end);
  free(time->w_stencil_comp_y_even_start);
  free(time->w_stencil_comp_y_even_end);
  free(time->w_stencil_comp_z_even_start);
  free(time->w_stencil_comp_z_even_end);
  free(time->w_stencil_comp_x_odd_start);
  free(time->w_stencil_comp_x_odd_end);
  free(time->w_stencil_comp_y_odd_start);
  free(time->w_stencil_comp_y_odd_end);
  free(time->w_stencil_comp_z_odd_start);
  free(time->w_stencil_comp_z_odd_end);

  free(time->w_rst_comp_x_start);
  free(time->w_rst_comp_x_end);
  free(time->w_rst_comp_y_start);
  free(time->w_rst_comp_y_end);
  free(time->w_rst_comp_z_start);
  free(time->w_rst_comp_z_end);

  free(time->hz_init_start);
  free(time->hz_init_end);
  free(time->hz_meta_start);
  free(time->hz_meta_end);
  free(time->hz_buffer_start);
  free(time->hz_buffer_end);
  free(time->hz_start);
  free(time->hz_end);
  free(time->hz_compress_start);
  free(time->hz_compress_end);
  free(time->hz_io_start);
  free(time->hz_io_end);
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

  free(time->agg_init_start);
  free(time->agg_init_end);
  free(time->agg_meta_start);
  free(time->agg_meta_end);
  free(time->agg_buf_start);
  free(time->agg_buf_end);
  free(time->agg_start);
  free(time->agg_end);
  free(time->agg_compress_start);
  free(time->agg_compress_end);
  free(time->agg_meta_cleanup_start);
  free(time->agg_meta_cleanup_end);

  free(time->io_start);
  free(time->io_end);
}
