/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */
#include "../../PIDX_inc.h"


void PIDX_init_timming_buffers1(PIDX_time time, int variable_count, int layout_count)
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

  time->hz_compress_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_compress_start, 0, sizeof(double) * variable_count);
  time->hz_compress_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_compress_end, 0, sizeof(double) * variable_count);

  time->hz_buffer_free_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_free_start, 0, sizeof(double) * variable_count);
  time->hz_buffer_free_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_buffer_free_end, 0, sizeof(double) * variable_count);

  time->hz_cleanup_start = malloc (sizeof(double) * variable_count);
  memset(time->hz_cleanup_start, 0, sizeof(double) * variable_count);
  time->hz_cleanup_end = malloc (sizeof(double) * variable_count);
  memset(time->hz_cleanup_end, 0, sizeof(double) * variable_count);

  time->hz_io_start = malloc (sizeof(double*) * variable_count);
  memset(time->hz_io_start, 0, sizeof(double*) * variable_count);
  time->hz_io_end = malloc (sizeof(double*) * variable_count);
  memset(time->hz_io_end, 0, sizeof(double*) * variable_count);



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



  // File io phase timings
  time->io_start = malloc (sizeof(double) * variable_count);
  memset(time->io_start, 0, sizeof(double) * variable_count);
  time->io_end = malloc (sizeof(double) * variable_count);
  memset(time->io_end, 0, sizeof(double) * variable_count);


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

  time->agg_compress_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_compress_start, 0, sizeof(double*) * variable_count);
  time->agg_compress_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_compress_end, 0, sizeof(double*) * variable_count);

  time->agg_meta_cleanup_start = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_cleanup_start, 0, sizeof(double*) * variable_count);
  time->agg_meta_cleanup_end = malloc (sizeof(double*) * variable_count);
  memset(time->agg_meta_cleanup_end, 0, sizeof(double*) * variable_count);


  // Aggregation phase timings
  for (int i = 0; i < variable_count; i++)
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

    time->agg_compress_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_compress_start[i], 0, sizeof(double) * layout_count);
    time->agg_compress_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_compress_end[i], 0, sizeof(double) * layout_count);

    time->agg_meta_cleanup_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_cleanup_start[i], 0, sizeof(double) * layout_count);
    time->agg_meta_cleanup_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->agg_meta_cleanup_end[i], 0, sizeof(double) * layout_count);

    time->hz_io_start[i] = malloc (sizeof(double) * layout_count);
    memset(time->hz_io_start[i], 0, sizeof(double) * layout_count);
    time->hz_io_end[i] = malloc (sizeof(double) * layout_count);
    memset(time->hz_io_end[i], 0, sizeof(double) * layout_count);
  }

}



void PIDX_delete_timming_buffers1(PIDX_time time, int variable_count)
{
  for (int i = 0; i < variable_count; i++)
  {
    free(time->agg_init_start[i]);
    free(time->agg_init_end[i]);
    free(time->agg_meta_start[i]);
    free(time->agg_meta_end[i]);
    free(time->agg_buf_start[i]);
    free(time->agg_buf_end[i]);
    free(time->agg_start[i]);
    free(time->agg_end[i]);
    free(time->agg_compress_start[i]);
    free(time->agg_compress_end[i]);
    free(time->agg_meta_cleanup_start[i]);
    free(time->agg_meta_cleanup_end[i]);
    free(time->hz_io_start[i]);
    free(time->hz_io_end[i]);
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

  free(time->io_start);
  free(time->io_end);

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

  return;
}
