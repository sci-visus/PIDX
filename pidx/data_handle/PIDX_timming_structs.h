/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2010-2018 ViSUS L.L.C.,
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


#ifndef __PIDX_TIMMING_STRUCTS_H
#define __PIDX_TIMMING_STRUCTS_H


struct PIDX_timming_struct
{
  double SX, EX;
  double sim_start, sim_end;

  double init_start, init_end;
  double set_reg_box_start, set_reg_box_end;
  double bit_string_start, bit_string_end;
  double layout_start, layout_end;
  double header_io_start, header_io_end;
  double group_cleanup_start, group_cleanup_end;
  double partition_start, partition_end;
  double partition_cleanup_start, partition_cleanup_end;

  double particle_meta_data_io_start, particle_meta_data_io_end;
  double particle_data_io_start, particle_data_io_end;

  double particle_reshuffling_start, particle_reshuffling_end;

  double *rst_init_start, *rst_init_end;
  double *rst_meta_data_create_start, *rst_meta_data_create_end;
  double *rst_meta_data_io_start, *rst_meta_data_io_end;
  double *rst_buffer_start, *rst_buffer_end;
  double *rst_write_read_start, *rst_write_read_end;
  double *rst_buff_agg_start, *rst_buff_agg_end;
  double *rst_buff_agg_free_start, *rst_buff_agg_free_end;
  double *rst_buff_agg_io_start, *rst_buff_agg_io_end;
  double *rst_cleanup_start, *rst_cleanup_end;

  double *hz_init_start, *hz_init_end;
  double *hz_meta_start, *hz_meta_end;
  double *hz_buffer_start, *hz_buffer_end;
  double *hz_start, *hz_end;
  double *hz_compress_start, *hz_compress_end;
  double *hz_buffer_free_start, *hz_buffer_free_end;
  double *hz_cleanup_start, *hz_cleanup_end;
  double **hz_io_start, **hz_io_end;

  double *chunk_init_start, *chunk_init_end;
  double *chunk_meta_start, *chunk_meta_end;
  double *chunk_buffer_start, *chunk_buffer_end;
  double *chunk_start, *chunk_end;
  double *chunk_buffer_free_start, *chunk_buffer_free_end;
  double *chunk_cleanup_start, *chunk_cleanup_end;

  double *compression_init_start, *compression_init_end;
  double *compression_start, *compression_end;

  double **agg_init_start, **agg_init_end;
  double **agg_meta_start, **agg_meta_end;
  double **agg_buf_start, **agg_buf_end;
  double **agg_start, **agg_end;
  double **agg_compress_start, **agg_compress_end;
  double **agg_meta_cleanup_start, **agg_meta_cleanup_end;

  double *io_start, *io_end;
};
typedef struct PIDX_timming_struct* PIDX_time;


#endif
