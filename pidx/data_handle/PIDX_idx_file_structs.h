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

#ifndef __PIDX_IDX_FILE_STRUCTS_H
#define __PIDX_IDX_FILE_STRUCTS_H


struct idx_file_struct
{
  int pidx_version;
  char metadata_version[8];


  enum PIDX_io_type io_type;                        /// I/O format and layout we want to use

  int current_time_step;                            /// The current timestep selected
  int first_tstep;                                  /// Index of the frist timestep
  int last_tstep;                                   /// Index of the last timestep


  int variable_pipe_length;                         /// pipes (combines) "variable_pipe_length" variables for io
  int variable_tracker[512];                        /// Which one of the 256 variables are present
  PIDX_variable variable[512];                      /// pointer to variable
  uint32_t variable_count;                          /// The number of variables contained in the dataset
  uint32_t particles_position_variable_index;       /// The index of the variable containing the particles position

  Agg_buffer **agg_buffer;                          /// aggregation related struct


  char filename[1024];                              /// The idx file path
  char filename_template[1024];
  char filename_partition[1024];
  char filename_template_partition[1024];

  int bits_per_block;                               /// Number of bits per block
  uint64_t samples_per_block;                       /// Number of samples in a block 2^bits_per_block
  int blocks_per_file;                              /// Number of blocks per file


  int max_file_count;                               /// maximum number of file (power_2(bounds[0]) * power_2(bounds[1]) * power_2(bounds[2])) / blocks_per_file * samples_per_block
  uint64_t bounds[PIDX_MAX_DIMENSIONS];             /// Logical bounds of the dataset
  uint64_t box_bounds[PIDX_MAX_DIMENSIONS];         /// Logical bounds of the box query
  double physical_bounds[PIDX_MAX_DIMENSIONS];      /// Physical bounds of the dataset (used in particle io)
  double physical_box_bounds[PIDX_MAX_DIMENSIONS];  /// Logical bounds of the box query (used in particle io)


  int maxh;                                         /// total number of hz levels
  char bitSequence[512];                            /// bitsequence used for HZ indexing, controls the layout
  char bitPattern[512];

  int compression_type;
  int compression_factor;
  float compression_bit_rate;
  uint64_t chunk_size[PIDX_MAX_DIMENSIONS];

  int particle_res_base;
  int particle_res_factor;

  int current_resolution;

  enum PIDX_endian_type endian;                     /// 1 for little endian and 0 for big endian
  int flip_endian;                                  /// 1 for flipping endianness required


  uint32_t partition_count[PIDX_MAX_DIMENSIONS];    /// number of partitions in the X, Y and Z dimesnions
  uint32_t partition_size[PIDX_MAX_DIMENSIONS];     /// size of partitions in each of the dimensions in voxels
  uint32_t partition_offset[PIDX_MAX_DIMENSIONS];   /// offset of each of the partition (n global index space)


  int cached_ts;                                    /// used for raw io, to cache meta data (1) or not (0)
};
typedef struct idx_file_struct* idx_dataset;


#endif
