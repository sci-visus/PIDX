/*
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

/**
 * \file PIDX_memory_layout_data_structs.h
 *
 * \author Sidharth Kumar
 *
 * Generic data structs-
 * -- N dimensional patch
 * -- Set of N dimension patches
 * -- HZ encoded buffer
 * -- Aggregation buffer
 */

#ifndef __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H
#define __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H


/// Struct to store the row/column major chunk of data given by application
struct PIDX_Ndim_empty_patch_struct
{
  int rank;
  int is_boundary_patch;
  off_t offset[PIDX_MAX_DIMENSIONS];       ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  size_t size[PIDX_MAX_DIMENSIONS];         ///< size (extents) in each of the dimensions for the data chunk
};
typedef struct PIDX_Ndim_empty_patch_struct* Ndim_empty_patch;



/// Struct to store the restructured grid
struct PIDX_grid_struct
{
  size_t patch_size[PIDX_MAX_DIMENSIONS];
  int total_patch_count[PIDX_MAX_DIMENSIONS];
  Ndim_empty_patch* patch;
};
typedef struct PIDX_grid_struct* PIDX_restructured_grid;


/// Struct to store the row/column major chunk of data given by application
struct PIDX_patch_struct
{
  int particle_count;
  off_t offset[PIDX_MAX_DIMENSIONS];       ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension) in the 3D global space
  size_t size[PIDX_MAX_DIMENSIONS];         ///< size (extents) in each of the dimensions for the data chunk
  unsigned char* buffer;                                ///< the data buffer
};
typedef struct PIDX_patch_struct* PIDX_patch;



struct PIDX_source_patch_index_struct
{
  int rank;
  int index;
};
typedef struct PIDX_source_patch_index_struct PIDX_source_patch_index;


/// Struct to store a group of patches also called a super patch
/// The super patch is a rectilinear 3D patch
struct PIDX_super_patch_struct
{
  uint8_t is_boundary_patch;                            ///< 1 for boundary patch 0 otherwise

  uint32_t patch_count;                                 ///< Number of patches in the super patch
  PIDX_patch *patch;                                    ///< Pointer to the patches that comprises the super patch

  PIDX_source_patch_index *source_patch;                ///< Rank and index of all the patches
  int max_patch_rank;                                   ///< Rank of the process that holds this super patch
  PIDX_patch restructured_patch;                        ///< Pointer to the restructured (super) patch
};
typedef struct PIDX_super_patch_struct* PIDX_super_patch;


/// Struct to store the HZ encoded data and meta-data
struct PIDX_HZ_buffer_struct
{
  // Flag set if the HZ buffer comes from a boundary patch
  int is_boundary_HZ_buffer;                            ///< 1 for boundary patch 0 otherwise

  // HZ related meta data
  int **nsamples_per_level;                             ///< number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  unsigned long long *start_hz_index;                   ///< Starting HZ index at of the data at all the HZ levels
  unsigned long long *end_hz_index;                     ///< Ending HZ index at of the data at all the HZ levels

  // HZ encoded data (for every level)
  unsigned char** buffer;                               ///< data buffer at all the HZ levels
};
typedef struct PIDX_HZ_buffer_struct* HZ_buffer;


/// Struct to store aggregated data and meta-data
struct PIDX_HZ_Agg_buffer_struct
{
  int file_number;                                      ///< Target file number for the aggregator
  int var_number;                                       ///< Target variable number for the aggregator
  int sample_number;                                    ///< Target sample index for the aggregator

  int no_of_aggregators;
  int agg_f;
  int aggregator_interval;

  int num_idx_blocks;
  unsigned long long buffer_size;                                 ///< Aggregator buffer size
  unsigned long long compressed_buffer_size;                      ///< Aggregator buffer size after compression
  unsigned long long *compressed_block_size;                      ///< Compressed size of each block in the aggregator
  unsigned char* buffer;                                ///< The actual aggregator buffer
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

#endif
