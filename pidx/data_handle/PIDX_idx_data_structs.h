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

/**
 * \file PIDX_idx_data_structs.h
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 *
 * IDX related data structs:
 * -- PIDX_variable
 * -- idx_dataset
 * -- idx_dataset_derived_metadata
 */

#ifndef __PIDX_IDX_DATA_STRUCTS_H
#define __PIDX_IDX_DATA_STRUCTS_H

#include "PIDX_memory_layout_data_structs.h"
#include "PIDX_timming_data_structs.h"
#include "PIDX_comm_structs.h"



struct PIDX_variable_struct
{  
  int is_particle;

  // General Info
  char var_name[1024];                                       ///< Variable name
  int vps;                                                   ///< values per sample, Vector(3), scalar(1), or n
  int bpv;                                                   ///< Number of bits each need
  PIDX_data_type type_name;                                  ///< Name of the type uint8, bob
  PIDX_data_layout data_layout;                              ///< Row major or column major

  // buffer (before, after HZ encoding phase)
  int sim_patch_count;                                       ///< The actual number of patches (application layout), most probably more than 1 in uintah
  PIDX_patch sim_patch[1024];                                ///< Pointer to the patches

  // buffer before aggregation
  int restructured_super_patch_count;                        ///< Number of groups of patches to be passed to aggregation phase
  PIDX_super_patch restructured_super_patch;                 ///< Pointer to the patch groups
  PIDX_super_patch chunked_super_patch;                      ///< Pointer to the patch group after block restructuring
  HZ_buffer hz_buffer;                                       ///< HZ encoded buffer of the patches

  int patch_group_count;
  PIDX_super_patch* rst_patch_group;
};
typedef struct PIDX_variable_struct* PIDX_variable;



struct PIDX_variable_group_struct
{
  int variable_index_tracker;
  int variable_count;
  int variable_tracker[256];
  PIDX_variable variable[256];

  int local_variable_index;
  int local_variable_count;

  //int agg_l_shared;
  int shared_start_layout_index;
  int shared_end_layout_index;
  int shared_layout_count;

  //int agg_l_nshared;
  int nshared_start_layout_index;
  int nshared_end_layout_index;
  int nshared_layout_count;

  int agg_level;
  PIDX_block_layout block_layout;
  PIDX_block_layout* block_layout_by_level;
};
typedef struct PIDX_variable_group_struct* PIDX_variable_group;






struct idx_metadata_cache_struct
{
  PIDX_metadata_cache meta_data_cache;
};
typedef struct idx_metadata_cache_struct* idx_metadata_cache;


/// idx_file
struct idx_file_struct
{
  enum PIDX_io_type io_type;              /// I/O format and layout we want to use
  
  int current_time_step;                  /// The current timestep selected

  int variable_count;                     /// The number of variables contained in the dataset
  int variable_group_count;
  int group_index_tracker;
  PIDX_variable_group variable_grp[16];
  
  char agg_list_filename[1024];

  char filename[1024];                    /// The idx file path
  char filename_partition[1024];
  char filename_template[1024];
  char filename_template_partition[1024];

  int first_tstep;                        /// Index of the frist timestep
  int last_tstep;                         /// Index of the last timestep

  int bits_per_block;                     /// Number of bits per block
  int blocks_per_file;                    /// Number of blocks per file
  size_t bounds[PIDX_MAX_DIMENSIONS];     /// Bounds of the dataset
  size_t box_bounds[PIDX_MAX_DIMENSIONS]; /// Bounds of the box query
  double physical_bounds[PIDX_MAX_DIMENSIONS];
  double physical_box_bounds[PIDX_MAX_DIMENSIONS];
  
  char bitSequence[512];
  char bitPattern[512];

  /// 0 No aggregation
  /// 1 Only aggregation
  int enable_agg;

  int compression_type;
  int compression_factor;
  float compression_bit_rate;
  size_t chunk_size[PIDX_MAX_DIMENSIONS];

  int file_zero_merge;

  /// 1 for little endian
  /// 0 for big endian
  enum PIDX_endian_type endian;                /// Endianess of the data

  /// 1 for flipping endian
  /// 0 for big endian
  int flip_endian;

  int agg_counter;

  int cached_ts;
};
typedef struct idx_file_struct* idx_dataset;


/// idx_dataset_derived_metadata
struct idx_dataset_derived_metadata_struct
{
  int pidx_version;
  char metadata_version[8];
  //int io_mode;


  PIDX_restructured_grid restructured_grid;

  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
  
  int fs_block_size;
  off_t start_fs_block;

  Agg_buffer **agg_buffer;

  int partition_count[PIDX_MAX_DIMENSIONS];
  int partition_size[PIDX_MAX_DIMENSIONS];
  int partition_offset[PIDX_MAX_DIMENSIONS];
  
  int start_layout_index;
  int end_layout_index;

  int layout_count;
  int reduced_res_from;
  int reduced_res_to;

  PIDX_time time;

  int raw_io_pipe_length;

  int aggregator_multiplier;

  //int shared_block_level;
  int total_partiton_level;

  int **block_bitmap;
  int ***block_offset_bitmap;

  int variable_pipe_length;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;



struct idx_debug_info_struct
{
  int simulate_rst_io;
  int simulate_rst;

  //FILE *rst_dump_fp;
  //int dump_rst_info;
  //char rst_dump_dir_name[512];

  //FILE *agg_dump_fp;
  //int dump_agg_info;
  //char agg_dump_dir_name[512];

  //FILE *io_dump_fp;
  //int dump_io_info;
  //char io_dump_dir_name[512];

  //FILE *process_size_and_offset_dump_fp;
  //int dump_process_state;
  //char process_state_dump_dir_name[512];

  FILE *local_dump_fp;
  FILE *mpi_dump_fp;
  int state_dump;

  int debug_rst;                               ///< Debug restructuring phase, works only on the test application
  int debug_hz;                                ///< Debug HZ encoding phase, works only on the test application

  /// Flags set by user
  int debug_do_rst;                            ///< User controlled flag to activate/deactivate restructuring phase
  int debug_do_chunk;                          ///< User controlled flag to activate/deactivate chunking phase
  int debug_do_compress;                       ///< User controlled flag to activate/deactivate compression
  int debug_do_hz;                             ///< User controlled flag to activate/deactivate hz encoding phase
  int debug_do_agg;                            ///< User controlled flag to activate/deactivate aggregation phase
  int debug_do_io;                             ///< User controlled flag to activate/deactivate I/O phase
};
typedef struct idx_debug_info_struct* idx_debug;

#endif
