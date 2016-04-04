/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

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


struct PIDX_timming_struct
{
   int header_counter;
   int variable_counter;
   double file_create_time;
   double partition_start_time;
   double partition_end_time;
   double populate_idx_start_time;
   double populate_idx_end_time;
   double sim_start, sim_end;
   double *write_init_start, *write_init_end;
   double *startup_start, *startup_end;
   double *init_start, *init_end;
   double *rst_start, *rst_end;
   double *rst_io_start, *rst_io_end;
   double *rst_meta_data_start_io, *rst_meta_data_end_io;
   double *hz_start, *hz_end;
   double *cleanup_start, *cleanup_end;
   double *finalize_start, *finalize_end;
   double *chunk_start, *chunk_end;
   double *buffer_start, *buffer_end;
   double *compression_start, *compression_end;
   double **agg_start, **agg_end;
   double **agg_buf_start, **agg_buf_end;

   double **windows_start, **windows_end;
   double **first_fence_start, **first_fence_end;
   double **second_fence_start, **second_fence_end;
   double **fence_free;
   double **io_per_process_start, **io_per_process_end;
   double **io_start, **io_end;
};
typedef struct PIDX_timming_struct* PIDX_time;


///
struct PIDX_variable_struct
{
  // Metadata
  int dump_meta_data_ON;                                                ///< Counter set if meta data dumping activated
  
  int io_state;

  // General Info
  char var_name[1024];                                                  ///< Variable name
  int values_per_sample;                                                ///< Vector(3), scalar(1), or n
  int bits_per_value;                                                   ///< Number of bits each need
  PIDX_data_type type_name;                                                 ///< Name of the type uint8, bob
  PIDX_data_layout data_layout;                                         ///< Row major or column major

  // Memory layout (before, after HZ encoding phase)
  int sim_patch_count;                                                  ///< The actual number of patches (application layout), most probably more than 1 in uintah
  Ndim_patch sim_patch[1024];                                           ///< Pointer to the patches
  HZ_buffer* hz_buffer;                                                 ///< HZ encoded buffer of the patches

  // Memory layout before aggregation 
  int patch_group_count;                                                ///< Number of groups of patches to be passed to aggregation phase
  Ndim_patch_group* rst_patch_group;                                    ///< Pointer to the patch groups
  Ndim_patch_group* chunk_patch_group;                                  ///< Pointer to the patch group after block restructuring

  // Block level layout
  PIDX_block_layout global_block_layout;                               ///< Block layout, specifically when variables might have different extents in the domain
  PIDX_block_layout* block_layout_by_level;                            ///< Block layout, specifically when variables might have different extents in the domain

  //Compression related
  int lossy_compressed_block_size;                                      ///< The expected size of the compressed buffer
};
typedef struct PIDX_variable_struct* PIDX_variable;


/// idx_file
struct idx_file_struct
{
  int current_time_step;                                                ///< Time step tracker
  
  int variable_count;
  int variable_index_tracker;
  PIDX_variable variable[1024];
  
  char filename[1024];
  int bits_per_block;
  int blocks_per_file;
  int64_t bounds[PIDX_MAX_DIMENSIONS];
  double transform[16];
  char bitSequence[512];
  char bitPattern[512];
  char filename_template[1024];                                         ///< Depends on the time step
  
  int64_t reg_patch_size[PIDX_MAX_DIMENSIONS];
  
  int compression_type;                                               ///< counter to enable/disable (1/0) compression
  int enable_rst;                                               ///< counter to enable/disable (1/0) compression
  /// 0 No aggregation
  /// 1 Hybrid aggregation
  /// 2 Only aggregation
  int enable_agg;                                               ///< counter to enable/disable (1/0) compression

  int compression_factor;
  int compression_bit_rate;
  int64_t chunk_size[PIDX_MAX_DIMENSIONS];                              ///< size of the block at which compression is applied eg. (4x4x4)
                                                                        ///< the current compression schemes only work in three dimensions
  int64_t chunked_bounds[PIDX_MAX_DIMENSIONS];                ///< Compressed global extents
};
typedef struct idx_file_struct* idx_dataset;


/// idx_dataset_derived_metadata
struct idx_dataset_derived_metadata_struct
{
  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
  
  int fs_block_size;
  off_t start_fs_block;

  Agg_buffer **agg_buffer;

  int dump_agg_info;
  char agg_dump_dir_name[512];
  int dump_io_info;
  char io_dump_dir_name[512];

  int color;
  int idx_count[PIDX_MAX_DIMENSIONS];          ///< Number of idx files in each dimensions

  int var_pipe_length;
  int parallel_mode;

  //extents of meta-data

  int64_t *rank_r_offset;                                                   ///< Offset of variables in each dimension
  int64_t *rank_r_count;                                                    ///< Count of variables in each dimension
  
  //int staged_aggregation;
  int agg_type;
  int start_layout_index;
  int end_layout_index;
  int perm_layout_count;
  int layout_count;
  int reduced_res_from;
  int reduced_res_to;

  PIDX_time time;

  int data_core_count;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;



struct idx_debug_info_struct
{
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
