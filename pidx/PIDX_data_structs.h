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
 * \file PIDX_data_structs.h
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 *
 * Basic data structures
 * -- N dimensional Box
 * -- Set of N dimension boxes
 * -- HZ encoded buffer
 * -- Aggregation buffer
 * -- IO_MODE
 * -- PIDX_variable;
 * -- idx_dataset;
 * -- idx_dataset_derived_metadata;
 *
 */

#ifndef __PIDX_DATA_STRUCTS_H
#define __PIDX_DATA_STRUCTS_H

/// Struct to store the row/column major chunk of data given by application
struct PIDX_Ndim_box_struct
{
  /// offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  int64 Ndim_box_offset[PIDX_MAX_DIMENSIONS];
  
  /// size (extents) in each of the dimensions for the data chunk
  int64 Ndim_box_size[PIDX_MAX_DIMENSIONS];
  
  /// the data buffer
  unsigned char* Ndim_box_buffer;
};
typedef struct PIDX_Ndim_box_struct* Ndim_box;


/// Struct to store a group of Ndim_buffer
struct PIDX_Ndim_box_group_struct
{
  /// decide the type of the group
  int box_group_type;
  
  /// how many Ndim_buffer are there in the group
  int box_count;
  
  /// Pointer to all the Ndim_buffer
  Ndim_box *box;
  
  ///
  int *source_box_rank;
  
  ///
  int max_box_rank;
  
  /// If restructuring used then this contains the offset of the power-two block
  int64 enclosing_box_offset[PIDX_MAX_DIMENSIONS];
  
  /// If restructuring used then this contains the extents of the power-two block
  int64 enclosing_box_size[PIDX_MAX_DIMENSIONS];
};
typedef struct PIDX_Ndim_box_group_struct* Ndim_box_group; 


/// Struct to store the HZ encoded data and meta-data
struct PIDX_HZ_buffer_struct
{
  /// starting HZ level
  int HZ_level_from;
  
  /// ending HZ level
  int HZ_level_to;
  
  /// number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  int64 *samples_per_level;
  
  /// Starting HZ index at of the data at all the HZ levels
  int64 *start_hz_index;
  
  /// Ending HZ index at of the data at all the HZ levels
  int64 *end_hz_index;
  
  /// HZ indices of the data (used only when no restructuring phsae is used)
  int64 *buffer_index;
  
  /// 
  int *missing_block_count_per_level;
  
  ///
  int **missing_block_index_per_level;
  
  /// data buffer at all the HZ levels
  unsigned char** buffer;
};
typedef struct PIDX_HZ_buffer_struct* HZ_buffer;


/// Struct to store aggregated data and meta-data
struct PIDX_HZ_Agg_buffer_struct
{
  /// Target file number for the aggregator
  int file_number;
  
  /// Target variable number for the aggregator
  int var_number;
  
  /// Target sample index for the aggregator
  int sample_number;
  
  /// Aggregator buffer size
  uint64 buffer_size;

  ///
  int ***rank_holder;
  
  /// The actual aggregator buffer
  unsigned char* buffer;
  
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

/// IO_MODE
/// \param PIDX_READ Read mode
/// \param PIDX_WRITE Write mode
///
enum IO_MODE { PIDX_READ, PIDX_WRITE};

/// PIDX_variable
struct PIDX_variable_struct
{
  int dump_meta_data_ON;
  
  char var_name[1024];
  int values_per_sample;
  int bits_per_value;
  char type_name[1024];
  PIDX_data_layout data_layout;

  int patch_count;
  Ndim_box patch[1024];
  HZ_buffer HZ_patch[1024];
  
  int patch_group_count;
  Ndim_box_group* patch_group_ptr;
  
  block_layout* VAR_global_block_layout;
  int *VAR_blocks_per_file;
  int VAR_existing_file_count;
  int *VAR_existing_file_index;
};
typedef struct PIDX_variable_struct* PIDX_variable;

/// idx_file
struct idx_file_struct
{
  int current_time_step;
  int variable_count;
  int variable_index_tracker;
  
  PIDX_variable variable[1024];
  
  char filename[1024];
  int bits_per_block;
  int blocks_per_file;
  int64* global_bounds;
  double transform[16];
  char bitSequence[512];
  char bitPattern[512];
  
  /// Depends on the time step
  char filename_template[1024];
};
typedef struct idx_file_struct* idx_dataset;

/// idx_dataset_derived_metadata
struct idx_dataset_derived_metadata_struct
{
  int *file_bitmap;
  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
  
  int fs_block_size;
  off_t start_fs_block;
  
  block_layout* global_block_layout;
  int *existing_blocks_index_per_file;
  int existing_file_count;
  int *existing_file_index;
  
  int aggregation_factor;
  Agg_buffer agg_buffer;
  int dump_agg_info;
  char agg_dump_dir_name[512];
  
  int color;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;

#endif
