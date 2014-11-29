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
 * \file Generic_data_structs.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Defining some of the very basic data struct
 * --N dimensional Box
 * --Set of N dimension boxes
 * --HZ encoded buffer
 * --Aggregation buffer
 *
 */

#ifndef __GENERIC_DATA_STRUCTS_H
#define __GENERIC_DATA_STRUCTS_H


/// Struct to store the row/column major chunk of data given by application
struct PIDX_Ndim_box_struct
{
  /// offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  long long Ndim_box_offset[PIDX_MAX_DIMENSIONS];
  
  /// size (extents) in each of the dimensions for the data chunk
  long long Ndim_box_size[PIDX_MAX_DIMENSIONS];
  
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
  long long enclosing_box_offset[PIDX_MAX_DIMENSIONS];
  
  /// If restructuring used then this contains the extents of the power-two block
  long long enclosing_box_size[PIDX_MAX_DIMENSIONS];
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
  long long *samples_per_level;
  
  /// Starting HZ index at of the data at all the HZ levels
  long long *start_hz_index;
  
  /// Ending HZ index at of the data at all the HZ levels
  long long *end_hz_index;
  
  /// HZ indices of the data (used only when no restructuring phsae is used)
  long long *buffer_index;
  
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
  int buffer_size;
  
  /// The actual aggregator buffer
  unsigned char* buffer;
  
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

#endif