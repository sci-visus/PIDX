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
 * \file PIDX_memory_layout_data_structs.h
 *
 * \author Sidharth Kumar
 *
 * Generic data structs-
 * -- N dimensional Box
 * -- Set of N dimension boxes
 * -- HZ encoded buffer
 * -- Aggregation buffer
 */

#ifndef __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H
#define __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H


/// Struct to store the row/column major chunk of data given by application
struct PIDX_Ndim_box_struct
{
  int64_t Ndim_box_offset[PIDX_MAX_DIMENSIONS];         ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  int64_t Ndim_box_size[PIDX_MAX_DIMENSIONS];           ///< size (extents) in each of the dimensions for the data chunk
  unsigned char* Ndim_box_buffer;                       ///< the data buffer
};
typedef struct PIDX_Ndim_box_struct* Ndim_box;


/// Struct to store a group of Ndim_buffer
struct PIDX_Ndim_box_group_struct
{
  int box_group_type;                                   ///< decide the type of the group 
  int box_count;                                        ///< how many Ndim_buffer are there in the group
  Ndim_box *box;                                        ///< Pointer to all the Ndim_buffer
  int *source_box_rank;                                 ///<
  int max_box_rank;                                     ///<
  int64_t enclosing_box_offset[PIDX_MAX_DIMENSIONS];    ///< If restructuring used then this contains the offset of the power-two block
  int64_t enclosing_box_size[PIDX_MAX_DIMENSIONS];      ///< If restructuring used then this contains the extents of the power-two block
};
typedef struct PIDX_Ndim_box_group_struct* Ndim_box_group; 


/// Struct to store the HZ encoded data and meta-data
struct PIDX_HZ_buffer_struct
{
  int HZ_level_from;                                    ///< starting HZ level
  int HZ_level_to;                                      ///< ending HZ level
  int64_t *samples_per_level;                           ///< number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  int64_t *start_hz_index;                              ///< Starting HZ index at of the data at all the HZ levels
  int64_t *end_hz_index;                                ///< Ending HZ index at of the data at all the HZ levels
  int64_t *buffer_index;                                ///< HZ indices of the data (used only when no restructuring phsae is used)
  int *missing_block_count_per_level;                   ///< 
  int **missing_block_index_per_level;                  ///<
  unsigned char** buffer;                               ///< data buffer at all the HZ levels
};
typedef struct PIDX_HZ_buffer_struct* HZ_buffer;


/// Struct to store aggregated data and meta-data
struct PIDX_HZ_Agg_buffer_struct
{
  int file_number;                                      ///< Target file number for the aggregator
  int var_number;                                       ///< Target variable number for the aggregator
  int sample_number;                                    ///< Target sample index for the aggregator
  uint64_t buffer_size;                                 ///< Aggregator buffer size
  uint64_t compressed_buffer_size;                      ///< Aggregator buffer size after compression
  uint64_t *compressed_block_size;                      ///< Compressed size of each block in the aggregator
  int ***rank_holder;                                   ///<
  unsigned char* buffer;                                ///< The actual aggregator buffer
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

#endif
