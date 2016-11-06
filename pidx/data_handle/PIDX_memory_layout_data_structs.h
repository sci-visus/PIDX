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
 * -- N dimensional patch
 * -- Set of N dimension patches
 * -- HZ encoded buffer
 * -- Aggregation buffer
 */

#ifndef __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H
#define __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H


/// Struct to store the row/column major chunk of data given by application
struct PIDX_Ndim_patch_struct
{
  unsigned long long offset[PIDX_MAX_DIMENSIONS];       ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  unsigned long long size[PIDX_MAX_DIMENSIONS];         ///< size (extents) in each of the dimensions for the data chunk
  unsigned char* buffer;                                ///< the data buffer
};
typedef struct PIDX_Ndim_patch_struct* Ndim_patch;


/// Struct to store a group of Ndim_buffer
struct PIDX_Ndim_patch_group_struct
{
  int data_source;
  int type;                                             ///< decide the type of the group
  int count;                                            ///< how many Ndim_buffer are there in the group
  Ndim_patch *patch;                                    ///< Pointer to all the Ndim_buffer
  int *source_patch_rank;                               ///<
  int max_patch_rank;                                   ///<
  Ndim_patch reg_patch;                                 ///< Pointer to the reg patch
};
typedef struct PIDX_Ndim_patch_group_struct* Ndim_patch_group;


struct PIDX_source_patch_index_struct{
  int rank;
  int index;
};
typedef struct PIDX_source_patch_index_struct PIDX_source_patch_index;


struct PIDX_Ndim_multi_patch_group_struct
{
  int data_source;
  int type;                                             ///< decide the type of the group
  int count;                                            ///< how many Ndim_buffer are there in the group
  Ndim_patch *patch;                                    ///< Pointer to all the Ndim_buffer
  PIDX_source_patch_index *source_patch;
  int max_patch_rank;                                   ///<
  Ndim_patch reg_patch;                                 ///< Pointer to the reg patch
};
typedef struct PIDX_Ndim_multi_patch_group_struct* Ndim_multi_patch_group;


/// Struct to store the HZ encoded data and meta-data
struct PIDX_HZ_buffer_struct
{
  int type;                                             ///< decide the type of the group
  int **nsamples_per_level;                             ///< number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  unsigned long long *samples_per_level;                           ///< number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  unsigned long long *start_hz_index;                              ///< Starting HZ index at of the data at all the HZ levels
  unsigned long long *end_hz_index;                                ///< Ending HZ index at of the data at all the HZ levels
  unsigned long long *buffer_index;                                ///< HZ indices of the data (used only when no restructuring phsae is used)
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
