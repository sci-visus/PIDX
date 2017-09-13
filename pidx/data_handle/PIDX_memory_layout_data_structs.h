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
struct PIDX_Ndim_empty_patch_struct
{
  int rank;
  int edge;
  unsigned long long offset[PIDX_MAX_DIMENSIONS];       ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  unsigned long long size[PIDX_MAX_DIMENSIONS];         ///< size (extents) in each of the dimensions for the data chunk
};
typedef struct PIDX_Ndim_empty_patch_struct* Ndim_empty_patch;



/// Struct to store the restructured grid
struct PIDX_grid_struct
{
  unsigned long long patch_size[PIDX_MAX_DIMENSIONS];
  int total_patch_count[PIDX_MAX_DIMENSIONS];
  Ndim_empty_patch* patch;
};
typedef struct PIDX_grid_struct* PIDX_restructured_grid;


/// Struct to store the row/column major chunk of data given by application
struct PIDX_patch_struct
{
  unsigned long long offset[PIDX_MAX_DIMENSIONS];       ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  unsigned long long size[PIDX_MAX_DIMENSIONS];         ///< size (extents) in each of the dimensions for the data chunk
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

  PIDX_source_patch_index *source_patch;                ///< Rank of all the patches
  int max_patch_rank;                                   ///< Rank of the process that holds this super patch
  PIDX_patch reg_patch;                                 ///< Pointer to the restructured (super) patch
};
typedef struct PIDX_super_patch_struct* PIDX_super_patch;


/// Struct to store the HZ encoded data and meta-data
struct PIDX_HZ_buffer_struct
{
  int is_boundary_HZ_buffer;                            ///< 1 for boundary patch 0 otherwise
  int **nsamples_per_level;                             ///< number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  unsigned long long *start_hz_index;                              ///< Starting HZ index at of the data at all the HZ levels
  unsigned long long *end_hz_index;                                ///< Ending HZ index at of the data at all the HZ levels
  unsigned long long *buffer_index;                                ///< HZ indices of the data (used only when no restructuring phsae is used)
  unsigned char** buffer;                               ///< data buffer at all the HZ levels
  unsigned long long *compressed_buffer_size;
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
