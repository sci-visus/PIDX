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
 
#ifndef __PIDX_BLOCKS_H
#define __PIDX_BLOCKS_H

extern const int PIDX_default_bits_per_block;
extern const int PIDX_default_blocks_per_file;

struct PIDX_block_layout_struct
{
  /// Total number of Levels
  int levels;           
  
  /// Number of filled blocks per level
  int *hz_block_count_array;
  
  /// Indices of filled blocks
  int ** hz_block_number_array;
};
typedef struct PIDX_block_layout_struct* PIDX_block_layout;


///
int PIDX_blocks_initialize_layout(PIDX_block_layout layout, int maxh, int bits_per_block);


///
int PIDX_blocks_create_layout(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, const char* bitPattern, PIDX_block_layout layout);


///
void PIDX_blocks_print_layout(PIDX_block_layout layout);


///
int PIDX_blocks_is_block_present(int block_number, PIDX_block_layout layout);


///
int PIDX_blocks_find_negative_offset(int blocks_per_file, int block_number, PIDX_block_layout layout );


///
void PIDX_blocks_free_layout(PIDX_block_layout layout);

#endif
