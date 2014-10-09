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

#include "Generic_data_structs.h"
#include "PIDX_utils.h"

extern const int PIDX_default_bits_per_block;
extern const int PIDX_default_blocks_per_file;

struct block_layout
{
  int    levels;                          // Total number of Levels
  int    *hz_block_count_array;           // Number of filled blocks per level
  int    ** hz_block_number_array;        // Indices of filled blocks
};
typedef struct block_layout block_layout;

int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, const char* bitPattern, block_layout* layout);

int initialize_block_layout(block_layout* layout, int maxh, int bits_per_block);

void print_block_layout(block_layout* layout);

int find_block_negative_offset(int blocks_per_file, int block_number, block_layout* layout );

int is_block_present(int block_number, block_layout* layout);

void destroyBlockBitmap(block_layout* layout);

#endif