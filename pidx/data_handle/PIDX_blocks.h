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


struct PIDX_block_layout_struct
{
  int resolution_from;
  int resolution_to;

  int efc;  ///existing file count

  int *file_bitmap;
  int *lbi;         /// index of last block in a file
  int *bcpf;        /// block count per file
  int *existing_file_index;
  int *inverse_existing_file_index;
  int *file_index;

  /// Indices of filled blocks
  int **hz_block_number_array;
};
typedef struct PIDX_block_layout_struct* PIDX_block_layout;


/// Initialize the block layout structure.
/// \param  layout create block layout
/// \param maxh  maximum HZ levels
/// \param bits_per_block  block size = pow (2, bits_per_block)
/// \return Error code
int PIDX_blocks_initialize_layout(PIDX_block_layout layout, int resolution_from, int resolution_to, int maxh, int bits_per_block);


///
int PIDX_blocks_create_layout(int bounding_box[2][5], int maxH, int bits_per_block, const char* bitPattern, PIDX_block_layout layout, int res_from, int res_to);


///
/// \brief PIDX_blocks_print_layout
/// \param layout
/// \param bits_per_block
/// \param blocks_per_file
///
void PIDX_blocks_print_layout(PIDX_block_layout layout, int bits_per_block);


///
/// \brief PIDX_blocks_is_block_present
/// \param block_number
/// \param bits_per_block
/// \param layout
/// \return
///
int PIDX_blocks_is_block_present(int block_number, int bits_per_block, PIDX_block_layout layout);


///
/// \brief PIDX_blocks_find_negative_offset
/// \param blocks_per_file
/// \param bits_per_block
/// \param block_number
/// \param layout
/// \return
///
int PIDX_blocks_find_negative_offset(int blocks_per_file, int bits_per_block, int block_number, PIDX_block_layout layout);



///
/// \brief PIDX_blocks_free_layout
/// \param layout
///
void PIDX_blocks_free_layout(int bits_per_block, int maxh, PIDX_block_layout layout);

void PIDX_free_layout(PIDX_block_layout layout);

#endif
