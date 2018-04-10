/*
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
