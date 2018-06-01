/*
 * BSD 3-Clause License
 * 
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

#ifndef __PIDX_HEADER_IO_H
#define __PIDX_HEADER_IO_H

struct PIDX_header_io_struct;
typedef struct PIDX_header_io_struct* PIDX_header_io_id;


///
/// \brief PIDX_header_io_init
/// \param idx_meta_data
/// \param idx_d
/// \param idx_c
/// \param first_index
/// \param last_index
/// \return
///
PIDX_header_io_id PIDX_header_io_init(idx_dataset idx_meta_data, idx_comm idx_c, idx_blocks idx_b, PIDX_restructured_grid restructured_grid, int fs_block_size, int first_index, int last_index);



///
/// \brief PIDX_header_io_write_idx
/// \param header_io
/// \param data_set_path
/// \param current_time_step
/// \return
///
PIDX_return_code PIDX_header_io_global_idx_write (PIDX_header_io_id header_io, char* data_set_path);




PIDX_return_code PIDX_header_io_partition_idx_write (PIDX_header_io_id header_io, char* data_set_path);



PIDX_return_code PIDX_header_io_raw_idx_write (PIDX_header_io_id header_io, char* data_set_path);


///
/// \brief PIDX_header_io_idx_file_create
/// \param header_io_id
/// \param block_layout
/// \param filename_template
/// \return
///
int PIDX_header_io_idx_file_create(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, char* filename_template);



///
/// \brief PIDX_header_io_idx_file_write
/// \param header_io_id
/// \param block_layout
/// \param file_name
/// \param file_name_template
/// \param mode
/// \return
///
PIDX_return_code PIDX_header_io_idx_file_write(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, char* file_name_template, int mode);




///
/// \brief PIDX_header_io_raw_dir_create
/// \param header_io_id
/// \param file_name
/// \return
///
int PIDX_header_io_raw_dir_create(PIDX_header_io_id header_io_id, char* file_name);



///
/// \brief PIDX_header_io_finalize
/// \param header_io
/// \return
///
int PIDX_header_io_finalize(PIDX_header_io_id header_io);

#endif
