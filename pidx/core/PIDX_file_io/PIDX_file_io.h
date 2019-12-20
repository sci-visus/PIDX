/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

#ifndef __PIDX_FILE_IO_H
#define __PIDX_FILE_IO_H 

/**
 * \file PIDX_file_io.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Move the data from the A aggregator processes to f files
 * 
 */

struct PIDX_file_io_struct
{
  idx_comm idx_c;

  idx_dataset idx;

  int fs_block_size;

  int first_index;
  int last_index;
};
typedef struct PIDX_file_io_struct* PIDX_file_io_id;


/// Creates the IO ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_file_io_id PIDX_file_io_init(idx_dataset idx_meta_data, idx_comm idx_c, int fs_block_size, int start_var_index, int end_var_index);



///
int PIDX_file_io_async_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* req, MPI_File *fp, char* filename_template);


PIDX_return_code PIDX_file_io_blocking_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, char* filename_template);


PIDX_return_code  PIDX_file_io_blocking_read(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, char* filename_template);

///
int PIDX_file_io_finalize(PIDX_file_io_id io_id);

#endif
