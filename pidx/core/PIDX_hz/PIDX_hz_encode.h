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

/**
 * \file PIDX_hz_encode.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Incore restructuring of data, transforming the mult-dimensional 
 * application layout into the hierarchial Z order layout of the
 * IDX format.
 * 
 */

#ifndef __PIDX_HZ_ENCODE_H
#define __PIDX_HZ_ENCODE_H

struct PIDX_hz_encode_struct
{
  idx_dataset idx;


  idx_comm idx_c;


  idx_debug idx_dbg;

  PIDX_metadata_cache meta_data_cache;
  //idx_metadata_cache cache;

  int fs_block_size;

  int** index;

  int first_index;
  int last_index;

  int resolution_from;
  int resolution_to;
};
typedef struct PIDX_hz_encode_struct* PIDX_hz_encode_id;


/// Creates the HZ encoding file ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_metadata_cache meta_data_cache, int fs_block_size, int start_var_index, int end_var_index);




///
/// \brief PIDX_hz_encode_meta_data_create
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_meta_data_destroy
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_buf_create
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id);




///
/// \brief PIDX_hz_encode_write_inverse
/// \param id
/// \param start_hz_index
/// \param end_hz_index
/// \return
///
PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index);


///
/// \brief PIDX_hz_encode_write
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_fast_write
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_fast_write(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_read
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id);




///
/// \brief PIDX_hz_encode_buf_destroy
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_finalize
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_set_resolution
/// \param id
/// \param resolution_from
/// \param resolution_to
/// \return
///
PIDX_return_code PIDX_hz_encode_set_resolution(PIDX_hz_encode_id id, int resolution_from, int resolution_to);





///
int PIDX_file_io_per_process(PIDX_hz_encode_id io_id, PIDX_block_layout block_layout, int MODE);



///
/// \brief HELPER_Hz_encode
/// \param id
/// \return
///
int HELPER_Hz_encode(PIDX_hz_encode_id id);

#endif
