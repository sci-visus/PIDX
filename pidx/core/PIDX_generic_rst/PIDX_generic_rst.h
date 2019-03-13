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

/**
 * \file PIDX_generic_rst.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Restructuring data from n cores to n' (n' <= n)
 * while keeping the data in mult-dimensional
 * application layout
 *
 */

#ifndef __PIDX_GENERIC_RST_NEW_H
#define __PIDX_GENERIC_RST_NEW_H


///
/// \brief The PIDX_generic_rst_struct struct is Struct for restructuring ID
///
struct PIDX_generic_rst_struct
{
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;


  idx_comm idx_c;


  idx_debug idx_dbg;

  int group_index;
  int first_index;
  int last_index;

  //dimension of the power-two volume imposed patch
  unsigned long long reg_patch_size[PIDX_MAX_DIMENSIONS];
  int reg_patch_grp_count;
  Ndim_patch_group* reg_patch_grp;

  unsigned long long *rank_r_offset;
  unsigned long long *rank_r_count;
};
typedef struct PIDX_generic_rst_struct* PIDX_generic_rst_id;



/*
 * Implementation in PIDX_generic_rst.c
 */
/// Creates an restructuring file ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_generic_rst_id The identifier associated with the task
PIDX_generic_rst_id PIDX_generic_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index);




///
/// \brief PIDX_generic_rst_finalize Tear down whatever was calculated for this particular combination of dimensions and bounds
/// \param id
/// \return
///
PIDX_return_code PIDX_generic_rst_finalize(PIDX_generic_rst_id id);




#if 1
/*
 * Implementation in PIDX_generic_rst_meta_data.c
 */
///
/// \brief PIDX_generic_rst_attach_restructuring_box
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_meta_data_create(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_meta_data_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_meta_data_write(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_meta_data_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_meta_data_destroy(PIDX_generic_rst_id rst_id);



/*
 * Implementation in PIDX_generic_rst_buffer.c
 */
///
/// \brief PIDX_generic_rst_buf_create Create the appropriate data structs to hold restructured output data
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_create(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_destroy(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_aggregate_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_aggregate_buf_create(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_aggregate_buf_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_aggregate_buf_destroy(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_buf_aggregate
/// \param rst_id
/// \param MODE
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_aggregate(PIDX_generic_rst_id rst_id, int MODE);



/*
 * Implementation in PIDX_generic_rst_io.c
 */
///
/// \brief PIDX_generic_rst_buf_aggregate_and_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_aggregate_and_write(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_buf_read_and_aggregate
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_read_and_aggregate(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_buf_aggregated_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_aggregated_write(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_buf_aggregated_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_buf_aggregated_read(PIDX_generic_rst_id rst_id);



///
/// \brief PIDX_generic_rst_forced_raw_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_forced_raw_read(PIDX_generic_rst_id rst_id);
/*
 * Implementation in PIDX_generic_rst_write.c
 */


///
/// \brief PIDX_generic_rst_staged_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_staged_write(PIDX_generic_rst_id rst_id);



/*
 * Implementation in PIDX_generic_rst_read.c
 */
///
/// \brief PIDX_generic_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_generic_rst_read(PIDX_generic_rst_id rst_id);



/*
 * Implementation in PIDX_generic_rst_verify.c
 */
///
/// \brief HELPER_rst
/// \param rst_id
/// \return
///
PIDX_return_code HELPER_generic_rst(PIDX_generic_rst_id rst_id);
#endif
#endif // __PIDX_generic_rst_NEW_H
