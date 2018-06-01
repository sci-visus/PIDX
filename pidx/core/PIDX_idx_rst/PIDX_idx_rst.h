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
 * \file PIDX_idx_rst.h
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Restructuring data from n cores to n' (n' <= n)
 * while keeping the data in mult-dimensional
 * application layout
 *
 */

#ifndef __PIDX_IDX_RST_NEW_H
#define __PIDX_IDX_RST_NEW_H


//Struct for restructuring ID
struct PIDX_idx_rst_struct
{
  idx_dataset idx_metadata;
  idx_debug idx_debug_metadata;
  idx_comm idx_c;

  PIDX_restructured_grid restructured_grid;

  int first_index;
  int last_index;

  int intersected_restructured_super_patch_count;
  PIDX_super_patch* intersected_restructured_super_patch;

  int sim_max_patch_count;
  uint64_t* sim_multi_patch_r_count;
  uint64_t* sim_multi_patch_r_offset;

  int maximum_neighbor_count;
};
typedef struct PIDX_idx_rst_struct* PIDX_idx_rst_id;


/*
 * Implementation in PIDX_rst.c
 */
///
/// \brief PIDX_idx_rst_init
/// \param idx_meta_data
/// \param idx_derived_ptr
/// \param var_start_index
/// \param var_end_index
/// \return
///
PIDX_idx_rst_id PIDX_idx_rst_init( idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_restructured_grid restructured_grid, int var_start_index, int var_end_index);




///
/// \brief PIDX_idx_rst_finalize
/// \param id
/// \return
///
PIDX_return_code PIDX_idx_rst_finalize(PIDX_idx_rst_id id);



/*
 * Implementation in PIDX_idx_rst_meta_data.c
 */
///
/// \brief PIDX_idx_rst_meta_data_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_meta_data_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_meta_data_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_destroy(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_buffer.c
 */
///
/// \brief PIDX_idx_rst_buf_create Create the appropriate data structs to hold restructured output data
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_destroy(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_aggregate_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_aggregate_buf_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_aggregate_buf_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_aggregate_buf_destroy(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_aggregate
/// \param rst_id
/// \param MODE
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregate(PIDX_idx_rst_id rst_id, int MODE);



/*
 * Implementation in PIDX_idx_rst_io.c
 */
///
/// \brief PIDX_idx_rst_buf_aggregated_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregated_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_aggregated_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregated_read(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_write.c
 */
///
/// \brief PIDX_idx_rst_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_staged_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_staged_write(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_read.c
 */
///
/// \brief PIDX_idx_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_read(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_forced_raw_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_forced_raw_read(PIDX_idx_rst_id rst_id);

/*
 * Implementation in PIDX_rst_verify.c
 */
///
/// \brief HELPER_rst
/// \param rst_id
/// \return
///
PIDX_return_code HELPER_idx_rst(PIDX_idx_rst_id rst_id);
#endif
