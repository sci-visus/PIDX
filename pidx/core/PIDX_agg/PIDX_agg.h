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

/**
 * \file PIDX_agg.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Move the hz encoded data sitting on n' cores
 * to A aggregator cores
 * 
 */

#ifndef __PIDX_AGG_H
#define __PIDX_AGG_H 

struct PIDX_agg_struct
{
  MPI_Win win;

  MPI_Win shard_block_win;

  idx_comm idx_c;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int gi;
  int fi;
  int li;
  //int lvi;

  int ***agg_r;
};

struct PIDX_agg_struct;
typedef struct PIDX_agg_struct* PIDX_agg_id;


/// Creates the Aggregation ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int fi, int li);



PIDX_return_code PIDX_agg_create_randomized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset, int var_offset, int file_status);

///
PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout);



///
PIDX_return_code PIDX_agg_buf_create_localized_aggregation(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int i1, int j1, int file_status);


PIDX_return_code PIDX_agg_buf_create_multiple_level(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset, int var_offset, int file_status);


PIDX_return_code PIDX_agg_buf_destroy(Agg_buffer agg_buffer);

///
PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout local_block_layout, int PIDX_MODE);


///
PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id agg_id, PIDX_block_layout local_block_layout);


PIDX_return_code PIDX_agg_create_global_partition_localized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset);


PIDX_return_code PIDX_agg_create_local_partition_localized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset);

///
PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id);

#endif //__PIDX_AGG_H
