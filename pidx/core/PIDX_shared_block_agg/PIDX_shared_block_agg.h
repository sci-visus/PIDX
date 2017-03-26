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

#ifndef __PIDX_SHARED_BLOCK_AGG_H
#define __PIDX_SHARED_BLOCK_AGG_H


struct PIDX_shared_block_agg_struct;
typedef struct PIDX_shared_block_agg_struct* PIDX_shared_block_agg_id;


/// Creates the Aggregation ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_shared_block_agg_id PIDX_shared_block_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int fi, int li);


PIDX_return_code PIDX_shared_block_agg_buf_create(PIDX_shared_block_agg_id agg_id, Agg_buffer ab);


PIDX_return_code PIDX_shared_block_agg_buf_destroy(PIDX_shared_block_agg_id id, Agg_buffer ab);

///
PIDX_return_code PIDX_shared_block_agg_global_and_local(PIDX_shared_block_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int PIDX_MODE);


///
PIDX_return_code PIDX_shared_block_agg_finalize(PIDX_shared_block_agg_id agg_id);

#endif
