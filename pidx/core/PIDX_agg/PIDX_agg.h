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

#ifndef __PIDX_AGG_H
#define __PIDX_AGG_H 

#include "../../PIDX_inc.h"

#include "./PIDX_local_agg/PIDX_local_agg.h"
//#include "./PIDX_global_agg/PIDX_global_agg.h"

struct PIDX_agg_struct;
typedef struct PIDX_agg_struct* PIDX_agg_id;


/// Creates the Aggregation ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int init_index, int start_var_index, int end_var_index);



#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id aggregator id
/// \param comm the communicator
/// \return error code
PIDX_return_code PIDX_agg_set_communicator(PIDX_agg_id io_id, MPI_Comm comm);


PIDX_return_code PIDX_agg_set_global_communicator(PIDX_agg_id agg_id, MPI_Comm comm);
#endif


PIDX_return_code PIDX_create_local_aggregation_comm(PIDX_agg_id agg_id);
PIDX_return_code PIDX_destroy_local_aggregation_comm(PIDX_agg_id agg_id);


///
PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout global_block_layout, PIDX_block_layout local_block_layout);



///
PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id agg_id, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout);


///
PIDX_return_code PIDX_agg_buf_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout, int i1, int j1);


///
PIDX_return_code PIDX_agg_buf_create_multiple_level(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout, int i1, int j1, int file_status);


///
PIDX_return_code PIDX_agg_buf_destroy(PIDX_agg_id agg_id, Agg_buffer agg_buffer);



///
PIDX_return_code PIDX_agg(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout local_block_layout, int PIDX_MODE, int v1, int bi);


///
PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout,  int PIDX_MODE, int vi, int bi);


///
PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id);

#endif //__PIDX_AGG_H
