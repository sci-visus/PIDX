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
 
#include "PIDX_data_structs.h"
#include "PIDX_utils.h"
#include "PIDX_error_codes.h"
#include "PIDX_point.h"

#ifndef __PIDX_AGG_H
#define __PIDX_AGG_H 

struct PIDX_agg_struct;
typedef struct PIDX_agg_struct* PIDX_agg_id;

PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, 
			  idx_dataset_derived_metadata idx_derived_ptr, 
			  MPI_Comm comm,
			  int start_var_index, 
			  int end_var_index);

int PIDX_agg_aggregate(PIDX_agg_id agg_id, Agg_buffer agg_buffer);

int PIDX_agg_aggregate_write_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer);

int PIDX_agg_buf_destroy(Agg_buffer agg_buffer);

int PIDX_agg_finalize(PIDX_agg_id agg_id);

#endif