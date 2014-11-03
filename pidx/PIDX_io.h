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

#ifndef __PIDX_IO_H
#define __PIDX_IO_H 

#include "PIDX_data_structs.h"

struct PIDX_io_struct;
typedef struct PIDX_io_struct* PIDX_io_id;


PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data,
			idx_dataset_derived_metadata idx_derived_ptr,
			int start_var_index, int end_var_index
 		      );

#if PIDX_HAVE_MPI
int PIDX_io_set_communicator(PIDX_io_id io_id, MPI_Comm comm);
#endif

int PIDX_io_file_create(PIDX_io_id io_id, 
			int time_step, 
			char* data_set_path, int MODE);

#if PIDX_HAVE_MPI
int PIDX_io_aggregated_IO(PIDX_io_id io_id, Agg_buffer agg_buffer, int MODE);
#endif

int PIDX_io_independent_IO_var(PIDX_io_id io_id, PIDX_variable* variable_ptr, int MODE);

int PIDX_io_finalize(PIDX_io_id io_id);

#endif