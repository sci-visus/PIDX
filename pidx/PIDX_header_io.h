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

#ifndef __PIDX_HEADER_IO_H
#define __PIDX_HEADER_IO_H

struct PIDX_header_io_struct;
typedef struct PIDX_header_io_struct* PIDX_header_io_id;

PIDX_header_io_id PIDX_header_io_init(idx_dataset idx_meta_data,
			idx_dataset_derived_metadata idx_derived_ptr,
			int start_var_index, int end_var_index );

#if PIDX_HAVE_MPI
int PIDX_header_io_set_communicator(PIDX_header_io_id header_io, MPI_Comm comm);
#endif

int PIDX_header_io_write_idx (PIDX_header_io_id header_io, char* data_set_path, int current_time_step);

int PIDX_header_io_file_create(PIDX_header_io_id header_io);

int PIDX_header_io_file_write(PIDX_header_io_id header_io_id);

int PIDX_header_io_finalize(PIDX_header_io_id header_io);

#endif
