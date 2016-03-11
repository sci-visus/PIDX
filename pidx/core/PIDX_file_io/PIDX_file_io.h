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

#ifndef __PIDX_FIlE_IO_H
#define __PIDX_FILE_IO_H 

/**
 * \file PIDX_io.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Move the data from the A aggregator processes to f files
 * 
 */

struct PIDX_io_struct;
typedef struct PIDX_io_struct* PIDX_io_id;


/// Creates the IO ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int init_index, int start_var_index, int end_var_index);



#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id IO id
/// \param comm the communicator
/// \return error code
int PIDX_io_set_communicator(PIDX_io_id io_id, MPI_Comm comm);
#endif



///
int PIDX_io_cached_data(uint32_t* cached_header);



///
int PIDX_io_file_create(PIDX_io_id io_id, int time_step, char* data_set_path, int MODE);



///
int PIDX_aggregated_io(PIDX_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, int MODE);



///
int PIDX_io_aggregated_read(PIDX_io_id io_id, Agg_buffer agg_buffer);



///
int PIDX_io_per_process(PIDX_io_id io_id, PIDX_block_layout block_layout, int MODE);



///
int PIDX_io_finalize(PIDX_io_id io_id);

#endif
