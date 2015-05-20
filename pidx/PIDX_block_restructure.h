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

#ifndef __PIDX_BOX_RST_H
#define __PIDX_BOX_RST_H


struct PIDX_chunk_id_struct;
typedef struct PIDX_chunk_id_struct* PIDX_chunk_id;


/// Creates a block restructuring ID
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_chunk_id The identifier associated with the task
PIDX_chunk_id PIDX_chunk_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int init_index, int start_var_index, int end_var_index);

#if PIDX_HAVE_MPI
int PIDX_chunk_set_communicator(PIDX_chunk_id id, MPI_Comm comm);
#endif

/// Restructures the voxels in (m x n x p) blocks
/// Example: The following 2D array
/// 1 2 3  4  5  6
/// 7 8 9 10 11 12
/// is usually stored in linear memory as 1 2 3 4 5 6 7 8 9 10 11 12
/// The restructured storage, in 3 x 2 blocks, is 1 2 3 7 8 9 4 5 6 10 11 12
PIDX_return_code PIDX_chunk_buf_create(PIDX_chunk_id id);


///
PIDX_return_code PIDX_chunk_write(PIDX_chunk_id id);



///
PIDX_return_code PIDX_chunk_read(PIDX_chunk_id id);



///
PIDX_return_code PIDX_chunk_buf_destroy(PIDX_chunk_id id);



///
PIDX_return_code PIDX_chunk_finalize(PIDX_chunk_id id);



///
int HELPER_chunking(PIDX_chunk_id id);

#endif
