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
 * \file PIDX_hz_encode.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Incore restructuring of data, transforming the mult-dimensional 
 * application layout into the hierarchial Z order layout of the
 * IDX format.
 * 
 */

#ifndef __PIDX_HZ_ENCODE_H
#define __PIDX_HZ_ENCODE_H


struct PIDX_hz_encode_struct;
typedef struct PIDX_hz_encode_struct* PIDX_hz_encode_id;


/// Creates the HZ encoding file ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int init_index, int start_var_index, int end_var_index);


#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id hz encoding id
/// \param comm the communicator
/// \return error code
int PIDX_hz_encode_set_communicator(PIDX_hz_encode_id id, MPI_Comm comm);
#endif



///
PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id);


///
PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id);


///
PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id);



///
PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id);



///
PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id);



///
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id);



///
PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id);



///
int HELPER_Hz_encode(PIDX_hz_encode_id id);

#endif
