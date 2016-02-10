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
 * \file PIDX_rst.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Restructuring data from n cores to n' (n' <= n)
 * while keeping the data in mult-dimensional 
 * application layout
 * 
 */

#ifndef __PIDX_RST_NEW_H
#define __PIDX_RST_NEW_H

/// The restructuring module is called only when there is MPI availiable, hence
/// the entire code is wrapped in PIDX_HAVE_MPI define


struct PIDX_rst_struct;
typedef struct PIDX_rst_struct* PIDX_rst_id;


/// Creates an restructuring file ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_rst_id The identifier associated with the task
PIDX_rst_id PIDX_rst_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int first_index, int var_start_index, int var_end_index);


#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id restructuring id
/// \param comm the communicator
/// \return error code
PIDX_return_code PIDX_rst_set_communicator(PIDX_rst_id id, MPI_Comm comm);
#endif



///
/// \brief PIDX_rst_attach_restructuring_box
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_meta_data_create(PIDX_rst_id rst_id);


PIDX_return_code PIDX_rst_meta_data_destroy(PIDX_rst_id rst_id);



/// Ceate the appropriate data structs to hold restructured output data

///
/// \brief PIDX_rst_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_create(PIDX_rst_id rst_id);



PIDX_return_code PIDX_rst_buf_aggregate(PIDX_rst_id rst_id);


///
/// \brief PIDX_rst_write Actually do the restructuring, using pre-calculated data associated with the id
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_write(PIDX_rst_id rst_id);


///
/// \brief PIDX_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_read(PIDX_rst_id rst_id);


///
/// \brief PIDX_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_destroy(PIDX_rst_id rst_id);


///
/// \brief PIDX_rst_finalize Tear down whatever was calculated for this particular combination of dimensions and bounds
/// \param id
/// \return
///
PIDX_return_code PIDX_rst_finalize(PIDX_rst_id id);


///
/// \brief HELPER_rst
/// \param rst_id
/// \return
///
PIDX_return_code HELPER_rst(PIDX_rst_id rst_id);

#endif // __PIDX_RST_NEW_H
