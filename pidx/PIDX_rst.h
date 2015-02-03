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

#if PIDX_HAVE_MPI

struct PIDX_rst_struct;
typedef struct PIDX_rst_struct* PIDX_rst_id;

PIDX_rst_id PIDX_rst_init( MPI_Comm comm, idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int var_start_index, int var_end_index);

int PIDX_rst_set_restructuring_box(PIDX_rst_id rst_id, int set_box_dim, int* box_dim);

/// actually do the restructuring, using pre-calculated data associated with the id
int PIDX_rst_restructure(PIDX_rst_id rst_id, PIDX_variable* variable);

int PIDX_rst_restructure_IO(PIDX_rst_id rst_id, PIDX_variable* variable, int MODE);
  
/// tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
int PIDX_rst_buf_destroy(PIDX_rst_id rst_id);

/// tear down whatever was calculated for this particular combination of dimensions and bounds
int PIDX_rst_finalize(PIDX_rst_id id);  

int64_t* PIDX_rst_get_box_dimension(PIDX_rst_id id);

int HELPER_rst(PIDX_rst_id rst_id, PIDX_variable* variable);

#endif // PIDX_HAVE_MPI
#endif // __PIDX_RST_NEW_H
