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
 * \file PIDX_multi_patch_rst.h
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Restructuring data from n cores to n' (n' <= n)
 * while keeping the data in mult-dimensional
 * application layout
 *
 */

#ifndef __PIDX_IDX_RST_NEW_H
#define __PIDX_IDX_RST_NEW_H


//Struct for restructuring ID
struct PIDX_idx_rst_struct
{
  idx_dataset idx_metadata;
  idx_dataset_derived_metadata idx_derived_metadata;
  idx_debug idx_debug_metadata;
  idx_comm idx_comm_metadata;

  int first_index;
  int last_index;
  int group_index;

  int intersected_restructured_super_patch_count;
  PIDX_super_patch* intersected_restructured_super_patch;

  unsigned long long sim_max_patch_count;
  unsigned long long* sim_multi_patch_r_count;
  unsigned long long* sim_multi_patch_r_offset;

  int maximum_neighbor_count;
};
typedef struct PIDX_idx_rst_struct* PIDX_idx_rst_id;


/*
 * Implementation in PIDX_rst.c
 */
///
/// \brief PIDX_idx_rst_init
/// \param idx_meta_data
/// \param idx_derived_ptr
/// \param var_start_index
/// \param var_end_index
/// \return
///
PIDX_idx_rst_id PIDX_idx_rst_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index);




///
/// \brief PIDX_idx_rst_finalize
/// \param id
/// \return
///
PIDX_return_code PIDX_idx_rst_finalize(PIDX_idx_rst_id id);



/*
 * Implementation in PIDX_idx_rst_meta_data.c
 */
///
/// \brief PIDX_idx_rst_meta_data_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_meta_data_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_meta_data_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_meta_data_destroy(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_buffer.c
 */
///
/// \brief PIDX_idx_rst_buf_create Create the appropriate data structs to hold restructured output data
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_destroy(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_aggregate_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_aggregate_buf_create(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_aggregate_buf_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_aggregate_buf_destroy(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_aggregate
/// \param rst_id
/// \param MODE
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregate(PIDX_idx_rst_id rst_id, int MODE);



/*
 * Implementation in PIDX_idx_rst_io.c
 */
///
/// \brief PIDX_idx_rst_buf_aggregated_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregated_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_buf_aggregated_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_buf_aggregated_read(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_write.c
 */
///
/// \brief PIDX_idx_rst_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_write(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_staged_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_staged_write(PIDX_idx_rst_id rst_id);



/*
 * Implementation in PIDX_idx_rst_read.c
 */
///
/// \brief PIDX_idx_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_read(PIDX_idx_rst_id rst_id);



///
/// \brief PIDX_idx_rst_forced_raw_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_idx_rst_forced_raw_read(PIDX_idx_rst_id rst_id);

/*
 * Implementation in PIDX_rst_verify.c
 */
///
/// \brief HELPER_rst
/// \param rst_id
/// \return
///
PIDX_return_code HELPER_idx_rst(PIDX_idx_rst_id rst_id);
#endif
