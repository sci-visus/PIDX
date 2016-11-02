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


///
/// \brief The PIDX_rst_struct struct is Struct for restructuring ID
///
struct PIDX_rst_struct
{
  //Passed by PIDX API
#if PIDX_HAVE_MPI
  MPI_Comm comm; //Communicator
#endif

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int init_index;
  int group_index;
  int first_index;
  int last_index;

  //dimension of the power-two volume imposed patch
  unsigned long long reg_patch_size[PIDX_MAX_DIMENSIONS];
  int reg_patch_grp_count;
  Ndim_patch_group* reg_patch_grp;

  unsigned long long *rank_r_offset;
  unsigned long long *rank_r_count;
};
typedef struct PIDX_rst_struct* PIDX_rst_id;



/*
 * Implementation in PIDX_rst.c
 */
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
/// \brief PIDX_rst_auto_set_reg_patch_size
/// \param rst_id
/// \param factor
/// \return
///
PIDX_return_code PIDX_rst_auto_set_reg_patch_size(PIDX_rst_id rst_id, int factor);



///
/// \brief PIDX_rst_set_reg_patch_size
/// \param rst_id
/// \param box_size
/// \return
///
PIDX_return_code PIDX_rst_set_reg_patch_size(PIDX_rst_id rst_id, PIDX_point box_size);



///
/// \brief PIDX_rst_set_reg_patch_size_from_bit_string
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_set_reg_patch_size_from_bit_string(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_finalize Tear down whatever was calculated for this particular combination of dimensions and bounds
/// \param id
/// \return
///
PIDX_return_code PIDX_rst_finalize(PIDX_rst_id id);



/*
 * Implementation in PIDX_rst_meta_data.c
 */
///
/// \brief PIDX_rst_attach_restructuring_box
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_meta_data_create(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_meta_data_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_meta_data_write(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_meta_data_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_meta_data_destroy(PIDX_rst_id rst_id);



/*
 * Implementation in PIDX_rst_buffer.c
 */
///
/// \brief PIDX_rst_buf_create Create the appropriate data structs to hold restructured output data
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_create(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_destroy(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_aggregate_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_aggregate_buf_create(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_aggregate_buf_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_aggregate_buf_destroy(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_buf_aggregate
/// \param rst_id
/// \param MODE
/// \return
///
PIDX_return_code PIDX_rst_buf_aggregate(PIDX_rst_id rst_id, int MODE);



/*
 * Implementation in PIDX_rst_io.c
 */
///
/// \brief PIDX_rst_buf_aggregate_and_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_aggregate_and_write(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_buf_read_and_aggregate
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_read_and_aggregate(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_buf_aggregated_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_aggregated_write(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_buf_aggregated_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_buf_aggregated_read(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_forced_raw_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_forced_raw_read(PIDX_rst_id rst_id);
/*
 * Implementation in PIDX_rst_write.c
 */
///
/// \brief PIDX_rst_write Actually do the restructuring, using pre-calculated data associated with the id
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_write(PIDX_rst_id rst_id);



///
/// \brief PIDX_rst_staged_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_staged_write(PIDX_rst_id rst_id);



/*
 * Implementation in PIDX_rst_read.c
 */
///
/// \brief PIDX_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_rst_read(PIDX_rst_id rst_id);



/*
 * Implementation in PIDX_rst_verify.c
 */
///
/// \brief HELPER_rst
/// \param rst_id
/// \return
///
PIDX_return_code HELPER_rst(PIDX_rst_id rst_id);

#endif // __PIDX_RST_NEW_H
