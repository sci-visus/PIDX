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

#ifndef __PIDX_MULTI_PATCH_RST_NEW_H
#define __PIDX_MULTI_PATCH_RST_NEW_H


//Struct for restructuring ID
struct PIDX_multi_patch_rst_struct
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
  idx_dataset_derived_metadata idx_derived;

  int init_index;
  int first_index;
  int last_index;
  int group_index;

  //dimension of the power-two volume imposed patch
  int64_t reg_patch_size[PIDX_MAX_DIMENSIONS];
  int reg_multi_patch_grp_count;
  Ndim_multi_patch_group* reg_multi_patch_grp;

  int64_t sim_max_patch_group_count;
  int64_t* sim_multi_patch_r_count;
  int64_t* sim_multi_patch_r_offset;
};
typedef struct PIDX_multi_patch_rst_struct* PIDX_multi_patch_rst_id;


/*
 * Implementation in PIDX_rst.c
 */
///
/// \brief PIDX_multi_patch_rst_init
/// \param idx_meta_data
/// \param idx_derived_ptr
/// \param first_index
/// \param var_start_index
/// \param var_end_index
/// \return
///
PIDX_multi_patch_rst_id PIDX_multi_patch_rst_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int first_index, int var_start_index, int var_end_index);



#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id restructuring id
/// \param comm the communicator
/// \return error code
PIDX_return_code PIDX_multi_patch_rst_set_communicator(PIDX_multi_patch_rst_id id, MPI_Comm comm);
#endif



///
/// \brief PIDX_multi_patch_rst_set_reg_patch_size
/// \param rst_id
/// \param factor
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_set_reg_patch_size(PIDX_multi_patch_rst_id rst_id, PIDX_point reg_box);




///
/// \brief PIDX_multi_patch_rst_auto_set_reg_patch_size
/// \param multi_patch_rst_id
/// \param factor
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_auto_set_reg_patch_size(PIDX_multi_patch_rst_id multi_patch_rst_id, int factor);



///
/// \brief PIDX_multi_patch_rst_finalize
/// \param id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_finalize(PIDX_multi_patch_rst_id id);



/*
 * Implementation in PIDX_multi_patch_rst_meta_data.c
 */
///
/// \brief PIDX_multi_patch_rst_meta_data_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_meta_data_create(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_meta_data_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_meta_data_write(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_meta_data_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_meta_data_destroy(PIDX_multi_patch_rst_id rst_id);



/*
 * Implementation in PIDX_multi_patch_rst_buffer.c
 */
///
/// \brief PIDX_multi_patch_rst_buf_create Create the appropriate data structs to hold restructured output data
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_create(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_buf_destroy Tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_destroy(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_aggregate_buf_create
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_aggregate_buf_create(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_aggregate_buf_destroy
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_aggregate_buf_destroy(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_buf_aggregate
/// \param rst_id
/// \param MODE
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_aggregate(PIDX_multi_patch_rst_id rst_id, int MODE);



/*
 * Implementation in PIDX_multi_patch_rst_io.c
 */
///
/// \brief PIDX_multi_patch_rst_buf_aggregate_and_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_and_write(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_buf_read_and_aggregate
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_read_and_aggregate(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_buf_aggregated_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_aggregated_write(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_buf_aggregated_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_buf_aggregated_read(PIDX_multi_patch_rst_id rst_id);



/*
 * Implementation in PIDX_multi_patch_rst_write.c
 */
///
/// \brief PIDX_multi_patch_rst_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_write(PIDX_multi_patch_rst_id rst_id);



///
/// \brief PIDX_multi_patch_rst_staged_write
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_staged_write(PIDX_multi_patch_rst_id rst_id);



/*
 * Implementation in PIDX_multi_patch_rst_read.c
 */
///
/// \brief PIDX_multi_patch_rst_read
/// \param rst_id
/// \return
///
PIDX_return_code PIDX_multi_patch_rst_read(PIDX_multi_patch_rst_id rst_id);



#endif // __PIDX_multi_patch_rst_NEW_H
