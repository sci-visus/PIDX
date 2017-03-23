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

struct PIDX_hz_encode_struct
{
  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;


  idx_comm idx_c;


  idx_debug idx_dbg;


  int** index;

  int group_index;

  int first_index;
  int last_index;

  int resolution_from;
  int resolution_to;
};
typedef struct PIDX_hz_encode_struct* PIDX_hz_encode_id;


/// Creates the HZ encoding file ID.
/// \param idx_meta_data All infor regarding the idx file passed from PIDX.c
/// \param idx_derived_ptr All derived idx related derived metadata passed from PIDX.c
/// \param start_var_index starting index of the variable on which the relevant operation is to be applied
/// \param end_var_index ending index of the variable on which the relevant operation is to be applied
/// \return PIDX_hz_encode_id The identifier associated with the task
PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, int start_var_index, int end_var_index);




///
/// \brief PIDX_hz_encode_meta_data_create
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_meta_data_destroy
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_buf_create
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id);




///
/// \brief PIDX_hz_encode_write_inverse
/// \param id
/// \param start_hz_index
/// \param end_hz_index
/// \return
///
PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index);


///
/// \brief PIDX_hz_encode_write
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_row_major_write
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_row_major_write(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_threshold_and_write
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_threshold_and_write(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_compress
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_compress(PIDX_hz_encode_id id);

///
/// \brief PIDX_hz_encode_read
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_buf_destroy
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id);



///
/// \brief PIDX_hz_encode_finalize
/// \param id
/// \return
///
PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id);


///
/// \brief PIDX_hz_encode_set_resolution
/// \param id
/// \param resolution_from
/// \param resolution_to
/// \return
///
PIDX_return_code PIDX_hz_encode_set_resolution(PIDX_hz_encode_id id, int resolution_from, int resolution_to);





///
int PIDX_file_io_per_process(PIDX_hz_encode_id io_id, PIDX_block_layout block_layout, int MODE);



///
/// \brief HELPER_Hz_encode
/// \param id
/// \return
///
int HELPER_Hz_encode(PIDX_hz_encode_id id);

#endif
