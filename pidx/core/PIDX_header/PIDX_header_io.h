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


///
/// \brief PIDX_header_io_init
/// \param idx_meta_data
/// \param idx_d
/// \param idx_c
/// \param first_index
/// \param last_index
/// \return
///
PIDX_header_io_id PIDX_header_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int first_index, int last_index);



///
/// \brief PIDX_header_io_write_idx
/// \param header_io
/// \param data_set_path
/// \param current_time_step
/// \return
///
PIDX_return_code PIDX_header_io_global_idx_write (PIDX_header_io_id header_io, char* data_set_path);




PIDX_return_code PIDX_header_io_partition_idx_write (PIDX_header_io_id header_io, char* data_set_path);



PIDX_return_code PIDX_header_io_raw_idx_write (PIDX_header_io_id header_io, char* data_set_path);


///
/// \brief PIDX_header_io_idx_file_create
/// \param header_io_id
/// \param block_layout
/// \param filename_template
/// \return
///
int PIDX_header_io_idx_file_create(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, char* filename_template);



///
/// \brief PIDX_header_io_idx_file_write
/// \param header_io_id
/// \param block_layout
/// \param file_name
/// \param file_name_template
/// \param mode
/// \return
///
PIDX_return_code PIDX_header_io_idx_file_write(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, char* file_name_template, int mode);




///
/// \brief PIDX_header_io_raw_dir_create
/// \param header_io_id
/// \param file_name
/// \return
///
int PIDX_header_io_raw_dir_create(PIDX_header_io_id header_io_id, char* file_name);



///
/// \brief PIDX_header_io_finalize
/// \param header_io
/// \return
///
int PIDX_header_io_finalize(PIDX_header_io_id header_io);

#endif
