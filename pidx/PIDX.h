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
 
#ifndef __PIDX_H
#define __PIDX_H

#include "PIDX_data_structs.h"
#include "PIDX_rst.h"
#include "PIDX_hz_encode.h"
#include "PIDX_agg.h"
#include "PIDX_io.h"

#include "PIDX_file_access_modes.h"
#include "PIDX_error_codes.h"
#include "PIDX_file_access_modes.h"
#include "PIDX_point.h"
#include "PIDX_data_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PIDX_row_major 0
#define PIDX_colume_major 1

struct PIDX_file_descriptor;
typedef struct PIDX_file_descriptor* PIDX_file;

typedef unsigned int PIDX_flags;

/// String describing a type. PIDX can support types like Float32 (e.g. for pressure) , 3*float64 (velocity vector field)
/// and more general structures. We use a string here so that a user who knows what he/she is doing can put an arbitrary
/// string, say "BOB". Obviously there will be need for a tool that understands the type "BOB".
typedef char* PIDX_type;
typedef int PIDX_data_layout;

extern const int PIDX_reccommended_bits_per_block;
extern const int PIDX_reccommended_blocks_per_file;


/*****************************************************************************************
 * Purpose:                                                                              *
 * Creates an IDX file.                                                                  *
 *                                                                                       *
 * Description:                                                                          *
 * PIDX_file_create is the primary function for creating IDX files;                      *
 * it creates a new IDX file with the specified name and mode                            *
 * specifying whether an existing file of same name should be overwritten.               *
 *                                                                                       *
 * Input parameters:                                                                     *
 * The name parameter specifies the name of the new file.                                *
 *                                                                                       *
 * The flags parameter specifies whether an existing file is to be overwritten.          *
 * It should be set to either PIDX_FILE_TRUNC to overwrite an existing file              *
 * or PIDX_FILE_EXCL, instructing the function to fail if the file already exists.       *
 * New files are always created in read-write mode, so the read-write and read-only      *
 * flags, PIDX_FILE_RDWR and PIDX_FILE_RDONLY, respectively, are not relevant in this    *
 * function. Further note that a specification of PIDX_FILE_RDONLY will be ignored;      *
 * the file will be created in read-write mode, regardless.                              *
 *                                                                                       *
 * Output parameters:                                                                    *
 * The flag parameter is the IDX file handler created by this function.                  *
 * The file handle returned, file, can be subsequently used to access the IDX            *
 * file until it is closed using PIDX_File_close.
 *                  
 * PIDX_return_code
 * PIDX_success if the task is completed correctly.
 * 
 *****************************************************************************************/
PIDX_return_code PIDX_file_create(const char* data_set_path, PIDX_flags flags, PIDX_file file);

/*****************************************************************************************
 * Purpose:                                                                               *
 * Opens an existing IDX file.                                                           *
 *                                                                                       *
 * Description:                                                                          *
 * PIDX_file_open is the primary function for accessing existing IDX files.              * 
 * This function opens the named file in the specified access mode.                      *
 *                                                                                       *
 * Note that PIDX_file_open does not create a file if it does not already exist;         *
 * see PIDX_file_create.                                                                 *
 *                                                                                       *
 * Input parameters:                                                                     *
 * The name parameter specifies the name of the file to be opened.                       *
 *                                                                                       *
 * The flags parameter specifies whether the file will be opened in                      *
 * read-write or read-only mode, PIDX_FILE_RDWR or PIDX_FILE_RDONLY, respectively.       *
 *                                                                                       *
 * Output parameters:                                                                    *
 * The return value is a file identifier for the open file; this file                    *
 * identifier should be closed by calling PIDX_File_close when it is no longer needed.   *
 * 
 * PIDX_return_code
 * PIDX_success if the task is completed correctly.
 *****************************************************************************************/
PIDX_return_code PIDX_file_open(const char* data_set_path, PIDX_flags  flags, PIDX_file file);

// Sets the bounding box for the IDX file associated with PIDX_file
PIDX_return_code PIDX_set_bounding_box_size(PIDX_file PIDX_file, PIDX_point box_size);

PIDX_return_code PIDX_get_bounding_box_size(PIDX_file PIDX_file, PIDX_point box_size);

PIDX_return_code PIDX_set_block_size(PIDX_file PIDX_file, const int block_size);

PIDX_return_code PIDX_get_block_size(PIDX_file PIDX_file, int* block_size);

PIDX_return_code PIDX_set_block_count(PIDX_file PIDX_file, const int block_count);

PIDX_return_code PIDX_get_block_count(PIDX_file PIDX_file, int* block_count);

PIDX_return_code PIDX_set_current_time_step(PIDX_file PIDX_file, const int time_step);

PIDX_return_code PIDX_get_current_time_step(PIDX_file PIDX_file, int* time_step);

/// Attach a communicator to a PIDX file descriptor
/// \return returns PIDX_OK if succeds or PIDX_FAIL otherwise
PIDX_return_code PIDX_set_communicator(PIDX_file file, MPI_Comm  comm);

PIDX_return_code PIDX_get_communicator(PIDX_file file, MPI_Comm* comm);

PIDX_return_code PIDX_variable_create(PIDX_file file, char* variable_name, unsigned int bits_per_sample, PIDX_type type_name, PIDX_variable variable);

PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable);

PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable);

PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool);

PIDX_return_code PIDX_get_bits_per_sample(PIDX_type type_name, unsigned int bits_per_sample);


/////////////////////////////////////////////////////////////////////////////////////////////
/// This is the only write function that a use should use for dumping data from a simulation.
PIDX_return_code PIDX_append_and_write_variable(PIDX_variable variable_ptr, 
			PIDX_point offset, 
			PIDX_point box_size, 
			const void* read_from_this_buffer, 
			PIDX_data_layout layout);

PIDX_return_code PIDX_read_next_variable(PIDX_variable variable_ptr, 
			PIDX_point offset, 
			PIDX_point box_size, 
			void* write_to_this_buffer, 
			PIDX_data_layout layout);

PIDX_return_code PIDX_get_box_count(PIDX_file file, int* box_coult);
PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point box_size);

PIDX_return_code PIDX_get_box_count_with_rank(PIDX_file file, int* box_coult, int MPI_rank);
PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, PIDX_point offset, PIDX_point box_size,  int MPI_rank);

PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_coult);
PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_coult);

PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index);
PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int  variable_index);

PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable variable);
PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable);

PIDX_return_code PIDX_read_variable(PIDX_variable variable_ptr, 
			PIDX_point offset, 
			PIDX_point box_size, 
			const void* read_from_this_buffer, 
			PIDX_data_layout layout);
			
PIDX_return_code PIDX_write_variable(PIDX_variable variable_ptr, 
			PIDX_point offset, 
			PIDX_point box_size, 
			const void* read_from_this_buffer, 
			PIDX_data_layout layout);
			
///Actually write the IDX file for all variables associated with file 
PIDX_return_code PIDX_flush(PIDX_file file);

///Perform all the necessary cleanups
PIDX_return_code PIDX_close(PIDX_file file);

#ifdef __cplusplus
}
#endif

#endif
