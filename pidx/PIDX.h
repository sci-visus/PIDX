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
 * \file PIDX.h
 *
 * \mainpage
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 * \author Giorgio Scorzelli
 * \author Valerio Pascucci
 * \date   10/09/14
 *
 * PIDX is an I/O library that enables HPC applications to write distributed 
 * multi-dimensional data directly into a hierarchical multi-resolution 
 * data format (IDX) with minimal overhead. 
 *
 */

#ifndef __PIDX_H
#define __PIDX_H

#include "PIDX_inc.h"

#ifdef __cplusplus
extern "C" {
#endif


struct PIDX_file_descriptor;
typedef struct PIDX_file_descriptor* PIDX_file;


///
double PIDX_get_time();


///
PIDX_return_code PIDX_time_step_caching_ON();


///
PIDX_return_code PIDX_time_step_caching_OFF();


/// Creates an IDX file.
/// PIDX_file_create is the primary function for creating IDX files;
/// it creates a new IDX file with the specified name and mode
/// specifying whether an existing file of same name should be overwritten.
/// \param  filename The filename parameter specifies the name of the new file.
/// \param flags  The flags parameter specifies whether an existing file is to be overwritten.
/// It should be set to either PIDX_FILE_TRUNC to overwrite an existing file
/// or PIDX_FILE_EXCL, instructing the function to fail if the file already exists.
/// New files are always created in read-write mode, so the read-write and read-only
/// flags, PIDX_FILE_RDWR and PIDX_FILE_RDONLY, respectively, are not relevant in this
/// function. Further note that a specification of PIDX_FILE_RDONLY will be ignored;
/// the file will be created in read-write mode, regardless.
/// \return file The file parameter is the IDX file handler created by this function.
/// The file handle returned, file, can be subsequently used to access the IDX
/// file until it is closed using PIDX_File_close.
/// \return PIDX_return_code The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file);


/// Opens an existing IDX file.
/// PIDX_file_open is the primary function for accessing existing IDX files.
/// This function opens the named file in the specified access mode.
/// Note that PIDX_file_open does not create a file if it does not already exist;
/// see PIDX_file_create.
/// \param filename The filename parameter specifies the name of the file to be opened.
/// \param flags The flags parameter specifies whether the file will be opened in
/// read-write or read-only mode, PIDX_FILE_RDWR or PIDX_FILE_RDONLY, respectively.
/// \return file The return value file, is a file identifier for the opened idx file; this 
/// identifier should be closed by calling PIDX_File_close when it is no longer needed.
/// \return PIDX_return_code The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file);


/// Get the PIDX_access associated with this file.
PIDX_return_code PIDX_get_access(PIDX_file file, PIDX_access *access);


/// Sets the dims of the IDX file.
/// \param file The IDX file handler.
/// \param dims Dimensions of the volume.
PIDX_return_code PIDX_set_dims(PIDX_file file, PIDX_point dims);


/// Gets the dims of the IDX file.
/// \param file The IDX file handler.
/// \param dims Dimensions of the volume will be returned here.
PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims);


/// Sets the transformation from logical dims to the physical coordinates (commonly used to visualize the volume)
/// \param file The IDX file handler.
/// \param transform standard 4x4 transformation matrix (right-handed / OpenGL style), see http://www.dirsig.org/docs/new/affine.html.
PIDX_return_code PIDX_set_transform(PIDX_file file, double transform[16]);


/// Gets the transformation from logical dims to the physical coordinates (commonly used to visualize the volume)
/// \param file The IDX file handler.
/// \param transform will return the standard 4x4 transformation matrix (right-handed / OpenGL style), see http://www.dirsig.org/docs/new/affine.html.
PIDX_return_code PIDX_get_transform(PIDX_file file, double transform[16]);


/// Gets the dims of the IDX file.
/// \param file The IDX file handler.
/// \param dims Dimensions of the volume will be returned here.
PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims);


/// Sets the block size of the IDX file.
/// \param file The IDX file handler.
/// \param block_size The block size.
PIDX_return_code PIDX_set_block_size(PIDX_file file, const int block_size);


/// Gets the block size of the IDX file.
/// \param file The IDX file handler.
/// \return block_size The block size.
PIDX_return_code PIDX_get_block_size(PIDX_file file, int* block_size);


/// Sets the number of blocks in an IDX file.
/// \param file The IDX file handler.
/// \param block_count Number of blocks per binary file.
PIDX_return_code PIDX_set_block_count(PIDX_file file, const int block_count);


/// Gets the number of blocks in an IDX file.
/// \param file The IDX file handler.
/// \return block_count Number of blocks per binary file.
PIDX_return_code PIDX_get_block_count(PIDX_file file, int* block_count);


/// Sets the current time step for IDX file.
/// The function should be used for time-varying dataset.
/// \param file The IDX file handler.
/// \param time_step The current time step.
PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int time_step);


/// Gets the current time step for IDX file.
/// \param file The IDX file handler.
/// \return time_step The current time step.
PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* time_step);


/// Creates a PIDX variable...
PIDX_return_code PIDX_variable_create(char* variable_name, unsigned int bits_per_sample, PIDX_type type_name, PIDX_variable* variable);


/// Gets the next PIDX variable...
PIDX_return_code PIDX_get_next_variable(PIDX_file file, PIDX_variable* variable);


/// Get the number of bits associated with a PIDX_type.
/// This function can be used to find the number of bits associated with a particular PIDX type.
/// \param PIDX_type 
/// \return Number of bits associated with type_name
PIDX_return_code PIDX_get_bits_per_sample(PIDX_type type_name, unsigned int bits_per_sample);



PIDX_return_code PIDX_variable_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout);


/// Write function used for dumping data from a simulation.
/// This function is used to write data in increasing order, typically suited for dumping data from a simulation.
/// \param variable The variable to be written.
/// \param offset The local offset of the data chunk associated with a process.
/// \param dims The box size of the chunk associated with the process.
/// \param src_buffer The data buffer that needs to be written by the process.
/// \param layout The current supported layouts are row major and column major.
PIDX_return_code PIDX_append_and_write_variable(PIDX_file file, PIDX_variable variable/*, PIDX_point offset, PIDX_point dims, const void* src_buffer, PIDX_data_layout layout*/);


/// Read function used for restarting a simulation from checkpoint dump.
/// This function is used to read data in the order in which the variables are laid out, typically suited for restarting a simulation from checkpoint.
/// \param variable The variale to be read.
/// \param offset The local offset of the data chunk associated with a process.
/// \param dims The box size of the chunk associated with the process.
/// \param dst_buffer The buffer to which data is read to.
/// \param layout The current supported layouts are row major and column major.
PIDX_return_code PIDX_read_next_variable(PIDX_variable variable, PIDX_point offset, PIDX_point dims, void* dst_buffer, PIDX_data_layout layout);


/// Enables dumping of meta-data.
/// This function enables writes of meta-data for any give variable in the form of box offset and sizes every process writes.
/// \param variable The variable handler.
PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable);


/// Disables dumping of meta-data.
/// This function disables writes of meta-data for any give variable in the form of box offset and sizes every process writes.
/// This is the default mode of the API.
/// \param variable The variable handler.
PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable);


/// Queries if the box meta-data is saved or not.
/// This function queries if meta-data associated with variable has been saved or not.
/// \param variable The variable handler.
/// \return 1 is meta-data is saved 0 otherwise.
PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool);


/// Gets the number of boxes.
/// This function gets the number of boxes written by all the processes.
/// \param file The file handler.
/// \return Number of boxes.
PIDX_return_code PIDX_get_box_count(PIDX_file file, int* box_count);


/// Gets the actual box meta-data associated with box_index.
/// \param file The file handler.
/// \param box_index The queried box index.
/// \return The offset of the queried box.
/// \return The length of the queried box.
PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point dims);


/// Gets the number of boxes associated with the given rank.
/// \param file The file handler.
/// \param MPI_rank  The rank for which one istrying to find the number of boxes.
/// \return Number of boxes.
PIDX_return_code PIDX_get_box_count_with_rank(PIDX_file file, int MPI_rank, int* box_count);


/// Gets the actual box meta-data associated with box_index and the MPI rank.
/// \param file The file handler.
/// \param box_index The queried box index.
/// \param MPI_rank The rank for which one istrying to find the number of boxes.
/// \return The offset of the queried box.
/// \return The length of the queried box.
PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, int MPI_rank, PIDX_point offset, PIDX_point dims);


/// Sets the number of variables in the IDX file.
/// \param file The IDX file handler.
/// \return The number of variables.
PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count);


/// Gets the number of variables in the IDX file.
/// \param file The IDX file handler.
/// \return The number of variables.
PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count);


/// Sets the index of the variable at which location it needs to be written at or read from in an IDX file.
/// \param file The IDX file handler.
/// \return Index of the variable.
PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int variable_index);


/// Gets the index of the variable in an IDX file.
/// starts with index 0, to (number of variable - 1) 
/// \param file The IDX file handler.
/// \return Index of the variable.
PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index);


///
PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable);


///
PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable);


///
PIDX_return_code PIDX_read_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void*  dst_buffer, PIDX_data_layout layout);


///
PIDX_return_code PIDX_write_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* src_buffer, PIDX_data_layout layout);


///Actually write the IDX file for all variables associated with file 
PIDX_return_code PIDX_flush(PIDX_file file);


///Perform all the necessary cleanups
PIDX_return_code PIDX_close(PIDX_file file);


///
PIDX_return_code PIDX_set_aggregation_factor(PIDX_file file, int agg_factor);


///
PIDX_return_code PIDX_get_aggregation_factor(PIDX_file file, int *agg_factor);


///
PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst);


///
PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz);


///
PIDX_return_code PIDX_dump_agg_info(PIDX_file file, int dump_agg_info);



///
PIDX_return_code PIDX_dump_io_info(PIDX_file file, int dump_io_info);


///
PIDX_return_code PIDX_set_resolution(PIDX_file file, int resolution_from, int resolution_to);


///
PIDX_return_code PIDX_get_resolution(PIDX_file file, int *resolution_from, int *resolution_to);


///
PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type);


///
PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type);


///
PIDX_return_code PIDX_set_lossy_compression_bit_rate(PIDX_file file, int compression_bit_rate);


///
PIDX_return_code PIDX_get_lossy_compression_bit_rate(PIDX_file file, int *compression_bit_rate);


///
PIDX_return_code PIDX_set_restructuring_box(PIDX_file file, PIDX_point restructured_box_size_point);


///
PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file);


///
PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file);


///
PIDX_return_code PIDX_debug_disable_compression(PIDX_file file);


///
PIDX_return_code PIDX_debug_disable_hz(PIDX_file file);


///
PIDX_return_code PIDX_debug_disable_agg(PIDX_file file);


///
PIDX_return_code PIDX_debug_disable_io(PIDX_file file);


///
PIDX_return_code PIDX_disable_rst(PIDX_file file);


///
PIDX_return_code PIDX_disable_agg(PIDX_file file);

#ifdef __cplusplus
}
#endif

#endif
