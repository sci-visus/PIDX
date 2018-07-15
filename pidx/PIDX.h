/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2010-2018 ViSUS L.L.C.,
 * Scientific Computing and Imaging Institute of the University of Utah
 *
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 *
 */

/*
 * Local partition read write
 * All idx mode read and write
 * group examples
 * documentation
 */

/**
 * \file PIDX.h
 *
 * \mainpage
 *
 * \author Sidharth Kumar
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


/*
 * Implementation in PIDX_file_create.c
 */
///
/// Creates an IDX file.
/// PIDX_file_create is the primary function for creating IDX files;
/// it creates a new IDX file with the specified name and mode
/// specifying whether an existing file of same name should be overwritten.
/// \param  filename The filename parameter specifies the name of the new file.
/// \param flags  The flags parameter specifies whether an existing file is to be overwritten.
/// It should be set to either PIDX_MODE_CREATE to overwrite an existing file
/// or PIDX_FILE_EXCL, instructing the function to fail if the file already exists.
/// New files are always created in read-write mode, so the read-write and read-only
/// flags, PIDX_FILE_RDWR and PIDX_FILE_RDONLY, respectively, are not relevant in this
/// function. Further note that a specification of PIDX_FILE_RDONLY will be ignored;
/// the file will be created in read-write mode, regardless.
/// \param access Used to manage file access between processes. Must be created prior to calling
/// this function using PIDX_create_access().
/// \param dims The dims parameter specifies the global domain dimensions
/// \param file Reference to PIDX_file object that will be initialized on successful completion of this function.
/// \return file The file parameter is the IDX file handler created by this function.
/// The file handle returned, file, can be subsequently used to access the IDX
/// file until it is closed using PIDX_File_close.
/// \return PIDX_return_code The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_point dims, PIDX_file* file);


/*
 * Implementation in PIDX_file_open.c
 */
///
/// Opens an existing IDX file.
/// PIDX_file_open is the primary function for accessing existing IDX files.
/// This function opens the named file in the specified access mode.
/// Note that PIDX_file_open does not create a file if it does not already exist;
/// see PIDX_file_create.
/// \param filename The filename parameter specifies the name of the file to be opened.
/// \param flags The flags parameter specifies whether the file will be opened in
/// read-write or read-only mode, PIDX_FILE_RDWR or PIDX_FILE_RDONLY, respectively.
/// \param access Used to manage file access between processes. Must be created prior to calling
/// this function using PIDX_create_access().
/// \param file Reference to PIDX_file object that will be initialized on successful completion of this function.
/// \return file: The return value file, is a file identifier for the opened idx file; this
/// identifier should be closed by calling PIDX_File_close when it is no longer needed.
/// \return PIDX_return_code: The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access, PIDX_point dims, PIDX_file* file);

/*
 *  Implementation in PIDX_metadata_parse.c
 */
// PIDX_metadata_parse parse the metadata .idx file (used by PIDX_file_open)
//
// There are different implementations for different versions of the metadata file
// PIDX_metadata_parse will use the corresponding file_open implementation
// for the specific version requested
//
PIDX_return_code PIDX_metadata_parse(FILE *fp, PIDX_file* file, char* version);
PIDX_return_code PIDX_metadata_parse_v6_0(FILE *fp, PIDX_file* file);
PIDX_return_code PIDX_metadata_parse_v6_1(FILE *fp, PIDX_file* file);

///
/// \brief PIDX_serial_file_open
/// This function reads an existing IDX file in serial
///
/// \param filename
/// \param flags
/// \param dims
/// \param file
/// \return
///
PIDX_return_code PIDX_serial_file_open(const char* filename, PIDX_flags flags, PIDX_point dims, PIDX_file* file);

// There are different implementations for different versions of the metadata file
// PIDX_file_open will use the corresponding file_open implementation for the specific version
PIDX_return_code PIDX_serial_file_open_v6(const char* filename, PIDX_flags flags, PIDX_point dims, PIDX_file* file);
PIDX_return_code PIDX_serial_file_open_v7(const char* filename, PIDX_flags flags, PIDX_point dims, PIDX_file* file);


///
/// \brief PIDX_query_box
/// \param file
/// \param box_dims
/// \return
///
PIDX_return_code PIDX_query_box(PIDX_file file, PIDX_point box_dims);



/*
 * Implementation in PIDX_close.c
 */
///
/// \brief PIDX_flush Actually write the IDX file for all variables associated with file
/// \param file
/// \return
///
PIDX_return_code PIDX_flush(PIDX_file file);



///
/// \brief PIDX_close Calls PIDX_flush and perform all the necessary cleanups
/// \param file
/// \return
///
PIDX_return_code PIDX_close(PIDX_file file);



/*
 * Implementation in PIDX_idx_set_get.c
 */
///
/// Sets the number of variables in the IDX file.
/// \param file The IDX file handler.
/// \return The number of variables.
///
PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count);



///
/// Gets the number of variables in the IDX file.
/// \param file The IDX file handler.
/// \return The number of variables.
///
PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count);



///
/// Sets the current time step for IDX file.
/// The function should be used for time-varying dataset.
/// \param file The IDX file handler.
/// \param time_step The current time step.
///
PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int time_step);



///
/// Gets the current time step for IDX file.
/// \param file The IDX file handler.
/// \return time_step The current time step.
///
PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* time_step);



///
/// Sets the block size of the IDX file.  An IDX binary file is composed of blocks.
/// \param file The IDX file handler.
/// \param block_size sets the block size of the IDX file.
/// \return PIDX_return_code The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_set_block_size(PIDX_file file, const int block_size);



///
/// Gets the block size of the IDX file. An IDX binary file is composed of blocks.
/// \param file The IDX file handler.
/// \return block_size: gets the block size of the IDX file.
/// \return PIDX_return_code: The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_get_block_size(PIDX_file file, int* block_size);



///
/// Sets the number of blocks in an IDX file. An IDX binary file is composed of blocks.
/// \param file The IDX file handler.
/// \param block_count Number of blocks per binary file.
/// \return PIDX_return_code: The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_set_block_count(PIDX_file file, const int block_count);



///
/// Gets the number of blocks in an IDX file.
/// \param file The IDX file handler.
/// \return block_count Number of blocks per binary file.
/// \return PIDX_return_code: The error code returned by the function.
/// It is PIDX_success if the task is completed correctly.
///
PIDX_return_code PIDX_get_block_count(PIDX_file file, int* block_count);



///
/// \brief PIDX_set_resolution
/// \param file
/// \param resolution_from
/// \param resolution_to
/// \return
///
PIDX_return_code PIDX_set_resolution(PIDX_file file, int resolution_to);



///
/// \brief PIDX_get_resolution
/// \param file
/// \param resolution_from
/// \param resolution_to
/// \return
///
PIDX_return_code PIDX_get_resolution(PIDX_file file, int *resolution_to);

  
///
/// \brief PIDX_get_box_for_resolution
/// \param file
/// \param resolution_to
/// \param offset of the box
/// \param size of the box
/// \return
///
PIDX_return_code PIDX_get_box_for_resolution(PIDX_file file, int resolution_to, PIDX_point offset, PIDX_point size, uint64_t* buffer_size);


///
/// \brief PIDX_set_partition_count
/// \param file
/// \param count_x
/// \param count_y
/// \param count_z
/// \return
///
PIDX_return_code PIDX_set_partition_count(PIDX_file file, int count_x, int count_y, int count_z);



///
/// \brief PIDX_get_partition_count
/// \param file
/// \param count_x
/// \param count_y
/// \param count_z
/// \return
///
PIDX_return_code PIDX_get_partition_count(PIDX_file file, int* count_x, int* count_y, int* count_z);



///
/// \brief PIDX_set_restructuring_box
/// \param file
/// \param restructured_box_size_point
/// \return
///
PIDX_return_code PIDX_set_restructuring_box(PIDX_file file, PIDX_point restructured_box_size_point);



///
/// \brief PIDX_get_restructuring_box
/// \param file
/// \param reg_patch_size
/// \return
///
PIDX_return_code PIDX_get_restructuring_box(PIDX_file file, PIDX_point reg_patch_size);



///
/// \brief PIDX_set_physical_dims
/// \param file
/// \param dims
/// \return
///
PIDX_return_code PIDX_set_physical_dims(PIDX_file file, PIDX_physical_point dims);



///
/// \brief PIDX_get_physical_dims
/// \param file
/// \param dims
/// \return
///
PIDX_return_code PIDX_get_physical_dims(PIDX_file file, PIDX_physical_point dims);



///
/// \brief PIDX_set_first_tstep
/// \param file
/// \param tstep
/// \return
///
PIDX_return_code PIDX_set_first_time_step(PIDX_file file, int tstep);



///
/// Gets the index of the first timestep in the IDX file.
/// \param file The IDX file handler.
/// \return The first timestep index.
///
PIDX_return_code PIDX_get_first_time_step(PIDX_file file, int* tstep);



///
/// \brief PIDX_set_last_tstep
/// \param file
/// \param tstep
/// \return
///
PIDX_return_code PIDX_set_last_time_step(PIDX_file file, int tstep);



///
/// Gets the index of the last timestep in the IDX file.
/// \param file The IDX file handler.
/// \return The last timestep index.
///
PIDX_return_code PIDX_get_last_time_step(PIDX_file file, int* tstep);



///
/// \brief PIDX_set_compression_type
/// \param file
/// \param compression_type
/// \return
///
PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type);



///
/// \brief PIDX_get_compression_type
/// \param file
/// \param compression_type
/// \return
///
PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type);





///
/// \brief PIDX_set_lossy_compression_bit_rate
/// \param file
/// \param var
/// \param compression_bit_rate
/// \return
///
PIDX_return_code PIDX_set_lossy_compression_bit_rate(PIDX_file file, PIDX_variable var, float compression_bit_rate);



///
/// \brief PIDX_set_average_compression_factor
/// \param file
/// \param compression_factor
/// \return
///
PIDX_return_code PIDX_set_average_compression_factor(PIDX_file file, int compression_factor, float bit_rate);



///
/// \brief PIDX_get_lossy_compression_bit_rate
/// \param file
/// \param compression_bit_rate
/// \return
///
PIDX_return_code PIDX_get_lossy_compression_bit_rate(PIDX_file file, int *compression_bit_rate);



///
/// \brief PIDX_set_io_mode
/// \param file
/// \param io_type
/// \return
///
PIDX_return_code PIDX_set_io_mode(PIDX_file file, enum PIDX_io_type io_type);



///
/// \brief PIDX_get_io_mode
/// \param file
/// \param io_type
/// \return
///
PIDX_return_code PIDX_get_io_mode(PIDX_file file, enum PIDX_io_type* io_type);



#if 0
///
/// \brief PIDX_set_wavelet_level
/// \param file
/// \param w_type
/// \return
///
PIDX_return_code PIDX_set_wavelet_level(PIDX_file file, int w_type);



///
/// \brief PIDX_set_wavelet_level
/// \param file
/// \param w_type
/// \return
///
PIDX_return_code PIDX_get_wavelet_level(PIDX_file file, int* w_type);



///
/// \brief PIDX_set_ROI_type
/// \param file
/// \param type
/// \return
///
PIDX_return_code PIDX_set_ROI_type(PIDX_file file, int type);



///
/// \brief PIDX_get_ROI_type
/// \param file
/// \param type
/// \return
///
PIDX_return_code PIDX_get_ROI_type(PIDX_file file, int* type);
#endif



///
/// \brief PIDX_set_variable_pile_length
/// \param file
/// \param var_pipe_length
/// \return
///
PIDX_return_code PIDX_set_variable_pile_length(PIDX_file file, int var_pipe_length);



///
/// \brief PIDX_get_variable_pile_length
/// \param file
/// \param var_pipe_length
/// \return
///
PIDX_return_code PIDX_get_variable_pile_length(PIDX_file file, int* var_pipe_length);



///
/// \brief PIDX_save_big_endian
/// \param file
/// \return
///
PIDX_return_code PIDX_save_big_endian(PIDX_file file);



///
/// \brief PIDX_save_little_endian
/// \param file
/// \return
///
PIDX_return_code PIDX_save_little_endian(PIDX_file file);



///
/// \brief PIDX_set_cache_time_step
/// \param file
/// \param ts
/// \return
///
PIDX_return_code PIDX_set_cache_time_step(PIDX_file file, int ts);



///
/// \brief PIDX_get_cache_time_step
/// \param file
/// \param ts
/// \return
///
PIDX_return_code PIDX_get_cache_time_step(PIDX_file file, int* ts);



#if 0
///
/// \brief PIDX_set_process_decomposition
/// \param file
/// \param np_x
/// \param np_y
/// \param np_z
/// \return
///
PIDX_return_code PIDX_set_process_decomposition(PIDX_file file, int np_x, int np_y, int np_z);



///
/// \brief PIDX_get_process_decomposition
/// \param file
/// \param np_x
/// \param np_y
/// \param np_z
/// \return
///
PIDX_return_code PIDX_get_process_decomposition(PIDX_file file, int* np_x, int* np_y, int* np_z);
#endif

///
/// \brief PIDX_set_comm
/// \param file
/// \param comm
/// \return
///
PIDX_return_code PIDX_set_comm(PIDX_file file, MPI_Comm comm);


///
/// \brief PIDX_get_comm
/// \param file
/// \param comm
/// \return
///
PIDX_return_code PIDX_get_comm(PIDX_file file, MPI_Comm *comm);
/*
 * Implementation in PIDX_variable.c
 */
///
/// \brief PIDX_variable_create Creates a PIDX variable.
/// \param variable_name
/// \param bits_per_sample
/// \param type_name
/// \param variable
/// \return
///
PIDX_return_code PIDX_variable_create(char* variable_name, unsigned int bits_per_sample, PIDX_data_type type_name, PIDX_variable* variable);



///
/// \brief PIDX_variable_write_data_layout
/// \param variable
/// \param offset The local offset of the data chunk associated with a process.
/// \param dims The box size of the chunk associated with the process.
/// \param read_from_this_buffer The data buffer that needs to be written by the process.
/// \param data_layout The current supported layouts are row major and column major.
/// \return
///
PIDX_return_code PIDX_variable_write_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout);



///
/// \brief PIDX_variable_write_particle_data_layout
/// \param variable
/// \param offset
/// \param dims
/// \param read_from_this_buffer
/// \param number_of_particles
/// \param data_layout
/// \return
///
PIDX_return_code PIDX_variable_write_particle_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, uint64_t number_of_particles, PIDX_data_layout data_layout);




///
/// \brief PIDX_variable_write_particle_data_physical_layout
/// \param variable
/// \param offset
/// \param dims
/// \param read_from_this_buffer
/// \param number_of_particles
/// \param data_layout
/// \return
///
PIDX_return_code PIDX_variable_write_particle_data_physical_layout(PIDX_variable variable, PIDX_physical_point offset, PIDX_physical_point dims, const void* read_from_this_buffer, uint64_t number_of_particles, PIDX_data_layout data_layout);



///
/// \brief PIDX_append_and_write_variable Write function used for dumping data from a simulation.
/// This function is used to write data in increasing order, typically suited for dumping data from a simulation.
/// \param file
/// \param variable The variable to be written.
/// \return
///
PIDX_return_code PIDX_append_and_write_variable(PIDX_file file, PIDX_variable variable);



///
/// \brief PIDX_get_next_variable Gets the next PIDX variable.
/// \param file
/// \param variable
/// \return
///
PIDX_return_code PIDX_get_next_variable(PIDX_file file, PIDX_variable* variable);



///
/// \brief PIDX_variable_read_data_layout
/// \param variable
/// \param offset The local offset of the data chunk associated with a process.
/// \param dims The box size of the chunk associated with the process.
/// \param read_from_this_buffer The buffer to which data is read to.
/// \param data_layout The current supported layouts are row major and column major.
/// \return
///
PIDX_return_code PIDX_variable_read_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, void* read_from_this_buffer, PIDX_data_layout data_layout);




///
/// \brief PIDX_variable_read_particle_data_layout
/// \param variable
/// \param offset
/// \param dims
/// \param write_to_this_buffer A buffer which will be allocated by PIDX to hold
///           sufficient storage for the particles being loaded.
/// \param number_of_particles The number of particles read
/// \param data_layout
/// \return
///
PIDX_return_code PIDX_variable_read_particle_data_layout(PIDX_variable variable, PIDX_physical_point offset, PIDX_physical_point dims, void** write_to_this_buffer, uint64_t* number_of_particles, PIDX_data_layout data_layout);




///
/// \brief PIDX_read_next_variable Read function used for restarting a simulation from checkpoint dump.
/// This function is used to read data in the order in which the variables are laid out, typically suited for restarting a simulation from checkpoint.
/// \param file
/// \param variable The variale to be read.
/// \return
///
PIDX_return_code PIDX_read_next_variable(PIDX_file file, PIDX_variable variable);



///
/// \brief PIDX_reset_variable_counter
/// \param file
/// \return
///
PIDX_return_code PIDX_reset_variable_counter(PIDX_file file);



///
/// Sets the index of the variable at which location it needs to be written at or read from in an IDX file.
/// \param file The IDX file handler.
/// \return Index of the variable.
///
PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int variable_index);



///
/// \brief PIDX_set_current_variable_by_name
/// \param file
/// \param variable_name
/// \return
///
PIDX_return_code PIDX_set_current_variable_by_name(PIDX_file file, const char* variable_name);



///
/// Gets the index of the variable in an IDX file.
/// starts with index 0, to (number of variable - 1)
/// \param file The IDX file handler.
/// \return Index of the variable.
///
PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index);



///
/// \brief PIDX_set_current_variable
/// \param file
/// \param variable
/// \return
///
PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable);



///
/// \brief PIDX_get_current_variable
/// \param file
/// \param variable
/// \return
///
PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable);



///
/// \brief PIDX_read_variable
/// \param file
/// \param variable
/// \param offset
/// \param dims
/// \param dst_buffer
/// \param layout
/// \return
///
PIDX_return_code PIDX_read_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void*  dst_buffer, PIDX_data_layout layout);



///
/// \brief PIDX_write_variable
/// \param file
/// \param variable
/// \param offset
/// \param dims
/// \param src_buffer
/// \param layout
/// \return
///
PIDX_return_code PIDX_write_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* src_buffer, PIDX_data_layout layout);


/*
 * Implementation in PIDX_meta_data.c
 */
///
/// \brief PIDX_set_meta_data_cache
/// \param file
/// \param cache
/// \return
///
PIDX_return_code PIDX_set_meta_data_cache(PIDX_file file, PIDX_metadata_cache cache);




///
/// \brief PIDX_set_meta_data_cache
/// \param file
/// \param cache
/// \return
///
PIDX_return_code PIDX_get_meta_data_cache(PIDX_file file, PIDX_metadata_cache* cache);




///
/// Gets the number of boxes.
/// This function gets the number of boxes written by all the processes.
/// \param file The file handler.
/// \return Number of boxes.
///
PIDX_return_code PIDX_get_box_count(PIDX_file file, int* box_count);



///
/// Gets the actual box meta-data associated with box_index.
/// \param file The file handler.
/// \param box_index The queried box index.
/// \return The offset of the queried box.
/// \return The length of the queried box.
///
PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point dims);



///
/// Gets the number of boxes associated with the given rank.
/// \param file The file handler.
/// \param MPI_rank  The rank for which one istrying to find the number of boxes.
/// \return Number of boxes.
///
PIDX_return_code PIDX_get_box_count_with_rank(PIDX_file file, int MPI_rank, int* box_count);



///
/// Gets the actual box meta-data associated with box_index and the MPI rank.
/// \param file The file handler.
/// \param box_index The queried box index.
/// \param MPI_rank The rank for which one istrying to find the number of boxes.
/// \return The offset of the queried box.
/// \return The length of the queried box.
///
PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, int MPI_rank, PIDX_point offset, PIDX_point dims);



///
/// Enables dumping of meta-data.
/// This function enables writes of meta-data for any give variable in the form of box offset and sizes every process writes.
/// \param variable The variable handler.
///
PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable);



///
/// Disables dumping of meta-data.
/// This function disables writes of meta-data for any give variable in the form of box offset and sizes every process writes.
/// This is the default mode of the API.
/// \param variable The variable handler.
///
PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable);



///
/// Queries if the box meta-data is saved or not.
/// This function queries if meta-data associated with variable has been saved or not.
/// \param variable The variable handler.
/// \return 1 is meta-data is saved 0 otherwise.
///
PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool);


/*
 * Implementation in PIDX_debug.c
 */
///
/// \brief PIDX_debug_disable_restructuring
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file);



///
/// \brief PIDX_debug_disable_chunking
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file);



///
/// \brief PIDX_debug_disable_compression
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_compression(PIDX_file file);



///
/// \brief PIDX_debug_disable_hz
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_hz(PIDX_file file);



///
/// \brief PIDX_debug_disable_agg
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_agg(PIDX_file file);



///
/// \brief PIDX_debug_disable_io
/// \param file
/// \return
///
PIDX_return_code PIDX_debug_disable_io(PIDX_file file);



///
/// \brief PIDX_debug_rst
/// \param file
/// \param debug_rst
/// \return
///
PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst);



///
/// \brief PIDX_debug_hz
/// \param file
/// \param debug_hz
/// \return
///
PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz);


///
/// \brief PIDX_disable_agg
/// \param file
/// \return
///
PIDX_return_code PIDX_disable_agg(PIDX_file file);



///
/// \brief PIDX_dump_state
/// \param file
/// \param process_state
/// \return
///
PIDX_return_code PIDX_dump_state(PIDX_file file, int process_current_state);
/*
 * Implementation in PIDX_debug.c
 */
///
/// \brief PIDX_default_bits_per_datatype Get the number of bits associated with a PIDX_type.
/// This function can be used to find the number of bits associated with a particular PIDX type.
/// \param type
/// \param bits Number of bits associated with type_name
/// \return
///
PIDX_return_code PIDX_default_bits_per_datatype(PIDX_data_type type, int* bits);



///
/// \brief PIDX_values_per_datatype
/// \param type
/// \param values
/// \param bits
/// \return
///
PIDX_return_code PIDX_values_per_datatype(PIDX_data_type type, int* values, int* bits);

#ifdef __cplusplus
}
#endif

#endif
