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

#include "PIDX_data_structs.h"
#ifndef __PIDX_RST_NEW_H
#define __PIDX_RST_NEW_H

struct PIDX_rst_struct;
typedef struct PIDX_rst_struct* PIDX_rst_id;

PIDX_rst_id PIDX_rst_init( MPI_Comm comm, idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int var_start_index, int var_end_index);

int PIDX_rst_set_restructuring_box(PIDX_rst_id rst_id, int set_box_dim, int* box_dim);

/* actually do the restructuring, using pre-calculated data associated with the id */
//int PIDX_rst_restructure(PIDX_rst_id rst_id, int samples_per_variable, MPI_Datatype datatype, Ndim_buffer* in_buf, Ndim_buffer* out_buf_array, int num_output_buffers);
int PIDX_rst_restructure(PIDX_rst_id rst_id, int samples_per_variable, MPI_Datatype datatype, Ndim_buffer* in_buf, Ndim_buffer_group* out_buf_array, int num_output_buffers);

int PIDX_rst_restructure_IO(PIDX_rst_id rst_id, int samples_per_variable, MPI_Datatype datatype, Ndim_buffer* in_buf, Ndim_buffer_group* out_buf_array, int num_output_buffers);
//int PIDX_rst_restructure_IO(PIDX_rst_id rst_id, Ndim_buffer* in_buf, Ndim_buffer* out_buf_array, int num_output_buffers);
  
/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
int PIDX_rst_buf_destroy(int count, Ndim_buffer_group* out_buf_array);

/* tear down whatever was calculated for this particular combination of dimensions and bounds */ 
int PIDX_rst_finalize(PIDX_rst_id id);  

int* PIDX_rst_get_box_dimension(PIDX_rst_id id);

void PIDX_rst_print_error(char *error_message, char* file, int line);

int HELPER_rst(Ndim_buffer_group* out_buf_array1, PIDX_rst_id rst_id, int num_output_buffers, int spv);
//int HELPER_rst(Ndim_buffer* out_buf_array1, PIDX_rst_id rst_id, int num_output_buffers, int spv);

#endif /* __PIDX_RST_NEW_H */
