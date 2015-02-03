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

#ifndef __PIDX_COMPRESSION_H
#define __PIDX_COMPRESSION_H 



struct PIDX_compression_id_struct;
typedef struct PIDX_compression_id_struct* PIDX_compression_id;


PIDX_compression_id PIDX_compression_init(idx_dataset idx_meta_data,
                        idx_dataset_derived_metadata idx_derived_ptr,
                        int start_var_index, int end_var_index
                      );

#if PIDX_HAVE_MPI
int PIDX_compression_set_communicator(PIDX_compression_id id, MPI_Comm comm);
#endif

int PIDX_compression_prepare(PIDX_compression_id id, PIDX_variable* variable);

int PIDX_compression_compress(PIDX_compression_id id, PIDX_variable* variable, int MODE);
  
int PIDX_compression_buf_destroy(PIDX_compression_id id);

int PIDX_compression_finalize(PIDX_compression_id id);  

#endif
