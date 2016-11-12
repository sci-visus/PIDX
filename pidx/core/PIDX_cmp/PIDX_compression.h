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


struct PIDX_comp_id_struct;
typedef struct PIDX_comp_id_struct* PIDX_comp_id;

PIDX_comp_id PIDX_compression_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, int start_var_index, int end_var_index );


///
PIDX_return_code PIDX_compression(PIDX_comp_id id);



///
PIDX_return_code PIDX_decompression(PIDX_comp_id id);



///
PIDX_return_code PIDX_compression_finalize(PIDX_comp_id id);
#endif
