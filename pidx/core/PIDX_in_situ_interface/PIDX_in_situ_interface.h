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
 * \file PIDX_in_situ_interface.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Interface to call external libraries
 *
 */

#ifndef __PIDX_IN_SITU_INTERFACE_H
#define __PIDX_IN_SITU_INTERFACE_H


//Struct for in situ interface ID
struct PIDX_in_situ_interface_struct
{
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;

  idx_debug idx_dbg;


  idx_comm idx_c;

  int first_index;
  int last_index;
  int group_index;

};
typedef struct PIDX_in_situ_interface_struct* PIDX_insitu_id;



PIDX_insitu_id PIDX_in_situ_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, int group_index, int var_start_index, int var_end_index);



PIDX_return_code PIDX_in_situ_finalize(PIDX_insitu_id id);



PIDX_return_code PIDX_in_situ_perform(PIDX_insitu_id insitu_id);


#endif
