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
 * \file PIDX_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_rst.h
 *
 */

#include "../../PIDX_inc.h"



PIDX_rst_id PIDX_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int var_start_index, int var_end_index)
{
  //Creating the restructuring ID
  PIDX_rst_id rst_id;
  rst_id = (PIDX_rst_id)malloc(sizeof (*rst_id));
  memset(rst_id, 0, sizeof (*rst_id));

  rst_id->idx = idx_meta_data;
  rst_id->idx_d = idx_d;
  rst_id->idx_c = idx_c;

  rst_id->group_index = 0;
  rst_id->first_index = var_start_index;
  rst_id->last_index = var_end_index;

  return (rst_id);
}



PIDX_return_code PIDX_rst_finalize(PIDX_rst_id rst_id)
{
  
  free(rst_id);
  rst_id = 0;
  
  return PIDX_success;
}
