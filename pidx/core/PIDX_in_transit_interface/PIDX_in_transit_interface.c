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
 * \file PIDX_in_transit_interface.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_in_transit_interface.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_in_transit_id PIDX_in_transit_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_in_transit_id in_transit_id;
  in_transit_id = (PIDX_in_transit_id)malloc(sizeof (*in_transit_id));
  memset(in_transit_id, 0, sizeof (*in_transit_id));

  in_transit_id->idx = idx_meta_data;
  in_transit_id->idx_derived = idx_derived;
  in_transit_id->idx_c = idx_c;
  in_transit_id->idx_dbg = idx_dbg;

  in_transit_id->group_index = 0;
  in_transit_id->first_index = var_start_index;
  in_transit_id->last_index = var_end_index;

  return in_transit_id;
}


PIDX_return_code PIDX_in_transit_finalize(PIDX_in_transit_id in_transit_id)
{
  free(in_transit_id);
  in_transit_id = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_in_transit_meta_data_create(PIDX_in_transit_id transit_id)
{
  return PIDX_success;
}



PIDX_return_code PIDX_in_transit_meta_data_destroy(PIDX_in_transit_id transit_id)
{
  return PIDX_success;
}



PIDX_return_code PIDX_in_transit_perform(PIDX_in_transit_id transit_id)
{
  return PIDX_success;
}
