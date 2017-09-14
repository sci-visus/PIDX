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
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_idx_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_idx_rst_id PIDX_idx_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_idx_rst_id idx_rst_id;
  idx_rst_id = (PIDX_idx_rst_id)malloc(sizeof (*idx_rst_id));
  memset(idx_rst_id, 0, sizeof (*idx_rst_id));

  idx_rst_id->idx_metadata = idx_meta_data;
  idx_rst_id->idx_derived_metadata = idx_derived;
  idx_rst_id->idx_comm_metadata = idx_c;
  idx_rst_id->idx_debug_metadata = idx_dbg;

  idx_rst_id->group_index = 0;
  idx_rst_id->first_index = var_start_index;
  idx_rst_id->last_index = var_end_index;

  idx_rst_id->maximum_neighbor_count = 256;

  idx_rst_id->intersected_restructured_super_patch_count = 0;
  idx_rst_id->sim_max_patch_count = 0;

  return (idx_rst_id);
}


PIDX_return_code PIDX_idx_rst_finalize(PIDX_idx_rst_id idx_rst_id)
{
  idx_rst_id->intersected_restructured_super_patch_count = 0;
  idx_rst_id->sim_max_patch_count = 0;
  free(idx_rst_id);
  idx_rst_id = 0;

  return PIDX_success;
}
