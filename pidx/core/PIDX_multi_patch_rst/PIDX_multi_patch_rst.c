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
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_multi_patch_rst_id PIDX_multi_patch_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_multi_patch_rst_id multi_patch_rst_id;
  multi_patch_rst_id = (PIDX_multi_patch_rst_id)malloc(sizeof (*multi_patch_rst_id));
  memset(multi_patch_rst_id, 0, sizeof (*multi_patch_rst_id));

  multi_patch_rst_id->idx = idx_meta_data;
  multi_patch_rst_id->idx_derived = idx_derived;
  multi_patch_rst_id->idx_c = idx_c;
  multi_patch_rst_id->idx_dbg = idx_dbg;

  multi_patch_rst_id->group_index = 0;
  multi_patch_rst_id->first_index = var_start_index;
  multi_patch_rst_id->last_index = var_end_index;

  multi_patch_rst_id->maximum_neighbor_count = 256;

  multi_patch_rst_id->reg_multi_patch_grp_count = 0;
  multi_patch_rst_id->sim_max_patch_group_count = 0;

  return (multi_patch_rst_id);
}


PIDX_return_code PIDX_multi_patch_rst_finalize(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  multi_patch_rst_id->reg_multi_patch_grp_count = 0;
  multi_patch_rst_id->sim_max_patch_group_count = 0;
  free(multi_patch_rst_id);
  multi_patch_rst_id = 0;

  return PIDX_success;
}
