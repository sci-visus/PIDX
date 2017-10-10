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
 * \file PIDX_raw_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_raw_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_raw_rst_id PIDX_raw_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_raw_rst_id raw_rst_id;
  raw_rst_id = (PIDX_raw_rst_id)malloc(sizeof (*raw_rst_id));
  memset(raw_rst_id, 0, sizeof (*raw_rst_id));

  raw_rst_id->idx = idx_meta_data;
  raw_rst_id->idx_derived = idx_derived;
  raw_rst_id->idx_c = idx_c;
  raw_rst_id->idx_dbg = idx_dbg;

  raw_rst_id->group_index = 0;
  raw_rst_id->first_index = var_start_index;
  raw_rst_id->last_index = var_end_index;

  raw_rst_id->maximum_neighbor_count = 256;

  raw_rst_id->reg_raw_grp_count = 0;
  raw_rst_id->sim_max_patch_group_count = 0;

  return (raw_rst_id);
}


PIDX_return_code PIDX_raw_rst_finalize(PIDX_raw_rst_id raw_rst_id)
{
  raw_rst_id->reg_raw_grp_count = 0;
  raw_rst_id->sim_max_patch_group_count = 0;
  free(raw_rst_id);
  raw_rst_id = 0;

  return PIDX_success;
}
