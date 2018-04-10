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

#include "../../PIDX_inc.h"



PIDX_return_code PIDX_hz_encode_set_resolution(PIDX_hz_encode_id id, int resolution_from, int resolution_to)
{
  if (resolution_from < 0 || resolution_from < 0)
    return PIDX_err_hz;

  id->resolution_from = resolution_from;
  id->resolution_to = resolution_to;

  return PIDX_success;
}


PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, idx_debug idx_dbg, idx_metadata_cache cache, int first_index, int last_index)
{
  PIDX_hz_encode_id hz_id;
  hz_id = (PIDX_hz_encode_id)malloc(sizeof (*hz_id));
  memset(hz_id, 0, sizeof (*hz_id));

  hz_id->idx = idx_meta_data;
  hz_id->idx_d = idx_d;
  hz_id->idx_c = idx_c;
  hz_id->idx_dbg = idx_dbg;
  hz_id->cache = cache;

  hz_id->group_index = 0;
  hz_id->first_index = first_index;
  hz_id->last_index = last_index;
  
  hz_id->resolution_from = 0;
  hz_id->resolution_to = 0;

  return hz_id;
}



PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id)
{
  
  free(id);
  id = 0;
  
  return PIDX_success;
}

