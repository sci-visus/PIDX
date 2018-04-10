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


PIDX_file_io_id PIDX_file_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int first_index, int last_index)
{
  PIDX_file_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_file_io_id)malloc(sizeof (*io_id));
  memset(io_id, 0, sizeof (*io_id));

  io_id->idx = idx_meta_data;
  io_id->idx_d = idx_d;
  io_id->idx_c = idx_c;

  io_id->group_index = 0;
  io_id->first_index = first_index;
  io_id->last_index = last_index;

  return io_id;
}



int PIDX_file_io_finalize(PIDX_file_io_id io_id)
{

  free(io_id);
  io_id = 0;

  return PIDX_success;
}
