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


PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  file->hz_from_shared = 0;
  file->hz_to_shared = idx->total_partiton_level;

  file->hz_from_non_shared = idx->total_partiton_level;
  file->hz_to_non_shared =  idx->maxh;

  if (file->hz_from_shared == file->hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
  }

  if (file->hz_from_non_shared == file->hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
  }

  if (file->hz_from_shared == file->hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
    var_grp->shared_layout_count = 0;
  }
  else
  {
    var_grp->shared_start_layout_index = (file->hz_from_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->shared_start_layout_index <= 0)
      var_grp->shared_start_layout_index = 0;

    var_grp->shared_end_layout_index = (file->hz_to_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->shared_end_layout_index <= 0)
      var_grp->shared_end_layout_index = 1;

    var_grp->shared_layout_count = var_grp->shared_end_layout_index - var_grp->shared_start_layout_index;
  }

  if (file->hz_from_non_shared == file->hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
    var_grp->nshared_layout_count = 0;
  }
  else
  {
    var_grp->nshared_start_layout_index = (file->hz_from_non_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->nshared_start_layout_index <= 0)
      var_grp->nshared_start_layout_index = 0;

    var_grp->nshared_end_layout_index = (file->hz_to_non_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->nshared_end_layout_index <= 0)
      var_grp->nshared_end_layout_index = 1;

    var_grp->nshared_layout_count = var_grp->nshared_end_layout_index - var_grp->nshared_start_layout_index;
  }

  return PIDX_success;
}




PIDX_return_code find_agg_level(PIDX_io file, int gi)
{
  int i = 0;
  int no_of_aggregators = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  int total_aggregator = 0;

  if (file->idx->enable_agg == 0)
    var_grp->agg_level = var_grp->shared_start_layout_index;
  else
  {
    for (i = 0; i < var_grp->shared_layout_count + var_grp->nshared_layout_count ; i++)
    {
      no_of_aggregators = var_grp->block_layout_by_level[i]->efc;
      total_aggregator = total_aggregator + no_of_aggregators;
      if (no_of_aggregators <= file->idx_c->lnprocs)
        var_grp->agg_level = i + 1;
    }
  }

  if (total_aggregator > file->idx_c->lnprocs)
    var_grp->agg_level = var_grp->shared_start_layout_index;

  return PIDX_success;
}
