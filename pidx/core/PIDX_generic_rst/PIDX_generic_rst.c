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
 * \file PIDX_generic_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_generic_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_generic_rst_id PIDX_generic_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_generic_rst_id generic_rst_id;
  generic_rst_id = (PIDX_generic_rst_id)malloc(sizeof (*generic_rst_id));
  memset(generic_rst_id, 0, sizeof (*generic_rst_id));

  generic_rst_id->idx = idx_meta_data;
  generic_rst_id->idx_d = idx_d;
  generic_rst_id->idx_c = idx_c;
  generic_rst_id->idx_dbg = idx_dbg;

  generic_rst_id->group_index = 0;
  generic_rst_id->first_index = var_start_index;
  generic_rst_id->last_index = var_end_index;

  return (generic_rst_id);
}


#if 0
PIDX_return_code PIDX_generic_rst_perform(PIDX_generic_rst_id id)
{
  int i = 0, v = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];

  for (i = 0; i < id->idx->regridded_process_count[0] * id->idx->regridded_process_count[1] * id->idx->regridded_process_count[2]; i++)
  {
    if (id->idx_c->grank = id->idx->regridded_patch[i]->rank)
    {
      id->reg_patch_grp_count = 1;

        for (v = id->first_index; v <= id->last_index; v++)
        {
          PIDX_variable var = var_grp->variable[v];
          Ndim_patch_group patch_group = var->rst_patch_group[0];
          patch_group->data_source = 0;

          for(j = 0; j < id->reg_patch_grp[i]->count; j++)
          {
            patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

            if (patch_group->patch[j]->buffer == NULL)
              return PIDX_err_rst;
            memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
          }
        }
    }
    else
      id->reg_patch_grp_count = 0;

  }

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_generic_rst_finalize(PIDX_generic_rst_id generic_rst_id)
{  
  free(generic_rst_id);
  generic_rst_id = 0;
  
  return PIDX_success;
}
