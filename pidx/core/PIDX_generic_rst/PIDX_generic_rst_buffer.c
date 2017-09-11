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

PIDX_return_code PIDX_generic_rst_buf_create(PIDX_generic_rst_id generic_rst_id)
{
  int j = 0, v = 0;
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  int cnt = 0, i = 0;
  for (v = generic_rst_id->first_index; v <= generic_rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    cnt = 0;
    for (i = 0; i < generic_rst_id->reg_patch_grp_count; i++)
    {
      if (generic_rst_id->idx_c->grank == generic_rst_id->reg_patch_grp[i]->max_patch_rank)
      {
        //fprintf(stderr, "Buffer create rank %d\n", generic_rst_id->idx_c->grank);
        Ndim_patch_group patch_group = var->rst_patch_group;
        for(j = 0; j < generic_rst_id->reg_patch_grp[i]->count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

          if (patch_group->patch[j]->buffer == NULL)
            return PIDX_err_rst;
          memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
        }
        cnt++;
        assert(cnt == 1);
      }
    }
  }

  return PIDX_success;
}




PIDX_return_code PIDX_generic_rst_buf_destroy(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  int j, v;
  for(v = generic_rst_id->first_index; v <= generic_rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for(j = 0; j < var_grp->variable[v]->rst_patch_group->count; j++)
    {
      free(var->rst_patch_group->patch[j]->buffer);
      var->rst_patch_group->patch[j]->buffer = 0;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_generic_rst_aggregate_buf_create(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
    return PIDX_success;

  int v = 0;
  for (v = generic_rst_id->first_index; v <= generic_rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    // copy the size and offset to output
    Ndim_patch out_patch = var->rst_patch_group->reg_patch;

    int nx = out_patch->size[0];
    int ny = out_patch->size[1];
    int nz = out_patch->size[2];

    var->rst_patch_group->reg_patch->buffer = malloc((uint32_t)(nx * ny * nz * (var->bpv/8) * var->vps));
    memset(var->rst_patch_group->reg_patch->buffer, 0,(uint32_t)(nx * ny * nz * (var->bpv/8) * var->vps));

    if (var->rst_patch_group->reg_patch->buffer == NULL)
      return PIDX_err_chunk;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_generic_rst_aggregate_buf_destroy(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = generic_rst_id->first_index; v <= generic_rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    free(var->rst_patch_group->reg_patch->buffer);
  }

  return PIDX_success;
}



PIDX_return_code PIDX_generic_rst_buf_aggregate(PIDX_generic_rst_id generic_rst_id, int mode)
{
  int v = 0;
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  for (v = generic_rst_id->first_index; v <= generic_rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    // copy the size and offset to output
    Ndim_patch_group patch_group = var->rst_patch_group;
    Ndim_patch out_patch = var->rst_patch_group->reg_patch;

    int nx = out_patch->size[0];
    int ny = out_patch->size[1];

    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
    for (r = 0; r < var->rst_patch_group->count; r++)
    {
      for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
      {
        for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
        {
          for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
          {
            index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
            send_o = index * var->vps * (var->bpv/8);
            send_c = (patch_group->patch[r]->size[0]);
            recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

            if (mode == PIDX_WRITE)
              memcpy(out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group->patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
            else
              memcpy(var->rst_patch_group->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
          }
        }
      }
    }
  }

  return PIDX_success;
}
