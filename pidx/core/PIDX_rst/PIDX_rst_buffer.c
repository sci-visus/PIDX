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


PIDX_return_code PIDX_rst_buf_create(PIDX_rst_id rst_id)
{
  int j = 0, v = 0, p = 0;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  if(rst_id->idx->enable_rst == 1)
  {
    int cnt = 0, i = 0;
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      cnt = 0;
      for (i = 0; i < rst_id->reg_patch_grp_count; i++)
      {
        if (rst_id->idx_c->grank == rst_id->reg_patch_grp[i]->max_patch_rank)
        {
          Ndim_patch_group patch_group = var->rst_patch_group[cnt];
          patch_group->data_source = 0;
          for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
          {
            patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

            if (patch_group->patch[j]->buffer == NULL)
              return PIDX_err_rst;
            memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
          }
          cnt++;
        }
      }
      if (cnt != var->patch_group_count)
        return PIDX_err_rst;
    }
  }
  else
  {
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p];
        for(j = 0; j < patch_group->count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps);

          memset(patch_group->patch[j]->buffer, 0, patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps);
        }
      }
    }
  }

  //assert(cnt == num_output_buffers);
  return PIDX_success;
}




PIDX_return_code PIDX_rst_buf_destroy(PIDX_rst_id rst_id)
{
  int i, j, v;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for(i = 0; i < var_grp->variable[v]->patch_group_count; i++)
    {
      for(j = 0; j < var_grp->variable[v]->rst_patch_group[i]->count; j++)
      {
        free(var->rst_patch_group[i]->patch[j]->buffer);
        var->rst_patch_group[i]->patch[j]->buffer = 0;
      }
    }
  }
  return PIDX_success;
}


PIDX_return_code PIDX_rst_aggregate_buf_create(PIDX_rst_id rst_id)
{
  int v = 0;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      Ndim_patch_group patch_group = var->rst_patch_group[g];
      Ndim_patch out_patch = var->rst_patch_group[g]->reg_patch;

      int nx = out_patch->size[0];
      int ny = out_patch->size[1];
      int nz = out_patch->size[2];

      patch_group->data_source = 1;
      var->rst_patch_group[g]->reg_patch->buffer = malloc(nx * ny * nz * (var->bpv/8) * var->vps);
      memset(var->rst_patch_group[g]->reg_patch->buffer, 0, nx * ny * nz * (var->bpv/8) * var->vps);

      if (var->rst_patch_group[g]->reg_patch->buffer == NULL)
        return PIDX_err_chunk;
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_rst_aggregate_buf_destroy(PIDX_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      free(var->rst_patch_group[g]->reg_patch->buffer);
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_rst_buf_aggregate(PIDX_rst_id rst_id, int mode)
{
#if 1
  int v = 0;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      Ndim_patch_group patch_group = var->rst_patch_group[g];
      Ndim_patch out_patch = var->rst_patch_group[g]->reg_patch;

      int nx = out_patch->size[0];
      int ny = out_patch->size[1];

      patch_group->data_source = 1;
      int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
      for (r = 0; r < var->rst_patch_group[g]->count; r++)
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
                memcpy(out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group[g]->patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
              else
                memcpy(var->rst_patch_group[g]->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
            }
          }
        }
      }
    }
  }
#else
  int v = 0;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  // loop through all groups
  int g = 0;
  for (g = 0; g < var0->patch_group_count; ++g)
  {
    // copy the size and offset to output
    Ndim_patch_group patch_group0 = var0->rst_patch_group[g];
    Ndim_patch out_patch0 = var0->rst_patch_group[g]->reg_patch;

    int nx = out_patch0->size[0];
    int ny = out_patch0->size[1];

    patch_group0->data_source = 1;
    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
    for (r = 0; r < var0->rst_patch_group[g]->count; r++)
    {
      for (k1 = patch_group0->patch[r]->offset[2]; k1 < patch_group0->patch[r]->offset[2] + patch_group0->patch[r]->size[2]; k1++)
      {
        for (j1 = patch_group0->patch[r]->offset[1]; j1 < patch_group0->patch[r]->offset[1] + patch_group0->patch[r]->size[1]; j1++)
        {
          for (i1 = patch_group0->patch[r]->offset[0]; i1 < patch_group0->patch[r]->offset[0] + patch_group0->patch[r]->size[0]; i1 = i1 + patch_group0->patch[r]->size[0])
          {
            index = ((patch_group0->patch[r]->size[0])* (patch_group0->patch[r]->size[1]) * (k1 - patch_group0->patch[r]->offset[2])) + ((patch_group0->patch[r]->size[0]) * (j1 - patch_group0->patch[r]->offset[1])) + (i1 - patch_group0->patch[r]->offset[0]);

            send_c = (patch_group0->patch[r]->size[0]);
            recv_o = (nx * ny * (k1 - out_patch0->offset[2])) + (nx * (j1 - out_patch0->offset[1])) + (i1 - out_patch0->offset[0]);

            for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
            {
              PIDX_variable var = var_grp->variable[v];
              Ndim_patch out_patch = var->rst_patch_group[g]->reg_patch;

              send_o = index * var->vps * (var->bpv/8);
              if (mode == PIDX_WRITE)
                memcpy(out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group[g]->patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
              else
                memcpy(var->rst_patch_group[g]->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
            }
          }
        }
      }
    }
  }

#endif

  return PIDX_success;
}
