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


// Creates the buffers for all the patches that constitutes the super patch
PIDX_return_code PIDX_idx_rst_buf_create(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];

  int j = 0, v = 0;
  int cnt = 0, i = 0;

  // Allocate buffer for restructured patches (super patch) for all variables
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    cnt = 0;
    PIDX_variable var = var_grp->variable[v];

    // Iterate through all the super patch a process intersects with
    for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      // If the process is the target rank of a super patch then allocate buffer for the patches of the super patch
      if (rst_id->idx_comm_metadata->grank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->restructured_super_patch;
        // Iterate through all the patches of the super patch and allocate buffer for them
        for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

          if (patch_group->patch[j]->buffer == NULL)
          {
            fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
            return PIDX_err_rst;
          }

          memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
        }

        cnt++;
        assert(cnt == 1);  // This is gauranteed because every process can hold at max only one suoer patch
      }
    }

    // This is gauranteed because every process can hold at max only one suoer patch
    if (cnt != var->restructured_super_patch_count)
    {
      fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
  }

  return PIDX_success;
}


// Free all the patches that constitute a super patch
PIDX_return_code PIDX_idx_rst_buf_destroy(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int j, v;

  // Iterate through all the variables
  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    // Iterate through all the patches of the super patch
    for(j = 0; j < var_grp->variable[v]->restructured_super_patch->patch_count; j++)
    {
      free(var->restructured_super_patch->patch[j]->buffer);
      var->restructured_super_patch->patch[j]->buffer = 0;
    }
  }

  return PIDX_success;
}


// Create the buffer that holds all the patches of a super patch into one single patch
PIDX_return_code PIDX_idx_rst_aggregate_buf_create(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    // The restructured patch is the sum of all the small patches
    PIDX_patch out_patch = var->restructured_super_patch->restructured_patch;

    var->restructured_super_patch->restructured_patch->buffer = malloc(out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var->bpv/8) * var->vps);
    memset(var->restructured_super_patch->restructured_patch->buffer, 0, out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var->bpv/8) * var->vps);

    //printf("[Restructuring] Buffer size of %d = %d\n", rst_id->idx_comm_metadata->grank, out_patch->size[0] * out_patch->size[1] * out_patch->size[2]);

    if (var->restructured_super_patch->restructured_patch->buffer == NULL)
    {
      fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
  }

  return PIDX_success;
}


// Free the restructured patch
PIDX_return_code PIDX_idx_rst_aggregate_buf_destroy(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    free(var->restructured_super_patch->restructured_patch->buffer);
  }

  return PIDX_success;
}



// Combine the small patches of a super patch into a restructured patch
PIDX_return_code PIDX_idx_rst_buf_aggregate(PIDX_idx_rst_id rst_id, int mode)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    PIDX_super_patch rst_super_patch = var->restructured_super_patch;
    PIDX_patch out_patch = rst_super_patch->restructured_patch;

    int nx = out_patch->size[0];
    int ny = out_patch->size[1];

    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;

    // Iterate through all the small patches of the super patch
    for (r = 0; r < var->restructured_super_patch->patch_count; r++)
    {
      PIDX_patch *patch = var->restructured_super_patch->patch;
      for (k1 = patch[r]->offset[2]; k1 < patch[r]->offset[2] + patch[r]->size[2]; k1++)
      {
        for (j1 = patch[r]->offset[1]; j1 < patch[r]->offset[1] + patch[r]->size[1]; j1++)
        {
          for (i1 = patch[r]->offset[0]; i1 < patch[r]->offset[0] + patch[r]->size[0]; i1 = i1 + patch[r]->size[0])
          {
            index = ((patch[r]->size[0])* (patch[r]->size[1]) * (k1 - patch[r]->offset[2])) + ((patch[r]->size[0]) * (j1 - patch[r]->offset[1])) + (i1 - patch[r]->offset[0]);
            send_o = index * var->vps * (var->bpv/8);
            send_c = (patch[r]->size[0]);
            recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

            if (mode == PIDX_WRITE)
              memcpy(rst_super_patch->restructured_patch->buffer + (recv_o * var->vps * (var->bpv/8)), patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
            else
              memcpy(patch[r]->buffer + send_o, rst_super_patch->restructured_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
          }
        }
      }
    }
  }

  return PIDX_success;
}
