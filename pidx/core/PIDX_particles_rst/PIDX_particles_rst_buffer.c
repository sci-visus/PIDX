/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2010-2018 ViSUS L.L.C.,
 * Scientific Computing and Imaging Institute of the University of Utah
 *
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 *
 */

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
PIDX_return_code PIDX_particles_rst_buf_create(PIDX_particles_rst_id rst_id)
{
  int j = 0, v = 0;
  int cnt = 0, i = 0;

  // Allocate buffer for restructured patches (super patch) for all variables
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    cnt = 0;
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    // Iterate through all the super patch a process intersects with
    for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      // If the process is the target rank of a super patch then allocate buffer for the patches of the super patch
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->restructured_super_patch;
        // Iterate through all the patches of the super patch and allocate buffer for them
        for (j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->particle_count * var->vps * var->bpv/8);

          if (patch_group->patch[j]->buffer == NULL)
          {
            fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
            return PIDX_err_rst;
          }

          memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->particle_count * var->vps * var->bpv/8));
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
PIDX_return_code PIDX_particles_rst_buf_destroy(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int j, v;

  // Iterate through all the variables
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    // Iterate through all the patches of the super patch
    for (j = 0; j < rst_id->idx_metadata->variable[v]->restructured_super_patch->patch_count; j++)
    {
      free(var->restructured_super_patch->patch[j]->buffer);
      var->restructured_super_patch->patch[j]->buffer = 0;
    }
  }

  return PIDX_success;
}


// Create the buffer that holds all the patches of a super patch into one single patch
PIDX_return_code PIDX_particles_rst_aggregate_buf_create(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  for (int v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    for (int i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      // If the process is the target rank of a super patch then allocate buffer for the patches of the super patch
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->restructured_super_patch;
        int total_particle_count = 0;
        // Iterate through all the patches of the super patch and allocate buffer for them
        for (int j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          total_particle_count = total_particle_count + patch_group->patch[j]->particle_count;
        }
        patch_group->restructured_patch->buffer = malloc(total_particle_count * var->vps * var->bpv/8);
        printf("[%d %d] my rank %d TPC %d\n", i, v, rst_id->idx_c->simulation_rank, total_particle_count);
        if (patch_group->restructured_patch->buffer == NULL)
        {
          fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        memset(patch_group->restructured_patch->buffer, 0, (total_particle_count * var->vps * var->bpv/8));
      }
    }
  }

  return PIDX_success;
}


// Free the restructured patch
PIDX_return_code PIDX_particles_rst_aggregate_buf_destroy(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    free(var->restructured_super_patch->restructured_patch->buffer);
  }

  return PIDX_success;
}



// Combine the small patches of a super patch into a restructured patch
PIDX_return_code PIDX_particles_rst_buf_aggregate(PIDX_particles_rst_id rst_id, int mode)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    PIDX_super_patch rst_super_patch = var->restructured_super_patch;
    PIDX_patch out_patch = rst_super_patch->restructured_patch;

    int nx = out_patch->size[0];
    int ny = out_patch->size[1];

    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;

    // Iterate through all the small patches of the super patch
    PIDX_patch *patch = var->restructured_super_patch->patch;
    for (r = 0; r < var->restructured_super_patch->patch_count; r++)
    {
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
