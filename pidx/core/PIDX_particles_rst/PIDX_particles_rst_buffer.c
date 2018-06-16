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
#if 1
  int cnt = 0;

  // Allocate buffer for restructured patches (super patch) for all variables
  for (int v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    cnt = 0;
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    // Iterate through all the super patch a process intersects with
    for (int i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      // If the process is the target rank of a super patch then allocate buffer for the patches of the super patch
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->restructured_super_patch;
        // Iterate through all the patches of the super patch and allocate buffer for them
        for (int j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->particle_count * var->vps * var->bpv / CHAR_BIT);

          if (patch_group->patch[j]->buffer == NULL)
          {
            fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
            return PIDX_err_rst;
          }

          memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->particle_count * var->vps * var->bpv / CHAR_BIT));
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

#endif
  return PIDX_success;
}


// Free all the patches that constitute a super patch
PIDX_return_code PIDX_particles_rst_buf_destroy(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  // Iterate through all the variables
  for (int v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    // Iterate through all the patches of the super patch
    for (int j = 0; j < rst_id->idx_metadata->variable[v]->restructured_super_patch->patch_count; j++)
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
          total_particle_count = total_particle_count + patch_group->patch[j]->particle_count;

        patch_group->restructured_patch->buffer = malloc(total_particle_count * var->vps * var->bpv / CHAR_BIT);
        patch_group->restructured_patch->particle_count = total_particle_count;
        //printf("var %d PC %d Buffer %d\n", v, patch_group->restructured_patch->particle_count, (total_particle_count * var->vps * var->bpv / CHAR_BIT));

        if (patch_group->restructured_patch->buffer == NULL)
        {
          fprintf(stderr, "[%s] [%d] malloc() failed.\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        memset(patch_group->restructured_patch->buffer, 0, (total_particle_count * var->vps * var->bpv / CHAR_BIT));
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

  for (int v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    free(var->restructured_super_patch->restructured_patch->buffer);
  }

  return PIDX_success;
}



// Combine the small patches of a super patch into a restructured patch
PIDX_return_code PIDX_particles_rst_buf_aggregate(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (uint32_t v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    int pcount = 0;
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    for (uint32_t r = 0; r < var->restructured_super_patch->patch_count; r++)
    {
      memcpy(var->restructured_super_patch->restructured_patch->buffer + (pcount * var->vps * (var->bpv / CHAR_BIT)), var->restructured_super_patch->patch[r]->buffer, var->restructured_super_patch->patch[r]->particle_count * var->vps * (var->bpv / CHAR_BIT));
      pcount = pcount + var->restructured_super_patch->patch[r]->particle_count;
    }
  }

#if 0
  char rank_filename[PATH_MAX];
  sprintf(rank_filename, "%s_b_%d", rst_id->idx_metadata->filename, rst_id->idx_c->simulation_rank);
  FILE *fp = fopen(rank_filename, "w");

  for (uint32_t p = 0; p < rst_id->idx_metadata->variable[0]->restructured_super_patch->restructured_patch->particle_count; p++)
  {
    //
    double px, py, pz;
    memcpy(&px, rst_id->idx_metadata->variable[0]->restructured_super_patch->restructured_patch->buffer + (p * 3 + 0) * sizeof(double), sizeof(double));
    memcpy(&py, rst_id->idx_metadata->variable[0]->restructured_super_patch->restructured_patch->buffer + (p * 3 + 1) * sizeof(double), sizeof(double));
    memcpy(&pz, rst_id->idx_metadata->variable[0]->restructured_super_patch->restructured_patch->buffer + (p * 3 + 2) * sizeof(double), sizeof(double));

    double d1;
    memcpy(&d1, rst_id->idx_metadata->variable[1]->restructured_super_patch->restructured_patch->buffer + (p) * sizeof(double), sizeof(double));

    double d2;
    memcpy(&d2, rst_id->idx_metadata->variable[2]->restructured_super_patch->restructured_patch->buffer + (p) * sizeof(double), sizeof(double));

    double d3;
    memcpy(&d3, rst_id->idx_metadata->variable[3]->restructured_super_patch->restructured_patch->buffer + (p) * sizeof(double), sizeof(double));

    int i1;
    memcpy(&i1, rst_id->idx_metadata->variable[4]->restructured_super_patch->restructured_patch->buffer + (p) * sizeof(int), sizeof(int));
    /*
    double e1, e2, e3, e4, e5, e6, e7, e8, e9;
    memcpy(&e1, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 0) * sizeof(double), sizeof(double));
    memcpy(&e2, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 1) * sizeof(double), sizeof(double));
    memcpy(&e3, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 2) * sizeof(double), sizeof(double));
    memcpy(&e4, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 3) * sizeof(double), sizeof(double));
    memcpy(&e5, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 4) * sizeof(double), sizeof(double));
    memcpy(&e6, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 5) * sizeof(double), sizeof(double));
    memcpy(&e7, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 6) * sizeof(double), sizeof(double));
    memcpy(&e8, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 7) * sizeof(double), sizeof(double));
    memcpy(&e9, rst_id->idx_metadata->variable[5]->restructured_super_patch->restructured_patch->buffer + (p * 9 + 8) * sizeof(double), sizeof(double));
    */
    //fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", e1, e2, e3, e4, e5, e6, e7, e8, e9);

    fprintf(fp, "%f %f %f %f %f %f %d\n", px, py, pz, d1, d2, d3, i1);
  }
  fclose(fp);
#endif


  return PIDX_success;
}
