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


#include "../../PIDX_inc.h"


/// Populates the following for every HZ level: nsamples_per_level, start_hz_index, end_hz_index
/// These are populated by calling visus function (align)
PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id)
{
  PIDX_variable var0 = id->idx->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  // Bounding box of the super patch the process is holding
  int **restructured_box;

  // For every HZ level what is the starting xyz coordinate
  int **start_xyz_per_hz_level;

  // For every HZ level what is the ending xyz coordinate
  int **end_xyz_per_hz_level;

  int maxH = id->idx->maxh;

  restructured_box = (int**) malloc(2 * sizeof (int*));
  memset(restructured_box, 0, 2 * sizeof (int*));
  restructured_box[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  restructured_box[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(restructured_box[0], 0, PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(restructured_box[1], 0, PIDX_MAX_DIMENSIONS * sizeof (int));

  if (maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH [%d] not set.\n", __FILE__, __LINE__, maxH);
    return PIDX_err_hz;
  }

  for (uint32_t v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    var->hz_buffer = malloc(sizeof(*(var->hz_buffer)));
    memset(var->hz_buffer, 0, sizeof(*(var->hz_buffer)));

    HZ_buffer hz_buf = var->hz_buffer;

    // if the patch is an edge patch or not
    hz_buf->is_boundary_HZ_buffer = var->chunked_super_patch->is_boundary_patch;

    hz_buf->start_hz_index = malloc(sizeof (uint64_t) * maxH);
    hz_buf->end_hz_index = malloc(sizeof (uint64_t) * maxH);
    memset(hz_buf->start_hz_index, 0, sizeof (uint64_t) * maxH);
    memset(hz_buf->end_hz_index, 0, sizeof (uint64_t) * maxH);

    hz_buf->nsamples_per_level = malloc(sizeof (int*) * maxH);
    memset(hz_buf->nsamples_per_level, 0, sizeof (int*) * maxH);

    start_xyz_per_hz_level = malloc(sizeof (int*) * maxH);
    end_xyz_per_hz_level = malloc(sizeof (int*) * maxH);
    memset(start_xyz_per_hz_level, 0, sizeof (int*) * maxH);
    memset(end_xyz_per_hz_level, 0, sizeof (int*) * maxH);

    for (uint32_t j = 0; j < maxH; j++)
    {
      start_xyz_per_hz_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(start_xyz_per_hz_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

      end_xyz_per_hz_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(end_xyz_per_hz_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

      hz_buf->nsamples_per_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(hz_buf->nsamples_per_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
    }

#if 0
    for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      restructured_box[0][d] = var->sim_patch[0]->offset[d];
      restructured_box[1][d] = var->sim_patch[0]->offset[d] + var->sim_patch[0]->size[d] - 1;
    }
    // In case we want to write a subset of the resolution and not all levels
    for (uint32_t j = 0; j < maxH - id->resolution_to; j++)
    {
      // Visus API call to compute start and end HZ for every HZ level
      Align((maxH - 1), j, id->idx->bitPattern, restructured_box, start_xyz_per_hz_level, end_xyz_per_hz_level, hz_buf->nsamples_per_level);

      Point3D startXYZ;
      startXYZ.x = start_xyz_per_hz_level[j][0];
      startXYZ.y = start_xyz_per_hz_level[j][1];
      startXYZ.z = start_xyz_per_hz_level[j][2];
      hz_buf->start_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, startXYZ);

      Point3D endXYZ;
      endXYZ.x = end_xyz_per_hz_level[j][0];
      endXYZ.y = end_xyz_per_hz_level[j][1];
      endXYZ.z = end_xyz_per_hz_level[j][2];

      hz_buf->end_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, endXYZ);
      if (id->idx_c->simulation_rank == 3)
        printf("Test %d : %lld - %lld -- %d %d %d\n", j, hz_buf->start_hz_index[j], hz_buf->end_hz_index[j], hz_buf->nsamples_per_level[j][0], hz_buf->nsamples_per_level[j][1], hz_buf->nsamples_per_level[j][2]);
    }
#endif

    for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      restructured_box[0][d] = var->chunked_super_patch->restructured_patch->offset[d] / id->idx->chunk_size[d];

      if (var->chunked_super_patch->restructured_patch->size[d] % id->idx->chunk_size[d] == 0)
        restructured_box[1][d] = (var->chunked_super_patch->restructured_patch->offset[d] / id->idx->chunk_size[d]) + (var->chunked_super_patch->restructured_patch->size[d] / id->idx->chunk_size[d]) - 1;
      else
        restructured_box[1][d] = (var->chunked_super_patch->restructured_patch->offset[d] / id->idx->chunk_size[d]) + ((var->chunked_super_patch->restructured_patch->size[d] / id->idx->chunk_size[d]) + 1) - 1;
    }

    // In case we want to write a subset of the resolution and not all levels
    for (uint32_t j = 0; j < maxH - id->resolution_to; j++)
    {
      // Visus API call to compute start and end HZ for every HZ level
      Align((maxH - 1), j, id->idx->bitPattern, restructured_box, start_xyz_per_hz_level, end_xyz_per_hz_level, hz_buf->nsamples_per_level);

      Point3D startXYZ;
      startXYZ.x = start_xyz_per_hz_level[j][0];
      startXYZ.y = start_xyz_per_hz_level[j][1];
      startXYZ.z = start_xyz_per_hz_level[j][2];
      hz_buf->start_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, startXYZ);

      Point3D endXYZ;
      endXYZ.x = end_xyz_per_hz_level[j][0];
      endXYZ.y = end_xyz_per_hz_level[j][1];
      endXYZ.z = end_xyz_per_hz_level[j][2];

      hz_buf->end_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, endXYZ);
    }

    for (uint32_t j = 0; j < maxH; j++)
    {
      free(start_xyz_per_hz_level[j]);
      free(end_xyz_per_hz_level[j]);
    }
    free(start_xyz_per_hz_level);
    free(end_xyz_per_hz_level);
  }

  free(restructured_box[0]);
  restructured_box[0] = 0;
  free(restructured_box[1]);
  restructured_box[1] = 0;
  free(restructured_box);
  restructured_box = 0;

  return PIDX_success;
}



/// free the HZ related meta data
PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id)
{
  PIDX_variable var0 = id->idx->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (uint32_t v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    free(var->hz_buffer->start_hz_index);
    free(var->hz_buffer->end_hz_index);

    for (int itr = 0; itr < id->idx->maxh; itr++)
      free(var->hz_buffer->nsamples_per_level[itr]);
    free(var->hz_buffer->nsamples_per_level);

    free(var->hz_buffer->buffer);
    var->hz_buffer->buffer = 0;

    free(var->hz_buffer);
    var->hz_buffer = 0;
  }

  return PIDX_success;
}
