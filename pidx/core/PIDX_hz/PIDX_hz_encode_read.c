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


PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  uint64_t z_order = 0, hz_order = 0, index = 0;
  int level = 0, cnt = 0, s = 0, number_levels = 0;
  uint64_t i = 0, j = 0, k = 0, l = 0;
  int bytes_for_datatype;
  uint64_t hz_index;
  uint64_t total_chunked_patch_size = 1;

  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable var0 = id->idx->variable[id->first_index];

  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};

  if (var0->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->count not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }

  if (maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }



  total_chunked_patch_size = 0;
  for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
  {
    chunked_patch_offset[l] = var0->chunked_super_patch->restructured_patch->offset[l] / id->idx->chunk_size[l];
    if (var0->chunked_super_patch->restructured_patch->size[l] % id->idx->chunk_size[l] == 0)
      chunked_patch_size[l] = var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l];
    else
      chunked_patch_size[l] = (var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l]) + 1;

    total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
  }

  number_levels = maxH - 1;
  Point3D xyzuv_Index;

  if (var0->data_layout == PIDX_row_major)
  {
    for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
      for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
        for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
        {
          index = (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
              + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
              + (i - chunked_patch_offset[0]);

          xyzuv_Index.x = i;
          xyzuv_Index.y = j;
          xyzuv_Index.z = k;

          z_order = 0;
          Point3D zero;
          zero.x = 0;
          zero.y = 0;
          zero.z = 0;
          memset(&zero, 0, sizeof (Point3D));

          for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
          {
            int bit = id->idx->bitPattern[number_levels];
            z_order |= ((uint64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
            PGET(xyzuv_Index, bit) >>= 1;
          }

          number_levels = maxH - 1;
          uint64_t lastbitmask = ((uint64_t) 1) << number_levels;
          z_order |= lastbitmask;
          while (!(1 & z_order)) z_order >>= 1;
          z_order >>= 1;

          hz_order = z_order;

          level = getLeveL(hz_order);

          hz_index = hz_order - var0->hz_buffer->start_hz_index[level];
          int v1;
          for (v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
            bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype), id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);

            /*
            if (id->idx_d->color == 0)
            {
              double x1, x2;
              memcpy(&x1, id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype), sizeof(double));
              memcpy(&x2, id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype) + sizeof(double), sizeof(double));
              fprintf(stderr, "xyz index - hz index :: %d - %d ----> %f %f\n", index, hz_index, x1, x2);
            }
            */

          }
        }
  }
  else
  {
    for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
      for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
        for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
        {
          xyzuv_Index.x = i;
          xyzuv_Index.y = j;
          xyzuv_Index.z = k;

          z_order = 0;
          Point3D zero;
          zero.x = 0;
          zero.y = 0;
          zero.z = 0;
          memset(&zero, 0, sizeof (Point3D));

          for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
          {
            int bit = id->idx->bitPattern[number_levels];
            z_order |= ((uint64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
            PGET(xyzuv_Index, bit) >>= 1;
          }

          number_levels = maxH - 1;
          uint64_t lastbitmask = ((uint64_t) 1) << number_levels;
          z_order |= lastbitmask;
          while (!(1 & z_order)) z_order >>= 1;
          z_order >>= 1;

          hz_order = z_order;
          level = getLeveL(hz_order);

          index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
              + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
              + (k - chunked_patch_offset[2]);

          int v1 = 0;
          for (v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
            bytes_for_datatype = id->idx->variable[v1]->bpv / 8;
            for (s = 0; s < id->idx->variable[v1]->vps; s++)
              memcpy(id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + ((index * id->idx->variable[v1]->vps) + s) * bytes_for_datatype * chunk_size, id->idx->variable[v1]->hz_buffer->buffer[level] + ((hz_index * id->idx->variable[v1]->vps + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);

          }
        }
  }

  return PIDX_success;
}

// Correct
PIDX_return_code PIDX_hz_encode_read_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index)
{
  uint64_t index = 0;
  int m = 0;
  uint64_t j = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  uint64_t hz_index;
  uint64_t total_chunked_patch_size = 1;
  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable var0 = id->idx->variable[id->first_index];

  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};

  if (var0->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->count not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }

  if (maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }

  total_chunked_patch_size = 0;
  for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
  {
    chunked_patch_offset[l] = var0->chunked_super_patch->restructured_patch->offset[l] / id->idx->chunk_size[l];
    if (var0->chunked_super_patch->restructured_patch->size[l] % id->idx->chunk_size[l] == 0)
      chunked_patch_size[l] = var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l];
    else
      chunked_patch_size[l] = (var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l]) + 1;

    total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
  }

  if (var0->data_layout == PIDX_row_major)
  {
    uint64_t hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (j = start_hz_index; j < end_hz_index; j++)
    {
      if (id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
      {
        hz_mins = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j];
        hz_maxes = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

        for (m = hz_mins; m < hz_maxes; m++)
        {
          mindex = m;
          maxH = id->idx->maxh - 1;
          uint64_t lastbitmask=((uint64_t)1)<<maxH;

          mindex <<= 1;
          mindex  |= 1;
          while ((lastbitmask & mindex) == 0) mindex <<= 1;
          mindex &= lastbitmask - 1;

          Point3D cnt;
          Point3D p  ;
          int n = 0;

          memset(&cnt,0,sizeof(Point3D));
          memset(&p  ,0,sizeof(Point3D));

          for (;mindex; mindex >>= 1,++n, maxH--)
          {
            int bit= id->idx->bitPattern[maxH];
            PGET(p,bit) |= (mindex & 1) << PGET(cnt,bit);
            ++PGET(cnt,bit);
          }

          if (p.x >= id->idx->box_bounds[0] || p.y >= id->idx->box_bounds[1] || p.z >= id->idx->box_bounds[2])
            continue;

          if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
            continue;

          if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
            continue;

          hz_index = m - hz_mins;
          index = (chunked_patch_size[0] * chunked_patch_size[1] * (p.z - chunked_patch_offset[2])) + (chunked_patch_size[0] * (p.y - chunked_patch_offset[1])) + (p.x - chunked_patch_offset[0]);

          for (v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   id->idx->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                 bytes_for_datatype);
          }
        }
      }
    }
  }
  else if (var0->data_layout == PIDX_column_major)
  {
    uint64_t hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (j = start_hz_index; j < end_hz_index; j++)
    {
      if (id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
      {
        hz_mins = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j];
        hz_maxes = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

        for (m = hz_mins; m < hz_maxes; m++)
        {
          mindex = m;
          maxH = id->idx->maxh - 1;
          uint64_t lastbitmask=((uint64_t)1)<<maxH;

          mindex <<= 1;
          mindex  |= 1;
          while ((lastbitmask & mindex) == 0) mindex <<= 1;
          mindex &= lastbitmask - 1;

          Point3D cnt;
          Point3D p  ;
          int n = 0;

          memset(&cnt,0,sizeof(Point3D));
          memset(&p  ,0,sizeof(Point3D));

          for (;mindex; mindex >>= 1,++n, maxH--)
          {
            int bit= id->idx->bitPattern[maxH];
            PGET(p,bit) |= (mindex & 1) << PGET(cnt,bit);
            ++PGET(cnt,bit);
          }

          if (p.x >= id->idx->box_bounds[0] || p.y >= id->idx->box_bounds[1] || p.z >= id->idx->box_bounds[2])
            continue;

          if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
            continue;

          if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
            continue;

          hz_index = m - hz_mins;
          index = (chunked_patch_size[2] * chunked_patch_size[1] * (p.x - chunked_patch_offset[0])) + (chunked_patch_size[2] * (p.y - chunked_patch_offset[1])) + (p.z - chunked_patch_offset[2]);

          for (v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   id->idx->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                   bytes_for_datatype);
          }
        }
      }
    }
  }

  return PIDX_success;
}
