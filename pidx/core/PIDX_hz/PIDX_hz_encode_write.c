/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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


PIDX_return_code PIDX_hz_encode_fast_write(PIDX_hz_encode_id id)
{
#if 0
  uint64_t z_order = 0, hz_order = 0, index = 0;
  int level = 0, cnt = 0, s = 0, number_levels = 0;
  uint64_t i = 0, j = 0, k = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  uint64_t hz_index;
  uint64_t total_chunked_patch_size = 1;
  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
  PIDX_variable var0 = id->idx->variable[id->first_index];
  int index_count = 0;

  // Basic checking
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

  // The process is not holding any patch after restructuring
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  // adjusted patch size and offset due to zfp chunking and compression
  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
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
    uint64_t b = 0;
    uint64_t hz_mins = 0, hz_maxes = 0;
    uint64_t block_min = 0, block_max = 0;
    for (j = 0; j < id->idx->maxh; j++)
    {
      if (id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
      {
        hz_mins = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j];
        hz_maxes = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

        block_min = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j] / id->idx->samples_per_block;
        block_max = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] / id->idx->samples_per_block;

        for (b = block_min; b < block_max; b++)
        {
          for (s = 0; s < id->idx->samples_per_block; s++)
          {
            uint64_t index = (chunked_patch_size[0] * chunked_patch_size[1] * (p.z - chunked_patch_offset[2])) + (chunked_patch_size[0] * (p.y - chunked_patch_offset[1])) + (p.x - chunked_patch_offset[0]);

            for (v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

              memcpy(id->idx->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                   id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);
            }
          }
        }
      }
    }
  }
#endif

  return PIDX_success;
}



// In this function we iterate through all the samples in the xyz order (application order), compute their HZ index and put them correctly in the hz buffer
PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  uint64_t z_order = 0, hz_order = 0, index = 0, hz_index = 0;
  int level = 0, number_levels = 0, bytes_for_datatype = 0, index_count = 0;

  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
  PIDX_variable var0 = id->idx->variable[id->first_index];

  // Basic checking
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

  // The process is not holding any patch after restructuring
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  // adjusted patch size and offset due to zfp chunking and compression
  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  for (int l = 0; l < PIDX_MAX_DIMENSIONS; l++)
  {
    chunked_patch_offset[l] = var0->chunked_super_patch->restructured_patch->offset[l] / id->idx->chunk_size[l];
    if (var0->chunked_super_patch->restructured_patch->size[l] % id->idx->chunk_size[l] == 0)
      chunked_patch_size[l] = var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l];
    else
      chunked_patch_size[l] = (var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l]) + 1;
  }

  number_levels = maxH - 1;
  Point3D xyzuv_Index;

  // This is for caching HZ indices
  // If caching is enabled then meta_data_cache will not be null
  PIDX_metadata_cache hz_cache = id->meta_data_cache;

  if (hz_cache != NULL)
  {
    // sets up the HZ cache (this is only done once)
    if (hz_cache->is_set == 0)
    {
      //if (id->idx_c->partition_rank == 0)
      //  fprintf(stderr, "Cache Setup\n");

      // The number of elements equals to the number of sample in the patch
      hz_cache->element_count = chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2];

      hz_cache->hz_level = malloc(hz_cache->element_count * sizeof(*hz_cache->hz_level));
      memset(hz_cache->hz_level, 0, hz_cache->element_count * sizeof(*hz_cache->hz_level));

      hz_cache->index_level = malloc(hz_cache->element_count * sizeof(*hz_cache->index_level));
      memset(hz_cache->index_level, 0, hz_cache->element_count * sizeof(*hz_cache->index_level));

      hz_cache->xyz_mapped_index = malloc(hz_cache->element_count * sizeof(*hz_cache->xyz_mapped_index));
      memset(hz_cache->xyz_mapped_index, 0, hz_cache->element_count * sizeof(*hz_cache->xyz_mapped_index));

      if (var0->data_layout == PIDX_row_major)
      {
        // For every sample in the patch, find the corresponding HZ level and HZ index
        // and copy the data to HZ encoded buffer (corresponding to a HZ level and a buffer for the HZ level)
        for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
          for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
            for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
            {
              index = (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                  + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                  + (i - chunked_patch_offset[0]);

              hz_cache->xyz_mapped_index[index_count] = index;

              xyzuv_Index.x = i;
              xyzuv_Index.y = j;
              xyzuv_Index.z = k;

              z_order = 0;
              Point3D zero;
              zero.x = 0;
              zero.y = 0;
              zero.z = 0;
              memset(&zero, 0, sizeof (Point3D));

              for (int cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
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

              // Global HZ order
              hz_order = z_order;

              // HZ level
              level = getLeveL(hz_order);
              hz_cache->hz_level[index_count] = level;

              if (level >= maxH - id->resolution_to)
                continue;

              for (int v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                // Local HZ order for every process
                hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
                hz_cache->index_level[index_count] = hz_index;

                bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;
                memcpy(id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }

              index_count++;
            }
      }
      else
      {
        for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
          for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
            for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
            {

              index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                      + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                      + (k - chunked_patch_offset[2]);

              hz_cache->xyz_mapped_index[index_count] = index;

              xyzuv_Index.x = i;
              xyzuv_Index.y = j;
              xyzuv_Index.z = k;

              z_order = 0;
              Point3D zero;
              zero.x = 0;
              zero.y = 0;
              zero.z = 0;
              memset(&zero, 0, sizeof (Point3D));

              for (int cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
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
              hz_cache->hz_level[index_count] = level;

              if (level >= maxH - id->resolution_to)
                continue;

              for (int v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
                hz_cache->index_level[index_count] = hz_index;

                bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;
                for (int s = 0; s < id->idx->variable[v1]->vps; s++)
                {
                  memcpy(id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                       id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                       bytes_for_datatype);
                }
              }

              index_count++;
            }
      }
      // The cache is setup
      hz_cache->is_set = 1;
    }

    // This condition is for using the cache
    else
    {
      //if (id->idx_c->partition_rank == 0)
      //  fprintf(stderr, "Cache Used\n");

      if (var0->data_layout == PIDX_row_major)
      {
        for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
          for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
            for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
            {
              index = (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                  + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                  + (i - chunked_patch_offset[0]);

              assert(index == hz_cache->xyz_mapped_index[index_count]);

              if (hz_cache->hz_level[index_count] >= maxH - id->resolution_to)
                continue;

              for (int v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
                bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

                memcpy(id->idx->variable[v1]->hz_buffer->buffer[hz_cache->hz_level[index_count]] + (hz_cache->index_level[index_count] * bytes_for_datatype),
                    id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (hz_cache->xyz_mapped_index[index_count] * bytes_for_datatype),
                    bytes_for_datatype);
              }

              index_count++;
            }

      }
      else
      {
        for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
          for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
            for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
            {

              if (hz_cache->hz_level[index_count] >= maxH - id->resolution_to)
                continue;

              index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                  + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                  + (k - chunked_patch_offset[2]);

              assert(index == hz_cache->xyz_mapped_index[index_count]);

              for (int v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;
                for (int s = 0; s < id->idx->variable[v1]->vps; s++)
                {
                  memcpy(id->idx->variable[v1]->hz_buffer->buffer[hz_cache->hz_level[index_count]] + (hz_cache->index_level[index_count] * bytes_for_datatype),
                      id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (hz_cache->xyz_mapped_index[index_count] * bytes_for_datatype),
                      bytes_for_datatype);
                }
              }
              index_count++;
            }
      }
    }
  }

  // If there is no caching enabled
  else
  {
    if (var0->data_layout == PIDX_row_major)
    {
      for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
        for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
          for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
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

            for (int cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
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

            if (level >= maxH - id->resolution_to)
              continue;

            for (int v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];

              bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

              memcpy(id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                   id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);
            }
          }
    }
    else
    {
      for (uint64_t k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
        for (uint64_t j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
          for (uint64_t i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
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

            for (int cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
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

            if (level >= maxH - id->resolution_to)
              continue;

            index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                + (k - chunked_patch_offset[2]);

            for (int v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              hz_index = hz_order - id->idx->variable[v1]->hz_buffer->start_hz_index[level];
              bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;
              for (int s = 0; s < id->idx->variable[v1]->vps; s++)
              {
                memcpy(id->idx->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }
            }
          }
    }
  }

  return PIDX_success;
}


// Alternate way of computing HZ ordering
// In this function we iterate through all the HZ level and then compute their corresponding xyz index and then mapping them
PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int res_from, int res_to)
{
  uint64_t index = 0, hz_index = 0;
  int bytes_for_datatype;

  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable var0 = id->idx->variable[id->first_index];
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

  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  for (int l = 0; l < PIDX_MAX_DIMENSIONS; l++)
  {
    chunked_patch_offset[l] = var0->chunked_super_patch->restructured_patch->offset[l] / id->idx->chunk_size[l];
    if (var0->chunked_super_patch->restructured_patch->size[l] % id->idx->chunk_size[l] == 0)
      chunked_patch_size[l] = var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l];
    else
      chunked_patch_size[l] = (var0->chunked_super_patch->restructured_patch->size[l] / id->idx->chunk_size[l]) + 1;
  }

  if (var0->data_layout == PIDX_row_major)
  {
    uint64_t hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (uint64_t j = res_from; j < res_to; j++)
    {
      if (id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
      {
        hz_mins = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j];
        hz_maxes = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

        for (uint64_t m = hz_mins; m < hz_maxes; m++)
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

          for (int v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(id->idx->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                 id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                 bytes_for_datatype);
          }
        }
      }
    }
  }
  else if (var0->data_layout == PIDX_column_major)
  {
    uint64_t hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (uint64_t j = res_from; j < res_to; j++)
    {
      if (id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
      {
        hz_mins = id->idx->variable[id->first_index]->hz_buffer->start_hz_index[j];
        hz_maxes = id->idx->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

        for (uint64_t m = hz_mins; m < hz_maxes; m++)
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

          for (int v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            bytes_for_datatype = ((id->idx->variable[v1]->bpv / 8) * chunk_size * id->idx->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(id->idx->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                 id->idx->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                 bytes_for_datatype);
          }
        }
      }
    }
  }

  return PIDX_success;
}
