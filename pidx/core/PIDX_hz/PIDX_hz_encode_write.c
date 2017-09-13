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


PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, s = 0, number_levels = 0;
  unsigned long long i = 0, j = 0, k = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  unsigned long long hz_index;
  unsigned long long total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  PIDX_variable var0 = var_grp->variable[id->first_index];

  if (var0->patch_group_count == 0)
    return PIDX_success;

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


  for (b = 0; b < var0->chunk_patch_group->patch_count; b++)
  {
    total_chunked_patch_size = 0;
    for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
    {
      chunked_patch_offset[l] = var0->chunk_patch_group->patch[b]->offset[l] / id->idx->chunk_size[l];
      chunked_patch_size[l] = var0->chunk_patch_group->patch[b]->size[l] / id->idx->chunk_size[l];
      total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
    }

    number_levels = maxH - 1;
    Point3D xyzuv_Index;

    if(var0->data_layout == PIDX_row_major)
    {
      //fprintf(stderr, "[%d] -> %d %d %d -- %d %d %d\n", id->idx_c->grank, chunked_patch_offset[0], chunked_patch_offset[1], chunked_patch_offset[2], chunked_patch_size[0], chunked_patch_size[1], chunked_patch_size[2]);
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
              z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            number_levels = maxH - 1;
            unsigned long long lastbitmask = ((unsigned long long) 1) << number_levels;
            z_order |= lastbitmask;
            while (!(1 & z_order)) z_order >>= 1;
            z_order >>= 1;

            hz_order = z_order;

            level = getLeveL(hz_order);

            if (level >= maxH - id->resolution_to)
              continue;

            for(v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              hz_index = hz_order - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
              bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
              //float x;
              //memcpy(&x, var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype), sizeof (float));
              //fprintf(stderr, "[%d] [%d] value %f\n", level, hz_index, x);

              memcpy(var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                   var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);
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
              z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            number_levels = maxH - 1;
            unsigned long long lastbitmask = ((unsigned long long) 1) << number_levels;
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

            for(v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              hz_index = hz_order - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
              bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
              for (s = 0; s < var_grp->variable[v1]->vps; s++)
              {
                memcpy(var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }
            }
          }
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_hz_encode_row_major_write(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, number_levels = 0;
  unsigned long long i = 0, j = 0, k = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  unsigned long long hz_index;
  unsigned long long offset_hz_index;
  unsigned long long total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
  const char* bit_string = id->idx->bitSequence + 1;
  int bs_len = strlen(bit_string);
  int bits_per_block = id->idx->bits_per_block;
  int block_index = 0;
  int total_number_of_blocks = (id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]) / id->idx_d->samples_per_block;
  Point3D inter_block_index;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  PIDX_variable var0 = var_grp->variable[id->first_index];

  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0};

  for (b = 0; b < var0->chunk_patch_group->patch_count; b++)
  {
    for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
    {
      chunked_patch_offset[l] = var0->chunk_patch_group->patch[b]->offset[l] / id->idx->chunk_size[l];
      chunked_patch_size[l] = var0->chunk_patch_group->patch[b]->size[l] / id->idx->chunk_size[l];
    }
  }


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

  Point3D* block_offset = malloc(sizeof(*block_offset) * total_number_of_blocks);
  memset(block_offset, 0, sizeof(*block_offset) * total_number_of_blocks);
  Point3D* block_nsamples = malloc(sizeof(*block_nsamples) * maxH);
  memset(block_nsamples, 0, sizeof(*block_nsamples) * maxH);
  Point3D* intra_stride = malloc(sizeof(*intra_stride) * maxH);
  memset(intra_stride, 0, sizeof(*intra_stride) * maxH);

#if 0
  Point3D* block_from = malloc(sizeof(*block_from) * maxH);
  memset(block_from, 0, sizeof(*block_from) * maxH);
  Point3D* block_to = malloc(sizeof(*block_to) * maxH);
  memset(block_to, 0, sizeof(*block_to) * maxH);
  Point3D* interleaved_block_nsamples = malloc(sizeof(*interleaved_block_nsamples) * maxH);
  memset(interleaved_block_nsamples, 0, sizeof(*interleaved_block_nsamples) * maxH);

  Point3D sub_vol_from = {chunked_patch_offset[0],
              chunked_patch_offset[1],
              chunked_patch_offset[2]};
  Point3D sub_vol_to = {(chunked_patch_offset[0] + chunked_patch_size[0] - 1),
              (chunked_patch_offset[1] + chunked_patch_size[1] - 1),
              (chunked_patch_offset[2] + chunked_patch_size[2] - 1)};

#endif

  int* samples_in_level = malloc(id->idx_d->maxh * sizeof(*samples_in_level));
  memset(samples_in_level, 0, (id->idx_d->maxh * sizeof(*samples_in_level)));
  for (j = 0; j < maxH; j++)
  {
    block_nsamples[j] = get_num_samples_per_block(bit_string, bs_len, j, bits_per_block);
    intra_stride[j] = get_intra_block_strides(bit_string, bs_len, j);
    samples_in_level[j] = var0->hz_buffer->nsamples_per_level[j][0] * var0->hz_buffer->nsamples_per_level[j][1] * var0->hz_buffer->nsamples_per_level[j][2];

#if 0
    if (j > id->idx->bits_per_block)
    {
      Point3D stride;
      get_grid( sub_vol_from, sub_vol_to, j, bit_string, bs_len, &(block_from[j]), &(block_to[j]), &stride);

      interleaved_block_nsamples[j].x = (block_to[j].x - block_from[j].x) / stride.x + 1;
      interleaved_block_nsamples[j].y = (block_to[j].y - block_from[j].y) / stride.y + 1;
      interleaved_block_nsamples[j].z = (block_to[j].z - block_from[j].z) / stride.z + 1;
    }
#endif
  }

  for (b = 0; b < var0->chunk_patch_group->patch_count; b++)
  {
    total_chunked_patch_size = 0;
    for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
    {
      chunked_patch_offset[l] = var0->chunk_patch_group->patch[b]->offset[l] / id->idx->chunk_size[l];
      chunked_patch_size[l] = var0->chunk_patch_group->patch[b]->size[l] / id->idx->chunk_size[l];
      total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
    }

    number_levels = maxH - 1;
    Point3D xyzuv_Index;

    if(var0->data_layout == PIDX_row_major)
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
              z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            number_levels = maxH - 1;
            unsigned long long lastbitmask = ((unsigned long long) 1) << number_levels;
            z_order |= lastbitmask;
            while (!(1 & z_order)) z_order >>= 1;
            z_order >>= 1;

            hz_order = z_order;
            level = getLeveL(hz_order);

            if (level >= maxH - id->resolution_to)
              continue;

            block_index = hz_order / id->idx_d->samples_per_block;

            if (samples_in_level[level] >= (int)pow(2, id->idx->bits_per_block))
            {
              if (hz_order % id->idx_d->samples_per_block == 0)
              {
                block_offset[block_index].x = i;
                block_offset[block_index].y = j;
                block_offset[block_index].z = k;

                inter_block_index.x = 0;
                inter_block_index.y = 0;
                inter_block_index.z = 0;
              }
              else
              {
                inter_block_index.x = (i - block_offset[block_index].x) / intra_stride[level].x;
                inter_block_index.y = (j - block_offset[block_index].y) / intra_stride[level].y;
                inter_block_index.z = (k - block_offset[block_index].z) / intra_stride[level].z;
              }

              offset_hz_index = (block_nsamples[level].x * block_nsamples[level].y * inter_block_index.z) +
                  (block_nsamples[level].x * inter_block_index.y) +
                  inter_block_index.x;

              for(v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                hz_index = block_index * id->idx_d->samples_per_block + offset_hz_index - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
                memcpy(var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }
            }
            else
            {
              for(v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                hz_index = hz_order - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
                memcpy(var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }
            }

#if 0
            if (var0->hz_buffer->nsamples_per_level[level][0] * var0->hz_buffer->nsamples_per_level[level][1] * var0->hz_buffer->nsamples_per_level[level][2] < (int)pow(2, id->idx->bits_per_block) && level > id->idx->bits_per_block)
            {
              inter_block_index.x = (i - block_from[level].x) / intra_stride[level].x;
              inter_block_index.y = (j - block_from[level].y) / intra_stride[level].y;
              inter_block_index.z = (k - block_from[level].z) / intra_stride[level].z;

              offset_hz_index = (interleaved_block_nsamples[level].x * interleaved_block_nsamples[level].y * inter_block_index.z) +
                  (interleaved_block_nsamples[level].x * inter_block_index.y) +
                  inter_block_index.x;

              //if (level == 3)
              //  fprintf(stderr, "XX [%d] [R %d] [%d %d %d] [B %d] F [%d %d %d] T [%d %d %d] S [%d %d %d] C [%d %d %d] Index %d\n", level, id->idx_c->grank, i, j, k, block_index, block_from[level].x, block_from[level].y, block_from[level].z, block_to[level].x, block_to[level].y, block_to[level].z, intra_stride[level].x, intra_stride[level].y, intra_stride[level].z, interleaved_block_nsamples[level].x, interleaved_block_nsamples[level].y, interleaved_block_nsamples[level].z, offset_hz_index);

              for(v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                hz_index = block_index;// * id->idx_d->samples_per_block + offset_hz_index;
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;

                memcpy(var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype),
                     var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                     bytes_for_datatype);
              }
              continue;
            }
#endif
          }
    }
  }


  free(samples_in_level);
  free(intra_stride);
  free(block_nsamples);
  free(block_offset);
#if 0
  free(block_from);
  free(block_to);
  free(interleaved_block_nsamples);
#endif
  return PIDX_success;
}




PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int res_from, int res_to)
{
  unsigned long long index = 0;
  int b = 0, m = 0;
  unsigned long long j = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  unsigned long long hz_index;
  unsigned long long total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  PIDX_variable var0 = var_grp->variable[id->first_index];

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


  for (b = 0; b < var0->chunk_patch_group->patch_count; b++)
  {
    total_chunked_patch_size = 0;

    for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
    {
      chunked_patch_offset[l] = var0->chunk_patch_group->patch[b]->offset[l] / id->idx->chunk_size[l];
      chunked_patch_size[l] = var0->chunk_patch_group->patch[b]->size[l] / id->idx->chunk_size[l];
      total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
    }

    if(var0->data_layout == PIDX_row_major)
    {
      unsigned long long hz_mins = 0, hz_maxes = 0, mindex = 0;
      for (j = res_from; j < res_to; j++)
      {
        if (var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
        {
          hz_mins = var_grp->variable[id->first_index]->hz_buffer->start_hz_index[j];
          hz_maxes = var_grp->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

          for (m = hz_mins; m < hz_maxes; m++)
          {
            mindex = m;
            maxH = id->idx_d->maxh - 1;
            unsigned long long lastbitmask=((unsigned long long)1)<<maxH;

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

            //if (p.x >= id->idx->bounds[0] || p.y >= id->idx->bounds[1] || p.z >= id->idx->bounds[2])
            //  continue;

            if (p.x >= id->idx->box_bounds[0] || p.y >= id->idx->box_bounds[1] || p.z >= id->idx->box_bounds[2])
              continue;

            if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
              continue;

            if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
              continue;

            hz_index = m - hz_mins;
            index = (chunked_patch_size[0] * chunked_patch_size[1] * (p.z - chunked_patch_offset[2])) + (chunked_patch_size[0] * (p.y - chunked_patch_offset[1])) + (p.x - chunked_patch_offset[0]);

            for(v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;

              memcpy(var_grp->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                   var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);
            }
          }
        }
      }
    }
    else if (var0->data_layout == PIDX_column_major)
    {
      unsigned long long hz_mins = 0, hz_maxes = 0, mindex = 0;
      for (j = res_from; j < res_to; j++)
      {
        if (var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][0] * var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][1] * var_grp->variable[id->first_index]->hz_buffer->nsamples_per_level[j][2] != 0)
        {
          hz_mins = var_grp->variable[id->first_index]->hz_buffer->start_hz_index[j];
          hz_maxes = var_grp->variable[id->first_index]->hz_buffer->end_hz_index[j] + 1;

          for (m = hz_mins; m < hz_maxes; m++)
          {
            mindex = m;
            maxH = id->idx_d->maxh - 1;
            unsigned long long lastbitmask=((unsigned long long)1)<<maxH;

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

            //if (p.x >= id->idx->bounds[0] || p.y >= id->idx->bounds[1] || p.z >= id->idx->bounds[2])
            //  continue;

            if (p.x >= id->idx->box_bounds[0] || p.y >= id->idx->box_bounds[1] || p.z >= id->idx->box_bounds[2])
              continue;

            if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
              continue;

            if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
              continue;

            hz_index = m - hz_mins;
            index = (chunked_patch_size[2] * chunked_patch_size[1] * (p.x - chunked_patch_offset[0])) + (chunked_patch_size[2] * (p.y - chunked_patch_offset[1])) + (p.z - chunked_patch_offset[2]);

            for(v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;

              memcpy(var_grp->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                   var_grp->variable[v1]->chunk_patch_group->patch[b]->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);
            }
          }
        }
      }
    }
  }


  return PIDX_success;
}

