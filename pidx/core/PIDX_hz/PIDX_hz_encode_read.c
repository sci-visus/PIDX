/*****************************************************
 **  PIDX Parallel I/O Library      **
 **  Copyright (c) 2010-2014 University of Utah   **
 **  Scientific Computing and Imaging Institute   **
 **  72 S Central Campus Drive, Room 3750     **
 **  Salt Lake City, UT 84112       **
 **             **
 **  PIDX is licensed under the Creative Commons  **
 **  Attribution-NonCommercial-NoDerivatives 4.0  **
 **  International License. See LICENSE.md.   **
 **             **
 **  For information about this project see:  **
 **  http://www.cedmav.com/pidx       **
 **  or contact: pascucci@sci.utah.edu    **
 **  For support: PIDX-support@visus.net    **
 **             **
 *****************************************************/


#include "../../PIDX_inc.h"


PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  int level = 0, cnt = 0, s = 0, number_levels = 0;
  unsigned long long i = 0, j = 0, k = 0, l = 0;
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

          hz_index = hz_order - var0->hz_buffer->start_hz_index[level];
          int v1;
          for(v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            hz_index = hz_order - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
            bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;

            memcpy(var_grp->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype), var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);

            /*
            if (id->idx_d->color == 0)
            {
              double x1, x2;
              memcpy(&x1, var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype), sizeof(double));
              memcpy(&x2, var_grp->variable[v1]->hz_buffer->buffer[level] + (hz_index * bytes_for_datatype) + sizeof(double), sizeof(double));
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

          index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
              + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
              + (k - chunked_patch_offset[2]);

          int v1 = 0;
          for(v1 = id->first_index; v1 <= id->last_index; v1++)
          {
            hz_index = hz_order - var_grp->variable[v1]->hz_buffer->start_hz_index[level];
            bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
            for (s = 0; s < var_grp->variable[v1]->vps; s++)
              memcpy(var_grp->variable[v1]->chunked_super_patch->restructured_patch->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype * chunk_size, var_grp->variable[v1]->hz_buffer->buffer[level] + ((hz_index * var_grp->variable[v1]->vps + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);

          }
        }
  }

  return PIDX_success;
}

// Correct
PIDX_return_code PIDX_hz_encode_read_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index)
{
  unsigned long long index = 0;
  int m = 0;
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

  if(var0->data_layout == PIDX_row_major)
  {
    unsigned long long hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (j = start_hz_index; j < end_hz_index; j++)
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

            memcpy(var_grp->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   var_grp->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                 bytes_for_datatype);
          }
        }
      }
    }
  }
  else if (var0->data_layout == PIDX_column_major)
  {
    unsigned long long hz_mins = 0, hz_maxes = 0, mindex = 0;
    for (j = start_hz_index; j < end_hz_index; j++)
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

            memcpy(var_grp->variable[v1]->chunked_super_patch->restructured_patch->buffer + (index * bytes_for_datatype),
                   var_grp->variable[v1]->hz_buffer->buffer[j] + (hz_index * bytes_for_datatype),
                   bytes_for_datatype);
          }
        }
      }
    }
  }

  return PIDX_success;
}
