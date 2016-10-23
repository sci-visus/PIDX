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

PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  //int n = 0, m = 0, d = 0, c = 0;
  //int index_count = 0;
  int b = 0, level = 0, cnt = 0, s = 0, y = 0, number_levels = 0;
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

  for (y = 0; y < var0->patch_group_count; y++)
  {
    /*
    if (var0->chunk_patch_group[y]->type == 2)
    {
      for (i = id->first_index; i <= id->last_index; i++)
      {
        PIDX_variable var = var_grp->variable[i];
        bytes_for_datatype = (var->bpv / 8) * var->vps;
        for (j = 0; j < maxH; j++)
        {
          if (var0->hz_buffer[y]->missing_block_count_per_level[j] != 0)
          {
            int samples_per_level = var0->hz_buffer[y]->end_hz_index[j] - var0->hz_buffer[y]->start_hz_index[j] + 1;

            int sample_start_index = 0;
            for (n = 0; n < var0->hz_buffer[y]->missing_block_count_per_level[j]; n++)
            {
              int sample = 0;
              for (sample = sample_start_index; sample < samples_per_level; sample++)
              {
                if (var0->hz_buffer[y]->start_hz_index[j] + sample == var0->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
                {
                  int src = (sample) * bytes_for_datatype;
                  int dest = (sample + id->idx_d->samples_per_block) * bytes_for_datatype;
                  int count = (samples_per_level - sample - id->idx_d->samples_per_block) * bytes_for_datatype;
                  //printf("Dest %d Source %d Count %d (%d)\n", src, dest, count, (var0->hz_buffer[y]->end_hz_index[j] - (var0->hz_buffer[y]->start_hz_index[j] + sample + id->idx_d->samples_per_block) + 1) * bytes_for_datatype);

                  memmove(var->hz_buffer[y]->buffer[j] + dest, var->hz_buffer[y]->buffer[j] + src, count);
                  //memset(var->hz_buffer[y]->buffer[j] + dest, 0, count);
                  sample_start_index = sample;
                  break;
                }
                sample++;
              }
            }
          }
        }
      }
    }
    */

    if (var0->chunk_patch_group[y]->type == 0)
    {
#if 0
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
#if 1
        total_chunked_patch_size = 1;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          chunked_patch_offset[d] = var0->chunk_patch_group[y]->patch[b]->offset[d] / id->idx->chunk_size[d];
          chunked_patch_size[d] = var0->chunk_patch_group[y]->patch[b]->size[d] / id->idx->chunk_size[d];
          total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[d];
        }

        index_count = 0;
        hz_tupple *tupple = malloc(total_chunked_patch_size * sizeof(hz_tupple));
        memset(tupple, 0, total_chunked_patch_size * sizeof(hz_tupple));

        number_levels = maxH - 1;
        PointND xyzuv_Index;

        if (var0->data_layout == PIDX_row_major)
        {
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                  {
                    xyzuv_Index.x = i;
                    xyzuv_Index.y = j;
                    xyzuv_Index.z = k;

                    z_order = 0;
                    PointND zero;
                    zero.x = 0;
                    zero.y = 0;
                    zero.z = 0;
                    memset(&zero, 0, sizeof (PointND));

                    for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
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

                    index = (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                        + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                        + (i - chunked_patch_offset[0]);

                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;

                    for(v = id->first_index; v <= id->last_index; v++)
                    {
                      PIDX_variable var = var_grp->variable[v];
                      tupple[index_count].value[v - id->first_index] = malloc(var->vps * sizeof(unsigned char*));
                      memset(tupple[index_count].value[v - id->first_index], 0, var->vps * sizeof(unsigned char*));

                      bytes_for_datatype = var->bpv / 8;

                      for (s = 0; s < var->vps; s++)
                      {
                        tupple[index_count].value[v - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memset(tupple[index_count].value[v - id->first_index][s], 0, bytes_for_datatype * chunk_size);

                        memcpy(tupple[index_count].value[v- id->first_index][s], var->chunk_patch_group[y]->patch[b]->buffer + ((index * var->vps) + s) * bytes_for_datatype * chunk_size, bytes_for_datatype * chunk_size);

                        tupple[index_count].index = hz_order;
                      }
                    }
                    index_count++;
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
                    xyzuv_Index.u = u;
                    xyzuv_Index.v = m;

                    z_order = 0;
                    PointND zero;
                    zero.x = 0;
                    zero.y = 0;
                    zero.z = 0;
                    zero.u = 0;
                    zero.v = 0;
                    memset(&zero, 0, sizeof (PointND));

                    for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
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

                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;

                    for(var = id->first_index; var <= id->last_index; var++)
                    {
                      tupple[index_count].value[var - id->first_index] = malloc(var_grp->variable[var]->vps * sizeof(unsigned char*));
                      memset(tupple[index_count].value[var - id->first_index], 0, var_grp->variable[var]->vps * sizeof(unsigned char*));

                      bytes_for_datatype = var_grp->variable[var]->bpv / 8;

                      for (s = 0; s < var_grp->variable[var]->vps; s++)
                      {
                        tupple[index_count].value[var - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memset(tupple[index_count].value[var - id->first_index][s], 0, bytes_for_datatype * chunk_size);

                        memcpy(tupple[index_count].value[var - id->first_index][s], var_grp->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[var]->vps) + s) * bytes_for_datatype * chunk_size, bytes_for_datatype * chunk_size);

                        tupple[index_count].index = hz_order;
                      }
                    }
                    index_count++;
                  }

        }
        qsort( tupple, total_chunked_patch_size, sizeof(hz_tupple), compare );

        for(var = id->first_index; var <= id->last_index; var++)
        {
          for(c = 0 ; c < maxH ; c++)
          {
            bytes_for_datatype = var_grp->variable[var]->bpv / 8;
            var_grp->variable[var]->hz_buffer[y]->buffer[c] = malloc(bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * var_grp->variable[var]->vps * chunk_size);
            memset(var_grp->variable[var]->hz_buffer[y]->buffer[c], 0, bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * var_grp->variable[var]->vps * chunk_size);
          }
        }

        cnt = 0;
        for(c = 0; c < maxH; c++)
        {
          for(s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
          {
            for(var = id->first_index; var <= id->last_index; var++)
            {
              for (i = 0; i < var_grp->variable[var]->vps; i++)
              {
                bytes_for_datatype = var_grp->variable[var]->bpv / 8;
                memcpy(var_grp->variable[var]->hz_buffer[y]->buffer[c] + ((s * var_grp->variable[var]->vps + i) * bytes_for_datatype * chunk_size), tupple[cnt].value[var - id->first_index][i], bytes_for_datatype * chunk_size);
                var_grp->variable[id->first_index]->hz_buffer[y]->buffer_index[cnt] = tupple[cnt].index;

                free(tupple[cnt].value[var - id->first_index][i]);
              }
              free(tupple[cnt].value[var - id->first_index]);
            }
            free(tupple[cnt].value);
            cnt++;
          }
        }
        free(tupple);
        tupple = 0;
#endif
      }
#endif
    }
    else
    {
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
        total_chunked_patch_size = 0;
        for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
        {
          chunked_patch_offset[l] = var0->chunk_patch_group[y]->patch[b]->offset[l] / id->idx->chunk_size[l];
          chunked_patch_size[l] = var0->chunk_patch_group[y]->patch[b]->size[l] / id->idx->chunk_size[l];
          total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
        }

        number_levels = maxH - 1;
        PointND xyzuv_Index;

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
                    PointND zero;
                    zero.x = 0;
                    zero.y = 0;
                    zero.z = 0;
                    memset(&zero, 0, sizeof (PointND));

                    for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
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

                    //if (level >= maxH - id->idx_d->resolution_to - 1)
                    //  continue;

                    hz_index = hz_order - var0->hz_buffer[y]->start_hz_index[level];
                    int v1;
                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;


  #if !SIMULATE_IO
                      /*
                      double x1, x2, x3, x4;
                      memcpy(&x1, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), sizeof(double));
                      memcpy(&x2, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + sizeof(double), sizeof(double));
                      memcpy(&x3, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + 2*sizeof(double), sizeof(double));
                      memcpy(&x4, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + 3*sizeof(double), sizeof(double));
                      printf("XXXXXXXXXXXXx x1 = %f %f %f %f\n", x1, x2, x3, x4);
                      */
                      memcpy(var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);
  #endif
                      /*
                      bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
                      for (s = 0; s < var_grp->variable[v1]->vps; s++)
                      {
                        memcpy(var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype * chunk_size, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * var_grp->variable[v1]->vps + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);
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
                    PointND zero;
                    zero.x = 0;
                    zero.y = 0;
                    zero.z = 0;
                    memset(&zero, 0, sizeof (PointND));

                    for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
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
                      hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
                      for (s = 0; s < var_grp->variable[v1]->vps; s++)
                      {
                        memcpy(var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype * chunk_size, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * var_grp->variable[v1]->vps + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);
                      }
                    }
                  }
        }
      }
    }
  }

  return PIDX_success;
}
