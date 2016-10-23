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

struct hz_tupple_double
{
  //pointer to variables -> pointer to samples per variable -> actual data
  unsigned char*** value;
  int index;
};
typedef struct hz_tupple_double hz_tupple;

static int compare( const void* a, const void* b);

PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, c = 0, s = 0, y = 0, number_levels = 0;
  unsigned long long i = 0, j = 0, k = 0, v = 0, l = 0, d = 0;
  int v1 = 0;
  int index_count = 0;
  int bytes_for_datatype;
  unsigned long long hz_index;
  unsigned long long total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
  int rank = 0;
  MPI_Comm_rank(id->comm, &rank);

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

  //printf("var0->patch_group_count = %d\n", var0->patch_group_count);
#if 1
  for (y = 0; y < var0->patch_group_count; y++)
  {
    if (var0->chunk_patch_group[y]->type == 0)
    {
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

                //if (level >= maxH - id->resolution_to)
                //if (level >= maxH)
                //  continue;

                index = (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                    + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                    + (i - chunked_patch_offset[0]);

                tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                for(v = id->first_index; v <= id->last_index; v++)
                {
                  PIDX_variable var = var_grp->variable[v];

                  if (level < maxH - id->resolution_to)
                    var->hz_buffer[y]->samples_per_level[level] = var->hz_buffer[y]->samples_per_level[level] + 1;

                  tupple[index_count].value[v - id->first_index] = malloc(var->vps * sizeof(unsigned char*));
                  memset(tupple[index_count].value[v - id->first_index], 0, var->vps * sizeof(unsigned char*));

                  //bytes_for_datatype = var->bpv / 8;
                  bytes_for_datatype = ((var->bpv / 8) * chunk_size) / id->idx->compression_factor;
                  for (s = 0; s < var->vps; s++)
                  {
                    tupple[index_count].value[v - id->first_index][s] = malloc(bytes_for_datatype);
                    memset(tupple[index_count].value[v - id->first_index][s], 0, bytes_for_datatype);

                    memcpy(tupple[index_count].value[v- id->first_index][s], var->chunk_patch_group[y]->patch[b]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

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

                //if (level >= maxH - id->resolution_to)
                //if (level >= maxH)
                //  continue;

                index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                    + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                    + (k - chunked_patch_offset[2]);

                tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                for(v1 = id->first_index; v1 <= id->last_index; v1++)
                {
                  if (level < maxH - id->resolution_to)
                    var_grp->variable[v1]->hz_buffer[y]->samples_per_level[level] = var_grp->variable[v1]->hz_buffer[y]->samples_per_level[level] + 1;

                  tupple[index_count].value[v1 - id->first_index] = malloc(var_grp->variable[v1]->vps * sizeof(unsigned char*));
                  memset(tupple[index_count].value[v1 - id->first_index], 0, var_grp->variable[v1]->vps * sizeof(unsigned char*));

                  //bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
                  bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size) / id->idx->compression_factor;
                  for (s = 0; s < var_grp->variable[v1]->vps; s++)
                  {
                    tupple[index_count].value[v1 - id->first_index][s] = malloc(bytes_for_datatype);
                    memset(tupple[index_count].value[v1 - id->first_index][s], 0, bytes_for_datatype);

                    memcpy(tupple[index_count].value[v1 - id->first_index][s], var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                    tupple[index_count].index = hz_order;
                  }
                }
                index_count++;
              }
        }
        qsort( tupple, total_chunked_patch_size, sizeof(hz_tupple), compare );

        for(v1 = id->first_index; v1 <= id->last_index; v1++)
        {
          for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
            //for(c = 0 ; c < maxH ; c++)
          {
            //bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
            bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size) / id->idx->compression_factor;

            var_grp->variable[v1]->hz_buffer[y]->buffer[c] = malloc(bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * var_grp->variable[v1]->vps);
            memset(var_grp->variable[v1]->hz_buffer[y]->buffer[c], 0, bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * var_grp->variable[v1]->vps);
          }
        }


        cnt = 0;
        for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
          //for (c = 0; c < maxH; c++)
        {
          for (s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
          {
            for (v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              for (i = 0; i < var_grp->variable[v1]->vps; i++)
              {
                //bytes_for_datatype = var_grp->variable[v1]->bpv / 8;
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size) / id->idx->compression_factor;
                memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[c] + ((s * var_grp->variable[v1]->vps + i) * bytes_for_datatype), tupple[cnt].value[v1 - id->first_index][i], bytes_for_datatype);
                var_grp->variable[/*id->first_index*/v1]->hz_buffer[y]->buffer_index[cnt] = tupple[cnt].index;
              }
            }
            cnt++;
          }
        }

        /*
    cnt = 0;
    for (c = 0; c < maxH; c++)
    {
      for (s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
      {
      for (v1 = id->first_index; v1 <= id->last_index; v1++)
      {
        for (i = 0; i < var_grp->variable[v1]->vps; i++)
        {
        free(tupple[cnt].value[v1 - id->first_index][i]);
        }
        free(tupple[cnt].value[v1 - id->first_index]);
      }
      free(tupple[cnt].value);
      cnt++;
      }
    }
    free(tupple);
    tupple = 0;
    */

        if (var0->data_layout == PIDX_row_major)
        {
          index_count = 0;
          for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
            for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
              for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
              {
                for(v = id->first_index; v <= id->last_index; v++)
                {
                  PIDX_variable var = var_grp->variable[v];
                  for (s = 0; s < var->vps; s++)
                  {
                    free(tupple[index_count].value[v - id->first_index][s]);
                  }
                  free(tupple[index_count].value[v - id->first_index]);
                }
                free(tupple[index_count].value);
                index_count++;
              }
        }
        free(tupple);



#endif
      }
    }
    else
    {
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
        index_count = 0;
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
#if 1
                unsigned long long lastbitmask = ((unsigned long long) 1) << number_levels;
                z_order |= lastbitmask;
                while (!(1 & z_order)) z_order >>= 1;
                z_order >>= 1;

                hz_order = z_order;

                level = getLeveL(hz_order);
                //printf("%lld %lld -> %lld %d\n", (long long)i, (long long)j, (long long)hz_order, level);

                //if (level >= maxH)
                if (level >= maxH - id->resolution_to)
                  continue;

                for(v1 = id->first_index; v1 <= id->last_index; v1++)
                {
                  hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                  bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
                  //printf("bytes_for_datatype = %d\n", bytes_for_datatype);
                  //for (s = 0; s < var_grp->variable[v1]->vps; s++)
                  //{
                  /*
                    int nprocs;
                    MPI_Comm_size(id->comm, &nprocs);
                    double x;
                    memcpy(&x, var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), bytes_for_datatype);
                    if (nprocs == 1)
                      printf("Dest %d %d Src %d Count %d byte size %d V %f\n", level, (int)((hz_index * var_grp->variable[v1]->vps + s) * chunk_size), (int)((index * var_grp->variable[v1]->vps) + s) * chunk_size, (int)chunk_size, var_grp->variable[v1]->bpv, x);
                  */

                  //memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * var_grp->variable[v1]->vps + s) * bytes_for_datatype), var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                  //if (rank == 1)
                  //{
                  memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype),
                       var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                       bytes_for_datatype);
                  //}
                  /*
            double x;
            memcpy(&x, var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);
            double y1;
            memcpy(&y1, var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), bytes_for_datatype);
            if (rank == 1)
            printf("R %d ----> HZ %d (%d) [V%d Y%d B%d] [%f %f] [Bytes %d]\n", index, hz_index, level, v1, y, b, x, y1, bytes_for_datatype);
            */
                  //}
                }
                #endif
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

                if (level >= maxH - id->resolution_to)
                  //if (level >= maxH)
                  continue;

                index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                    + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                    + (k - chunked_patch_offset[2]);

                for(v1 = id->first_index; v1 <= id->last_index; v1++)
                {
                  hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                  bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size) / (var_grp->variable[v1]->bpv / id->idx->compression_bit_rate);
                  for (s = 0; s < var_grp->variable[v1]->vps; s++)
                  {
                    memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * var_grp->variable[v1]->vps + s) * bytes_for_datatype),
                         var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * var_grp->variable[v1]->vps) + s) * bytes_for_datatype,
                         bytes_for_datatype);
                  }
                }
              }
        }
      }
    }

    /*
  int rank = 0;
  MPI_Comm_rank(id->comm, &rank);
  int missing_sample_count = 0;
  if (var0->chunk_patch_group[y]->type == 2)
  {
    for (i = id->first_index; i <= id->last_index; i++)
    {
    PIDX_variable var = var_grp->variable[i];
    bytes_for_datatype = (var->bpv / 8) * var->vps;
    for (j = id->idx_d->res_from; j < maxH - id->idx_d->res_to; j++)
    {
      int skipped_blocks = 0;
      int skipped_samples = 0;
      if (var0->hz_buffer[y]->missing_block_count_per_level[j] != 0)
      {
      int samples_per_level = var0->hz_buffer[y]->end_hz_index[j] - var0->hz_buffer[y]->start_hz_index[j] + 1;
      missing_sample_count = var0->hz_buffer[y]->missing_block_count_per_level[j] * id->idx_d->samples_per_block;
      int adjusted_buffer_size = (samples_per_level- missing_sample_count) * bytes_for_datatype;

      int sample_start_index = 0;
      for (n = 0; n < var0->hz_buffer[y]->missing_block_count_per_level[j]; n++)
      {
        int sample = 0;
        for (sample = sample_start_index; sample < samples_per_level; sample++)
        {
        if (var0->hz_buffer[y]->start_hz_index[j] + sample == var0->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
        {
          skipped_samples = skipped_blocks * id->idx_d->samples_per_block;
          int dest = (sample - skipped_samples) * bytes_for_datatype;
          int src = (sample + id->idx_d->samples_per_block - skipped_samples) * bytes_for_datatype;
          int count = (samples_per_level - sample - id->idx_d->samples_per_block) * bytes_for_datatype;
          //printf("[%d] src dest count %d %d %d (%d - %d - %d)\n", rank, src, dest, count, samples_per_level, sample, id->idx_d->samples_per_block);
          if (count < 0)
          break;

#if !SIMULATE_IO
          memmove(var->hz_buffer[y]->buffer[j] + dest, var->hz_buffer[y]->buffer[j] + src, count);
#endif
          skipped_blocks++;
          sample_start_index = sample;
          break;
        }
        sample++;
        }
      }
#if !SIMULATE_IO
      //printf("adjusted_buffer_size = %d\n", adjusted_buffer_size);
      unsigned char* temp_buffer = realloc(var->hz_buffer[y]->buffer[j], adjusted_buffer_size);
      if (temp_buffer == NULL)
      {
        fprintf(stderr, "realloc error File %s Line no %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
      else
        var->hz_buffer[y]->buffer[j] = temp_buffer;
#endif
      }
    }
    }
  }
  */
  }
#endif
  return PIDX_success;
}



PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index)
{
  unsigned long long index = 0;
  int b = 0, y = 0, m = 0;
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


  //
  for (y = 0; y < var0->patch_group_count; y++)
  {
#if !SIMULATE_IO
    for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
    {
      total_chunked_patch_size = 0;

      for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
      {
        chunked_patch_offset[l] = var0->chunk_patch_group[y]->patch[b]->offset[l] / id->idx->chunk_size[l];
        chunked_patch_size[l] = var0->chunk_patch_group[y]->patch[b]->size[l] / id->idx->chunk_size[l];
        total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
      }

      if(var0->data_layout == PIDX_row_major)
      {
        unsigned long long hz_mins = 0, hz_maxes = 0, mindex = 0;
        //for (j = 0; j < id->idx_d->maxh - id->resolution_to; j++)
        for (j = start_hz_index; j < end_hz_index; j++)
        {

          //if (rank == 1)
          //  printf("At level %d\n", j);
          if (var_grp->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][0] * var_grp->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][1] * var_grp->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][2] != 0)
          {
            hz_mins = var_grp->variable[id->first_index]->hz_buffer[y]->start_hz_index[j];// out_buf_array[buffer_count]->allign_start_hz[j];
            hz_maxes = var_grp->variable[id->first_index]->hz_buffer[y]->end_hz_index[j] + 1;// out_buf_array[buffer_count]->allign_end_hz[j] + 1;

            //if (rank == 1)
            //  printf("[%d (%d - %d)] ----> %d %d\n", j, maxH, id->resolution_to, hz_mins, hz_maxes);

            for (m = hz_mins; m < hz_maxes; m++)
            {
              mindex = m;
              maxH = id->idx_d->maxh - 1;
              unsigned long long lastbitmask=((unsigned long long)1)<<maxH;

              mindex <<= 1;
              mindex  |= 1;
              while ((lastbitmask & mindex) == 0) mindex <<= 1;
                mindex &= lastbitmask - 1;

              PointND cnt;
              PointND p  ;
              int n = 0;

              memset(&cnt,0,sizeof(PointND));
              memset(&p  ,0,sizeof(PointND));

              for (;mindex; mindex >>= 1,++n, maxH--)
              {
                int bit= id->idx->bitPattern[maxH];
                PGET(p,bit) |= (mindex & 1) << PGET(cnt,bit);
                ++PGET(cnt,bit);
              }

              if (p.x >= id->idx->bounds[0] || p.y >= id->idx->bounds[1] || p.z >= id->idx->bounds[2])
              {
                //printf("A [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
              {
                //printf("B [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
              {
                //printf("C [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              hz_index = m - hz_mins;
              index = (chunked_patch_size[0] * chunked_patch_size[1] * (p.z - chunked_patch_offset[2])) + (chunked_patch_size[0] * (p.y - chunked_patch_offset[1])) + (p.x - chunked_patch_offset[0]);

#if 1
              for(v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                //hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
#if !SIMULATE_IO
                memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[j] + (hz_index * bytes_for_datatype),
                        var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                        bytes_for_datatype);

                /*
                if (rank == 1)
                {
                  double xx;
                  memcpy(&xx, var_grp->variable[v1]->hz_buffer[y]->buffer[j] + (hz_index * bytes_for_datatype), sizeof(double));
                  printf("[%d] value %d %d - %f\n", j, index, hz_index, xx);
                }
                */
#endif
              }
#endif

            }
          }

        }
      }
    }
#endif
  }
  //
  return PIDX_success;
}


static int compare( const void* a, const void* b)
{
  unsigned long long int_a = ((const hz_tupple*)a)->index;
  unsigned long long int_b = ((const hz_tupple*)b)->index;

  if ( int_a == int_b ) return 0;
  else if ( int_a < int_b ) return -1;
  else return 1;
}
