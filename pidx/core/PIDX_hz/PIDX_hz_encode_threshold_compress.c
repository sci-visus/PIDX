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
#include <zfp.h>

PIDX_return_code PIDX_hz_encode_threshold_and_write(PIDX_hz_encode_id id)
{
  unsigned long long z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, s = 0, y = 0, number_levels = 0;
  unsigned long long i = 0, j = 0, k = 0, l = 0;
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

  for (y = 0; y < var0->patch_group_count; y++)
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
      Point3D xyzuv_Index;

      if(var0->data_layout == PIDX_row_major)
      {
        //printf("[%d] -> %d %d %d -- %d %d %d\n", id->idx_c->grank, chunked_patch_offset[0], chunked_patch_offset[1], chunked_patch_offset[2], chunked_patch_size[0], chunked_patch_size[1], chunked_patch_size[2]);
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
                hz_index = hz_order - var_grp->variable[v1]->hz_buffer[y]->start_hz_index[level];
                bytes_for_datatype = ((var_grp->variable[v1]->bpv / 8) * chunk_size * var_grp->variable[v1]->vps) / id->idx->compression_factor;
                //float x;
                //memcpy(&x, var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), sizeof (float));
                //printf("[%d] [%d] value %f\n", level, hz_index, x);

                if ((var_grp->variable[v1]->bpv / 8) == sizeof (float) && strcmp(var_grp->variable[v1]->type_name, FLOAT32) == 0)
                {
                  float temp;
                  memcpy(&temp, var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                         bytes_for_datatype);
                  if (temp < 0.5)
                    temp = 0;

                  memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype),
                       &temp,
                       bytes_for_datatype);
                }

                //memcpy(var_grp->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype),
                //     var_grp->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                //     bytes_for_datatype);
              }
            }
      }
    }

  }
  return PIDX_success;
}


PIDX_return_code PIDX_hz_encode_compress(PIDX_hz_encode_id id)
{
  int p = 0, c = 0, v = 0, bytes_for_datatype = 0;
  int maxH = id->idx_d->maxh;

  // Allocate actual HZ buffer for the variables
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      //var->hz_buffer[p]->buffer = (unsigned char**)malloc( maxH * sizeof (unsigned char*));
      //memset(var->hz_buffer[p]->buffer, 0,  maxH * sizeof (unsigned char*));
      bytes_for_datatype = var->bpv / 8;
      for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
      {
        void* buf = var->hz_buffer[p]->buffer[c];
        int dim_x = var->hz_buffer[p]->nsamples_per_level[c][0];
        int dim_y = var->hz_buffer[p]->nsamples_per_level[c][1];
        int dim_z = var->hz_buffer[p]->nsamples_per_level[c][2];

        // only compress from the second block onwards
        if (c <= id->idx->bits_per_block)
        {
          var->hz_buffer[p]->compressed_buffer_size[c] = dim_x * dim_y * dim_z * bytes_for_datatype;
          continue;
        }
        zfp_type type = (bytes_for_datatype == 4) ? zfp_type_float : zfp_type_double;
        zfp_field* field = zfp_field_3d(buf, type, dim_x, dim_y, dim_z);
        zfp_stream* zfp = zfp_stream_open(NULL);
        zfp_stream_set_accuracy(zfp, 0, type);
        size_t max_compressed_bytes = zfp_stream_maximum_size(zfp, field);
        unsigned char* output = (unsigned char*)malloc(max_compressed_bytes + 16);
        bitstream* stream = stream_open(output + 16, max_compressed_bytes);
        zfp_stream_set_bit_stream(zfp, stream);
        //printf("[%d] [Dim %d %d %d] [BD %d] MCB %d\n", c, dim_x, dim_y, dim_z, type, max_compressed_bytes);
        size_t compressed_bytes = zfp_compress(zfp, field);
        if (compressed_bytes == 0)
          puts("ERROR: Something wrong happened during compression\n");
        var->hz_buffer[p]->compressed_buffer_size[c] = compressed_bytes;
        size_t original_bytes = dim_x * dim_y * dim_z * bytes_for_datatype;
        if (compressed_bytes + 16 > original_bytes)
          puts("WARNING: compressed size > original size\n");
        if (compressed_bytes > 0x7FFFFFFF)
          puts("WARNING: compressed size does not fit in an int");
        *((unsigned int*)output) = (unsigned int)compressed_bytes; // first 4 bytes = size of compressed stream
        ((int*)output)[1] = dim_x; // next 4 bytes = dim x
        ((int*)output)[2] = dim_y; // next 4 bytes = dim y
        ((int*)output)[3] = dim_z; // next 4 bytes = dim z
        free(buf);
        var->hz_buffer[p]->buffer[c] = output;
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
      }
    }
  }

  return PIDX_success;
}
