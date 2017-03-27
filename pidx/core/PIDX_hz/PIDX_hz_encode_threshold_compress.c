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


PIDX_return_code PIDX_hz_encode_block_wise_compress(PIDX_hz_encode_id id)
{
  int block_offset = 0;
  int b = 0, p = 0, c = 0, v = 0, bytes_for_datatype = 0;
  int maxH = id->idx_d->maxh;
  int block_per_level = 0;
  int samples_per_level = 0;
  const char* bit_string = id->idx->bitSequence + 1;
  int bs_len = strlen(bit_string);
  int bits_per_block = id->idx->bits_per_block;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      bytes_for_datatype = var->bpv / 8;
      for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
      {
        int dim_x = var->hz_buffer[p]->nsamples_per_level[c][0];
        int dim_y = var->hz_buffer[p]->nsamples_per_level[c][1];
        int dim_z = var->hz_buffer[p]->nsamples_per_level[c][2];
        samples_per_level = dim_x * dim_y * dim_z;
        //if (id->idx_c->grank == 0)
        //  printf("[%d] -----> %d %d %d\n", c, dim_x, dim_y, dim_z);

        // only compress from the second block onwards
        if (samples_per_level >= id->idx_d->samples_per_block)
          goto bcast;
        else
          var->hz_buffer[p]->compressed_buffer_size[c] = samples_per_level * bytes_for_datatype;
      }
    }
  }

  bcast:
  id->idx->compression_start_level = c;
  //MPI_Allreduce(&c, &(id->idx->compression_start_level), 1, MPI_INT, MPI_MAX, id->idx_c->global_comm);

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      bytes_for_datatype = var->bpv / 8;
      for (c = id->idx->compression_start_level; c < maxH - id->resolution_to; c++)
      {
        int dim_x = var->hz_buffer[p]->nsamples_per_level[c][0];
        int dim_y = var->hz_buffer[p]->nsamples_per_level[c][1];
        int dim_z = var->hz_buffer[p]->nsamples_per_level[c][2];
        samples_per_level = dim_x * dim_y * dim_z;
        block_per_level = samples_per_level /  id->idx_d->samples_per_block;

        Point3D block_nsamples = get_num_samples_per_block(bit_string, bs_len, c, bits_per_block);

        block_offset = 24;
        var->hz_buffer[p]->compressed_buffer_size[c] = 24;
        for (b = 0; b < block_per_level; b++)
        {
          void* buf = var->hz_buffer[p]->buffer[c] + (b * id->idx_d->samples_per_block) * bytes_for_datatype;
          zfp_type type = (bytes_for_datatype == 4) ? zfp_type_float : zfp_type_double;
          zfp_field* field = zfp_field_3d(buf, type, block_nsamples.x, block_nsamples.y, block_nsamples.z);
          zfp_stream* zfp = zfp_stream_open(NULL);
          zfp_stream_set_accuracy(zfp, /*id->idx->zfp_precisison*/0.5, type);
          size_t max_compressed_bytes = zfp_stream_maximum_size(zfp, field);
          unsigned char* output = (unsigned char*)malloc(max_compressed_bytes + 4);
          bitstream* stream = stream_open(output + 4, max_compressed_bytes);
          zfp_stream_set_bit_stream(zfp, stream);

          size_t compressed_bytes = zfp_compress(zfp, field);
          if (compressed_bytes == 0)
            puts("ERROR: Something wrong happened during compression\n");
          var->hz_buffer[p]->compressed_buffer_size[c] = var->hz_buffer[p]->compressed_buffer_size[c] + compressed_bytes + 4;
          size_t original_bytes = block_nsamples.x * block_nsamples.y * block_nsamples.z * bytes_for_datatype;
          if (compressed_bytes + 4 > original_bytes)
            puts("WARNING: compressed size > original size\n");
          if (compressed_bytes > 0x7FFFFFFF)
            puts("WARNING: compressed size does not fit in an int");
          *((unsigned int*)output) = (unsigned int)compressed_bytes; // first 4 bytes = size of compressed stream

          //printf("[%d] size %d %d %d CMP Offset %d CMP Bytes %d\n", b, block_nsamples.x, block_nsamples.y, block_nsamples.z, block_offset, compressed_bytes);
          memcpy(var->hz_buffer[p]->buffer[c] + block_offset, output,  (compressed_bytes + 4));
          free(output);
          zfp_field_free(field);
          zfp_stream_close(zfp);
          stream_close(stream);
          block_offset = block_offset + (compressed_bytes + 4);
        }

        ((int*)var->hz_buffer[p]->buffer[c])[0] = dim_x;
        ((int*)var->hz_buffer[p]->buffer[c])[1] = dim_y;
        ((int*)var->hz_buffer[p]->buffer[c])[2] = dim_z;
        ((int*)var->hz_buffer[p]->buffer[c])[3] = block_nsamples.x;
        ((int*)var->hz_buffer[p]->buffer[c])[4] = block_nsamples.y;
        ((int*)var->hz_buffer[p]->buffer[c])[5] = block_nsamples.z;
      }
    }
  }

#if 0
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      bytes_for_datatype = var->bpv / 8;
      for (c = id->idx->compression_start_level; c < maxH - id->resolution_to; c++)
      {
        int dim_x = var->hz_buffer[p]->nsamples_per_level[c][0];
        int dim_y = var->hz_buffer[p]->nsamples_per_level[c][1];
        int dim_z = var->hz_buffer[p]->nsamples_per_level[c][2];
        samples_per_level = dim_x * dim_y * dim_z;
        block_per_level = samples_per_level /  id->idx_d->samples_per_block;

        block_offset = 24;
        int ldimx = ((int*)var->hz_buffer[p]->buffer[c])[0];
        int ldimy = ((int*)var->hz_buffer[p]->buffer[c])[1];
        int ldimz = ((int*)var->hz_buffer[p]->buffer[c])[2];

        int bdimx = ((int*)var->hz_buffer[p]->buffer[c])[3];
        int bdimy = ((int*)var->hz_buffer[p]->buffer[c])[4];
        int bdimz = ((int*)var->hz_buffer[p]->buffer[c])[5];

        printf("[%d] ldim: %d %d %d bdim %d %d %d\n", c, ldimx, ldimy, ldimz, bdimx, bdimy, bdimz);

        int bc = 0;
        int offset = 0;
        while(bc != (ldimx/bdimx) * (ldimy/bdimy) * (ldimz/bdimz))
        {
          unsigned int compressed_bytes = *(unsigned int*)(var->hz_buffer[p]->buffer[c] + 24 + offset);
          printf("[Block %d] Compressed size %d\n", bc, compressed_bytes);
          offset = offset + compressed_bytes + 4;
          bc++;
        }
      }
    }
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_hz_encode_compress(PIDX_hz_encode_id id)
{
  int p = 0, c = 0, v = 0, bytes_for_datatype = 0;
  int maxH = id->idx_d->maxh;
  int samples_per_level = 0;

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
        samples_per_level = dim_x * dim_y * dim_z;

        // only compress from the second block onwards
        if (samples_per_level < id->idx_d->samples_per_block)
        {
          var->hz_buffer[p]->compressed_buffer_size[c] = samples_per_level * bytes_for_datatype;
          continue;
        }

        zfp_type type = (bytes_for_datatype == 4) ? zfp_type_float : zfp_type_double;
        zfp_field* field = zfp_field_3d(buf, type, dim_x, dim_y, dim_z);
        zfp_stream* zfp = zfp_stream_open(NULL);
        zfp_stream_set_accuracy(zfp, 1.0, type);
        size_t max_compressed_bytes = zfp_stream_maximum_size(zfp, field);
        unsigned char* output = (unsigned char*)malloc(max_compressed_bytes + 16);
        bitstream* stream = stream_open(output + 16, max_compressed_bytes);
        zfp_stream_set_bit_stream(zfp, stream);
        //printf("[%d] [Dim %d %d %d] [BD %d] MCB %d\n", c, dim_x, dim_y, dim_z, type, max_compressed_bytes);
        size_t compressed_bytes = zfp_compress(zfp, field);
        if (compressed_bytes == 0)
          puts("ERROR: Something wrong happened during compression\n");
        var->hz_buffer[p]->compressed_buffer_size[c] = compressed_bytes + 16;
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


        //if (c == maxH - id->resolution_to - 1)
        //  printf("Compressed %d Original %d\n", var->hz_buffer[p]->compressed_buffer_size[c], dim_x * dim_y * dim_z * bytes_for_datatype);
      }
    }
  }

  return PIDX_success;
}
