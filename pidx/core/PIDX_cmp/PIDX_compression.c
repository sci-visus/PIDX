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
 * \file PIDX_cmp.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_compression.h
 *
 */

#include "../../PIDX_inc.h"

#include <zfp.h>

int compress_buffer(PIDX_comp_id comp_id, unsigned char* buffer, int nx, int ny, int nz, char* base_type, int bps, int vps, float bit_rate);
int decompress_buffer(PIDX_comp_id comp_id, unsigned char** buffer, int nx, int ny, int nz, char* base_type, int bps, int vps, float bit_rate);

///Struct for restructuring ID
struct PIDX_comp_id_struct
{
  idx_dataset idx;

  idx_comm idx_c;

  int first_index;
  int last_index;

};


int compress_buffer(PIDX_comp_id comp_id, unsigned char* buffer, int nx, int ny, int nz, char* base_type, int bps, int vps, float bit_rate)
{
  uint64_t total_bytes = 0;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    uint64_t* chunk_dim = comp_id->idx->chunk_size;
    assert(chunk_dim[0] == 4 && chunk_dim[1] == 4 && chunk_dim[2] == 4);
    uint64_t total_chunk_dim = (uint64_t)chunk_dim[0] * (uint64_t)chunk_dim[1] * (uint64_t)chunk_dim[2];
    //uint64_t chunk_bytes = total_chunk_dim * bytes_per_sample;
    //zfp_type type = (bytes_per_sample == 4) ? zfp_type_float : zfp_type_double;

    uint64_t chunk_bytes = total_chunk_dim * bps;
    zfp_type type;
    if (strcmp(base_type, "int") == 0) {
      type = bps == 32 ? zfp_type_int32 : zfp_type_int64;
    }
    else if (strcmp(base_type, "float") == 0) {
      type = bps == 32 ? zfp_type_float : zfp_type_double;
    }
    else {
      assert(0);
    }

    zfp_field* field = zfp_field_3d(NULL, type, nx, ny, nz);
    zfp_stream* zfp = zfp_stream_open(NULL);
    zfp_stream_set_rate(zfp, bit_rate, type, 3, 0);
    uint64_t bytes_max = zfp_stream_maximum_size(zfp, field);
    bytes_max = bytes_max * vps;

    unsigned char* output = malloc(bytes_max);
    bitstream* stream = stream_open(output, bytes_max);
    zfp_stream_set_bit_stream(zfp, stream);
    uint64_t i = 0;
    uint64_t length = (uint64_t)nx * (uint64_t)ny * (uint64_t)nz;

    for (i = 0; i < length * bps * vps; i += chunk_bytes)
    {
      uint64_t bits = 0;
      switch (type) {
      case zfp_type_float:
        bits = zfp_encode_block_float_3(zfp, (float*)(buffer + i));
        break;
      case zfp_type_double:
        bits = zfp_encode_block_double_3(zfp, (double*)(buffer + i));
        break;
      case zfp_type_int32:
        bits = zfp_encode_block_int32_3(zfp, (const int32*)(buffer + i));
        break;
      case zfp_type_int64:
        bits = zfp_encode_block_int64_3(zfp, (const int64*)(buffer + i));
        break;
      default:
        assert(0);
        break;
      }

      assert(bits % CHAR_BIT == 0);
      total_bytes += bits / CHAR_BIT;
    }
    memcpy(buffer, output, total_bytes);
    free(output);
    zfp_stream_close(zfp);
    stream_close(stream);
    zfp_field_free(field);
  }

  return total_bytes;
}


int decompress_buffer(PIDX_comp_id comp_id, unsigned char** buffer, int nx, int ny, int nz, char* base_type, int bps, int vps, float bit_rate)
{
   uint64_t total_bytes = 0;

   if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
   {
     uint64_t* chunk_dim = comp_id->idx->chunk_size;
     assert(chunk_dim[0] == 4 && chunk_dim[1] == 4 && chunk_dim[2] == 4);
     uint64_t total_chunk_dim = (uint64_t)chunk_dim[0] * (uint64_t)chunk_dim[1] * (uint64_t)chunk_dim[2];

     //uint64_t chunk_bytes = total_chunk_dim * bytes_per_sample;
     //zfp_type type = (bytes_per_sample == 4) ? zfp_type_float : zfp_type_double;

     uint64_t chunk_bytes = total_chunk_dim * bps;
     zfp_type type;
     if (strcmp(base_type, "int") == 0) {
       type = bps == 32 ? zfp_type_int32 : zfp_type_int64;
     }
     else if (strcmp(base_type, "float") == 0) {
       type = bps == 32 ? zfp_type_float : zfp_type_double;
     }
     else {
       assert(0);
     }

     zfp_stream* zfp = zfp_stream_open(NULL);
     zfp_stream_set_rate(zfp, bit_rate, type, 3, 0);

     unsigned char* temp_buffer = malloc(nx * ny * nz * bps * vps);
     int compression_factor = bps*8/comp_id->idx->compression_bit_rate;
     bitstream* stream = stream_open(*buffer, (nx * ny * nz * bps * vps) / compression_factor);
     zfp_stream_set_bit_stream(zfp, stream);
     uint64_t i = 0;
     uint64_t length = (uint64_t)nx * (uint64_t)ny * (uint64_t)nz;

     for (i = 0; i < length * bps * vps; i += chunk_bytes)
     {
       uint64_t bits = 0;
       switch (type) {
       case zfp_type_float:
         bits = zfp_decode_block_float_3(zfp, (const float*)(temp_buffer + i));
         break;
       case zfp_type_double:
         bits = zfp_decode_block_double_3(zfp, (const double*)(temp_buffer + i));
         break;
       case zfp_type_int32:
         bits = zfp_decode_block_int32_3(zfp, (const int32*)(temp_buffer + i));
         break;
       case zfp_type_int64:
         bits = zfp_decode_block_int64_3(zfp, (const int64*)(temp_buffer + i));
         break;
       default:
         assert(0);
         break;
       }

       assert(bits % CHAR_BIT == 0);
       total_bytes += bits / CHAR_BIT;
     }

     unsigned char *temp_buffer2 = realloc(*buffer, nx * ny * nz * bps * vps);
     if (temp_buffer2 == NULL)
     {
       fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
       return PIDX_err_rst;
     }
     else
       *buffer = temp_buffer2;
     memcpy(*buffer, temp_buffer, nx * ny * nz * bps * vps);

     free(temp_buffer);
     zfp_stream_close(zfp);
     stream_close(stream);
   }

   return total_bytes;
}

PIDX_comp_id PIDX_compression_init(idx_dataset idx_meta_data,
                                   idx_comm idx_c, int start_var_index, int end_var_index)
{
  PIDX_comp_id comp_id;

  comp_id = malloc(sizeof (*comp_id));
  memset(comp_id, 0, sizeof (*comp_id));

  comp_id->idx = idx_meta_data;
  comp_id->idx_c = idx_c;

  comp_id->first_index = start_var_index;
  comp_id->last_index = end_var_index;

  return comp_id;
}

PIDX_return_code PIDX_compression(PIDX_comp_id comp_id)
{
  if (comp_id->idx->compression_type == PIDX_NO_COMPRESSION || comp_id->idx->compression_type == PIDX_CHUNKING_ONLY)
    return PIDX_success;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    int v;

    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = comp_id->idx->variable[v];
      PIDX_patch patch = var->chunked_super_patch->restructured_patch;
      unsigned char* buffer = patch->buffer;
      int nx = patch->size[0];
      int ny = patch->size[1];
      int nz = patch->size[2];
      float bit_rate = comp_id->idx->compression_bit_rate;


      int ncomps = 0;
      int bits = 0;
      char base_type[10];
      PIDX_decompose_type(var->type_name, base_type, &ncomps, &bits);
      //PIDX_get_datatype_details(var->type_name, &values, &bits);


      size_t compressed_bytes = compress_buffer(comp_id, buffer, nx, ny, nz, bits/CHAR_BIT, base_type, ncomps, bit_rate);
      unsigned char* temp_buffer = realloc(patch->buffer, compressed_bytes);
      if (temp_buffer == NULL)
        return PIDX_err_compress;
      else
        patch->buffer = temp_buffer;
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_decompression(PIDX_comp_id comp_id)
{
  if (comp_id->idx->compression_type == PIDX_NO_COMPRESSION || comp_id->idx->compression_type == PIDX_CHUNKING_ONLY)
    return PIDX_success;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    int v, ret = 0;
    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = comp_id->idx->variable[v];

      PIDX_patch patch = var->chunked_super_patch->restructured_patch;
      int nx = patch->size[0];
      int ny = patch->size[1];
      int nz = patch->size[2];
      float bit_rate = comp_id->idx->compression_bit_rate;

      //PIDX_get_datatype_details(var->type_name, &values, &bits);
      int ncomps = 0;
      int bits = 0;
      char base_type[10];
      PIDX_decompose_type(var->type_name, base_type, &ncomps, &bits);

      ret = decompress_buffer(comp_id, &patch->buffer, nx, ny, nz, bits/CHAR_BIT, base_type, ncomps, bit_rate);
      if (ret == -1)
        return PIDX_err_compress;
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_compression_finalize(PIDX_comp_id comp_id)
{
  free(comp_id);
  comp_id = 0;
  return PIDX_success;
}
