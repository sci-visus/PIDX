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


///Struct for restructuring ID
struct PIDX_comp_id_struct
{
  /// Passed by PIDX API
  MPI_Comm comm;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;


  idx_comm idx_c;

  int first_index;
  int last_index;

};


int compress_buffer(PIDX_comp_id comp_id, unsigned char* buffer,
                    int nx, int ny, int nz, int bytes_per_sample, float bit_rate)
{
  size_t total_bytes = 0;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    unsigned long long* chunk_dim = comp_id->idx->chunk_size;
    assert(chunk_dim[0] == 4 && chunk_dim[1] == 4 && chunk_dim[2] == 4);
    size_t total_chunk_dim = (size_t)chunk_dim[0] * (size_t)chunk_dim[1] * (size_t)chunk_dim[2];
    size_t chunk_bytes = total_chunk_dim * bytes_per_sample;
    zfp_type type = (bytes_per_sample == 4) ? zfp_type_float : zfp_type_double;
    zfp_field* field = zfp_field_3d(NULL, type, nx, ny, nz);
    zfp_stream* zfp = zfp_stream_open(NULL);
    zfp_stream_set_rate(zfp, bit_rate, type, 3, 0);
    size_t bytes_max = zfp_stream_maximum_size(zfp, field);

    unsigned char* output = malloc(bytes_max);
    bitstream* stream = stream_open(output, bytes_max);
    zfp_stream_set_bit_stream(zfp, stream);
    size_t i = 0;
    size_t length = (size_t)nx * (size_t)ny * (size_t)nz;
    for (i = 0; i < length * bytes_per_sample; i += chunk_bytes)
    {
      size_t bits = 0;
      if (type == zfp_type_float)
        bits = zfp_encode_block_float_3(zfp, (float*)(buffer + i));
      else if (type == zfp_type_double)
        bits = zfp_encode_block_double_3(zfp, (double*)(buffer + i));

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


int decompress_buffer(PIDX_comp_id comp_id, unsigned char* buffer, int nx, int ny, int nz,
                      int bytes_per_sample, float bit_rate)
{
   size_t total_bytes = 0;

   if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
   {
     unsigned long long* chunk_dim = comp_id->idx->chunk_size;
     assert(chunk_dim[0] == 4 && chunk_dim[1] == 4 && chunk_dim[2] == 4);
     size_t total_chunk_dim = (size_t)chunk_dim[0] * (size_t)chunk_dim[1] * (size_t)chunk_dim[2];
     size_t chunk_bytes = total_chunk_dim * bytes_per_sample;
     zfp_type type = (bytes_per_sample == 4) ? zfp_type_float : zfp_type_double;
     zfp_field* field = zfp_field_3d(NULL, type, nx, ny, nz);
     zfp_stream* zfp = zfp_stream_open(NULL);
     zfp_stream_set_rate(zfp, bit_rate, type, 3, 0);

     unsigned char* temp_buffer = malloc(nx * ny * nz * bytes_per_sample);
     bitstream* stream = stream_open(buffer, (nx * ny * nz * bytes_per_sample) / comp_id->idx->compression_factor);
     zfp_stream_set_bit_stream(zfp, stream);
     size_t i = 0;
     size_t length = (size_t)nx * (size_t)ny * (size_t)nz;

     for (i = 0; i < length * bytes_per_sample; i += chunk_bytes)
     {
       size_t bits = 0;
       if (type == zfp_type_float)
         bits = zfp_decode_block_float_3(zfp, (float*)(temp_buffer + i));
       else if (type == zfp_type_double)
         bits = zfp_decode_block_double_3(zfp, (double*)(temp_buffer + i));

       assert(bits % CHAR_BIT == 0);
       total_bytes += bits / CHAR_BIT;
     }

     unsigned char *temp_buffer2 = realloc(buffer, nx * ny * nz * bytes_per_sample);
     if (temp_buffer2 == NULL)
     {
       fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
       return PIDX_err_rst;
     }
     else
       buffer = temp_buffer2;

     memcpy(buffer, temp_buffer, nx * ny * nz * bytes_per_sample);

     free(temp_buffer);
     zfp_stream_close(zfp);
     stream_close(stream);
     zfp_field_free(field);
   }

   return total_bytes;
}

PIDX_comp_id PIDX_compression_init(idx_dataset idx_meta_data,
                                   idx_dataset_derived_metadata idx_derived,
                                   idx_comm idx_c, int start_var_index, int end_var_index)
{
  PIDX_comp_id comp_id;

  comp_id = malloc(sizeof (*comp_id));
  memset(comp_id, 0, sizeof (*comp_id));

  comp_id->idx = idx_meta_data;
  comp_id->idx_derived = idx_derived;
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
    PIDX_variable_group var_grp = comp_id->idx->variable_grp[0];
    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      PIDX_patch patch = var->chunked_super_patch->restructured_patch;
      unsigned char* buffer = patch->buffer;
      int nx = patch->size[0];
      int ny = patch->size[1];
      int nz = patch->size[2];
      float bit_rate = comp_id->idx->compression_bit_rate;
      int compressed_bytes = compress_buffer(comp_id, buffer, nx, ny, nz, var->bpv/8, bit_rate);
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
    PIDX_variable_group var_grp = comp_id->idx->variable_grp[0];
    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];

      PIDX_patch patch = var->chunked_super_patch->restructured_patch;
      unsigned char* buffer = patch->buffer;
      int nx = patch->size[0];
      int ny = patch->size[1];
      int nz = patch->size[2];
      float bit_rate = comp_id->idx->compression_bit_rate;

      ret = decompress_buffer(comp_id, buffer, nx, ny, nz, var->bpv/8, bit_rate);
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
