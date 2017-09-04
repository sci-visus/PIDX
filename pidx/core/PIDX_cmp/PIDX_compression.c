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
#if PIDX_HAVE_MPI
  /// Passed by PIDX API
  MPI_Comm comm;
#endif

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
#if PIDX_HAVE_ZFP
  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP || comp_id->idx->compression_type == PIDX_CHUNKING_ZFP_63_COEFFICIENT)
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
    //printf("BM %d\n", bytes_max);
    unsigned char* output = malloc(bytes_max);
    bitstream* stream = stream_open(output, bytes_max);
    zfp_stream_set_bit_stream(zfp, stream);
    size_t i = 0;
    size_t length = (size_t)nx * (size_t)ny * (size_t)nz;
    for (i = 0; i < length * bytes_per_sample; i += chunk_bytes)
    {
      size_t bits = 0;
      if (type == zfp_type_float)
        bits = zfp_encode_block2_float_3(zfp, (float*)(buffer + i));
      else if (type == zfp_type_double)
        bits = zfp_encode_block2_double_3(zfp, (double*)(buffer + i));
      memcpy(buffer + total_bytes, output, bits / CHAR_BIT);
      assert(bits % CHAR_BIT == 0);
      total_bytes += bits / CHAR_BIT;
    }
    free(output);
    zfp_stream_close(zfp);
    stream_close(stream);
    zfp_field_free(field);
  }
#endif
  return total_bytes;
}


int decompress_buffer(PIDX_comp_id comp_id, unsigned char* buffer, int nx, int ny, int nz,
                      int bytes_per_sample, float bit_rate)
{
#if 0
   int i = 0;
   unsigned long long total_chunk_size = comp_id->idx->chunk_size[0] * comp_id->idx->chunk_size[1] * comp_id->idx->chunk_size[2];

   if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
   {
     size_t typesize;
     typesize = bytes_per_sample;
     zfp_params params;
     if (bytes_per_sample == sizeof(double))
       params.type = ZFP_TYPE_DOUBLE;
     else if (bytes_per_sample == sizeof(float))
        params.type = ZFP_TYPE_FLOAT;

     params.nx = comp_id->idx->chunk_size[0];
     params.ny = comp_id->idx->chunk_size[1];
     params.nz = comp_id->idx->chunk_size[2];
     zfp_set_rate(&params, bit_rate);

     unsigned char* temp_buffer = malloc(length * bytes_per_sample);
     memset(temp_buffer, 0,length * bytes_per_sample);

     for ( i = 0; i < length * typesize; i = i + total_chunk_size * typesize)
     {
       int success = zfp_decompress(&params, temp_buffer + i , buffer + (i / ((typesize * 8) / bit_rate)), total_chunk_size * (bit_rate / 8));
       if (success == 0) // failure
       {
         free(temp_buffer);
         return -1;
       }
     }

     unsigned char* temp_buffer1 = realloc(buffer, length * bytes_per_sample);
     if (temp_buffer1 == NULL)
       return PIDX_err_compress;
     else
       memcpy(buffer, temp_buffer, length * bytes_per_sample);

  }
#endif

   size_t total_bytes = 0;
 #if PIDX_HAVE_ZFP
   if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP || comp_id->idx->compression_type == PIDX_CHUNKING_ZFP_63_COEFFICIENT)
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
     //printf("BM %d\n", bytes_max);
     unsigned char* output = malloc(bytes_max);
     bitstream* stream = stream_open(output, bytes_max);
     zfp_stream_set_bit_stream(zfp, stream);
     size_t i = 0;
     size_t length = (size_t)nx * (size_t)ny * (size_t)nz;
     for (i = 0; i < length * bytes_per_sample; i += chunk_bytes)
     {
       size_t bits = 0;
       if (type == zfp_type_float)
         bits = zfp_encode_block2_float_3(zfp, (float*)(buffer + i));
       else if (type == zfp_type_double)
         bits = zfp_encode_block2_double_3(zfp, (double*)(buffer + i));
       memcpy(buffer + total_bytes, output, bits / CHAR_BIT);
       assert(bits % CHAR_BIT == 0);
       total_bytes += bits / CHAR_BIT;
     }
     free(output);
     zfp_stream_close(zfp);
     stream_close(stream);
     zfp_field_free(field);
   }
 #endif
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
  if (comp_id->idx->compression_type == PIDX_NO_COMPRESSION ||
      comp_id->idx->compression_type == PIDX_CHUNKING_ONLY)
  {
    return PIDX_success;
  }

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP || comp_id->idx->compression_type == PIDX_CHUNKING_ZFP_63_COEFFICIENT)
  {

    int v, p, b;
    PIDX_variable_group var_grp = comp_id->idx->variable_grp[0];
    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        for (b = 0; b < var->chunk_patch_group[p]->count; b++)
        {
          Ndim_patch patch = var->chunk_patch_group[p]->patch[b];
          unsigned char* buffer = patch->buffer;
          int nx = patch->size[0];
          int ny = patch->size[1];
          int nz = patch->size[2];
          float bit_rate = comp_id->idx->compression_bit_rate;
          int compressed_bytes = compress_buffer(comp_id, buffer, nx, ny, nz, var->bpv/8, bit_rate);
          //fprintf(stderr, "%d %d %d %d %f CMP %d\n", nx, ny, nz, var->bpv/8, bit_rate, compressed_bytes);
          unsigned char* temp_buffer = realloc(patch->buffer, compressed_bytes);
          if (temp_buffer == NULL)
            return PIDX_err_compress;
          else
            patch->buffer = temp_buffer;
        }
      }
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_decompression(PIDX_comp_id comp_id)
{
#if 1
  if (comp_id->idx->compression_type == PIDX_NO_COMPRESSION || comp_id->idx->compression_type == PIDX_CHUNKING_ONLY)
    return PIDX_success;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
#if PIDX_HAVE_ZFP
    int v, p, b, ret = 0;

    PIDX_variable_group var_grp = comp_id->idx->variable_grp[0];
    PIDX_variable var0 = var_grp->variable[comp_id->first_index];

    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var0->patch_group_count; p++)
      {
        for (b = 0; b < var0->chunk_patch_group[p]->count; b++)
        {
          Ndim_patch patch = var->chunk_patch_group[p]->patch[b];
          unsigned char* buffer = patch->buffer;
          int nx = patch->size[0];
          int ny = patch->size[1];
          int nz = patch->size[2];
          float bit_rate = comp_id->idx->compression_bit_rate;

          //if (rank == 0)
          //fprintf(stderr, "Before [%d] element count %d byte size %d bit rate %d\n", rank, element_count*var->bpv/8, var->bpv/8, comp_id->idx->compression_bit_rate);

          ret = decompress_buffer(comp_id, buffer, nx, ny, nz, var->bpv/8, bit_rate);
          if (ret == -1)
            return PIDX_err_compress;

          //if (rank == 0)
          //fprintf(stderr, "After [%d] Compressed element count = %d\n", rank, compressed_element_count);

          /*
          if (compressed_element_count <= 0)
            return PIDX_err_compress;

          unsigned char* temp_buffer = realloc(patch->buffer, compressed_element_count);
          if (temp_buffer == NULL)
            return PIDX_err_compress;
          else
            patch->buffer = temp_buffer;
          */
        }
      }
    }
#else
    //fprintf(stderr, "Compression Library not found.\n");
    return PIDX_err_compress;
#endif
  }
#endif
  return PIDX_success;
}

PIDX_return_code PIDX_compression_finalize(PIDX_comp_id comp_id)
{
  free(comp_id);
  comp_id = 0;
  return PIDX_success;
}
