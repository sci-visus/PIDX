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
 * \file PIDX_compression.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_compression.h
 *
 */

#include "PIDX_inc.h"

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
  
  int init_index;
  int first_index;
  int last_index;

};

#if PIDX_HAVE_ZFP
int compress_buffer(PIDX_comp_id comp_id, unsigned char* buffer, int length)
{
  int total_size = 0;
  int i = 0;
  int64_t total_chunk_size = comp_id->idx->chunk_size[0] * comp_id->idx->chunk_size[1] * comp_id->idx->chunk_size[2] * comp_id->idx->chunk_size[3] * comp_id->idx->chunk_size[4];

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    size_t typesize, outsize;
    typesize = sizeof(double);
    zfp_params params;
    params.type = ZFP_TYPE_DOUBLE;
    params.nx = comp_id->idx->chunk_size[0];
    params.ny = comp_id->idx->chunk_size[1];
    params.nz = comp_id->idx->chunk_size[2];
    total_size = 0;
    zfp_set_rate(&params, comp_id->idx->compression_bit_rate);
    outsize = zfp_estimate_compressed_size(&params);

    unsigned char* zip = malloc(outsize);
    memset(zip, 0, outsize);

    for ( i = 0; i < length; i = i + total_chunk_size * typesize)
    {
      memset(zip, 0, outsize);
      outsize = zfp_compress(&params, buffer + i , zip, outsize);
      if (outsize == 0)
      {
        fprintf(stderr, "compression failed\n");
        return -1;
      }
      memcpy(buffer + total_size, zip, outsize);
      total_size = total_size + outsize;
    }
    free(zip);
  }
  return total_size;
}

int decompress_buffer(PIDX_comp_id comp_id, void* buffer, int length)
{
  int64_t compression_block_num_elems = comp_id->idx->chunk_size[0] *
                                        comp_id->idx->chunk_size[1] *
                                        comp_id->idx->chunk_size[2];
  unsigned int type_size = sizeof(double); // DUONG_HARDCODE
  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP) // zfp lossy compression
  {
    zfp_params zfp;
    zfp.type = ZFP_TYPE_DOUBLE; // DUONG_HARDCODE
    zfp.nx = comp_id->idx->chunk_size[0];
    zfp.ny = comp_id->idx->chunk_size[1];
    zfp.nz = comp_id->idx->chunk_size[2];
    unsigned int bit_rate = (unsigned int) comp_id->idx->compression_bit_rate;
    zfp_set_rate(&zfp, bit_rate);
    unsigned int input_offset = 0;
    unsigned int output_offset = 0;
    unsigned int input_buffer_size = (unsigned int) length;
    unsigned int output_buffer_size = (length / (bit_rate / 8.0)) * type_size; // DUONG_HARDCODE
    unsigned char* srcPtr = buffer;
    unsigned char* dstPtr = (unsigned char*) malloc(output_buffer_size);
    memset(dstPtr, 0, output_buffer_size);
    unsigned int compressed_block_size = compression_block_num_elems * (bit_rate / 8.0);
    unsigned int uncompressed_block_size = compression_block_num_elems * type_size;
    while (input_offset < input_buffer_size) // decompress one compression block at a time
    {
      int success = zfp_decompress(&zfp, dstPtr + output_offset, srcPtr + input_offset, compressed_block_size);
      if (success == 0) // failure
      {
        free(dstPtr);
        return -1;
      }
      input_offset += compressed_block_size;
      output_offset += uncompressed_block_size;
      assert(output_offset <= output_buffer_size);
    }
    memcpy(buffer, dstPtr, output_buffer_size);
    free(dstPtr);
    return 0;
  }
  return 0;
}

#endif

PIDX_comp_id PIDX_compression_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, int init_index, int start_var_index, int end_var_index)
{
  PIDX_comp_id comp_id;

  comp_id = malloc(sizeof (*comp_id));
  memset(comp_id, 0, sizeof (*comp_id));

  comp_id->idx = idx_meta_data;
  comp_id->idx_derived = idx_derived;

  comp_id->init_index = init_index;
  comp_id->first_index = start_var_index;
  comp_id->last_index = end_var_index;
  
  return comp_id;
}

#if PIDX_HAVE_MPI
int PIDX_compression_set_communicator(PIDX_comp_id comp_id, MPI_Comm comm)
{
  if (comp_id == NULL)
    return PIDX_err_id;

  comp_id->comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_compression(PIDX_comp_id comp_id)
{
  if (comp_id->idx->compression_type == PIDX_NO_COMPRESSION || comp_id->idx->compression_type == PIDX_CHUNKING_ONLY)
    return PIDX_success;

  if (comp_id->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
#if PIDX_HAVE_ZFP
    int v, p, b;
    PIDX_variable var0 = comp_id->idx->variable[comp_id->first_index];

    for (v = comp_id->first_index; v <= comp_id->last_index; v++)
    {
      PIDX_variable var = comp_id->idx->variable[v];
      for (p = 0; p < var0->patch_group_count; p++)
      {
        for (b = 0; b < var0->chunk_patch_group[p]->count; b++)
        {
          Ndim_patch patch = var->chunk_patch_group[p]->patch[b];
          unsigned char* buffer = patch->buffer;
          int element_count = patch->size[0] * patch->size[1] * patch->size[2] * patch->size[3] * patch->size[4];
          int compressed_element_count;

          compressed_element_count = compress_buffer(comp_id, buffer, element_count);
          //printf("Compressed element count = %d\n", compressed_element_count);
          if (compressed_element_count < 0)
            return PIDX_err_compress;

          unsigned char* temp_buffer = realloc(patch->buffer, compressed_element_count);

          if (temp_buffer == NULL)
            return PIDX_err_compress;
          else
            patch->buffer = temp_buffer;
        }
      }
    }
#else
    //printf("Compression Library not found.\n");
    return PIDX_err_compress;
#endif
  }

  return PIDX_success;
}

PIDX_return_code PIDX_decompression(PIDX_comp_id comp_id)
{
  return -1;
}

PIDX_return_code PIDX_compression_finalize(PIDX_comp_id comp_id)
{
  free(comp_id);
  comp_id = 0;
  return PIDX_success;
}
