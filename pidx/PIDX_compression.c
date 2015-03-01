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
struct PIDX_compression_id_struct 
{
#if PIDX_HAVE_MPI
  /// Passed by PIDX API
  MPI_Comm comm;
#endif
  
  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_variable_index;
  int end_variable_index;

#if PIDX_HAVE_LOSSY_ZFP
  zfp_params *params;
#endif
};

PIDX_compression_id PIDX_compression_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  PIDX_compression_id compression_id;

  compression_id = malloc(sizeof (*compression_id));
  memset(compression_id, 0, sizeof (*compression_id));

  compression_id->idx_ptr = idx_meta_data;
  compression_id->idx_derived_ptr = idx_derived_ptr;
  compression_id->start_variable_index = start_var_index;
  compression_id->end_variable_index = end_var_index;

#if PIDX_HAVE_LOSSY_ZFP
  compression_id->params = malloc((end_var_index - start_var_index + 1) * sizeof(*compression_id->params));
  memset(compression_id->params, 0, (end_var_index - start_var_index + 1) * sizeof(*compression_id->params));
#endif
  
  return compression_id;
  
  /*
  size_t count = (size_t)8 * 8 * 8 * 8;
  size_t size = (type == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
  size_t inbytes = count * size;
  size_t bufbytes = 1024 + inbytes;
  void* buffer = malloc(bufbytes);
  FPZ* fpz = fpzip_write_to_buffer(buffer, bufbytes);
  */  
}

#if PIDX_HAVE_MPI
int PIDX_compression_set_communicator(PIDX_compression_id compression_id, MPI_Comm comm)
{
  compression_id->comm = comm;
  //MPI_Comm_dup(comm, &compression_id->comm);
  return 0;
}
#endif

int PIDX_compression_prepare(PIDX_compression_id compression_id)
{
#if PIDX_HAVE_LOSSY_ZFP
  if (compression_id->idx_ptr->compression_type == 1)           ///< Lossy
  {
    int var = 0;
    int dp;
    double rate = 16;
    uint mx, my, mz;
    size_t typesize;
    size_t insize;
    
    for(var = compression_id->start_variable_index; var <= compression_id->end_variable_index; var++)
    {
      if (compression_id->idx_ptr->variable[var]->bits_per_value / 8 == sizeof(float))
      {
        compression_id->params[var].type = ZFP_TYPE_FLOAT;
        dp = 0;
      }
      else if (compression_id->idx_ptr->variable[var]->bits_per_value / 8 == sizeof(double))
      {
        compression_id->params[var].type = ZFP_TYPE_DOUBLE;
        dp = 1;
      }
      
      /// size of floating-point type in bytes
      typesize = dp ? sizeof(double) : sizeof(float);
      
      compression_id->params[var].nx = compression_id->idx_ptr->compression_block_size[0];
      compression_id->params[var].ny = compression_id->idx_ptr->compression_block_size[1];
      compression_id->params[var].nz = compression_id->idx_ptr->compression_block_size[2];
      
      /// set compression mode
#ifdef PIDX_MORE_COMPRESSION_MODES
      double tolerance = 0;
      uint precision = 0;
      uint minbits = 0;
      uint maxbits = 0;
      uint maxprec = 0;
      int minexp = INT_MIN;
      
      zfp_set_accuracy(&compression_id->params[var], tolerance);
      zfp_set_precision(&compression_id->params[var], precision);
      compression_id->params[var].minbits = minbits;
      compression_id->params[var].maxbits = maxbits;
      compression_id->params[var].maxprec = maxprec;
      compression_id->params[var].minexp = minexp;
#endif
      
      zfp_set_rate(&compression_id->params[var], rate);
      
      /// effective array dimensions
      mx = max(compression_id->idx_ptr->compression_block_size[0], 1u);
      my = max(compression_id->idx_ptr->compression_block_size[1], 1u);
      mz = max(compression_id->idx_ptr->compression_block_size[2], 1u);

      /// size of floating-point type in bytes 
      typesize = dp ? sizeof(double) : sizeof(float);

      /// allocate space for uncompressed and compressed fields
      insize = mx * my * mz * typesize;
      compression_id->idx_ptr->variable[var]->lossy_compressed_block_size = zfp_estimate_compressed_size(&compression_id->params[var]);
      //printf("before zip : after zip %d (%d x %d x %d x %d) %d\n", (int)insize, (int)compression_id->idx_ptr->compression_block_size[0], (int)compression_id->idx_ptr->compression_block_size[1], (int)compression_id->idx_ptr->compression_block_size[2], (int)typesize, (int)compression_id->idx_ptr->variable[var]->lossy_compressed_block_size);
    }
  }
  else if (compression_id->idx_ptr->compression_type == 0)      /// Lossless
  {
    
  }
#endif
  return 0;
}

int PIDX_compression_compress(PIDX_compression_id compression_id)
{
  /*
  size_t outsize;
  for(var = compression_id->start_variable_index; var <= compression_id->end_variable_index; var++)
  {
    
    outsize = zfp_compress(&compression_id->params[var], f, zip, outsize);
    if (outsize == 0) 
    {
      fprintf(stderr, "compression failed\n");
      return EXIT_FAILURE;
    }
    rate = CHAR_BIT * (double)outsize / (mx * my * mz);
  }
  */
  return -1;
}
  
int PIDX_compression_buf_destroy(PIDX_compression_id compression_id)
{
  return -1;
}

int PIDX_compression_finalize(PIDX_compression_id compression_id)
{
#if PIDX_HAVE_LOSSY_ZFP
  free(compression_id->params);
  compression_id->params = 0;
#endif
  
  free(compression_id);
  compression_id = 0;
  
  return -1;
}