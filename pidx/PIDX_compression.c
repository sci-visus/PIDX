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
  /// Passed by PIDX API
  MPI_Comm comm;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_variable_index;
  int end_variable_index;
};

PIDX_compression_id PIDX_compression_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  /*
  zfp_params params;
  size_t outsize;
  
  char mode = 0;
  mode = 'r';
  double rate = 16;
  uint precision = 0;
  double tolerance = 0;
  uint minbits = 0;
  uint maxbits = 0;
  uint maxprec = 0;
  int minexp = INT_MIN;
  
  uint type = ZFP_TYPE_DOUBLE;
  /// set array type and size
  params.type = type;
  params.nx = 8;
  params.ny = 8;
  params.nz = 8;
  
  switch (mode) {
    case 'a':
      zfp_set_accuracy(&params, tolerance);
      break;
    case 'p':
      zfp_set_precision(&params, precision);
      break;
    case 'r':
      zfp_set_rate(&params, rate);
      break;
    case 'c':
      params.minbits = minbits;
      params.maxbits = maxbits;
      params.maxprec = maxprec;
      params.minexp = minexp;
      break;
    default:
      fprintf(stderr, "must specify compression parameters via -a, -c, -p, or -r\n");
      return NULL;
  }
  
  outsize = zfp_estimate_compressed_size(&params);
  printf("outsize = %ld\n", outsize);

  size_t count = (size_t)8 * 8 * 8 * 8;
  size_t size = (type == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
  size_t inbytes = count * size;
  size_t bufbytes = 1024 + inbytes;
  void* buffer = malloc(bufbytes);
  FPZ* fpz = fpzip_write_to_buffer(buffer, bufbytes);
  */
  return NULL;
}

#if PIDX_HAVE_MPI
int PIDX_compression_set_communicator(PIDX_compression_id id, MPI_Comm comm)
{
  return -1;
}
#endif

int PIDX_compression_prepare(PIDX_compression_id id, PIDX_variable* variable)
{
  return -1;
}

int PIDX_compression_compress(PIDX_compression_id id, PIDX_variable* variable, int MODE)
{
  return -1;
}
  
int PIDX_compression_buf_destroy(PIDX_compression_id id)
{
  return -1;
}

int PIDX_compression_finalize(PIDX_compression_id id)
{
  return -1;
}