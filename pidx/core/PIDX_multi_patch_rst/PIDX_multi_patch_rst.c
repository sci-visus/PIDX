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
 * \file PIDX_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"

//Struct for restructuring ID
struct PIDX_multi_patch_rst_struct
{
  //Passed by PIDX API
#if PIDX_HAVE_MPI
  MPI_Comm comm; //Communicator
#endif

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;
  
  int init_index;
  int first_index;
  int last_index;
  
  //int if_perform_rst;

  //dimension of the power-two volume imposed patch
  unsigned long long reg_patch_size[PIDX_MAX_DIMENSIONS];
  int reg_patch_grp_count;
  Ndim_patch_group* reg_patch_grp;
};

#if 0
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);
static int getPowerOftwo(int x);
#endif


#if 0
/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];
  
  return !(check_bit);
}


/// Function to find the power of 2 of an integer value (example 5->8)
static int getPowerOftwo(int x)
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}
#endif


PIDX_multi_patch_rst_id PIDX_multi_patch_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, int first_index, int var_start_index, int var_end_index)
{
  //Creating the restructuring ID
  PIDX_multi_patch_rst_id multi_patch_rst_id;
  multi_patch_rst_id = (PIDX_multi_patch_rst_id)malloc(sizeof (*multi_patch_rst_id));
  memset(multi_patch_rst_id, 0, sizeof (*multi_patch_rst_id));

  multi_patch_rst_id->idx = idx_meta_data;
  multi_patch_rst_id->idx_derived = idx_derived;

  multi_patch_rst_id->init_index = first_index;
  multi_patch_rst_id->first_index = var_start_index;
  multi_patch_rst_id->last_index = var_end_index;

  return (multi_patch_rst_id);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_multi_patch_rst_set_communicator(PIDX_multi_patch_rst_id multi_patch_rst_id, MPI_Comm comm)
{
  if (multi_patch_rst_id == NULL)
    return PIDX_err_id;

  multi_patch_rst_id->comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_multi_patch_rst_meta_data_create(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_buf_create(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_write(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_read(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_multi_patch_rst_buf_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_read(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_multi_patch_rst_aggregate_buf_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_write(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_meta_data_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_multi_patch_rst_finalize(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  free(multi_patch_rst_id);
  multi_patch_rst_id = 0;

  return PIDX_success;
}
