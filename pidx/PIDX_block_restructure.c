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
 * \file PIDX_block_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_block_rst.h
 *
 */

#include "PIDX_inc.h"


///Struct for restructuring ID
struct PIDX_block_rst_id_struct 
{
#if PIDX_HAVE_MPI
  /// Passed by PIDX API
  MPI_Comm comm;
#endif
  
  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, block, file name template and more
  idx_dataset idx_ptr;
  
  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_variable_index;
  int end_variable_index;
};

PIDX_block_rst_id PIDX_block_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  PIDX_block_rst_id block_rst_id;

  //Creating the block rst ID
  block_rst_id = (PIDX_block_rst_id)malloc(sizeof (*block_rst_id)); 
  memset(block_rst_id, 0, sizeof (*block_rst_id));
  
  block_rst_id->idx_ptr = idx_meta_data;
  block_rst_id->idx_derived_ptr = idx_derived_ptr;

  block_rst_id->start_variable_index = start_var_index;
  block_rst_id->end_variable_index = end_var_index;
      
  return block_rst_id;
}

#if PIDX_HAVE_MPI
int PIDX_block_rst_set_communicator(PIDX_block_rst_id block_rst_id, MPI_Comm comm)
{
  block_rst_id->comm = comm;
  //MPI_Comm_dup(comm, &block_rst_id->comm);
  return 0;
}
#endif

int PIDX_block_rst_prepare(PIDX_block_rst_id block_rst_id)
{
  int var = 0, y = 0, b = 0;

  //printf("var %d %d\n", block_rst_id->start_variable_index, block_rst_id->end_variable_index);
  for (var = block_rst_id->start_variable_index; var <= block_rst_id->end_variable_index; var++)
  {
    //printf("y = %d\n", block_rst_id->idx_ptr->variable[block_rst_id->start_variable_index]->patch_group_count);
    for (y = 0; y < block_rst_id->idx_ptr->variable[block_rst_id->start_variable_index]->patch_group_count; y++)
    {
      for (b = 0; b < block_rst_id->idx_ptr->variable[block_rst_id->start_variable_index]->patch_group_ptr[y]->box_count; b++) 
      {
#if 0
        for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
        {
          compressed_patch_offset[l] = id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[l] / id->idx_ptr->compression_block_size[l];
          compressed_patch_size[l] = id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_size[l] / id->idx_ptr->compression_block_size[l];  
          total_compressed_patch_size = total_compressed_patch_size * compressed_patch_size[l];
        }
#endif
      }
    }
  }
  return -1;
}

int PIDX_block_rst_compress(PIDX_block_rst_id block_rst_id)
{
  return -1;
}

int PIDX_block_rst_buf_destroy(PIDX_block_rst_id block_rst_id)
{
  return -1;
}

int PIDX_block_rst_finalize(PIDX_block_rst_id block_rst_id)
{
  return -1;
}