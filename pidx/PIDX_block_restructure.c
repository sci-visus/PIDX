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

// one restructured box per group
// TODO: detect wrong input
int PIDX_block_rst_prepare(PIDX_block_rst_id block_rst_id)
{
  // compute the intra compression block strides
  int64_t *compression_block_size = block_rst_id->idx_ptr->compression_block_size;
  int64_t intra_cblock_stride[PIDX_MAX_DIMENSIONS]; // stride inside a compression block
  intra_cblock_stride[0] = 1;
  int d = 0;
  for (d = 1; d < PIDX_MAX_DIMENSIONS; ++d)
  {
    intra_cblock_stride[d] = intra_cblock_stride[d - 1] * compression_block_size[d - 1];
  }

  // loop through all variables
  int v = 0;
  for (v = block_rst_id->start_variable_index; v <= block_rst_id->end_variable_index; ++v)
  {
    PIDX_variable var = block_rst_id->idx_ptr->variable[v];
    int bytes_per_value = var->bits_per_value / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      Ndim_box_group box_group = var->patch_group_ptr[g];
      Ndim_box_group out_box = var->post_rst_block[g];
      int64_t *group_size = box_group->enclosing_box_size;
      memcpy(&out_box->box[0]->Ndim_box_size, &box_group->enclosing_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
      memcpy(&out_box->box[0]->Ndim_box_offset, &box_group->enclosing_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

      // compute the strides of the group
      int64_t group_stride[PIDX_MAX_DIMENSIONS]; // stride inside a group
      group_stride[0] = 1;
      for (d = 1; d < PIDX_MAX_DIMENSIONS; ++d)
      {
        group_stride[d] = group_stride[d - 1] * group_size[d - 1];
      }

      // compute the number of elements in the group
      int64_t num_elems_group = 1; // number of elements in the group
      for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
      {
        if (group_size[d] > 0)
        {
          num_elems_group *= group_size[d];
        }
      }

      // malloc the storage for all elements in the output array
      out_box->box[0]->Ndim_box_buffer = malloc(bytes_per_value * num_elems_group);

      int64_t *group_offset = box_group->enclosing_box_offset;
      // loop through all boxes
      int b = 0;
      for (b = 0; b < box_group->box_count; ++b)
      {
        Ndim_box box = box_group->box[b];
        int64_t *box_size = box->Ndim_box_size;
        int64_t *box_offset = box->Ndim_box_offset; // global offset of the box

        // compute the number of elements in the box, and the number of elements before the box in the same group
        int64_t num_elems_box = 1; // number of elements in the box
        int64_t num_prev_elems = 0; // number of elements that precede this box in the group
        for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
        {
          if (box_size[d] > 0)
          {
            int64_t local_offset = box_offset[d] - group_offset[d];
            num_prev_elems += local_offset * group_stride[d];
            num_elems_box *= box_size[d];
          }
        }

        // compute strides of the box in all dimensions
        // stride[i] = the number of elements between two consecutive indices in the ith dimension
        int64_t box_stride[PIDX_MAX_DIMENSIONS];
        box_stride[0] = 1; // assume x (or dimension [0]) is the fastest varying dimension
        for (d = 1; d < PIDX_MAX_DIMENSIONS; ++d)
        {
          box_stride[d] = box_stride[d - 1] * box_size[d - 1];
        }

        // loop through the elements to find their new positions
        int64_t i = 0;
        for (i = 0; i < num_elems_box; i += compression_block_size[0]) // we copy k consecutive bytes at a time
        {
          // compute the 5D index of the ith element
          int64_t j = i + num_prev_elems;
          int64_t index[PIDX_MAX_DIMENSIONS];
          // the following calculation is based on the formula
          // linear index = index[0] * stride[0] + index[1] * stride[1] + ... + index[n] * stride[n] (n == PIDX_MAX_DIMENSIONS)
          // in which stride[1] is a multiple of stride[0], stride[2] is a multiple of stride[1] and so on
          // and index[i] * stride[i] < stride[i + 1]
          for (d = 0; d + 1 < PIDX_MAX_DIMENSIONS; ++d)
          {
            if (j > 0 && box_stride[d + 1] > 0 && box_stride[d] > 0)
            {
              int64_t k = j % box_stride[d + 1];
              index[d] = k / box_stride[d];
              j -= k;
            }
            else // if (j == 0)
            {
              index[d] = 0;
            }
          }
          if (box_stride[PIDX_MAX_DIMENSIONS - 1] > 0 && j >= 0)
          {
            index[PIDX_MAX_DIMENSIONS - 1] = j / box_stride[PIDX_MAX_DIMENSIONS - 1];
          }
          else
          {
            index[PIDX_MAX_DIMENSIONS - 1] = 0;
          }

          // compute the output linear index
          j = 0; // compression block-level linear index
          box_stride[0] = 1; // re-use this, but now means the stride at compression block level
          for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
          {
            int64_t cbz = compression_block_size[d];
            if (cbz > 0)
            {
              if (d > 0)
              { // reduce the stride based on the new blocking scheme
                box_stride[d] = box_stride[d - 1] * (group_size[d - 1] / compression_block_size[d - 1]);
              }
              if (index[d] > 0)
              {
                j += (index[d] / cbz) * box_stride[d];
                j += (index[d] % cbz) * intra_cblock_stride[d];
              }
            }
          }

          // copy k consecutive elements at once
          memcpy(&out_box->box[0]->Ndim_box_buffer[j * bytes_per_value], &box->Ndim_box_buffer[i * bytes_per_value], bytes_per_value * compression_block_size[0]);
        }
      }
    }
  }
  return 0;
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