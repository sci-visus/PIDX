/*****************************************************
 **  PIDX Parallel I/O Library            **
 **  Copyright (c) 2010-2014 University of Utah   **
 **  Scientific Computing and Imaging Institute   **
 **  72 S Central Campus Drive, Room 3750       **
 **  Salt Lake City, UT 84112             **
 **                         **
 **  PIDX is licensed under the Creative Commons  **
 **  Attribution-NonCommercial-NoDerivatives 4.0  **
 **  International License. See LICENSE.md.     **
 **                         **
 **  For information about this project see:    **
 **  http://www.cedmav.com/pidx           **
 **  or contact: pascucci@sci.utah.edu        **
 **  For support: PIDX-support@visus.net      **
 **                         **
 *****************************************************/

/**
 * \file PIDX_chunk.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_chunk.h
 *
 */

#include "../../PIDX_inc.h"

///Struct for restructuring ID
struct PIDX_chunk_id_struct
{
  idx_comm idx_c;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, block, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;

  //int if_perform_chunk;

  int group_index;
  int first_index;
  int last_index;
};

PIDX_chunk_id PIDX_chunk_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, int start_var_index, int end_var_index)
{
  PIDX_chunk_id chunk_id;

  //Creating the block rst ID
  chunk_id = (PIDX_chunk_id)malloc(sizeof (*chunk_id));
  memset(chunk_id, 0, sizeof (*chunk_id));

  chunk_id->idx = idx_meta_data;
  chunk_id->idx_derived = idx_derived;
  chunk_id->idx_c = idx_c;

  chunk_id->group_index = 0;

  chunk_id->first_index = start_var_index;
  chunk_id->last_index = end_var_index;

  return chunk_id;
}



PIDX_return_code PIDX_chunk_meta_data_create(PIDX_chunk_id chunk_id)
{
  int v = 0, j = 0;
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];
  if (var0->patch_group_count == 0)
    return PIDX_success;

  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    var->chunk_patch_group = malloc(sizeof(*(var->chunk_patch_group)));
    memset(var->chunk_patch_group, 0, sizeof(*(var->chunk_patch_group)));

    PIDX_super_patch out_patch = var->chunk_patch_group;
    PIDX_super_patch in_patch = var->rst_patch_group;

    out_patch->patch_count = 1;

    out_patch->is_boundary_patch = in_patch->is_boundary_patch;

    out_patch->patch = malloc(sizeof(*(out_patch->patch)) * out_patch->patch_count);
    memset(out_patch->patch, 0, sizeof(*(out_patch->patch)) * out_patch->patch_count);

    out_patch->reg_patch = malloc(sizeof(*(out_patch->reg_patch)));
    memset(out_patch->reg_patch, 0, sizeof(*(out_patch->reg_patch)));

    for(j = 0; j < out_patch->patch_count; j++)
    {
      out_patch->patch[j] = malloc(sizeof(*(out_patch->patch[j])));
      memset(out_patch->patch[j], 0, sizeof(*(out_patch->patch[j])));

      memcpy(out_patch->patch[j]->size, in_patch->reg_patch->size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
      memcpy(out_patch->patch[j]->offset, in_patch->reg_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

    }
    memcpy(out_patch->reg_patch->offset, in_patch->reg_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
    memcpy(out_patch->reg_patch->size, in_patch->reg_patch->size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  }

  return PIDX_success;
}

// one restructured patch per group
// TODO: detect wrong input
PIDX_return_code PIDX_chunk_buf_create(PIDX_chunk_id chunk_id)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->patch_group_count == 0)
    return PIDX_success;

  int v = 0, j = 0;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    int bytes_per_value = var->bpv / 8;


    PIDX_super_patch out_patch = var->chunk_patch_group;
    PIDX_super_patch in_patch = var->rst_patch_group;

    unsigned long long *group_size = in_patch->reg_patch->size;

    for(j = 0; j < out_patch->patch_count; j++)
    {

      //if (chunk_id->idx->compression_type == PIDX_CHUNKING_ONLY || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP_WAVELET)
      //{
      // compute the number of elements in the group
      unsigned long long num_elems_group = 1; // number of elements in the group
      int d;
      for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
      {
        if (out_patch->patch[j]->size[d] % chunk_id->idx->chunk_size[d] == 0)
          num_elems_group *= group_size[d];
        else
          num_elems_group *= (((group_size[d]/chunk_id->idx->chunk_size[d]) + 1) * chunk_id->idx->chunk_size[d]);
      }

      // malloc the storage for all elements in the output array
      out_patch->patch[j]->buffer = malloc(bytes_per_value * num_elems_group);
      memset(out_patch->patch[j]->buffer, 0, bytes_per_value * num_elems_group);
      //}
      //else if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
      //{
      //  out_patch->patch[j]->buffer = malloc(out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * bytes_per_value * var->vps);
      //  memset(out_patch->patch[j]->buffer, 0, out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * bytes_per_value * var->vps);
      //}
    }

  }
  return PIDX_success;
}


PIDX_return_code PIDX_chunk(PIDX_chunk_id chunk_id, int MODE)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->patch_group_count == 0)
    return PIDX_success;

  unsigned long long v, j;
  if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION || chunk_id->idx->compression_type == PIDX_ZFP_COMPRESSION)
  {
    for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      int bytes_per_value = var->bpv / 8;


      PIDX_super_patch out_patch = var->chunk_patch_group;
      PIDX_super_patch in_patch = var->rst_patch_group;

      for(j = 0; j < out_patch->patch_count; j++)
      {
        if (MODE == PIDX_WRITE)
          memcpy(out_patch->patch[j]->buffer, in_patch->reg_patch->buffer, out_patch->reg_patch->size[0] * out_patch->reg_patch->size[1] * out_patch->reg_patch->size[2] * bytes_per_value * var->vps);
        else
          memcpy(in_patch->reg_patch->buffer, out_patch->patch[j]->buffer, out_patch->reg_patch->size[0] * out_patch->reg_patch->size[1] * out_patch->reg_patch->size[2] * bytes_per_value * var->vps);
      }

    }
    return PIDX_success;
  }

  //if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION || chunk_id->idx->compression_type == PIDX_ZFP_COMPRESSION)
  //  return PIDX_success;

  // compute the intra compression block strides
  unsigned long long *chunk_size = chunk_id->idx->chunk_size;
  unsigned long long  cbz = 1;
  int d;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    cbz = cbz * chunk_size[d];

  // loop through all variables
  v = 0;

  for (v = chunk_id->first_index; v <= chunk_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;


    // copy the size and offset to output
    PIDX_super_patch in_patch = var->rst_patch_group;
    PIDX_super_patch out_patch = var->chunk_patch_group;

    int nx=(int)out_patch->patch[0]->size[0];
    if (out_patch->patch[0]->size[0] % chunk_id->idx->chunk_size[0] != 0)
      nx = ((out_patch->patch[0]->size[0] / chunk_id->idx->chunk_size[0]) + 1) * chunk_id->idx->chunk_size[0];

    int ny=(int)out_patch->patch[0]->size[1];
    if (out_patch->patch[0]->size[1] % chunk_id->idx->chunk_size[1] != 0)
      ny = ((out_patch->patch[0]->size[1] / chunk_id->idx->chunk_size[1]) + 1) * chunk_id->idx->chunk_size[1];

    int nz=(int)out_patch->patch[0]->size[2];
    if (out_patch->patch[0]->size[2] % chunk_id->idx->chunk_size[2] != 0)
      nz = ((out_patch->patch[0]->size[2] / chunk_id->idx->chunk_size[2]) + 1) * chunk_id->idx->chunk_size[2];

    unsigned char* temp_buffer = malloc(nx * ny * nz * var->bpv/8 * var->vps);
    if (temp_buffer == NULL)
      return PIDX_err_chunk;

    if (MODE == PIDX_WRITE)
      memcpy(temp_buffer, in_patch->reg_patch->buffer, nx * ny * nz * var->bpv/8 * var->vps);

    /*
    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
    for (r = 0; r < var->rst_patch_group->count; r++)
    {
    for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
    {
      for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
      {
      for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
      {
        index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
        send_o = index * var->vps * (var->bpv/8);
        send_c = (patch_group->patch[r]->size[0]);
        recv_o = ((nx) * (ny) * (k1 - out_patch->patch[0]->offset[2])) + ((nx)* (j1 - out_patch->patch[0]->offset[1])) + (i1 - out_patch->patch[0]->offset[0]);
        if (MODE == PIDX_WRITE)
        memcpy(temp_buffer + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group->patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
      }
      }
    }
    }
    */

    double* p=(double*)temp_buffer;
    int dz=4*nx*(ny-(ny/4)*4);
    int dy=4*(nx-(nx/4)*4);
    int dx=4;

    int s1 = 0;
    int z, y, x;
    int zz, yy, xx;
    for (s1 = 0; s1 < var->vps; s1++)
    {
      for (z=0;z<nz;z+=4)
      {
        double* s = (double*)(out_patch->patch[0]->buffer+((z/4)*((ny+3)/4)*((nx+3)/4))*(var->bpv/8)*cbz) + nx*ny*nz*s1;
        for (y=0;y<ny;y+=4)
        {
          for (x=0;x<nx;x+=4)
          {
            unsigned long long diff=(z/4)*(dz+4*(ny/4)*(dy+4*(nx/4)*dx))+(y/4)*(dy+4*(nx/4)*dx)+(x/4)*dx;
            double* q=p+var->vps*diff+s1;

            for (zz = 0; zz < 4; ++zz)
            {
              for (yy = 0; yy < 4; ++yy)
              {
                for (xx = 0; xx < 4; ++xx)
                {
                  int i = xx + yy * 4 + zz * 4 * 4;
                  int j = xx + yy * nx + zz * nx * ny;

                  if (MODE == PIDX_WRITE)
                    s[nx*ny*nz*s1 + i*var->vps+s1] = q[j*var->vps + s1];
                  else
                    q[j*var->vps + s1] = s[nx*ny*nz*s1 + i*var->vps+s1];
                }
              }
            }
            s += cbz;
          }
        }
      }
    }

    if (MODE == PIDX_READ)
      memcpy(in_patch->reg_patch->buffer, temp_buffer, nx * ny * nz * var->bpv/8 * var->vps);

#if 0
    for (r = 0; r < var->rst_patch_group->count; r++)
    {
      for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
      {
        for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
        {
          for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
          {
            index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
            send_o = index * var->vps * (var->bpv/8);
            send_c = (patch_group->patch[r]->size[0]);
            recv_o = ((nx) * (ny) * (k1 - out_patch->patch[0]->offset[2])) + ((nx)* (j1 - out_patch->patch[0]->offset[1])) + (i1 - out_patch->patch[0]->offset[0]);
            if (MODE == PIDX_READ)
              memcpy(var->rst_patch_group->patch[r]->buffer + send_o, temp_buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
          }
        }
      }
    }
#endif

    free(temp_buffer);

  }

  return PIDX_success;
}



PIDX_return_code PIDX_chunk_meta_data_destroy(PIDX_chunk_id chunk_id)
{
  int j, var;
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->patch_group_count == 0)
    return PIDX_success;

  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {

    for(j = 0; j < var_grp->variable[var]->chunk_patch_group->patch_count; j++)
    {
      free(var_grp->variable[var]->chunk_patch_group->patch[j]);
      var_grp->variable[var]->chunk_patch_group->patch[j] = 0;
    }

    free(var_grp->variable[var]->chunk_patch_group->patch);
    var_grp->variable[var]->chunk_patch_group->patch = 0;

    free(var_grp->variable[var]->chunk_patch_group->reg_patch);
    var_grp->variable[var]->chunk_patch_group->reg_patch = 0;

    free(var_grp->variable[var]->chunk_patch_group);
    var_grp->variable[var]->chunk_patch_group = 0;
  }

  return PIDX_success;
}

PIDX_return_code PIDX_chunk_buf_destroy(PIDX_chunk_id chunk_id)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->patch_group_count == 0)
    return PIDX_success;

  int j = 0, var = 0;
  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {

    for(j = 0; j < var_grp->variable[var]->chunk_patch_group->patch_count; j++)
    {
      free(var_grp->variable[var]->chunk_patch_group->patch[j]->buffer);
      var_grp->variable[var]->chunk_patch_group->patch[j]->buffer = 0;
    }

  }

  return PIDX_success;
}


PIDX_return_code PIDX_chunk_finalize(PIDX_chunk_id chunk_id)
{
  free(chunk_id);
  chunk_id = 0;

  return PIDX_success;
}


int HELPER_chunking(PIDX_chunk_id id)
{
  return PIDX_success;
}
