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
  int v = 0;
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    var->chunked_super_patch = malloc(sizeof(*(var->chunked_super_patch)));
    memset(var->chunked_super_patch, 0, sizeof(*(var->chunked_super_patch)));

    PIDX_super_patch out_patch = var->chunked_super_patch;
    PIDX_super_patch in_patch = var->restructured_super_patch;

    out_patch->patch_count = 1;
    out_patch->is_boundary_patch = in_patch->is_boundary_patch;

    out_patch->restructured_patch = malloc(sizeof(*(out_patch->restructured_patch)));
    memset(out_patch->restructured_patch, 0, sizeof(*(out_patch->restructured_patch)));

    memcpy(out_patch->restructured_patch->offset, in_patch->restructured_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
    memcpy(out_patch->restructured_patch->size, in_patch->restructured_patch->size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  }

  return PIDX_success;
}

// one restructured patch per group
PIDX_return_code PIDX_chunk_buf_create(PIDX_chunk_id chunk_id)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int v = 0;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    int bytes_per_value = var->bpv / CHAR_BIT;

    PIDX_super_patch out_patch = var->chunked_super_patch;
    unsigned long long *group_size = out_patch->restructured_patch->size;
    unsigned long long num_elems_group = 1;
    int d;
    for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    {
      if (out_patch->restructured_patch->size[d] % chunk_id->idx->chunk_size[d] == 0)
        num_elems_group *= group_size[d];
      else
        num_elems_group *= (((group_size[d]/chunk_id->idx->chunk_size[d]) + 1) * chunk_id->idx->chunk_size[d]);
    }

    // malloc the storage for all elements in the output array
    // printf("[Chunking] Buffer size of %d = %d\n", chunk_id->idx_c->grank, num_elems_group);
    out_patch->restructured_patch->buffer = malloc(bytes_per_value * num_elems_group);
    memset(out_patch->restructured_patch->buffer, 0, bytes_per_value * num_elems_group);
  }

  return PIDX_success;
}


PIDX_return_code PIDX_chunk(PIDX_chunk_id chunk_id, int MODE)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  unsigned long long v;
  if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      int bytes_per_value = var->bpv / CHAR_BIT;

      PIDX_super_patch out_patch = var->chunked_super_patch;
      PIDX_super_patch in_patch = var->restructured_super_patch;

      //printf("Buffer size %d %d %d - %d %d\n", out_patch->restructured_patch->size[0], out_patch->restructured_patch->size[1], out_patch->restructured_patch->size[2], bytes_per_value, var->vps);

      if (MODE == PIDX_WRITE)
        memcpy(out_patch->restructured_patch->buffer, in_patch->restructured_patch->buffer, out_patch->restructured_patch->size[0] * out_patch->restructured_patch->size[1] * out_patch->restructured_patch->size[2] * bytes_per_value * var->vps);
      else
        memcpy(in_patch->restructured_patch->buffer, out_patch->restructured_patch->buffer, out_patch->restructured_patch->size[0] * out_patch->restructured_patch->size[1] * out_patch->restructured_patch->size[2] * bytes_per_value * var->vps);

      //if (chunk_id->idx_derived->color == 1)
      //{
      //double t1;
      //memcpy(&t1, in_patch->restructured_patch->buffer, sizeof(double));
      //printf("Value %f\n", t1);
      //}
    }

    return PIDX_success;
  }

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

    // copy the size and offset to output
    PIDX_super_patch in_patch = var->restructured_super_patch;
    PIDX_super_patch out_patch = var->chunked_super_patch;

    int nx=(int)out_patch->restructured_patch->size[0];
    if (out_patch->restructured_patch->size[0] % chunk_id->idx->chunk_size[0] != 0)
      nx = ((out_patch->restructured_patch->size[0] / chunk_id->idx->chunk_size[0]) + 1) * chunk_id->idx->chunk_size[0];

    int ny=(int)out_patch->restructured_patch->size[1];
    if (out_patch->restructured_patch->size[1] % chunk_id->idx->chunk_size[1] != 0)
      ny = ((out_patch->restructured_patch->size[1] / chunk_id->idx->chunk_size[1]) + 1) * chunk_id->idx->chunk_size[1];

    int nz=(int)out_patch->restructured_patch->size[2];
    if (out_patch->restructured_patch->size[2] % chunk_id->idx->chunk_size[2] != 0)
      nz = ((out_patch->restructured_patch->size[2] / chunk_id->idx->chunk_size[2]) + 1) * chunk_id->idx->chunk_size[2];

    if (strcmp(var->type_name, FLOAT64) == 0)
    {
    double* p = (double*)in_patch->restructured_patch->buffer;
    int dz = 4 * nx*(ny-(ny/4)*4);
    int dy = 4 * (nx-(nx/4)*4);
    int dx = 4;

    int s1 = 0;
    int z, y, x;
    int zz, yy, xx;
    for (s1 = 0; s1 < var->vps; s1++)
    {
      for (z=0;z<nz;z+=4)
      {
        double* s = (double*)(out_patch->restructured_patch->buffer + ((z/4)*((ny+3)/4)*((nx+3)/4))*(var->bpv/CHAR_BIT)*cbz) + nx*ny*nz*s1;
        for (y=0;y<ny;y+=4)
        {
          for (x=0;x<nx;x+=4)
          {
            unsigned long long diff = (z/4) * (dz + 4 * (ny/4) * (dy+4 * (nx/4) * dx)) + (y/4) * (dy+4*(nx/4)*dx) + (x/4)*dx;
            double* q = p + var->vps * diff + s1;

            for (zz = 0; zz < 4; ++zz)
            {
              for (yy = 0; yy < 4; ++yy)
              {
                for (xx = 0; xx < 4; ++xx)
                {
                  int i = xx + yy * 4 + zz * 4 * 4;
                  int j = xx + yy * nx + zz * nx * ny;

                  if (MODE == PIDX_WRITE)
                  {
                    if (j * var->vps + s1 + var->vps * diff + s1 >= var->restructured_super_patch->restructured_patch->size[0] * var->restructured_super_patch->restructured_patch->size[1] * var->restructured_super_patch->restructured_patch->size[2])
                      continue;

                    s[nx*ny*nz*s1 + i*var->vps+s1] = q[j*var->vps + s1];
                  }
                  else
                  {
                    if (j * var->vps + s1 + var->vps * diff + s1 >= var->restructured_super_patch->restructured_patch->size[0] * var->restructured_super_patch->restructured_patch->size[1] * var->restructured_super_patch->restructured_patch->size[2])
                      continue;

                    q[j*var->vps + s1] = s[nx*ny*nz*s1 + i*var->vps+s1];
                  }
                }
              }
            }
            s += cbz;
          }
        }
      }
    }
    }
    else
    {
        float* p = (float*)in_patch->restructured_patch->buffer;
        int dz = 4 * nx*(ny-(ny/4)*4);
        int dy = 4 * (nx-(nx/4)*4);
        int dx = 4;

        int s1 = 0;
        int z, y, x;
        int zz, yy, xx;
        for (s1 = 0; s1 < var->vps; s1++)
        {
          for (z=0;z<nz;z+=4)
          {
            float* s = (float*)(out_patch->restructured_patch->buffer + ((z/4)*((ny+3)/4)*((nx+3)/4))*(var->bpv/CHAR_BIT)*cbz) + nx*ny*nz*s1;
            for (y=0;y<ny;y+=4)
            {
              for (x=0;x<nx;x+=4)
              {
                unsigned long long diff = (z/4) * (dz + 4 * (ny/4) * (dy+4 * (nx/4) * dx)) + (y/4) * (dy+4*(nx/4)*dx) + (x/4)*dx;
                float* q = p + var->vps * diff + s1;

                for (zz = 0; zz < 4; ++zz)
                {
                  for (yy = 0; yy < 4; ++yy)
                  {
                    for (xx = 0; xx < 4; ++xx)
                    {
                      int i = xx + yy * 4 + zz * 4 * 4;
                      int j = xx + yy * nx + zz * nx * ny;

                      if (MODE == PIDX_WRITE)
                      {
                        if (j * var->vps + s1 + var->vps * diff + s1 >= var->restructured_super_patch->restructured_patch->size[0] * var->restructured_super_patch->restructured_patch->size[1] * var->restructured_super_patch->restructured_patch->size[2])
                          continue;

                        s[nx*ny*nz*s1 + i*var->vps+s1] = q[j*var->vps + s1];
                      }
                      else
                      {
                        if (j * var->vps + s1 + var->vps * diff + s1 >= var->restructured_super_patch->restructured_patch->size[0] * var->restructured_super_patch->restructured_patch->size[1] * var->restructured_super_patch->restructured_patch->size[2])
                          continue;

                        q[j*var->vps + s1] = s[nx*ny*nz*s1 + i*var->vps+s1];
                      }
                    }
                  }
                }
                s += cbz;
              }
            }
          }
        }
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_chunk_meta_data_destroy(PIDX_chunk_id chunk_id)
{
  int var;
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {
    free(var_grp->variable[var]->chunked_super_patch->restructured_patch);
    free(var_grp->variable[var]->chunked_super_patch);
  }

  return PIDX_success;
}

PIDX_return_code PIDX_chunk_buf_destroy(PIDX_chunk_id chunk_id)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int var = 0;
  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
    free(var_grp->variable[var]->chunked_super_patch->restructured_patch->buffer);

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
