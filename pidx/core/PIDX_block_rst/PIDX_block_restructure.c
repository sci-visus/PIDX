/*
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */

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

//Struct for restructuring ID
struct PIDX_chunk_id_struct
{
  idx_comm idx_c;

  // Contains all relevant IDX file info
  // Blocks per file, samples per block, bitmask, block, file name template and more
  idx_dataset idx;

  // Contains all derieved IDX file info
  // number of files, files that are ging to be populated
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
    int bytes_per_value = (var->bpv / CHAR_BIT) * var->vps;

    PIDX_super_patch out_patch = var->chunked_super_patch;
    size_t *group_size = out_patch->restructured_patch->size;
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
    //printf("\n \n Out Patch is before chunking: \n");
    //printf("%s\n", out_patch->restructured_patch->buffer);
  }

  return PIDX_success;
}


void returnbuffer(unsigned char * q, unsigned char * s, int cbz, int nx, int ny, int nz, int compval, int bits, int mode)
{
  int dz = 4 * nx*(ny-(ny/4)*4);
  int dy = 4 * (nx-(nx/4)*4);
  int dx = 4;

  unsigned char * temp = q;
  unsigned char * temp2 = s;
  int x,y,z, zz,yy,xx;
  for (z=0;z<nz;z+=4)
  {
    s = temp2 + ((z/4) * ((ny+3)/4) * ((nx+3)/4)) * cbz * (bits/CHAR_BIT);
    for (y=0;y<ny;y+=4)
    {
      for (x=0;x<nx;x+=4)
      {
        unsigned long long diff = (z/4) * (dz + 4 * (ny/4) * (dy+4 * (nx/4) * dx)) + (y/4) * (dy+4*(nx/4)*dx) + (x/4)*dx;
        q = temp + diff * (bits/CHAR_BIT);
        for (zz = 0; zz < 4; zz++)
        {
          for (yy = 0; yy < 4; yy++)
          {
            for (xx = 0; xx < 4; xx++)
            {
              int i = xx + yy * 4 + zz * 4 * 4;
              int j = xx + yy * nx + zz * nx * ny;

              if (j + diff >= compval)
                continue;
              if (mode == PIDX_WRITE)
                memcpy(s + i * (bits/CHAR_BIT), q + j * (bits/CHAR_BIT), (bits/CHAR_BIT));
              else
                memcpy(q + j * (bits/CHAR_BIT), s + i * (bits/CHAR_BIT), (bits/CHAR_BIT));
            }
          }
        }
        s += cbz * (bits/CHAR_BIT);
      }
    }
  }
}



PIDX_return_code PIDX_chunk(PIDX_chunk_id chunk_id, int MODE)
{
  PIDX_variable_group var_grp = chunk_id->idx->variable_grp[chunk_id->group_index];
  PIDX_variable var0 = var_grp->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int i = 0, j = 0, d = 0;
  unsigned long long v;
  if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      int bytes_per_value = (var->bpv / CHAR_BIT) * var->vps;

      PIDX_super_patch out_patch = var->chunked_super_patch;
      PIDX_super_patch in_patch = var->restructured_super_patch;

      if (MODE == PIDX_WRITE)
        memcpy(out_patch->restructured_patch->buffer, in_patch->restructured_patch->buffer, out_patch->restructured_patch->size[0] * out_patch->restructured_patch->size[1] * out_patch->restructured_patch->size[2] * bytes_per_value);
      else
        memcpy(in_patch->restructured_patch->buffer, out_patch->restructured_patch->buffer, out_patch->restructured_patch->size[0] * out_patch->restructured_patch->size[1] * out_patch->restructured_patch->size[2] * bytes_per_value);

    }

    return PIDX_success;
  }

  // compute the intra compression block strides
  size_t *chunk_size = chunk_id->idx->chunk_size;
  //printf("Chunk size is %lld\n", *chunk_size);
  size_t  cbz = 1;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d){
    cbz = cbz * chunk_size[d];
    //printf("some terms are %lld, %lld, %d\n", cbz, chunk_size[d], PIDX_MAX_DIMENSIONS);
  }

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

    //printf("value of nx is %d\n", nx);

    int ny=(int)out_patch->restructured_patch->size[1];
    if (out_patch->restructured_patch->size[1] % chunk_id->idx->chunk_size[1] != 0)
      ny = ((out_patch->restructured_patch->size[1] / chunk_id->idx->chunk_size[1]) + 1) * chunk_id->idx->chunk_size[1];
    //printf("value of ny is %d\n", ny);
    int nz=(int)out_patch->restructured_patch->size[2];
    if (out_patch->restructured_patch->size[2] % chunk_id->idx->chunk_size[2] != 0)
      nz = ((out_patch->restructured_patch->size[2] / chunk_id->idx->chunk_size[2]) + 1) * chunk_id->idx->chunk_size[2];

    int values = 0;
    int bits = 0;
    if (strcmp(var->type_name, FLOAT32) == 0)
    {
      values = 1;
      bits = 32;
    }
    else if (strcmp(var->type_name, FLOAT32_GA) == 0)
    {
      values = 2;
      bits = 32;
    }
    else if (strcmp(var->type_name, FLOAT32_RGB) == 0)
    {
      values = 3;
      bits = 32;
    }
    else if (strcmp(var->type_name, FLOAT32_RGBA) == 0)
    {
      values = 4;
      bits = 32;
    }

    else if (strcmp(var->type_name, FLOAT64) == 0)
    {
      values = 1;
      bits = 64;
    }
    else if (strcmp(var->type_name, FLOAT64_GA) == 0)
    {
      values = 2;
      bits = 64;
    }
    else if (strcmp(var->type_name, FLOAT64_RGB) == 0)
    {
      values = 3;
      bits = 64;
    }
    else if (strcmp(var->type_name, FLOAT64_RGBA) == 0)
    {
      values = 4;
      bits = 64;
    }
    else if (strcmp(var->type_name, FLOAT64_7STENCIL) == 0)
    {
      values = 7;
      bits = 64;
    }

    int dim = values;
    int sz = dim*nx*ny*nz;

    unsigned char * s = out_patch->restructured_patch->buffer;
    unsigned char * q = in_patch->restructured_patch->buffer;

    unsigned char ** sdim = malloc(dim*(bits/CHAR_BIT));
    for (i = 0; i < dim; i++)
      sdim[i] = malloc((sz/dim) * (bits/CHAR_BIT));

    unsigned char ** op = malloc(dim*(bits/CHAR_BIT));
    for (i = 0; i < dim; i++)
      op[i] = malloc((sz/dim) * (bits/CHAR_BIT));

    unsigned char ** qp = malloc(dim*(bits/CHAR_BIT));
    for (i = 0; i < dim; i++)
      qp[i] = malloc((sz/dim) * (bits/CHAR_BIT));

    int compval = in_patch->restructured_patch->size[0] * in_patch->restructured_patch->size[1] * in_patch->restructured_patch->size[2];
    if (MODE == PIDX_WRITE)
    {
      int k = 0;
      for (i = 0; i < dim; i++)
      {
        k = i;
        for (j = 0; j < (sz/dim); j++)
        {
          if(k < sz)
          {
            memcpy(sdim[i] + j*(bits/CHAR_BIT), in_patch->restructured_patch->buffer + k* (bits/CHAR_BIT), (bits/CHAR_BIT));
            k+=dim;
          }
        }
        k = 0;
      }
      int d = 0;
      for (d = 0; d < dim; d++)
        returnbuffer(sdim[d], op[d], cbz, nx, ny, nz, compval, bits, PIDX_WRITE);

      int ctr = 0, row = 0;
      for (i = 0; i < sz; i++)
      {
        if(ctr == (sz/dim) ){
          row+=1;
          ctr = 0;
        }
        memcpy(s + i*(bits/CHAR_BIT), op[row] + ctr*(bits/CHAR_BIT), (bits/CHAR_BIT));
        ctr++;
      }
    }
    else
    {
      int ctr = 0, row = 0;
      for (i = 0; i < sz; i++)
      {
        if(ctr == (sz/dim) ){
          row+=1;
          ctr = 0;
        }
        memcpy(op[row] + ctr*(bits/CHAR_BIT), s + i*(bits/CHAR_BIT), (bits/CHAR_BIT));
        ctr++;
      }

      for (d = 0; d < dim; d++)
        returnbuffer(qp[d], op[d], cbz, nx, ny, nz, compval, bits, PIDX_READ);

      //copy data to the input buffer
      ctr = 0;
      int pos = 0;
      for (i = 0; i < dim; i++)
      {
        for (j = 0; j < (sz/dim); j++)
        {
          memcpy(q + (pos + ctr)*(bits/CHAR_BIT), qp[i] + j*(bits/CHAR_BIT), (bits/CHAR_BIT));
          ctr+=dim;
        }
        pos+=1;
        ctr = 0;
      }
    }

    for (i = 0; i < dim; i++)
      free(sdim[i]);
    free(sdim);
    for (i = 0; i < dim; i++)
      free(op[i]);
    free(op);
    for (i = 0; i < dim; i++)
      free(qp[i]);
    free(qp);

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
