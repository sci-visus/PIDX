/*
 * BSD 3-Clause License
 * 
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

static void returnbuffer(unsigned char * q, unsigned char * s, int cbz, int nx, int ny, int nz, uint64_t compval, int bits, int mode);

//Struct for restructuring ID
struct PIDX_chunk_id_struct
{
  idx_comm idx_c;

  // Contains all relevant IDX file info
  // Blocks per file, samples per block, bitmask, block, file name template and more
  idx_dataset idx;


  int first_index;
  int last_index;
};

PIDX_chunk_id PIDX_chunk_init(idx_dataset idx_meta_data, idx_comm idx_c, int start_var_index, int end_var_index)
{
  PIDX_chunk_id chunk_id;

  //Creating the block rst ID
  chunk_id = (PIDX_chunk_id)malloc(sizeof (*chunk_id));
  memset(chunk_id, 0, sizeof (*chunk_id));

  chunk_id->idx = idx_meta_data;
  chunk_id->idx_c = idx_c;

  chunk_id->first_index = start_var_index;
  chunk_id->last_index = end_var_index;

  return chunk_id;
}



PIDX_return_code PIDX_chunk_meta_data_create(PIDX_chunk_id chunk_id)
{
  int v = 0;
  PIDX_variable var0 = chunk_id->idx->variable[chunk_id->first_index];
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = chunk_id->idx->variable[v];

    var->chunked_super_patch = malloc(sizeof(*(var->chunked_super_patch)));
    memset(var->chunked_super_patch, 0, sizeof(*(var->chunked_super_patch)));

    PIDX_super_patch out_patch = var->chunked_super_patch;
    PIDX_super_patch in_patch = var->restructured_super_patch;

    out_patch->patch_count = 1;
    out_patch->is_boundary_patch = in_patch->is_boundary_patch;

    out_patch->restructured_patch = malloc(sizeof(*(out_patch->restructured_patch)));
    memset(out_patch->restructured_patch, 0, sizeof(*(out_patch->restructured_patch)));

    memcpy(out_patch->restructured_patch->offset, in_patch->restructured_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
    memcpy(out_patch->restructured_patch->size, in_patch->restructured_patch->size, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
  }

  return PIDX_success;
}



// one restructured patch per group
PIDX_return_code PIDX_chunk_buf_create(PIDX_chunk_id chunk_id)
{
  PIDX_variable var0 = chunk_id->idx->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int v = 0;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = chunk_id->idx->variable[v];
    int bytes_per_value = (var->bpv / CHAR_BIT) * var->vps;

    PIDX_super_patch out_patch = var->chunked_super_patch;
    uint64_t *group_size = out_patch->restructured_patch->size;
    uint64_t num_elems_group = 1;
    int d;
    for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    {
      if (out_patch->restructured_patch->size[d] % chunk_id->idx->chunk_size[d] == 0)
        num_elems_group *= group_size[d];
      else
        num_elems_group *= (((group_size[d]/chunk_id->idx->chunk_size[d]) + 1) * chunk_id->idx->chunk_size[d]);
    }

    // malloc the storage for all elements in the output array
    // printf("[Chunking] Buffer size of %d = %d\n", chunk_id->idx_c->simulation_rank, num_elems_group);
    out_patch->restructured_patch->buffer = malloc(bytes_per_value * num_elems_group);
    memset(out_patch->restructured_patch->buffer, 0, bytes_per_value * num_elems_group);
    //printf("\n \n Out Patch is before chunking: \n");
    //printf("%s\n", out_patch->restructured_patch->buffer);
  }

  return PIDX_success;
}


static void returnbuffer(unsigned char * q, unsigned char * s, int cbz, int nx, int ny, int nz, uint64_t compval, int bits, int mode)
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
        uint64_t diff = (z/4) * (dz + 4 * (ny/4) * (dy+4 * (nx/4) * dx)) + (y/4) * (dy+4*(nx/4)*dx) + (x/4)*dx;
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
  PIDX_variable var0 = chunk_id->idx->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  // when no compression is used then chunking only involves copying the input buffer to the chunked buffer
  if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    for (uint32_t v = chunk_id->first_index; v <= chunk_id->last_index; v++)
    {
      PIDX_variable var = chunk_id->idx->variable[v];
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
  uint64_t *chunk_size = chunk_id->idx->chunk_size;
  uint64_t  cbz = chunk_size[0] * chunk_size[1] * chunk_size[2];


  // loop through all variables
  for (uint32_t v = chunk_id->first_index; v <= chunk_id->last_index; ++v)
  {
    PIDX_variable var = chunk_id->idx->variable[v];

    // copy the size and offset to output
    PIDX_super_patch in_patch = var->restructured_super_patch;
    PIDX_super_patch out_patch = var->chunked_super_patch;

    // output buffers have to multiples of fours
    int nx=(int)out_patch->restructured_patch->size[0];
    if (out_patch->restructured_patch->size[0] % chunk_id->idx->chunk_size[0] != 0)
      nx = ((out_patch->restructured_patch->size[0] / chunk_id->idx->chunk_size[0]) + 1) * chunk_id->idx->chunk_size[0];

    int ny=(int)out_patch->restructured_patch->size[1];
    if (out_patch->restructured_patch->size[1] % chunk_id->idx->chunk_size[1] != 0)
      ny = ((out_patch->restructured_patch->size[1] / chunk_id->idx->chunk_size[1]) + 1) * chunk_id->idx->chunk_size[1];

    int nz=(int)out_patch->restructured_patch->size[2];
    if (out_patch->restructured_patch->size[2] % chunk_id->idx->chunk_size[2] != 0)
      nz = ((out_patch->restructured_patch->size[2] / chunk_id->idx->chunk_size[2]) + 1) * chunk_id->idx->chunk_size[2];

    int values_per_sample = 0;
    int bits = 0;
    PIDX_get_datatype_details(var->type_name, &values_per_sample, &bits);

    unsigned char *sdim = malloc(nx*ny*nz * (bits/CHAR_BIT));
    unsigned char *op = malloc(nx*ny*nz * (bits/CHAR_BIT));

    uint32_t k = 0;
    uint64_t compval = in_patch->restructured_patch->size[0] * in_patch->restructured_patch->size[1] * in_patch->restructured_patch->size[2];
    if (MODE == PIDX_WRITE)
    {
      for (uint32_t i = 0; i < values_per_sample; i++)
      {
        //memset(sdim, 0, nx*ny*nz * (bits/CHAR_BIT));
        //memset(op, 0, nx*ny*nz * (bits/CHAR_BIT));

        k = i;
        for (uint32_t j = 0; j < nx*ny*nz; j++)
        {
          memcpy(sdim + j*(bits/CHAR_BIT), in_patch->restructured_patch->buffer + k* (bits/CHAR_BIT), (bits/CHAR_BIT));
          k+=values_per_sample;
        }

        returnbuffer(sdim, op, cbz, nx, ny, nz, compval, bits, PIDX_WRITE);
        memcpy(out_patch->restructured_patch->buffer + i * nx*ny*nz * (bits/CHAR_BIT), op, nx*ny*nz * (bits/CHAR_BIT));
      }
    }
    else
    {
      k = 0;
      for (uint32_t i = 0; i < values_per_sample; i++)
      {
        //memset(sdim, 0, nx*ny*nz * (bits/CHAR_BIT));
        //memset(op, 0, nx*ny*nz * (bits/CHAR_BIT));

        memcpy(op, out_patch->restructured_patch->buffer + i * nx*ny*nz * (bits/CHAR_BIT), nx*ny*nz * (bits/CHAR_BIT));
        returnbuffer(sdim, op, cbz, nx, ny, nz, compval, bits, PIDX_READ);

        k = i;
        for (uint32_t j = 0; j < nx*ny*nz; j++)
        {
          //double x;
          //memcpy(&x, sdim + j*(bits/CHAR_BIT), sizeof(double));
          //printf("vlue = %f\n", x);

          memcpy(in_patch->restructured_patch->buffer + (k)*(bits/CHAR_BIT), sdim + j*(bits/CHAR_BIT), (bits/CHAR_BIT));
          k+=values_per_sample;
        }
      }
    }

    free(sdim);
    free(op);
  }

  return PIDX_success;
}



PIDX_return_code PIDX_chunk_meta_data_destroy(PIDX_chunk_id chunk_id)
{
  int var;
  PIDX_variable var0 = chunk_id->idx->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {
    free(chunk_id->idx->variable[var]->chunked_super_patch->restructured_patch);
    free(chunk_id->idx->variable[var]->chunked_super_patch);
  }

  return PIDX_success;
}

PIDX_return_code PIDX_chunk_buf_destroy(PIDX_chunk_id chunk_id)
{
  PIDX_variable var0 = chunk_id->idx->variable[chunk_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int var = 0;
  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
    free(chunk_id->idx->variable[var]->chunked_super_patch->restructured_patch->buffer);

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
