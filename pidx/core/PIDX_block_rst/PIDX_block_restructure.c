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
#if PIDX_HAVE_MPI
  /// Passed by PIDX API
  MPI_Comm comm;
#endif

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, block, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;

  //int if_perform_chunk;

  int init_index;
  int first_index;
  int last_index;
};

PIDX_chunk_id PIDX_chunk_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, int init_index, int start_var_index, int end_var_index)
{
  PIDX_chunk_id chunk_id;

  //Creating the block rst ID
  chunk_id = (PIDX_chunk_id)malloc(sizeof (*chunk_id));
  memset(chunk_id, 0, sizeof (*chunk_id));

  chunk_id->idx = idx_meta_data;
  chunk_id->idx_derived = idx_derived;

  chunk_id->init_index = init_index;
  chunk_id->first_index = start_var_index;
  chunk_id->last_index = end_var_index;

  /*
  if (chunk_id->idx->enable_compression == 1)
  {
    // No Chunking and compression without restructuring
    if (chunk_id->idx->enable_rst != 1)
    {
      chunk_id->idx->enable_compression = 0;
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        chunk_id->idx->chunk_size[d] = 1;
      chunk_id->idx->compression_bit_rate = 64;
    }
    else
    {
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        if (chunk_id->idx->bounds[d] % chunk_id->idx->chunk_size[d] != 0)
        {
          chunk_id->idx->enable_compression = 0;
          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            chunk_id->idx->chunk_size[d] = 1;
          chunk_id->idx->compression_bit_rate = 64;

          break;
        }
      }
    }
  }
  */

  return chunk_id;
}


#if PIDX_HAVE_MPI
int PIDX_chunk_set_communicator(PIDX_chunk_id chunk_id, MPI_Comm comm)
{
  if (chunk_id == NULL)
    return PIDX_err_id;

  chunk_id->comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_chunk_meta_data_create(PIDX_chunk_id chunk_id)
{
  int v = 0, p = 0, j = 0;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = chunk_id->idx->variable[v];

    var->chunk_patch_group = malloc(sizeof(*var->chunk_patch_group) * var->patch_group_count);
    memset(var->chunk_patch_group, 0, sizeof(*var->chunk_patch_group) * var->patch_group_count);

    for (p = 0; p < var->patch_group_count; p++)
    {
      var->chunk_patch_group[p] = malloc(sizeof(*(var->chunk_patch_group[p])));
      memset(var->chunk_patch_group[p], 0, sizeof(*(var->chunk_patch_group[p])));

      Ndim_patch_group out_patch = var->chunk_patch_group[p];
      Ndim_patch_group in_patch = var->rst_patch_group[p];

      if (chunk_id->idx->compression_type == PIDX_CHUNKING_ONLY || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP)
        out_patch->count = 1;
      else
        out_patch->count = in_patch->count;

      if (in_patch->data_source == 1)
        out_patch->count = 1;

      out_patch->type = in_patch->type;

      out_patch->patch = malloc(sizeof(*(out_patch->patch)) * out_patch->count);
      memset(out_patch->patch, 0, sizeof(*(out_patch->patch)) * out_patch->count);

      out_patch->reg_patch = malloc(sizeof(*(out_patch->reg_patch)));
      memset(out_patch->reg_patch, 0, sizeof(*(out_patch->reg_patch)));

      for(j = 0; j < out_patch->count; j++)
      {
        out_patch->patch[j] = malloc(sizeof(*(out_patch->patch[j])));
        memset(out_patch->patch[j], 0, sizeof(*(out_patch->patch[j])));

        if (chunk_id->idx->compression_type == PIDX_CHUNKING_ONLY || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP || in_patch->data_source == 1)
        {
          memcpy(out_patch->patch[j]->size, in_patch->reg_patch->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(out_patch->patch[j]->offset, in_patch->reg_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
        else if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
        {
          memcpy(out_patch->patch[j]->offset, in_patch->patch[j]->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(out_patch->patch[j]->size, in_patch->patch[j]->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
      }
      memcpy(out_patch->reg_patch->offset, in_patch->reg_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
      memcpy(out_patch->reg_patch->size, in_patch->reg_patch->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
    }
  }

  return PIDX_success;
}

// one restructured patch per group
// TODO: detect wrong input
PIDX_return_code PIDX_chunk_buf_create(PIDX_chunk_id chunk_id)
{
#if !SIMULATE_IO
  int v = 0, p = 0, j = 0;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
  {
    PIDX_variable var = chunk_id->idx->variable[v];
    int bytes_per_value = var->bits_per_value / 8;

    for (p = 0; p < var->patch_group_count; p++)
    {
      Ndim_patch_group out_patch = var->chunk_patch_group[p];
      Ndim_patch_group in_patch = var->rst_patch_group[p];

      int64_t *group_size = in_patch->reg_patch->size;

      for(j = 0; j < out_patch->count; j++)
      {

        if (chunk_id->idx->compression_type == PIDX_CHUNKING_ONLY || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP || in_patch->data_source == 1)
        {
          // compute the number of elements in the group
          int64_t num_elems_group = 1; // number of elements in the group
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
        }
        else if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
        {
          out_patch->patch[j]->buffer = malloc(out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * out_patch->patch[j]->size[3] * out_patch->patch[j]->size[4] * bytes_per_value * var->values_per_sample);
          memset(out_patch->patch[j]->buffer, 0, out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * out_patch->patch[j]->size[3] * out_patch->patch[j]->size[4] * bytes_per_value * var->values_per_sample);
        }
      }

      /*
      for(j = 0; j < in_patch->count; j++)
      {
        if (chunk_id->idx->compression_type == PIDX_CHUNKING_ONLY || chunk_id->idx->compression_type == PIDX_CHUNKING_ZFP)
        {
          char filename[100];
          sprintf(filename, "%d_%d_%d",p, j, rank);
          FILE *fp = fopen (filename, "w");
          double dv;
          int i;

          fprintf(fp, "[%d] [%d] Offset Count %d %d %d :: %d %d %d\n", p, j, (int)in_patch->patch[j]->offset[0], (int)in_patch->patch[j]->offset[1], (int)in_patch->patch[j]->offset[2], (int)in_patch->patch[j]->size[0], (int)in_patch->patch[j]->size[1], (int)in_patch->patch[j]->size[2]);
          for (i = 0; i < in_patch->patch[j]->size[0] * in_patch->patch[j]->size[1] * in_patch->patch[j]->size[2]; i++)
          {
            //printf("R%d [%d] %d %d : %d\n", rank, i, p, j, var->bits_per_value/8);
            memcpy(&dv, in_patch->patch[j]->buffer + (i*var->bits_per_value / 8), var->bits_per_value/8);
            fprintf(fp, "%f\n", dv);
          }
          fclose(fp);
        }
      }
      */
    }
  }
#endif
  return PIDX_success;
}


PIDX_return_code PIDX_chunk(PIDX_chunk_id chunk_id, int MODE)
{
  int64_t v,p,j;
  if (chunk_id->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    for (v = chunk_id->first_index; v <= chunk_id->last_index; v++)
    {
      PIDX_variable var = chunk_id->idx->variable[v];
      int bytes_per_value = var->bits_per_value / 8;

      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group out_patch = var->chunk_patch_group[p];
        Ndim_patch_group in_patch = var->rst_patch_group[p];

        for(j = 0; j < out_patch->count; j++)
        {
          if (in_patch->data_source == 0)
          {
#if !SIMULATE_IO
            if (MODE == PIDX_WRITE)
              memcpy(out_patch->patch[j]->buffer, in_patch->patch[j]->buffer, out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * out_patch->patch[j]->size[3] * out_patch->patch[j]->size[4] * bytes_per_value * var->values_per_sample);
            else
              memcpy(in_patch->patch[j]->buffer, out_patch->patch[j]->buffer, out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * out_patch->patch[j]->size[3] * out_patch->patch[j]->size[4] * bytes_per_value * var->values_per_sample);
#endif
          }
          else
          {
            if (MODE == PIDX_WRITE)
              memcpy(out_patch->patch[j]->buffer, in_patch->reg_patch->buffer, out_patch->reg_patch->size[0] * out_patch->reg_patch->size[1] * out_patch->reg_patch->size[2] * out_patch->reg_patch->size[3] * out_patch->reg_patch->size[4] * bytes_per_value * var->values_per_sample);
            else
              memcpy(in_patch->patch[j]->buffer, out_patch->patch[j]->buffer, out_patch->patch[j]->size[0] * out_patch->patch[j]->size[1] * out_patch->patch[j]->size[2] * out_patch->patch[j]->size[3] * out_patch->patch[j]->size[4] * bytes_per_value * var->values_per_sample);
          }
        }
      }
    }
    return PIDX_success;
  }


  // compute the intra compression block strides
  int64_t *chunk_size = chunk_id->idx->chunk_size;
  int64_t  cbz = 1;
  int d;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    cbz = cbz * chunk_size[d];

  // loop through all variables
  v = 0;

  for (v = chunk_id->first_index; v <= chunk_id->last_index; ++v)
  {
    PIDX_variable var = chunk_id->idx->variable[v];
    //int bytes_per_value = var->bits_per_value / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      Ndim_patch_group patch_group = var->rst_patch_group[g];
      Ndim_patch_group out_patch = var->chunk_patch_group[g];

      int nx=(int)out_patch->patch[0]->size[0];
      if (out_patch->patch[0]->size[0] % chunk_id->idx->chunk_size[0] != 0)
        nx = ((out_patch->patch[0]->size[0] / chunk_id->idx->chunk_size[0]) + 1) * chunk_id->idx->chunk_size[0];

      int ny=(int)out_patch->patch[0]->size[1];
      if (out_patch->patch[0]->size[1] % chunk_id->idx->chunk_size[1] != 0)
        ny = ((out_patch->patch[0]->size[1] / chunk_id->idx->chunk_size[1]) + 1) * chunk_id->idx->chunk_size[1];

      int nz=(int)out_patch->patch[0]->size[2];
      if (out_patch->patch[0]->size[2] % chunk_id->idx->chunk_size[2] != 0)
        nz = ((out_patch->patch[0]->size[2] / chunk_id->idx->chunk_size[2]) + 1) * chunk_id->idx->chunk_size[2];

      unsigned char* temp_buffer = malloc(nx * ny * nz * var->bits_per_value/8 * var->values_per_sample);
      if (temp_buffer == NULL)
        return PIDX_err_chunk;


      int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
      for (r = 0; r < var->rst_patch_group[g]->count; r++)
      {
        for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
        {
          for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
          {
            for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
            {
              index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
              send_o = index * var->values_per_sample * (var->bits_per_value/8);
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = ((nx) * (ny) * (k1 - out_patch->patch[0]->offset[2])) + ((nx)* (j1 - out_patch->patch[0]->offset[1])) + (i1 - out_patch->patch[0]->offset[0]);
#if !SIMULATE_IO
              if (MODE == PIDX_WRITE)
                memcpy(temp_buffer + (recv_o * var->values_per_sample * (var->bits_per_value/8)), var->rst_patch_group[g]->patch[r]->buffer + send_o, send_c * var->values_per_sample * (var->bits_per_value/8));
#endif
            }
          }
        }
      }

      double* p=(double*)temp_buffer;
      int dz=4*nx*(ny-(ny/4)*4);
      int dy=4*(nx-(nx/4)*4);
      int dx=4;      

      int s1 = 0;
      int z, y, x;
      int zz, yy, xx;      
      for (s1 = 0; s1 < var->values_per_sample; s1++)
      {
        for (z=0;z<nz;z+=4)
        {
#if !SIMULATE_IO
          double* s = (double*)(out_patch->patch[0]->buffer+((z/4)*((ny+3)/4)*((nx+3)/4))*(var->bits_per_value/8)*cbz) + nx*ny*nz*s1;
#endif
          for (y=0;y<ny;y+=4)
          {
            for (x=0;x<nx;x+=4)
            {
              int64_t diff=(z/4)*(dz+4*(ny/4)*(dy+4*(nx/4)*dx))+(y/4)*(dy+4*(nx/4)*dx)+(x/4)*dx;
              double* q=p+var->values_per_sample*diff+s1;

              for (zz = 0; zz < 4; ++zz)
              {
                for (yy = 0; yy < 4; ++yy)
                {
                  for (xx = 0; xx < 4; ++xx)
                  {
                    int i = xx + yy * 4 + zz * 4 * 4;
                    int j = xx + yy * nx + zz * nx * ny;
                    //printf("%d %d %f\n", i, j, s[i]);
#if !SIMULATE_IO
                    if (MODE == PIDX_WRITE)
                      s[nx*ny*nz*s1 + i*var->values_per_sample+s1] = q[j*var->values_per_sample + s1];
                    else
                    {
                      //printf("VAL %f\n", s[nx*ny*nz*s1 + i*var->values_per_sample+s1]);
                      q[j*var->values_per_sample + s1] = s[nx*ny*nz*s1 + i*var->values_per_sample+s1];
                    }
#endif
                  }
                }
              }
#if !SIMULATE_IO
              s += cbz;
#endif
            }
          }
        }
      }

      for (r = 0; r < var->rst_patch_group[g]->count; r++)
      {
        for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
        {
          for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
          {
            for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
            {
              index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
              send_o = index * var->values_per_sample * (var->bits_per_value/8);
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = ((nx) * (ny) * (k1 - out_patch->patch[0]->offset[2])) + ((nx)* (j1 - out_patch->patch[0]->offset[1])) + (i1 - out_patch->patch[0]->offset[0]);
#if !SIMULATE_IO
              if (MODE == PIDX_READ)
                memcpy(var->rst_patch_group[g]->patch[r]->buffer + send_o, temp_buffer + (recv_o * var->values_per_sample * (var->bits_per_value/8)), send_c * var->values_per_sample * (var->bits_per_value/8));
#endif
            }
          }
        }
      }

      free(temp_buffer);
    }
  }

  /*
  Remove the comment to print debug info
  double dv;
  int i;
  for (v = chunk_id->first_index; v <= chunk_id->last_index; ++v)
  {
    PIDX_variable var = chunk_id->idx->variable[v];

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      char filename[100];
      sprintf(filename, "%d_%d",g, rank);
      FILE *fp = fopen (filename, "w");

      Ndim_patch_group out_patch = var->chunk_patch_group[g];
      fprintf(fp, "Offset Count %d %d %d :: %d %d %d\n", (int)out_patch->reg_patch_offset[0], (int)out_patch->reg_patch_offset[1], (int)out_patch->reg_patch_offset[2], (int)out_patch->reg_patch_size[0], (int)out_patch->reg_patch_size[1], (int)out_patch->reg_patch_size[2]);
      for (i = 0; i < out_patch->reg_patch_size[0] * out_patch->reg_patch_size[1] * out_patch->reg_patch_size[2]; i++)
      {
        memcpy(&dv, out_patch->patch[0]->buffer + (i*var->bits_per_value/8), var->bits_per_value/8);
        fprintf(fp, "%f\n", dv);
      }
      fclose(fp);
    }
  }
  */

  return PIDX_success;
}



PIDX_return_code PIDX_chunk_meta_data_destroy(PIDX_chunk_id chunk_id)
{
  int j, p, var;

  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {
    for (p = 0; p < chunk_id->idx->variable[var]->patch_group_count; p++)
    {
      for(j = 0; j < chunk_id->idx->variable[var]->chunk_patch_group[p]->count; j++)
      {
        free(chunk_id->idx->variable[var]->chunk_patch_group[p]->patch[j]);
        chunk_id->idx->variable[var]->chunk_patch_group[p]->patch[j] = 0;
      }

      free(chunk_id->idx->variable[var]->chunk_patch_group[p]->patch);
      chunk_id->idx->variable[var]->chunk_patch_group[p]->patch = 0;

      free(chunk_id->idx->variable[var]->chunk_patch_group[p]->reg_patch);
      chunk_id->idx->variable[var]->chunk_patch_group[p]->reg_patch = 0;

      free(chunk_id->idx->variable[var]->chunk_patch_group[p]);
      chunk_id->idx->variable[var]->chunk_patch_group[p] = 0;
    }
    free(chunk_id->idx->variable[var]->chunk_patch_group);
    chunk_id->idx->variable[var]->chunk_patch_group = 0;
  }

  return PIDX_success;
}

PIDX_return_code PIDX_chunk_buf_destroy(PIDX_chunk_id chunk_id)
{
#if !SIMULATE_IO
  int j = 0, p = 0, var = 0;
  for (var = chunk_id->first_index; var <= chunk_id->last_index; var++)
  {
    for (p = 0; p < chunk_id->idx->variable[var]->patch_group_count; p++)
    {
      for(j = 0; j < chunk_id->idx->variable[var]->chunk_patch_group[p]->count; j++)
      {
        free(chunk_id->idx->variable[var]->chunk_patch_group[p]->patch[j]->buffer);
        chunk_id->idx->variable[var]->chunk_patch_group[p]->patch[j]->buffer = 0;
      }
    }
  }
#endif
  return PIDX_success;
}


PIDX_return_code PIDX_chunk_finalize(PIDX_chunk_id chunk_id)
{
  //TODO?
  free(chunk_id);
  chunk_id = 0;

  return PIDX_success;
}


int HELPER_chunking(PIDX_chunk_id id)
{
  return PIDX_success;
}
