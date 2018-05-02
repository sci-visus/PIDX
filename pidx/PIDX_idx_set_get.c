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

#include "PIDX_file_handler.h"


///
/// \brief PIDX_validate
/// \param file
/// \return
///
static PIDX_return_code PIDX_validate(PIDX_file file);

PIDX_return_code PIDX_set_meta_data_cache(PIDX_file file, PIDX_metadata_cache cache)
{
  if (!file)
    return PIDX_err_file;

  if (!cache)
    return PIDX_err_file;

  file->meta_data_cache = cache;
  return PIDX_success;
}



PIDX_return_code PIDX_get_meta_data_cache(PIDX_file file, PIDX_metadata_cache* cache)
{
  if (!file)
    return PIDX_err_file;

  if (!cache)
    return PIDX_err_file;

  *cache = file->meta_data_cache;
  return PIDX_success;
}


PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count)
{
  if (!file)
    return PIDX_err_file;

  if (variable_count <= 0)
    return PIDX_err_count;

  file->idx->variable_count = variable_count;
  file->idx->variable_pipe_length = file->idx->variable_count;

  return PIDX_success;
}



PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count)
{
  if (!file)
    return PIDX_err_file;

  *variable_count = file->idx->variable_count;

  return PIDX_success;
}



PIDX_return_code PIDX_set_physical_dims(PIDX_file file, PIDX_physical_point dims)
{
  if (!file)
    return PIDX_err_file;

  memcpy(file->idx->physical_bounds, dims, (sizeof(double) * PIDX_MAX_DIMENSIONS));
  memcpy(file->idx->physical_box_bounds, dims, (sizeof(double) * PIDX_MAX_DIMENSIONS));

  return PIDX_success;
}



PIDX_return_code PIDX_get_physical_dims(PIDX_file file, PIDX_physical_point dims)
{
  if (!file)
    return PIDX_err_file;

  memcpy(dims, file->idx->physical_bounds, (sizeof(double) * PIDX_MAX_DIMENSIONS));

  return PIDX_success;
}




PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int current_time_step)
{
  if (!file)
    return PIDX_err_file;

  if (current_time_step < 0)
    return PIDX_err_time;

  file->idx->current_time_step = current_time_step;

  return PIDX_success;
}



PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* current_time_step)
{
  if (!file)
    return PIDX_err_file;

  *current_time_step = file->idx->current_time_step;

  return PIDX_success;
}



PIDX_return_code PIDX_set_block_size(PIDX_file file, const int bits_per_block)
{
  if (!file)
    return PIDX_err_file;

  if (bits_per_block <= 0)
    return PIDX_err_block;

  file->idx->bits_per_block = bits_per_block;
  file->idx->samples_per_block = (int)pow(2, bits_per_block);

  return PIDX_validate(file);
}



PIDX_return_code PIDX_get_block_size(PIDX_file file, int* bits_per_block)
{
  if (!file)
    return PIDX_err_file;

  *bits_per_block = file->idx->bits_per_block;

  return PIDX_success;
}



PIDX_return_code PIDX_set_block_count(PIDX_file file, const int blocks_per_file)
{
  if (!file)
    return PIDX_err_file;

  if (blocks_per_file <= 0)
    return PIDX_err_block;

  file->idx->blocks_per_file = blocks_per_file;

  return PIDX_success;
}



PIDX_return_code PIDX_get_block_count(PIDX_file file, int* blocks_per_file)
{
  if (!file)
    return PIDX_err_file;

  *blocks_per_file = file->idx->blocks_per_file;

  return PIDX_success;
}



PIDX_return_code PIDX_set_resolution(PIDX_file file, int hz_from, int hz_to)
{
  if (file == NULL)
    return PIDX_err_file;

  file->idx_b->reduced_res_from = hz_from;
  file->idx_b->reduced_res_to = hz_to;

  return PIDX_success;
}



PIDX_return_code PIDX_get_resolution(PIDX_file file, int *hz_from, int *hz_to)
{
  if (file == NULL)
    return PIDX_err_file;

  *hz_from = file->idx_b->reduced_res_from;
  *hz_to = file->idx_b->reduced_res_to;

  return PIDX_success;
}



PIDX_return_code PIDX_set_partition_count(PIDX_file file, int count_x, int count_y, int count_z)
{
  if (count_x < 0 || count_y < 0 || count_z < 0)
    return PIDX_err_box;

  if (file == NULL)
    return PIDX_err_file;

  file->idx->partition_count[0] = count_x;
  file->idx->partition_count[1] = count_y;
  file->idx->partition_count[2] = count_z;

  return PIDX_success;
}



PIDX_return_code PIDX_get_partition_count(PIDX_file file, int* count_x, int* count_y, int* count_z)
{
  if (file == NULL)
    return PIDX_err_file;

  *count_x = file->idx->partition_count[0];
  *count_y = file->idx->partition_count[1];
  *count_z = file->idx->partition_count[2];

  return PIDX_success;
}



PIDX_return_code PIDX_set_restructuring_box(PIDX_file file, PIDX_point reg_patch_size)
{
  if (file == NULL)
    return PIDX_err_file;

  if (((reg_patch_size[0] & (reg_patch_size[0] - 1)) != 0) && ((reg_patch_size[1] & (reg_patch_size[1] - 1)) != 0) && ((reg_patch_size[2] & (reg_patch_size[2] - 1)) != 0))
  {
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "Warning in %s %d restructuring box needs to be power of two in size %d %d %d\n", __FILE__, __LINE__, (int)reg_patch_size[0], (int)reg_patch_size[1], (int)reg_patch_size[2]);
    //return PIDX_err_box;
  }

  memcpy(file->restructured_grid->patch_size, reg_patch_size, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  return PIDX_success;
}



PIDX_return_code PIDX_get_restructuring_box(PIDX_file file, PIDX_point reg_patch_size)
{
  if (file == NULL)
    return PIDX_err_file;

  memcpy(reg_patch_size, file->restructured_grid->patch_size, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  return PIDX_success;
}



PIDX_return_code PIDX_set_first_tstep(PIDX_file file, int tstep)
{
  if (!file)
    return PIDX_err_file;

  file->idx->first_tstep = tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_get_first_tstep(PIDX_file file, int* tstep)
{
  if (!file)
    return PIDX_err_file;

  *tstep = file->idx->first_tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_set_last_tstep(PIDX_file file, int tstep)
{
  if (!file)
    return PIDX_err_file;

  file->idx->last_tstep = tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_get_last_tstep(PIDX_file file, int* tstep)
{
  if (!file)
    return PIDX_err_file;

  *tstep = file->idx->last_tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type)
{
  if (!file)
    return PIDX_err_file;

  if (compression_type != PIDX_NO_COMPRESSION && compression_type != PIDX_CHUNKING_ONLY && compression_type != PIDX_CHUNKING_ZFP)
    return PIDX_err_unsupported_compression_type;

  file->idx->compression_type = compression_type;

  if (file->idx->compression_type == PIDX_NO_COMPRESSION)
    return PIDX_success;
  else if (file->idx->compression_type == PIDX_CHUNKING_ONLY || file->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    file->idx->chunk_size[0] = 4;
    file->idx->chunk_size[1] = 4;
    file->idx->chunk_size[2] = 4;

    int reduce_by_sample = 1;
    uint64_t total_chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
    if (reduce_by_sample == 1)
    {
      file->idx->bits_per_block = file->idx->bits_per_block - (int)log2(total_chunk_size);
      file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

      if (file->idx->bits_per_block <= 0)
      {
        file->idx->bits_per_block = 0;
        file->idx->samples_per_block = 1;
      }
    }
    else
      file->idx->blocks_per_file = file->idx->blocks_per_file / total_chunk_size;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type)
{
  if (!file)
    return PIDX_err_file;

   *compression_type = file->idx->compression_type;

  return PIDX_success;
}


PIDX_return_code PIDX_set_average_compression_factor(PIDX_file file, int compression_factor, float bit_rate)
{
  if (!file)
    return PIDX_err_file;

  file->idx->compression_factor = compression_factor;
  file->idx->compression_bit_rate = bit_rate;

  file->idx->bits_per_block = file->idx->bits_per_block + 6;
  file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

  return PIDX_success;
}


PIDX_return_code PIDX_set_lossy_compression_bit_rate(PIDX_file file, float compression_bit_rate)
{
  if (!file)
    return PIDX_err_file;

  if (file->idx->compression_type == PIDX_CHUNKING_ONLY)
    return PIDX_success;

  //if (file->idx->compression_type != PIDX_CHUNKING_ZFP)
  //  return PIDX_err_unsupported_compression_type;

  file->idx->compression_bit_rate = compression_bit_rate;

  if (file->idx->compression_bit_rate == 64)
  {
    file->idx->bits_per_block = file->idx->bits_per_block;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);
    file->idx->compression_factor = 1;
  }
  if (file->idx->compression_bit_rate == 32)
  {
    file->idx->bits_per_block = file->idx->bits_per_block;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);
    file->idx->compression_factor = 2;
  }
  if (file->idx->compression_bit_rate == 16)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 1;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 4;
  }
  if (file->idx->compression_bit_rate == 8)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 2;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 8;
  }
  if (file->idx->compression_bit_rate == 4)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 3;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 16;
  }
  if (file->idx->compression_bit_rate == 2)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 4;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 32;
  }

  // TODO : major hack here
  if (file->idx->compression_bit_rate == 1)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 5;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 64;
  }

  if (file->idx->compression_bit_rate == 0.5)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 6;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 128;
  }

  if (file->idx->compression_bit_rate == 0.25)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 7;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 256;
  }

  if (file->idx->compression_bit_rate == 0.125)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 8;
    file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 512;
  }

  if (file->idx->bits_per_block <= 0)
  {
    file->idx->bits_per_block = 0;
    file->idx->samples_per_block = 1;
  }

  return PIDX_validate(file);
}



PIDX_return_code PIDX_get_lossy_compression_bit_rate(PIDX_file file, int *compression_bit_rate)
{
  if (!file)
    return PIDX_err_file;

  *compression_bit_rate = file->idx->compression_bit_rate;

  return PIDX_success;
}



PIDX_return_code PIDX_set_io_mode(PIDX_file file, enum PIDX_io_type io_type)
{
  if (file == NULL)
    return PIDX_err_file;

  file->idx->io_type = io_type;

  return PIDX_success;
}



PIDX_return_code PIDX_get_io_mode(PIDX_file file, enum PIDX_io_type* io_type)
{
  if (file == NULL)
    return PIDX_err_file;

  *io_type = file->idx->io_type;

  return PIDX_success;
}



#if 0
PIDX_return_code PIDX_set_wavelet_level(PIDX_file file, int w_level)
{
  if (file == NULL)
    return PIDX_err_file;

  file->idx_d->wavelet_levels = w_level;

  return PIDX_success;
}



PIDX_return_code PIDX_get_wavelet_level(PIDX_file file, int* w_level)
{
  if (file == NULL)
    return PIDX_err_file;

  *w_level = file->idx_d->wavelet_levels;

  return PIDX_success;
}



PIDX_return_code PIDX_set_ROI_type(PIDX_file file, int type)
{
  if (file == NULL)
    return PIDX_err_file;

  file->ROI_writes = 1;

  return PIDX_success;
}



PIDX_return_code PIDX_get_ROI_type(PIDX_file file, int* type)
{
  if (file == NULL)
    return PIDX_err_file;

  *type = 0;
  //TODO

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_set_variable_pile_length(PIDX_file file, int var_pipe_length)
{
  if (!file)
    return PIDX_err_file;

  if (var_pipe_length < 0)
    return PIDX_err_size;

  file->idx->variable_pipe_length = var_pipe_length;

  return PIDX_success;
}



PIDX_return_code PIDX_get_variable_pile_length(PIDX_file file, int* var_pipe_length)
{
  if (!file)
    return PIDX_err_file;

  if (var_pipe_length < 0)
    return PIDX_err_size;

  *var_pipe_length = file->idx->variable_pipe_length;

  return PIDX_success;
}



PIDX_return_code PIDX_save_big_endian(PIDX_file file)
{
  file->idx->endian = 0;

  unsigned int endian = 1;
  char *c = (char*)&endian;
  if (*c)
    file->idx->flip_endian = 1;

  return PIDX_success;
}



PIDX_return_code PIDX_save_little_endian(PIDX_file file)
{
  file->idx->endian = PIDX_LITTLE_ENDIAN;

  unsigned int endian = 1;
  char *c = (char*)&endian;
  if (!*c)
    file->idx->flip_endian = 1;

  return PIDX_success;
}



PIDX_return_code PIDX_set_cache_time_step(PIDX_file file, int ts)
{
  if (file == NULL)
    return PIDX_err_file;

  file->idx->cached_ts = ts;

  return PIDX_success;
}



PIDX_return_code PIDX_get_cache_time_step(PIDX_file file, int* ts)
{
  if (file == NULL)
    return PIDX_err_file;

  *ts = file->idx->cached_ts;

  return PIDX_success;
}


/*
PIDX_return_code PIDX_set_process_decomposition(PIDX_file file, int np_x, int np_y, int np_z)
{
  if (file == NULL)
    return PIDX_err_file;

  file->idx_c->gnproc_x = np_x;
  file->idx_c->gnproc_y = np_y;
  file->idx_c->gnproc_z = np_z;

  return PIDX_success;
}




PIDX_return_code PIDX_get_process_decomposition(PIDX_file file, int* np_x, int* np_y, int* np_z)
{
  if (file == NULL)
    return PIDX_err_file;

  *np_x = file->idx_c->gnproc_x;
  *np_y = file->idx_c->gnproc_y;
  *np_z = file->idx_c->gnproc_z;

  return PIDX_success;
}
*/


PIDX_return_code PIDX_get_comm(PIDX_file file, MPI_Comm *comm)
{
  if (!file)
    return PIDX_err_file;

   *comm = file->idx_c->simulation_comm;

  return PIDX_success;
}


PIDX_return_code PIDX_set_comm(PIDX_file file, MPI_Comm comm)
{
  if (!file)
    return PIDX_err_file;

   file->idx_c->simulation_comm = comm;
   file->idx_c->partition_comm = comm;

   MPI_Comm_rank(file->idx_c->simulation_comm, &(file->idx_c->simulation_rank));
   MPI_Comm_size(file->idx_c->simulation_comm, &(file->idx_c->simulation_nprocs));
   MPI_Comm_rank(file->idx_c->partition_comm, &(file->idx_c->partition_rank));
   MPI_Comm_size(file->idx_c->partition_comm, &(file->idx_c->partition_nprocs));

  return PIDX_success;
}



static PIDX_return_code PIDX_validate(PIDX_file file)
{
  uint64_t dims;
  uint64_t adjusted_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_bounds[0] = file->idx->bounds[0] / file->idx->chunk_size[0];
  adjusted_bounds[1] = file->idx->bounds[1] / file->idx->chunk_size[1];
  adjusted_bounds[2] = file->idx->bounds[2] / file->idx->chunk_size[2];

  if (PIDX_inner_product(&dims, adjusted_bounds))
    return PIDX_err_size;

  //fprintf(stderr, "dims %d spb %d\n", dims, file->idx->samples_per_block);

  if (dims < file->idx->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx->bits_per_block = getNumBits(file->idx->samples_per_block) - 1;
  }

  if (file->idx->bits_per_block == 0)
  {
    file->idx->bits_per_block = 0;
    file->idx->samples_per_block = 1;
  }

  // other validations...
  // TODO

  return PIDX_success;
}

