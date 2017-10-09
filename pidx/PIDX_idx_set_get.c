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

#include "PIDX_file_handler.h"


///
/// \brief PIDX_validate
/// \param file
/// \return
///
static PIDX_return_code PIDX_validate(PIDX_file file);

PIDX_return_code PIDX_set_meta_data_cache(PIDX_file file, PIDX_meta_data_cache cache)
{
  if(!file)
    return PIDX_err_file;

  if(!cache)
    return PIDX_err_file;

  file->idx_cache->meta_data_cache = cache;
  return PIDX_success;
}



PIDX_return_code PIDX_get_meta_data_cache(PIDX_file file, PIDX_meta_data_cache* cache)
{
  if(!file)
    return PIDX_err_file;

  if(!cache)
    return PIDX_err_file;

  *cache = file->idx_cache->meta_data_cache;
  return PIDX_success;
}


PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count)
{
  if(!file)
    return PIDX_err_file;

  if(variable_count <= 0)
    return PIDX_err_count;

  file->idx->variable_count = variable_count;
  file->idx->variable_grp[0]->variable_count = variable_count;
  file->idx_d->var_pipe_length = file->idx->variable_count - 1 > 1 ? file->idx->variable_count - 1 : 1;

  int total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
  if (total_header_size % file->idx_d->fs_block_size)
    file->idx_d->start_fs_block++;

  return PIDX_success;
}



PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count)
{
  if(!file)
    return PIDX_err_file;

  *variable_count = file->idx->variable_count;

  return PIDX_success;
}



PIDX_return_code PIDX_set_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;

  memcpy(file->idx->transform, transform, (sizeof(double) * 16));

  return PIDX_success;
}



PIDX_return_code PIDX_get_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;

  memcpy(transform, file->idx->transform, (sizeof(double) * 16));

  return PIDX_success;
}



PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int current_time_step)
{
  if(!file)
    return PIDX_err_file;

  if(current_time_step < 0)
    return PIDX_err_time;

  file->idx->current_time_step = current_time_step;

  if (file->idx->cached_ts == -1)
    file->idx->cached_ts = current_time_step;

  return PIDX_success;
}



PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* current_time_step)
{
  if(!file)
    return PIDX_err_file;

  *current_time_step = file->idx->current_time_step;

  return PIDX_success;
}



PIDX_return_code PIDX_set_block_size(PIDX_file file, const int bits_per_block)
{
  if(!file)
    return PIDX_err_file;

  if(bits_per_block <= 0)
    return PIDX_err_block;

  file->idx->bits_per_block = bits_per_block;
  file->idx_d->samples_per_block = (int)pow(2, bits_per_block);

  return PIDX_validate(file);
}



PIDX_return_code PIDX_get_block_size(PIDX_file file, int* bits_per_block)
{
  if(!file)
    return PIDX_err_file;

  *bits_per_block = file->idx->bits_per_block;

  return PIDX_success;
}



PIDX_return_code PIDX_set_block_count(PIDX_file file, const int blocks_per_file)
{
  if(!file)
    return PIDX_err_file;

  if(blocks_per_file <= 0)
    return PIDX_err_block;

  file->idx->blocks_per_file = blocks_per_file;

  int total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
  if (total_header_size % file->idx_d->fs_block_size)
    file->idx_d->start_fs_block++;

  return PIDX_success;
}



PIDX_return_code PIDX_get_block_count(PIDX_file file, int* blocks_per_file)
{
  if(!file)
    return PIDX_err_file;

  *blocks_per_file = file->idx->blocks_per_file;

  return PIDX_success;
}



PIDX_return_code PIDX_set_resolution(PIDX_file file, int hz_from, int hz_to)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->reduced_res_from = hz_from;
  file->idx_d->reduced_res_to = hz_to;

  return PIDX_success;
}



PIDX_return_code PIDX_get_resolution(PIDX_file file, int *hz_from, int *hz_to)
{
  if(file == NULL)
    return PIDX_err_file;

  *hz_from = file->idx_d->reduced_res_from;
  *hz_to = file->idx_d->reduced_res_to;

  return PIDX_success;
}



PIDX_return_code PIDX_set_partition_count(PIDX_file file, int count_x, int count_y, int count_z)
{
  if(count_x < 0 || count_y < 0 || count_z < 0)
    return PIDX_err_box;

  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->partition_count[0] = count_x;
  file->idx_d->partition_count[1] = count_y;
  file->idx_d->partition_count[2] = count_z;

  return PIDX_success;
}



PIDX_return_code PIDX_get_partition_count(PIDX_file file, int* count_x, int* count_y, int* count_z)
{
  if(file == NULL)
    return PIDX_err_file;

  *count_x = file->idx_d->partition_count[0];
  *count_y = file->idx_d->partition_count[1];
  *count_z = file->idx_d->partition_count[2];

  return PIDX_success;
}



PIDX_return_code PIDX_set_restructuring_box(PIDX_file file, PIDX_point reg_patch_size)
{
  if (file == NULL)
    return PIDX_err_file;

  if (((reg_patch_size[0] & (reg_patch_size[0] - 1)) != 0) && ((reg_patch_size[1] & (reg_patch_size[1] - 1)) != 0) && ((reg_patch_size[2] & (reg_patch_size[2] - 1)) != 0))
  {
    if (file->idx_c->grank == 0)
      fprintf(stderr, "Warning in %s %d restructuring box needs to be power of two in size %d %d %d (will be over written!!!!)\n", __FILE__, __LINE__, (int)reg_patch_size[0], (int)reg_patch_size[1], (int)reg_patch_size[2]);

    return PIDX_err_box;
  }

  memcpy(file->idx_d->restructured_grid->patch_size, reg_patch_size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  return PIDX_success;
}



PIDX_return_code PIDX_get_restructuring_box(PIDX_file file, PIDX_point reg_patch_size)
{
  if (file == NULL)
    return PIDX_err_file;

  memcpy(reg_patch_size, file->idx_d->restructured_grid->patch_size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  return PIDX_success;
}



PIDX_return_code PIDX_set_first_tstep(PIDX_file file, int tstep)
{
  if(!file)
    return PIDX_err_file;

  file->idx->first_tstep = tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_get_first_tstep(PIDX_file file, int* tstep)
{
  if(!file)
    return PIDX_err_file;

  *tstep = file->idx->first_tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_set_last_tstep(PIDX_file file, int tstep)
{
  if(!file)
    return PIDX_err_file;

  file->idx->last_tstep = tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_get_last_tstep(PIDX_file file, int* tstep)
{
  if(!file)
    return PIDX_err_file;

  *tstep = file->idx->last_tstep;

  return PIDX_success;
}



PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type)
{
  if(!file)
    return PIDX_err_file;

  if (compression_type != PIDX_NO_COMPRESSION && compression_type != PIDX_CHUNKING_ONLY && compression_type != PIDX_CHUNKING_ZFP&& compression_type != PIDX_ZFP_COMPRESSION && compression_type != PIDX_CHUNKING_ZFP_63_COEFFICIENT && compression_type != PIDX_CHUNKING_AVERAGE)
    return PIDX_err_unsupported_compression_type;

  file->idx->compression_type = compression_type;

  if (file->idx->compression_type == PIDX_NO_COMPRESSION)
    return PIDX_success;
  else if (file->idx->compression_type == PIDX_CHUNKING_ONLY || file->idx->compression_type == PIDX_CHUNKING_ZFP || file->idx->compression_type == PIDX_CHUNKING_ZFP_63_COEFFICIENT || file->idx->compression_type == PIDX_CHUNKING_AVERAGE)
  {
    PIDX_point chunk_size;

    chunk_size[0] = 4;
    chunk_size[1] = 4;
    chunk_size[2] = 4;

    memcpy(file->idx->chunk_size, chunk_size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

    int reduce_by_sample = 1;
    unsigned long long total_chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
    if (reduce_by_sample == 1)
    {
      file->idx->bits_per_block = file->idx->bits_per_block - (int)log2(total_chunk_size);
      file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

      if (file->idx->bits_per_block <= 0)
      {
        file->idx->bits_per_block = 1;
        file->idx_d->samples_per_block = 2;
      }
    }
    else
      file->idx->blocks_per_file = file->idx->blocks_per_file / total_chunk_size;
  }

  return PIDX_success;
}

PIDX_return_code PIDX_set_zfp_precisison(PIDX_file file, float precisison)
{
  if(!file)
    return PIDX_err_file;

  if (file->idx->compression_type != PIDX_ZFP_COMPRESSION)
    return PIDX_err_file;

  file->idx->zfp_precisison = precisison;

  return PIDX_success;
}


PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type)
{
  if(!file)
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
  file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

  return PIDX_success;
}


PIDX_return_code PIDX_set_lossy_compression_bit_rate(PIDX_file file, float compression_bit_rate)
{
  if (!file)
    return PIDX_err_file;

  //if (file->idx->compression_type != PIDX_CHUNKING_ZFP)
  //  return PIDX_err_unsupported_compression_type;

  file->idx->compression_bit_rate = compression_bit_rate;

  if (file->idx->compression_bit_rate == 64)
  {
    file->idx->bits_per_block = file->idx->bits_per_block;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
    file->idx->compression_factor = 1;
  }
  if (file->idx->compression_bit_rate == 32)
  {
    file->idx->bits_per_block = file->idx->bits_per_block;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
    file->idx->compression_factor = 2;
  }
  if (file->idx->compression_bit_rate == 16)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 1;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 4;
  }
  if (file->idx->compression_bit_rate == 8)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 2;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 8;
  }
  if (file->idx->compression_bit_rate == 4)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 3;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 16;
  }
  if (file->idx->compression_bit_rate == 2)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 4;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 32;
  }

  // TODO : major hack here
  if (file->idx->compression_bit_rate == 1)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 5;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 64;
  }

  if (file->idx->compression_bit_rate == 0.5)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 6;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 128;
  }

  if (file->idx->compression_bit_rate == 0.25)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 7;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 256;
  }

  if (file->idx->compression_bit_rate == 0.125)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 8;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 512;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_get_lossy_compression_bit_rate(PIDX_file file, int *compression_bit_rate)
{
  if(!file)
    return PIDX_err_file;

  *compression_bit_rate = file->idx->compression_bit_rate;

  return PIDX_success;
}



PIDX_return_code PIDX_set_io_mode(PIDX_file file, int io_type)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx->io_type = io_type;

  return PIDX_success;
}



PIDX_return_code PIDX_get_io_mode(PIDX_file file, int* io_type)
{
  if(file == NULL)
    return PIDX_err_file;

  *io_type = file->idx->io_type;

  return PIDX_success;
}




PIDX_return_code PIDX_set_wavelet_level(PIDX_file file, int w_level)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->wavelet_levels = w_level;

  return PIDX_success;
}



PIDX_return_code PIDX_get_wavelet_level(PIDX_file file, int* w_level)
{
  if(file == NULL)
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



PIDX_return_code PIDX_set_raw_io_pipe_length(PIDX_file file, int pipe_length)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->raw_io_pipe_length = pipe_length;

  return PIDX_success;
}



PIDX_return_code PIDX_get_raw_io_pipe_length(PIDX_file file, int *pipe_length)
{
  if(file == NULL)
    return PIDX_err_file;

  *pipe_length = file->idx_d->raw_io_pipe_length;

  return PIDX_success;
}



PIDX_return_code PIDX_set_variable_pile_length(PIDX_file file, int var_pipe_length)
{
  if(!file)
    return PIDX_err_file;

  if (var_pipe_length < 0)
    return PIDX_err_size;

  file->idx_d->var_pipe_length = var_pipe_length;

  return PIDX_success;
}



PIDX_return_code PIDX_get_variable_pile_length(PIDX_file file, int* var_pipe_length)
{
  if(!file)
    return PIDX_err_file;

  if (var_pipe_length < 0)
    return PIDX_err_size;

  *var_pipe_length = file->idx_d->var_pipe_length;

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
  file->idx->endian = 1;

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



PIDX_return_code PIDX_get_comm(PIDX_file file, MPI_Comm *comm)
{
  if(!file)
    return PIDX_err_file;

   *comm = file->idx_c->global_comm;

  return PIDX_success;
}


PIDX_return_code PIDX_set_comm(PIDX_file file, MPI_Comm comm)
{
  if(!file)
    return PIDX_err_file;

   file->idx_c->global_comm = comm;
   file->idx_c->local_comm = comm;

   MPI_Comm_rank(file->idx_c->global_comm, &(file->idx_c->grank));
   MPI_Comm_size(file->idx_c->global_comm, &(file->idx_c->gnprocs));
   MPI_Comm_rank(file->idx_c->local_comm, &(file->idx_c->lrank));
   MPI_Comm_size(file->idx_c->local_comm, &(file->idx_c->lnprocs));

  return PIDX_success;
}



static PIDX_return_code PIDX_validate(PIDX_file file)
{
  unsigned long long dims;
  unsigned long long adjusted_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_bounds[0] = file->idx->bounds[0] / file->idx->chunk_size[0];
  adjusted_bounds[1] = file->idx->bounds[1] / file->idx->chunk_size[1];
  adjusted_bounds[2] = file->idx->bounds[2] / file->idx->chunk_size[2];

  if (PIDX_inner_product(&dims, adjusted_bounds))
    return PIDX_err_size;

  if (dims < file->idx_d->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx_d->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx->bits_per_block = getNumBits(file->idx_d->samples_per_block) - 1;
  }

  if (file->idx->bits_per_block == 0)
  {
    file->idx->bits_per_block = 1;
    file->idx_d->samples_per_block = 2;
  }

  // other validations...
  // TODO

  return PIDX_success;
}

