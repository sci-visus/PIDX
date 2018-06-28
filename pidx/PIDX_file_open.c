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

/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_point dims, PIDX_file* file)
{
  int i;
  char file_name_skeleton[1024];

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;

  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  (*file)->flags = flags;

  (*file)->idx = (idx_dataset)malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));

  (*file)->idx_c = malloc(sizeof (*((*file)->idx_c)));
  memset((*file)->idx_c, 0, sizeof (*((*file)->idx_c)));

  (*file)->idx_b = malloc(sizeof (*((*file)->idx_b)));
  memset((*file)->idx_b, 0, sizeof (*((*file)->idx_b)));

  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));

  (*file)->time = malloc(sizeof (*((*file)->time)));
  memset((*file)->time, 0, sizeof (*((*file)->time)));
  (*file)->time->sim_start = PIDX_get_time();

  (*file)->meta_data_cache = malloc(sizeof (*((*file)->meta_data_cache)));
  memset((*file)->meta_data_cache, 0, sizeof (*((*file)->meta_data_cache)));

  (*file)->restructured_grid = malloc(sizeof(*(*file)->restructured_grid ));
  memset((*file)->restructured_grid , 0, sizeof(*(*file)->restructured_grid));

  (*file)->idx_c->simulation_comm = access_type->comm;
  (*file)->idx_c->partition_comm = access_type->comm;
  MPI_Comm_rank((*file)->idx_c->simulation_comm, &((*file)->idx_c->simulation_rank));
  MPI_Comm_size((*file)->idx_c->simulation_comm, &((*file)->idx_c->simulation_nprocs));
  MPI_Comm_rank((*file)->idx_c->partition_comm, &((*file)->idx_c->partition_rank));
  MPI_Comm_size((*file)->idx_c->partition_comm, &((*file)->idx_c->partition_nprocs));

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    (*file)->idx->partition_count[i] = 1;
    (*file)->idx->partition_offset[i] = 0;
  }

  (*file)->idx_dbg->debug_do_rst = 1;
  (*file)->idx_dbg->debug_do_chunk = 1;
  (*file)->idx_dbg->debug_do_compress = 1;
  (*file)->idx_dbg->debug_do_hz = 1;
  (*file)->idx_dbg->debug_do_agg = 1;
  (*file)->idx_dbg->debug_do_io = 1;
  (*file)->idx_dbg->debug_rst = 0;
  (*file)->idx_dbg->debug_hz = 0;

  (*file)->idx_b->reduced_resolution_factor = 0;

  (*file)->idx_c->color = 0;

  (*file)->idx->io_type = PIDX_IDX_IO;
  //(*file)->enable_raw_dump = 0;

  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx_dbg->enable_agg = 1;
  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';

  sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_partition, "%s_0.idx", file_name_skeleton);

#if 0
  if ((*file)->idx->partition_count[0] == 1 && (*file)->idx->partition_count[1] == 1 && (*file)->idx->partition_count[2] == 1)
    sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  else
    sprintf((*file)->idx->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_c->color);
#endif

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);

  (*file)->idx->compression_bit_rate = 64;
  (*file)->idx->compression_factor = 1;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);
  (*file)->idx->maxh = 0;
  (*file)->idx->max_file_count = 0;
  (*file)->fs_block_size = 0;
  (*file)->idx->pidx_version = 1;
  (*file)->idx_dbg->debug_file_output_state = PIDX_NO_META_DATA_DUMP;
  (*file)->idx->endian = PIDX_LITTLE_ENDIAN;

  int var = 0;

  if ((*file)->idx_c->simulation_rank == 0)
  {
    FILE *fp = fopen((*file)->idx->filename, "r");
    if (fp == NULL)
    {
      fprintf(stderr, "Error Opening %s\n", (*file)->idx->filename);
      return PIDX_err_file;
    }

    char line [512];

    while (fgets(line, sizeof (line), fp) != NULL)
    {
      line[strcspn(line, "\r\n")] = 0;

      // find the version number in the file
      if (strcmp(line, "(version)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        strncpy((*file)->idx->metadata_version, line, 8);
        break;
      }
    }

    // Parse the metadata file
    if ((*file)->idx->metadata_version > 0)
    {
      if (PIDX_metadata_parse(fp, file, (*file)->idx->metadata_version) != PIDX_success)
        return PIDX_err_metadata;
    }
    else
      return PIDX_err_metadata;

    fclose(fp);
  }
  
  MPI_Bcast((*file)->idx->bounds, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->box_bounds, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->restructured_grid->patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->chunk_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->endian), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->pidx_version), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->metadata_version), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->blocks_per_file), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->bits_per_block), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->variable_count), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->variable_count), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->bitSequence, 512, MPI_CHAR, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->partition_count, PIDX_MAX_DIMENSIONS, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->partition_size, PIDX_MAX_DIMENSIONS, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast((*file)->idx->partition_offset, PIDX_MAX_DIMENSIONS, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->compression_bit_rate), 1, MPI_FLOAT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->compression_type), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->idx->io_type), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  MPI_Bcast(&((*file)->fs_block_size), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
  
  //printf("reading version %d\n",(*file)->idx->metadata_version);
  if ((*file)->idx->io_type == PIDX_IDX_IO)
  {
    (*file)->idx->maxh = strlen((*file)->idx->bitSequence);
    for (uint32_t i = 0; i <= (*file)->idx->maxh; i++)
      (*file)->idx->bitPattern[i] = RegExBitmaskBit((*file)->idx->bitSequence, i);
  }

  //printf("reading version %d\n",(*file)->idx_d->metadata_version);

  if ((*file)->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    if ((*file)->idx->compression_bit_rate == 64)
      (*file)->idx->compression_factor = 1;
    if ((*file)->idx->compression_bit_rate == 32)
      (*file)->idx->compression_factor = 2;
    if ((*file)->idx->compression_bit_rate == 16)
      (*file)->idx->compression_factor = 4;
    if ((*file)->idx->compression_bit_rate == 8)
      (*file)->idx->compression_factor = 8;
    if ((*file)->idx->compression_bit_rate == 4)
      (*file)->idx->compression_factor = 16;
    if ((*file)->idx->compression_bit_rate == 2)
      (*file)->idx->compression_factor = 32;
    if ((*file)->idx->compression_bit_rate == 1)
      (*file)->idx->compression_factor = 64;
  }

  if ((*file)->idx->io_type != PIDX_RAW_IO)
    (*file)->idx->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);

  if ((*file)->idx_c->simulation_rank != 0)
  {
    for (var = 0; var < (*file)->idx->variable_count; var++)
    {
      (*file)->idx->variable[var] = malloc(sizeof (*((*file)->idx->variable[var])));
      memset((*file)->idx->variable[var], 0, sizeof (*((*file)->idx->variable[var])));
    }
  }

  for (var = 0; var < (*file)->idx->variable_count; var++)
  {
    MPI_Bcast(&((*file)->idx->variable[var]->bpv), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
    MPI_Bcast(&((*file)->idx->variable[var]->vps), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
    MPI_Bcast((*file)->idx->variable[var]->var_name, 512, MPI_CHAR, 0, (*file)->idx_c->simulation_comm);
    MPI_Bcast((*file)->idx->variable[var]->type_name, 512, MPI_CHAR, 0, (*file)->idx_c->simulation_comm);

    (*file)->idx->variable[var]->sim_patch_count = 0;
  }

#if 0
  if ((*file)->idx_c->simulation_rank == 0)
  {
    int ret;
    struct stat stat_buf;
    ret = stat((*file)->idx->filename, &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    (*file)->fs_block_size = stat_buf.st_blksize;
  }

  MPI_Bcast(&((*file)->fs_block_size), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);
#endif
  (*file)->idx->flip_endian = 0;

  unsigned int endian = 1;
  int current_endian = 0;
  char *c = (char*)&endian;
  if (*c)
    current_endian = 1;
  else
    current_endian = 0;

  if (current_endian == (*file)->idx->endian)
    (*file)->idx->flip_endian = 0;
  else
    (*file)->idx->flip_endian = 1;

  if (dims != NULL)
    memcpy(dims, (*file)->idx->bounds, (sizeof(uint64_t) * PIDX_MAX_DIMENSIONS));

  //if (physical_dims != NULL)
  //  memcpy(physical_dims, (*file)->idx->physical_bounds, (sizeof(double) * PIDX_MAX_DIMENSIONS));

  return PIDX_success;
}

// TODO merge serial_file_open and file_open

PIDX_return_code PIDX_serial_file_open(const char* filename, PIDX_flags flags, PIDX_point dims, PIDX_file* file)
{
  int i;
  char file_name_skeleton[1024];
  
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
  
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );
  
  (*file)->flags = flags;
  
  (*file)->idx = (idx_dataset)malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));
  
  (*file)->idx_c = malloc(sizeof (*((*file)->idx_c)));
  memset((*file)->idx_c, 0, sizeof (*((*file)->idx_c)));
  
  (*file)->idx_b = malloc(sizeof (*((*file)->idx_b)));
  memset((*file)->idx_b, 0, sizeof (*((*file)->idx_b)));
  
  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));
  
  (*file)->time = malloc(sizeof (*((*file)->time)));
  memset((*file)->time, 0, sizeof (*((*file)->time)));
  (*file)->time->sim_start = 0;
  
  (*file)->meta_data_cache = malloc(sizeof (*((*file)->meta_data_cache)));
  memset((*file)->meta_data_cache, 0, sizeof (*((*file)->meta_data_cache)));
  
  (*file)->restructured_grid = malloc(sizeof(*(*file)->restructured_grid ));
  memset((*file)->restructured_grid , 0, sizeof(*(*file)->restructured_grid));
  
//  (*file)->idx_c->simulation_comm = access_type->comm;
//  (*file)->idx_c->partition_comm = access_type->comm;
  
  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    (*file)->idx->partition_count[i] = 1;
    (*file)->idx->partition_offset[i] = 0;
  }
  
  (*file)->idx_dbg->debug_do_rst = 1;
  (*file)->idx_dbg->debug_do_chunk = 1;
  (*file)->idx_dbg->debug_do_compress = 1;
  (*file)->idx_dbg->debug_do_hz = 1;
  (*file)->idx_dbg->debug_do_agg = 1;
  (*file)->idx_dbg->debug_do_io = 1;
  (*file)->idx_dbg->debug_rst = 0;
  (*file)->idx_dbg->debug_hz = 0;
  
  (*file)->idx_b->reduced_resolution_factor = 0;
  
  (*file)->idx_c->color = 0;
  
  (*file)->idx->io_type = PIDX_IDX_IO;
  //(*file)->enable_raw_dump = 0;
  
  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx_dbg->enable_agg = 1;
  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;
  
  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';

  sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_partition, "%s_0.idx", file_name_skeleton);

#if 0
  if ((*file)->idx->partition_count[0] == 1 && (*file)->idx->partition_count[1] == 1 && (*file)->idx->partition_count[2] == 1)
    sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  else
    sprintf((*file)->idx->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_c->color);
#endif

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);

  (*file)->idx->compression_bit_rate = 64;
  (*file)->idx->compression_factor = 1;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);
  (*file)->idx->maxh = 0;
  (*file)->idx->max_file_count = 0;
  (*file)->fs_block_size = 0;

  (*file)->idx->pidx_version = 1;

  (*file)->idx_dbg->debug_file_output_state = PIDX_NO_META_DATA_DUMP;

  (*file)->idx->endian = PIDX_LITTLE_ENDIAN;

  int var = 0;
  char line [512];

  FILE *fp = fopen((*file)->idx->filename, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error Opening %s\n", (*file)->idx->filename);
    return PIDX_err_file;
  }

  while (fgets(line, sizeof (line), fp) != NULL)
  {
    line[strcspn(line, "\r\n")] = 0;

    // find the version number in the file
    if (strcmp(line, "(version)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      strncpy((*file)->idx->metadata_version, line, 8);

      break;
    }
  }

  // Parse the metadata file
  if ((*file)->idx->metadata_version > 0)
  {
    if (PIDX_metadata_parse(fp, file, (*file)->idx->metadata_version) != PIDX_success)
      return PIDX_err_metadata;
  }
  else
    return PIDX_err_metadata;

  fclose(fp);

  (*file)->idx->variable_count = (*file)->idx->variable_count;

  if ((*file)->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    if ((*file)->idx->compression_bit_rate == 64)
      (*file)->idx->compression_factor = 1;
    if ((*file)->idx->compression_bit_rate == 32)
      (*file)->idx->compression_factor = 2;
    if ((*file)->idx->compression_bit_rate == 16)
      (*file)->idx->compression_factor = 4;
    if ((*file)->idx->compression_bit_rate == 8)
      (*file)->idx->compression_factor = 8;
    if ((*file)->idx->compression_bit_rate == 4)
      (*file)->idx->compression_factor = 16;
    if ((*file)->idx->compression_bit_rate == 2)
      (*file)->idx->compression_factor = 32;
    if ((*file)->idx->compression_bit_rate == 1)
      (*file)->idx->compression_factor = 64;
  }


  if ((*file)->idx->io_type != PIDX_RAW_IO)
    (*file)->idx->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);

  for (var = 0; var < (*file)->idx->variable_count; var++)
    (*file)->idx->variable[var]->sim_patch_count = 0;

  (*file)->idx->flip_endian = 0;

  unsigned int endian = 1;
  int current_endian = 0;
  char *c = (char*)&endian;
  if (*c)
    current_endian = 1;
  else
    current_endian = 0;

  if (current_endian == (*file)->idx->endian)
    (*file)->idx->flip_endian = 0;
  else
    (*file)->idx->flip_endian = 1;

  memcpy(dims, (*file)->idx->bounds, (sizeof(uint64_t) * PIDX_MAX_DIMENSIONS));

  return PIDX_success;
}




PIDX_return_code PIDX_query_box(PIDX_file file, PIDX_point box_dims)
{
  if (!file)
    return PIDX_err_file;

  memcpy(file->idx->box_bounds, box_dims, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  return PIDX_success;
}
