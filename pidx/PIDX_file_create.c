/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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


/// Function to create IDX file descriptor (based on flags and access)
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_point dims, PIDX_file* file)
{
  if (flags != PIDX_MODE_CREATE && flags != PIDX_MODE_EXCL){
    fprintf(stderr,"[%s] [%d]\n", __FILE__, __LINE__);
    return PIDX_err_unsupported_flags;
  }

  if (flags == PIDX_MODE_EXCL)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0){
      fprintf(stderr,"[%s] [%d]\n", __FILE__, __LINE__);
      return PIDX_err_file_exists;
    }
  }

  uint64_t i = 0;
  int ret;
  char file_name_skeleton[PIDX_FILE_PATH_LENGTH];

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename){
    fprintf(stderr,"[%s] [%d]\n", __FILE__, __LINE__);
    return PIDX_err_name;
  }

  *file = malloc(sizeof (*(*file)));
  memset(*file, 0, sizeof (*(*file)));

  (*file)->idx = malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));

  (*file)->idx_b = malloc(sizeof (*((*file)->idx_b)));
  memset((*file)->idx_b, 0, sizeof (*((*file)->idx_b)));

  (*file)->time = malloc(sizeof (*((*file)->time)));
  memset((*file)->time, 0, sizeof (*((*file)->time)));

  (*file)->time->sim_start = PIDX_get_time();

  (*file)->restructured_grid = malloc(sizeof(*(*file)->restructured_grid ));
  memset((*file)->restructured_grid , 0, sizeof(*(*file)->restructured_grid));

  (*file)->idx_c = malloc(sizeof (*((*file)->idx_c)));
  memset((*file)->idx_c, 0, sizeof (*((*file)->idx_c)));

  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));

  (*file)->flags = flags;

  if (dims != NULL)
  {
    memcpy((*file)->idx->bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
    memcpy((*file)->idx->box_bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
  }

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);

  if (dims != NULL)
  {
    if (dims[0] * dims[1] * dims[2] < (*file)->idx->samples_per_block)
    {
      // ensure blocksize is a subset of the total volume.
      (*file)->idx->samples_per_block = getPowerOf2(dims[0] * dims[1] * dims[2]) >> 1;
      (*file)->idx->bits_per_block = getNumBits((*file)->idx->samples_per_block) - 1;
    }
  }
  //fprintf(stderr, "BPB %d SPB %d D %d\n", (*file)->idx->bits_per_block, (*file)->idx->samples_per_block, getPowerOf2(dims[0] * dims[1] * dims[2]));

  if ((*file)->idx->bits_per_block == 0)
  {
    (*file)->idx->bits_per_block = 0;
    (*file)->idx->samples_per_block = 1;
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    (*file)->idx->partition_count[i] = 1;
    if (dims != NULL)
      (*file)->idx->partition_size[i] = getPowerOf2(dims[i]);
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

  (*file)->idx_c->color = 0;
  (*file)->idx->io_type = PIDX_IDX_IO;

  (*file)->idx_b->reduced_resolution_factor = 0;

  (*file)->idx_c->simulation_comm = access_type->comm;
  (*file)->idx_c->partition_comm = access_type->comm;
  MPI_Comm_rank((*file)->idx_c->simulation_comm, &((*file)->idx_c->simulation_rank));
  MPI_Comm_size((*file)->idx_c->simulation_comm, &((*file)->idx_c->simulation_nprocs));
  MPI_Comm_rank((*file)->idx_c->partition_comm, &((*file)->idx_c->partition_rank));
  MPI_Comm_size((*file)->idx_c->partition_comm, &((*file)->idx_c->partition_nprocs));


  (*file)->idx_dbg->enable_agg = 1;

  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;

  (*file)->idx->pidx_version = 1;

  (*file)->idx->restructuring_factor[0] = 1;
  (*file)->idx->restructuring_factor[1] = 1;
  (*file)->idx->restructuring_factor[2] = 1;

  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';

  sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_partition, "%s_0.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_time_template, "time%%09d/");

  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);

  (*file)->idx->particle_regridding_factor = 1;
  (*file)->restructured_grid->patch_size[0] = -1;
  (*file)->restructured_grid->patch_size[1] = -1;
  (*file)->restructured_grid->patch_size[2] = -1;

  (*file)->idx->compression_factor = 1;
  (*file)->idx->compression_bit_rate = 64;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx->particles_position_variable_index = 0;
  (*file)->idx->particle_res_base = 32;
  (*file)->idx->particle_res_factor = 2;
  (*file)->idx->particle_number = 0;
  
  (*file)->idx->maxh = 0;
  (*file)->idx->max_file_count = 0;
  (*file)->fs_block_size = 0;

  (*file)->idx_dbg->debug_file_output_state = PIDX_NO_META_DATA_DUMP;

  if ((*file)->idx_c->simulation_rank == 0)
  {
    //TODO: close and delete the file (there is a way to do this automatically by fopen...)
    struct stat stat_buf;
    FILE *dummy = fopen(".dummy.txt", "w");
    fclose(dummy);
    ret = stat(".dummy.txt", &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if defined _MSC_VER
    (*file)->fs_block_size = 512; // TODO: double check
#else
    (*file)->fs_block_size = stat_buf.st_blksize;
#endif
  }

  MPI_Bcast(&((*file)->fs_block_size), 1, MPI_INT, 0, (*file)->idx_c->simulation_comm);

  (*file)->idx->flip_endian = 0;

  unsigned int endian = 1;
  char *c = (char*)&endian;
  if (*c)
    (*file)->idx->endian = PIDX_LITTLE_ENDIAN;
  else
    (*file)->idx->endian = PIDX_BIG_ENDIAN;

  return PIDX_success;
}
