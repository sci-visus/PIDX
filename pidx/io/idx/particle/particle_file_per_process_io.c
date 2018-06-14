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

#include "../../../PIDX_inc.h"
static PIDX_return_code group_meta_data_init(PIDX_io file, int svi, int evi);
static PIDX_return_code PIDX_meta_data_write(PIDX_io file, int svi);
static PIDX_return_code PIDX_meta_data_read(PIDX_io file, int svi);

PIDX_return_code PIDX_particle_file_per_process_write(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;

  time->header_io_start = PIDX_get_time();
  if (group_meta_data_init(file, svi, evi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();

  time->particle_meta_data_io_start = MPI_Wtime();
  if (PIDX_meta_data_write(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->particle_meta_data_io_end = MPI_Wtime();

  time->particle_data_io_start = MPI_Wtime();
  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  PIDX_variable var0 = file->idx->variable[svi];
  for (int p = 0; p < var0->sim_patch_count; p++)
  {
    char file_name[PATH_MAX];
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));
    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, file->idx_c->simulation_rank, p);

    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    uint64_t data_offset = 0;
    for (int si = svi; si < evi; si++)
    {
      PIDX_variable var = file->idx->variable[si];

      int sample_count;
      int bits_per_sample;
      PIDX_get_datatype_details(var->type_name, &sample_count, &bits_per_sample);

      file->idx->variable_tracker[si] = 1;

#if 0
      fprintf(stderr, "[Inside PIDX] Rank %d has %d particles [%lld %lld %lld - %lld %lld %lld]\n", file->idx_c->simulation_rank, var->sim_patch[p]->particle_count, var->sim_patch[p]->offset[0], var->sim_patch[p]->offset[1], var->sim_patch[p]->offset[2], var->sim_patch[p]->size[0], var->sim_patch[p]->size[1], var->sim_patch[p]->size[2]);
      if (si == 0)
      {
        double particle_x = 0, particle_y, particle_z = 0;
        for (int pt = 0; pt < var->sim_patch[p]->particle_count; pt++)
        {
          memcpy(&particle_x, var->sim_patch[p]->buffer + (pt * sample_count + 0) * sizeof(double), sizeof (double));
          memcpy(&particle_y, var->sim_patch[p]->buffer + (pt * sample_count + 1) * sizeof(double), sizeof (double));
          memcpy(&particle_z, var->sim_patch[p]->buffer + (pt * sample_count + 2) * sizeof(double), sizeof (double));
          printf("%f\t%f\t%f\n", particle_x, particle_y, particle_z);
        }
      }
      printf("[%s] [%d] = %d %d %d\n", var->type_name, var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT), var->sim_patch[p]->particle_count, sample_count, (bits_per_sample/CHAR_BIT));
#endif

      uint64_t buffer_size = var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT);
      uint64_t write_count = pwrite(fp, var->sim_patch[p]->buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      data_offset = data_offset + buffer_size;
    }

    close(fp);
  }

  free (directory_path);
  time->particle_data_io_end = MPI_Wtime();

  return PIDX_success;
}



PIDX_return_code PIDX_particle_file_per_process_read(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;

  //if (group_meta_data_init(file, svi, evi) != PIDX_success)
  //{
  //  fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
  //  return PIDX_err_file;
  //}

  time->particle_meta_data_io_start = MPI_Wtime();
  if (PIDX_meta_data_read(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->particle_meta_data_io_end = MPI_Wtime();

  time->particle_data_io_start = MPI_Wtime();
  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  PIDX_variable var0 = file->idx->variable[svi];
  for (int p = 0; p < var0->sim_patch_count; p++)
  {
    char file_name[PATH_MAX];
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));
    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, file->idx_c->simulation_rank, p);

    int fp = open(file_name, O_RDONLY);

    uint64_t data_offset = 0;
    for (int si = svi; si < evi; si++)
    {
      PIDX_variable var = file->idx->variable[si];

      int sample_count;
      int bits_per_sample;
      PIDX_get_datatype_details(var->type_name, &sample_count, &bits_per_sample);

      file->idx->variable_tracker[si] = 1;

      uint64_t buffer_size = var0->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT);

      *(var->sim_patch[p]->read_particle_buffer) = malloc(buffer_size);
      uint64_t write_count = pread(fp, *var->sim_patch[p]->read_particle_buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      data_offset = data_offset + buffer_size;
    }

    close(fp);
  }

  free (directory_path);
  time->particle_data_io_end = MPI_Wtime();

  return PIDX_success;
}



/* Will:
 * The metadata file structure is as follows:
 * Header:
 *   global number of procs: f64
 *   max patch count on a rank: f64
 *
 * List of patch data, with nprocs * max_count entries
 * Each entry consists of:
 *   local patch count: f64
 *   box.lower: 3*f64
 *   box.upper: 3*f64
 *   local particle count: f64
 *
 * TODO: Wouldn't this be a lot easier to construct and write if we had
 * some structs for the types? Then we also wouldn't need to store the
 * counts as f64 (they should be u64).
 */
static PIDX_return_code PIDX_meta_data_write(PIDX_io file, int svi)
{
  double *global_patch;
  PIDX_variable var0 = file->idx->variable[svi];
  int max_patch_count;
  int patch_count = var0->sim_patch_count;
  // TODO WILL: Why the max patch count across all ranks? Ranks could have differing numbers of patches
  // right? So would we leave some unused space in the file here?
  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);

  double *local_patch = malloc(sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));
  memset(local_patch, 0, sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));

  local_patch[0] = (double)patch_count;
  for (int i = 0; i < patch_count; i++)
  {
    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * (2 * PIDX_MAX_DIMENSIONS + 1) + d + 1] = var0->sim_patch[i]->physical_offset[d];

    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * (2 * PIDX_MAX_DIMENSIONS + 1) + PIDX_MAX_DIMENSIONS + d + 1] = var0->sim_patch[i]->physical_size[d];

    local_patch[i * (2 * PIDX_MAX_DIMENSIONS + 1) + 2*PIDX_MAX_DIMENSIONS + 1] = var0->sim_patch[i]->particle_count;
  }

  //for (int i = 0; i < (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1); i++)
  //  printf("Local [%d] [np %d] ----> %f\n", i, file->idx_c->simulation_nprocs, local_patch[i]);

  global_patch = malloc((file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));
  memset(global_patch, 0,(file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));

  MPI_Allgather(local_patch, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, global_patch + 2, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, file->idx_c->simulation_comm);

  global_patch[0] = (double)file->idx_c->simulation_nprocs;
  global_patch[1] = (double)max_patch_count;

  char file_path[PATH_MAX];

  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(file_path, "%s_OFFSET_SIZE", directory_path);
  free(directory_path);
  if (file->idx_c->simulation_rank == 1 || file->idx_c->simulation_nprocs == 1)
  {
    int fp = open(file_path, O_CREAT | O_WRONLY, 0664);

    //for (int i = 0; i < (file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2); i++)
    //  printf("[%d] [np %d] ----> %f\n", i, file->idx_c->simulation_nprocs, global_patch[i]);

    uint64_t write_count = pwrite(fp, global_patch, (file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double), 0);
    if (write_count != (file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);
  }

  free(local_patch);
  free(global_patch);

  return PIDX_success;
}



static PIDX_return_code PIDX_meta_data_read(PIDX_io file, int svi)
{
  double *global_patch;
  PIDX_variable var0 = file->idx->variable[svi];

  int max_patch_count;
  int patch_count = var0->sim_patch_count;
  // TODO WILL: Why the max patch count across all ranks? Ranks could have differing numbers of patches
  // right? So would we leave some unused space in the file here?
  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);

  double *local_patch = malloc(sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));
  memset(local_patch, 0, sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));

  global_patch = malloc((file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));
  memset(global_patch, 0,(file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));

  char file_path[PATH_MAX];
  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(file_path, "%s_OFFSET_SIZE", directory_path);
  free(directory_path);
  if (file->idx_c->simulation_rank == 0 || file->idx_c->simulation_nprocs == 1)
  {
    int fp = open(file_path, O_RDONLY);
    uint64_t write_count = pread(fp, global_patch, (file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double), 0);
    if (write_count != (file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double))
    {
      fprintf(stderr, "[%s] [%d] pread() failed.  %d != %d\n", __FILE__, __LINE__, (int)write_count, (int)((file->idx_c->simulation_nprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double)) );
      return PIDX_err_io;
    }
    close(fp);
  }


  MPI_Scatter(global_patch + 2, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, local_patch, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, 0, file->idx_c->simulation_comm);

  global_patch[0] = (double)file->idx_c->simulation_nprocs;
  global_patch[1] = (double)max_patch_count;

  int pcounter = 0;
  local_patch[0] = (double)patch_count;
  for (int i = 0; i < patch_count; i++)
  {
    //for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    //  local_patch[i * PIDX_MAX_DIMENSIONS + d + 1] = var0->sim_patch[i]->physical_offset[d];

    //for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    //  local_patch[i * PIDX_MAX_DIMENSIONS + PIDX_MAX_DIMENSIONS + d + 1] = var0->sim_patch[i]->physical_size[d];

    var0->sim_patch[i]->particle_count = local_patch[i * PIDX_MAX_DIMENSIONS + 2*PIDX_MAX_DIMENSIONS + 1];
    //printf("Particle count = %d\n", var0->sim_patch[i]->particle_count);

    pcounter++;
  }

  free(local_patch);
  free(global_patch);

  return PIDX_success;
}



static PIDX_return_code group_meta_data_init(PIDX_io file, int svi, int evi)
{
  // Creates the file heirarchy and writes the header info for all binary files
  if (init_raw_headers_layout(file, svi, evi, file->idx->filename) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}
