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


#define INVERT_ENDIANESS 1
#include "../../PIDX_inc.h"
static int maximum_neighbor_count = 256;
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);


static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);


PIDX_return_code PIDX_rst_buf_aggregate_and_write(PIDX_rst_id rst_id)
{
  int v;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  int g = 0;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  for (g = 0; g < var0->patch_group_count; ++g)
  {
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->rank, g);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0, v_end = 0;
    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;
    for (v_start = svi; v_start < evi; v_start = v_start + (rst_id->idx_d->raw_io_pipe_length + 1))
    {
      v_end = ((v_start + rst_id->idx_d->raw_io_pipe_length) >= (evi)) ? (evi - 1) : (v_start + rst_id->idx_d->raw_io_pipe_length);

      // copy the size and offset to output
      PIDX_variable var_start = var_grp->variable[v_start];
      Ndim_patch_group patch_group = var_start->rst_patch_group[g];
      Ndim_patch out_patch = var_start->rst_patch_group[g]->reg_patch;

      int nx = out_patch->size[0];
      int ny = out_patch->size[1];
      int nz = out_patch->size[2];

      int bits = 0;
      for (v = v_start; v <= v_end; v++)
      {
        PIDX_variable var = var_grp->variable[v];
        bits = bits + (var->bpv/8) * var->vps;
      }

      unsigned char* reg_patch_buffer = malloc(nx * ny * nz * bits);
      memset(reg_patch_buffer, 0, nx * ny * nz * bits);
      if (reg_patch_buffer == NULL)
        return PIDX_err_chunk;

      int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
      for (r = 0; r < var_start->rst_patch_group[g]->count; r++)
      {
        for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
        {
          for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
          {
            for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
            {
              index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
              send_o = index;
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

              for (v = v_start; v <= v_end; v++)
              {
                int v1 = 0;
                int data_offset = 0;
                for (v1 = v_start; v1 < v; v1++)
                {
                  data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));
                }
                PIDX_variable var = var_grp->variable[v];
                memcpy(reg_patch_buffer + data_offset + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group[g]->patch[r]->buffer + send_o * var->vps * (var->bpv/8), send_c * var->vps * (var->bpv/8));
              }
            }
          }
        }
      }

      int data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

      int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
      ssize_t write_count = pwrite(fp, reg_patch_buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      free(reg_patch_buffer);
      reg_patch_buffer = 0;
    }
    close(fp);
    free(file_name);
  }

  free(directory_path);

  return PIDX_success;
}



PIDX_return_code PIDX_rst_buf_read_and_aggregate(PIDX_rst_id rst_id)
{
  int v;
  MPI_File fh;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, rst_id->idx->current_time_step);

  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      Ndim_patch_group patch_group = var->rst_patch_group[g];
      Ndim_patch out_patch = var->rst_patch_group[g]->reg_patch;

      int nx = out_patch->size[0];
      int ny = out_patch->size[1];
      int nz = out_patch->size[2];

      var->rst_patch_group[g]->reg_patch->buffer = malloc(nx * ny * nz * (var->bpv/8) * var->vps);
      memset(var->rst_patch_group[g]->reg_patch->buffer, 0, nx * ny * nz * (var->bpv/8) * var->vps);

      if (var->rst_patch_group[g]->reg_patch->buffer == NULL)
        return PIDX_err_chunk;

      int data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

      int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var->vps * (var->bpv/8));
      char *file_name;
      file_name = malloc(PATH_MAX * sizeof(*file_name));
      memset(file_name, 0, PATH_MAX * sizeof(*file_name));

      sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->rank, g);

      MPI_Status status;
      int ret = 0;
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stdout, "Line %d File %s File opening %s\n", __LINE__, __FILE__, file_name);
        return PIDX_err_rst;
      }

      ret = MPI_File_read_at(fh, data_offset, out_patch->buffer, (buffer_size), MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }

      ret = MPI_File_close(&fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }

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
              send_o = index * var->vps * (var->bpv/8);
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

              memcpy(var->rst_patch_group[g]->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
            }
          }
        }
      }

      free(var->rst_patch_group[g]->reg_patch->buffer);
      var->rst_patch_group[g]->reg_patch->buffer = 0;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_rst_buf_aggregated_write(PIDX_rst_id rst_id)
{
  int g = 0;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  for (g = 0; g < var0->patch_group_count; ++g)
  {
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->rank, g);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0;
    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;
    for (v_start = svi; v_start < evi; v_start = v_start + 1)
    {
      // copy the size and offset to output
      PIDX_variable var_start = var_grp->variable[v_start];
      Ndim_patch out_patch = var_start->rst_patch_group[g]->reg_patch;

      int bits = 0;
      PIDX_variable var = var_grp->variable[v_start];
      bits = (var->bpv/8) * var->vps;

      int data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

      int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
      ssize_t write_count = pwrite(fp, var_start->rst_patch_group[g]->reg_patch->buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
    close(fp);
    free(file_name);
  }

  free(directory_path);

  return PIDX_success;
}


PIDX_return_code PIDX_rst_buf_aggregated_read(PIDX_rst_id rst_id)
{
  int g = 0;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  for (g = 0; g < var0->patch_group_count; ++g)
  {
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->rank, g);
    int fp = open(file_name, O_RDONLY, 0664);
    if (fp == -1)
    {
      fprintf(stderr, "[%s] [%d] open() failed while trying to open %s\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    int v_start = 0;
    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;
    for (v_start = svi; v_start < evi; v_start = v_start + 1)
    {
      // copy the size and offset to output
      PIDX_variable var_start = var_grp->variable[v_start];
      Ndim_patch out_patch = var_start->rst_patch_group[g]->reg_patch;

      int bits = 0;
      PIDX_variable var = var_grp->variable[v_start];
      bits = (var->bpv/8) * var->vps;

      int data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

      int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
      ssize_t read_count = pread(fp, var_start->rst_patch_group[g]->reg_patch->buffer, buffer_size, data_offset);
      if (read_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() %s offset %d failed %lld %lld\n", __FILE__, __LINE__, file_name, data_offset, (unsigned long long)read_count, (unsigned long long)buffer_size);
        return PIDX_err_io;
      }
    }
    close(fp);
    free(file_name);
  }

  free(directory_path);

  return PIDX_success;
}




PIDX_return_code PIDX_rst_forced_raw_read(PIDX_rst_id rst_id)
{
  rst_id->idx_d->var_pipe_length = rst_id->idx->variable_count - 1;
  if (rst_id->idx_d->var_pipe_length == 0)
    rst_id->idx_d->var_pipe_length = 1;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  int svi = rst_id->first_index;
  int evi = rst_id->last_index;


  char *idx_directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  idx_directory_path = malloc(sizeof(*idx_directory_path) * PATH_MAX);
  memset(idx_directory_path, 0, sizeof(*idx_directory_path) * PATH_MAX);
  strncpy(idx_directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", idx_directory_path);
  sprintf(size_path, "%s_SIZE", idx_directory_path);
  free(idx_directory_path);

  uint32_t number_cores = 0;
  int fp = open(size_path, O_RDONLY);
  ssize_t write_count = pread(fp, &number_cores, sizeof(uint32_t), 0);
  if (write_count != sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

#if INVERT_ENDIANESS
  uint32_t temp_number_cores = 0;
  if (rst_id->idx->flip_endian == 1)
  {
    bit32_reverse_endian((unsigned char*)&number_cores, (unsigned char*)&temp_number_cores);
    number_cores = temp_number_cores;
  }
#endif

  uint32_t max_patch_count = 0;
  write_count = pread(fp, &max_patch_count, sizeof(uint32_t), sizeof(uint32_t));
  if (write_count != sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

#if INVERT_ENDIANESS
  uint32_t temp_max_patch_count = 0;
  if (rst_id->idx->flip_endian == 1)
  {
    bit32_reverse_endian((unsigned char*)&max_patch_count, (unsigned char*)&temp_max_patch_count);
    max_patch_count = temp_max_patch_count;
  }
#endif

  uint32_t *size_buffer = malloc((number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t));
  memset(size_buffer, 0, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t));

  write_count = pread(fp, size_buffer, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t), 2 * sizeof(uint32_t));
  if (write_count != (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  close(fp);

#if INVERT_ENDIANESS
  int buff_i;
  for(buff_i = 0; buff_i < (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)); buff_i++){
    uint32_t temp_value = 0;
    if (rst_id->idx->flip_endian == 1)
    {
      bit32_reverse_endian((unsigned char*)&size_buffer[buff_i], (unsigned char*)&temp_value);
      size_buffer[buff_i] = temp_value;
    }
  }
#endif

  uint32_t *offset_buffer = malloc((number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t));
  memset(offset_buffer, 0, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t));

  int fp1 = open(offset_path, O_RDONLY);
  write_count = pread(fp1, offset_buffer, (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t), 2 * sizeof(uint32_t));
  if (write_count != (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)) * sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }
  close(fp1);

#if INVERT_ENDIANESS
  for(buff_i = 0; buff_i<  (number_cores * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)); buff_i++){
    uint32_t temp_value = 0;
    if (rst_id->idx->flip_endian == 1)
    {
      bit32_reverse_endian((unsigned char*)&offset_buffer[buff_i], (unsigned char*)&temp_value);
      offset_buffer[buff_i] = temp_value;
    }
  }
#endif

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, rst_id->idx->current_time_step);

  int pc1 = 0;
  for (pc1 = 0; pc1 < var_grp->variable[svi]->sim_patch_count; pc1++)
  {

    int n = 0, m = 0, d = 0;
    Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = var_grp->variable[svi]->sim_patch[pc1]->offset[d];
      local_proc_patch->size[d] = var_grp->variable[svi]->sim_patch[pc1]->size[d];
    }

    // PC - - - - - - - - - -   PC - - - - - - - - - -      PC - - - - -
    // 0  1 2 3 4 5 6 7 8 9 10   11 12 13 14 15 16 17 18 19 20 21  22
    Ndim_patch n_proc_patch = (Ndim_patch)malloc(sizeof (*n_proc_patch));
    memset(n_proc_patch, 0, sizeof (*n_proc_patch));
    int p_counter = 1;
    int pc = 0, pc_index = 0;

    Ndim_patch_group patch_grp;
    patch_grp = malloc(sizeof(*(patch_grp)));
    memset(patch_grp, 0, sizeof(*(patch_grp)));

    int *source_patch_id = malloc(sizeof(int) * maximum_neighbor_count);
    memset(source_patch_id, 0, sizeof(int) * maximum_neighbor_count);

    patch_grp->source_patch_rank = malloc(sizeof(int) * maximum_neighbor_count);
    patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
    patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
    memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
    memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
    memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

    int patch_count = 0;
    for (n = 0; n < number_cores; n++)
    {
      pc = (int)offset_buffer[n * (max_patch_count * PIDX_MAX_DIMENSIONS + 1)];
      pc_index = n * (max_patch_count * PIDX_MAX_DIMENSIONS + 1);
      for (m = 0; m < pc; m++)
      {
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          n_proc_patch->offset[d] = (unsigned long long)offset_buffer[pc_index + m * PIDX_MAX_DIMENSIONS + d + 1];
          n_proc_patch->size[d] = (unsigned long long)size_buffer[pc_index + m * PIDX_MAX_DIMENSIONS + d + 1];
        }

        if (intersectNDChunk(local_proc_patch, n_proc_patch))
        {
          //sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, n, m);
          patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
          memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));
          for (d = 0; d < 3; d++)
          {
            //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
            if (n_proc_patch->offset[d] <= local_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) - local_proc_patch->offset[d] + 1;
            }
            else if (local_proc_patch->offset[d] <= n_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) >= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) - n_proc_patch->offset[d] + 1;
            }
            else if (( local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) <= (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) && local_proc_patch->offset[d] >= n_proc_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = local_proc_patch->size[d];
            }
            else if (( n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) && n_proc_patch->offset[d] >= local_proc_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = n_proc_patch->size[d];
            }
          }
          patch_grp->source_patch_rank[patch_count] = n;
          source_patch_id[patch_count] = m;
          //printf("patch_count %d -> %d %d %d\n", patch_count, patch_grp->patch[patch_count]->size[0], patch_grp->patch[patch_count]->size[1], patch_grp->patch[patch_count]->size[2]);

          patch_count++;
          if (patch_count >= maximum_neighbor_count)
          {
            maximum_neighbor_count = maximum_neighbor_count * 2;
            //printf("patch_count = %d maximum_neighbor_count = %d\n", patch_count, maximum_neighbor_count);

            int *temp_buffer1 = realloc(source_patch_id, sizeof(int) * maximum_neighbor_count);
            if (temp_buffer1 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              source_patch_id = temp_buffer1;

            int *temp_buffer2 = realloc(patch_grp->source_patch_rank, maximum_neighbor_count * sizeof(int));
            if (temp_buffer2 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->source_patch_rank = temp_buffer2;

            Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
            if (temp_buffer3 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->patch = temp_buffer3;
          }
        }
        p_counter++;
      }
    }

    free(local_proc_patch);
    local_proc_patch = 0;
    free(n_proc_patch);
    n_proc_patch = 0;

    unsigned long long k1 = 0, j1 = 0, i1 = 0, i = 0, j = 0;
    int count1 = 0, send_o = 0, send_c = 0, index = 0;
    unsigned long long sim_patch_offsetx[PIDX_MAX_DIMENSIONS];
    unsigned long long sim_patch_countx[PIDX_MAX_DIMENSIONS];

    unsigned char ***temp_patch_buffer;
    temp_patch_buffer = malloc(sizeof(*temp_patch_buffer) * (evi - svi + 1));
    memset(temp_patch_buffer, 0, sizeof(*temp_patch_buffer) * (evi - svi + 1));
    for (i = 0; i <= (evi - svi); i++)
    {
      PIDX_variable var = var_grp->variable[i + svi];
      temp_patch_buffer[i] = malloc(sizeof(*(temp_patch_buffer[i])) * patch_count);
      memset(temp_patch_buffer[i], 0, sizeof(*(temp_patch_buffer[i])) * patch_count);

      for (j = 0; j < patch_count; j++)
      {
        temp_patch_buffer[i][j] = malloc(sizeof(*(temp_patch_buffer[i][j])) * patch_grp->patch[j]->size[0] * patch_grp->patch[j]->size[1] * patch_grp->patch[j]->size[2] * var->bpv/8 * var->vps);
        memset(temp_patch_buffer[i][j], 0, sizeof(*(temp_patch_buffer[i][j])) * patch_grp->patch[j]->size[0] * patch_grp->patch[j]->size[1] * patch_grp->patch[j]->size[2] * var->bpv/8 * var->vps);
      }
    }

    for (i = 0; i < patch_count; i++)
    {
      sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, patch_grp->source_patch_rank[i], source_patch_id[i]);
      int fpx = open(file_name, O_RDONLY);

      pc_index = patch_grp->source_patch_rank[i] * (max_patch_count * PIDX_MAX_DIMENSIONS + 1);
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        sim_patch_offsetx[d] = (unsigned long long)offset_buffer[pc_index + source_patch_id[i] * PIDX_MAX_DIMENSIONS + d + 1];
        sim_patch_countx[d] = (unsigned long long)size_buffer[pc_index + source_patch_id[i] * PIDX_MAX_DIMENSIONS + d + 1];
      }

      count1 = 0;

      for (k1 = patch_grp->patch[i]->offset[2]; k1 < patch_grp->patch[i]->offset[2] + patch_grp->patch[i]->size[2]; k1++)
      {
        for (j1 = patch_grp->patch[i]->offset[1]; j1 < patch_grp->patch[i]->offset[1] + patch_grp->patch[i]->size[1]; j1++)
        {
          for (i1 = patch_grp->patch[i]->offset[0]; i1 < patch_grp->patch[i]->offset[0] + patch_grp->patch[i]->size[0]; i1 = i1 + patch_grp->patch[i]->size[0])
          {
            unsigned long long *sim_patch_offset = sim_patch_offsetx;// var_grp->variable[svi]->sim_patch[pc1]->offset;
            unsigned long long *sim_patch_count = sim_patch_countx;// var_grp->variable[svi]->sim_patch[pc1]->size;

            index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                (i1 - sim_patch_offset[0]);

            int start_index = 0, other_offset = 0, v1 = 0;
            for (start_index = svi; start_index <= evi; start_index = start_index + 1)
            {
              other_offset = 0;
              for (v1 = 0; v1 < start_index; v1++)
              {
                PIDX_variable var1 = var_grp->variable[v1];
                other_offset = other_offset + ((var1->bpv/8) * var1->vps * sim_patch_countx[0] * sim_patch_countx[1] * sim_patch_countx[2]);
              }

              PIDX_variable var = var_grp->variable[start_index];
              send_o = index * var->vps;
              send_c = patch_grp->patch[i]->size[0] * var->vps;

              size_t preadc = pread(fp, temp_patch_buffer[start_index - svi][i] + (count1 * send_c * var->bpv/8), send_c * var->bpv/8, other_offset + (send_o * var->bpv/8));

              if (preadc != send_c * var->bpv/8)
              {
                fprintf(stderr, "[%s] [%d] Error in pread [%d %d]\n", __FILE__, __LINE__, (int)preadc, (int)send_c * var->bpv/8);
                return PIDX_err_rst;
              }
            }
            count1++;
          }
        }
      }
      close(fpx);
    }

    int r, recv_o;
    for (r = 0; r < patch_count; r++)
    {
      for (k1 = patch_grp->patch[r]->offset[2]; k1 < patch_grp->patch[r]->offset[2] + patch_grp->patch[r]->size[2]; k1++)
      {
        for (j1 = patch_grp->patch[r]->offset[1]; j1 < patch_grp->patch[r]->offset[1] + patch_grp->patch[r]->size[1]; j1++)
        {
          for (i1 = patch_grp->patch[r]->offset[0]; i1 < patch_grp->patch[r]->offset[0] + patch_grp->patch[r]->size[0]; i1 = i1 + patch_grp->patch[r]->size[0])
          {
            index = ((patch_grp->patch[r]->size[0])* (patch_grp->patch[r]->size[1]) * (k1 - patch_grp->patch[r]->offset[2])) + ((patch_grp->patch[r]->size[0]) * (j1 - patch_grp->patch[r]->offset[1])) + (i1 - patch_grp->patch[r]->offset[0]);

            int start_index;
            for (start_index = svi; start_index <= evi; start_index = start_index + 1)
            {
              PIDX_variable var = var_grp->variable[start_index];

              send_o = index * var->vps * (var->bpv/8);
              send_c = (patch_grp->patch[r]->size[0]);
              recv_o = (var_grp->variable[start_index]->sim_patch[pc1]->size[0] * var_grp->variable[start_index]->sim_patch[pc1]->size[1] * (k1 - var_grp->variable[start_index]->sim_patch[pc1]->offset[2])) + (var_grp->variable[start_index]->sim_patch[pc1]->size[0] * (j1 - var_grp->variable[start_index]->sim_patch[pc1]->offset[1])) + (i1 - var_grp->variable[start_index]->sim_patch[pc1]->offset[0]);

              memcpy(var_grp->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)), temp_patch_buffer[start_index - svi][r] + send_o, send_c * var->vps * (var->bpv/8));

#if INVERT_ENDIANESS
              if (rst_id->idx->flip_endian == 1)
              {
                if (var->bpv/8 == 4 || var->bpv/8 == 12)
                {
                  int y = 0;
                  float temp;
                  float temp2;

                  for (y = 0; y < send_c * var->vps * (var->bpv/8) / sizeof(float); y++)
                  {
                    memcpy(&temp, var_grp->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(float)), sizeof(float));
                    bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
                    memcpy(var_grp->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(float)), &temp2, sizeof(float));
                  }
                }
                else if (var->bpv/8 == 8 || var->bpv/8 == 24)
                {
                  int y = 0;
                  double temp;
                  double temp2;

                  for (y = 0; y < send_c * var->vps * (var->bpv/8) / sizeof(double); y++)
                  {
                    memcpy(&temp, var_grp->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(double)), sizeof(double));
                    bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
                    memcpy(var_grp->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(double)), &temp2, sizeof(double));
                  }
                }
              }
#endif
            }
          }
        }
      }
      free(patch_grp->patch[r]);
    }

    for (i = 0; i <= (evi - svi); i++)
    {
      for (j = 0; j < patch_count; j++)
        free(temp_patch_buffer[i][j]);

      free(temp_patch_buffer[i]);
    }
    free(temp_patch_buffer);

    free(patch_grp->source_patch_rank);
    free(patch_grp->patch);
    free(patch_grp->reg_patch);
    free(patch_grp);
    free(source_patch_id);
  }

  free(file_name);
  free(data_set_path);
  free(directory_path);

  free(offset_buffer);
  free(size_buffer);

  return PIDX_success;
}




static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
 {
   int d = 0, check_bit = 0;
   for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
     check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

   return !(check_bit);
 }



#if INVERT_ENDIANESS

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
    unsigned char *data = ((unsigned char *)val) + 3;
    unsigned char *out = (unsigned char *)outbuf;

    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out = *data;

    return;
}

static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
    unsigned char *data = ((unsigned char *)val) + 7;
    unsigned char *out = (unsigned char *)outbuf;

    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out++ = *data--;
    *out = *data;

    return;
}

#endif
