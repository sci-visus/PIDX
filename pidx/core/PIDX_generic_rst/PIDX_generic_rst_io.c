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


#define INVERT_ENDIANESS 1
#include "../../PIDX_inc.h"

#if 1

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);


PIDX_return_code PIDX_generic_rst_buf_aggregate_and_write(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  int v;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, generic_rst_id->idx->filename, strlen(generic_rst_id->idx->filename) - 4);

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char time_template[512];
  sprintf(time_template, "%%s/%s/%%d_0", file->idx->filename_time_template);
  sprintf(file_name, time_template, directory_path, generic_rst_id->idx->current_time_step, generic_rst_id->idx_c->grank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  int v_start = 0, v_end = 0;
  int svi = generic_rst_id->first_index;
  int evi = generic_rst_id->last_index + 1;
  for (v_start = svi; v_start < evi; v_start = v_start + (generic_rst_id->idx_d->raw_io_pipe_length + 1))
  {
    v_end = ((v_start + generic_rst_id->idx_d->raw_io_pipe_length) >= (evi)) ? (evi - 1) : (v_start + generic_rst_id->idx_d->raw_io_pipe_length);

    // copy the size and offset to output
    PIDX_variable var_start = var_grp->variable[v_start];
    Ndim_patch_group patch_group = var_start->rst_patch_group;
    Ndim_patch out_patch = var_start->rst_patch_group->reg_patch;

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
    for (r = 0; r < var_start->rst_patch_group->count; r++)
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
              memcpy(reg_patch_buffer + data_offset + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group->patch[r]->buffer + send_o * var->vps * (var->bpv/8), send_c * var->vps * (var->bpv/8));
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


  free(directory_path);

  return PIDX_success;
}



PIDX_return_code PIDX_generic_rst_buf_read_and_aggregate(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  int v;
  MPI_File fh;
  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, generic_rst_id->idx->filename, strlen(generic_rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, generic_rst_id->idx->current_time_step);

  for (v = generic_rst_id->first_index; v <= generic_rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // copy the size and offset to output
    Ndim_patch_group patch_group = var->rst_patch_group;
    Ndim_patch out_patch = var->rst_patch_group->reg_patch;

    int nx = out_patch->size[0];
    int ny = out_patch->size[1];
    int nz = out_patch->size[2];

    var->rst_patch_group->reg_patch->buffer = malloc(nx * ny * nz * (var->bpv/8) * var->vps);
    memset(var->rst_patch_group->reg_patch->buffer, 0, nx * ny * nz * (var->bpv/8) * var->vps);

    if (var->rst_patch_group->reg_patch->buffer == NULL)
      return PIDX_err_chunk;

    int data_offset = 0, v1 = 0;
    for (v1 = 0; v1 < v; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var->vps * (var->bpv/8));
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_0", directory_path, generic_rst_id->idx->current_time_step, generic_rst_id->idx_c->grank);

    MPI_Status status;
    int ret = 0;
    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Line %d File %s File opening %s\n", __LINE__, __FILE__, file_name);
      return PIDX_err_rst;
    }

    ret = MPI_File_read_at(fh, data_offset, out_patch->buffer, (buffer_size), MPI_BYTE, &status);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    ret = MPI_File_close(&fh);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
    for (r = 0; r < var->rst_patch_group->count; r++)
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

            memcpy(var->rst_patch_group->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
          }
        }
      }
    }

    free(var->rst_patch_group->reg_patch->buffer);
    var->rst_patch_group->reg_patch->buffer = 0;

  }

  return PIDX_success;
}


PIDX_return_code PIDX_generic_rst_buf_aggregated_write(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count == 0)
      return PIDX_success;

  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, generic_rst_id->idx->filename, strlen(generic_rst_id->idx->filename) - 4);

  // loop through all groups
  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  sprintf(file_name, "%s/time%09d/%d_0", directory_path, generic_rst_id->idx->current_time_step, generic_rst_id->idx_c->grank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  int v_start = 0;
  int svi = generic_rst_id->first_index;
  int evi = generic_rst_id->last_index + 1;
  for (v_start = svi; v_start < evi; v_start = v_start + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = var_grp->variable[v_start];
    Ndim_patch out_patch = var_start->rst_patch_group->reg_patch;

    int bits = 0;
    PIDX_variable var = var_grp->variable[v_start];
    bits = (var->bpv/8) * var->vps;

    int data_offset = 0, v1 = 0;
    for (v1 = 0; v1 < v_start; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;

    if (generic_rst_id->idx->flip_endian == 1)
    {
      PIDX_variable curr_var = var_grp->variable[v_start];
      if (curr_var->bpv/8 == 4 || curr_var->bpv/8 == 12)
      {
        int y = 0;
        float temp;
        float temp2;

        for (y = 0; y < buffer_size / sizeof(float); y++)
        {
          memcpy(&temp, var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(float)), sizeof(float));
          bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(float)), &temp2, sizeof(float));
        }
      }
      else if (curr_var->bpv/8 == 8 || curr_var->bpv/8 == 24)
      {
        int y = 0;
        double temp;
        double temp2;

        for (y = 0; y < buffer_size / sizeof(double); y++)
        {
          memcpy(&temp, var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(double)), sizeof(double));
          bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(double)), &temp2, sizeof(double));
        }
      }
    }


    ssize_t write_count = pwrite(fp, var_start->rst_patch_group->reg_patch->buffer, buffer_size, data_offset);
    if (write_count != buffer_size)
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }
  close(fp);
  free(file_name);


  free(directory_path);

  return PIDX_success;
}


PIDX_return_code PIDX_generic_rst_buf_aggregated_read(PIDX_generic_rst_id generic_rst_id)
{
  PIDX_variable_group var_grp = generic_rst_id->idx->variable_grp[generic_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[generic_rst_id->first_index];

  // This process does not have any patch to process (after restructuring)
  if (var0->patch_group_count)
    return PIDX_success;

  int g = 0;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, generic_rst_id->idx->filename, strlen(generic_rst_id->idx->filename) - 4);

  // loop through all groups
  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  sprintf(file_name, "%s/time%09d/%d_%d", directory_path, generic_rst_id->idx->current_time_step, generic_rst_id->idx_c->grank, g);
  int fp = open(file_name, O_RDONLY, 0664);
  if (fp == -1)
  {
    fprintf(stderr, "[%s] [%d] open() failed while trying to open %s\n", __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }

  int v_start = 0;
  int svi = generic_rst_id->first_index;
  int evi = generic_rst_id->last_index + 1;
  for (v_start = svi; v_start < evi; v_start = v_start + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = var_grp->variable[v_start];
    Ndim_patch out_patch = var_start->rst_patch_group->reg_patch;

    int bits = 0;
    PIDX_variable var = var_grp->variable[v_start];
    bits = (var->bpv/8) * var->vps;

    int data_offset = 0, v1 = 0;
    for (v1 = 0; v1 < v_start; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
    ssize_t read_count = pread(fp, var_start->rst_patch_group->reg_patch->buffer, buffer_size, data_offset);
    if (read_count != buffer_size)
    {
      fprintf(stderr, "[%s] [%d] pwrite() %s offset %d failed %lld %lld\n", __FILE__, __LINE__, file_name, data_offset, (unsigned long long)read_count, (unsigned long long)buffer_size);
      return PIDX_err_io;
    }

    if (generic_rst_id->idx->flip_endian == 1)
    {
      PIDX_variable curr_var = var_grp->variable[v_start];
      if (curr_var->bpv/8 == 4 || curr_var->bpv/8 == 12)
      {
        int y = 0;
        float temp;
        float temp2;

        for (y = 0; y < buffer_size / sizeof(float); y++)
        {
          memcpy(&temp, var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(float)), sizeof(float));
          bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(float)), &temp2, sizeof(float));
        }
      }
      else if (curr_var->bpv/8 == 8 || curr_var->bpv/8 == 24)
      {
        int y = 0;
        double temp;
        double temp2;

        for (y = 0; y < buffer_size / sizeof(double); y++)
        {
          memcpy(&temp, var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(double)), sizeof(double));
          bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(var_start->rst_patch_group->reg_patch->buffer + (y * sizeof(double)), &temp2, sizeof(double));
        }
      }
    }
  }
  close(fp);
  free(file_name);


  free(directory_path);

  return PIDX_success;
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
#endif
