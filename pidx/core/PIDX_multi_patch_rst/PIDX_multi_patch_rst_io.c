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
 * \file PIDX_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"



PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_and_write(PIDX_multi_patch_rst_id rst_id)
{


  int v;
  char *directory_path;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  int g = 0;
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  for (g = 0; g < var0->patch_group_count; ++g)
  {
    //int bytes_per_value = var->bpv / 8;
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->grank, g);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0, v_end = 0;
    int start_var_index = rst_id->first_index;
    int end_var_index = rst_id->last_index + 1;
    for (v_start = start_var_index; v_start < end_var_index; v_start = v_start + (rst_id->idx_derived->raw_io_pipe_length + 1))
    {
      v_end = ((v_start + rst_id->idx_derived->raw_io_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (v_start + rst_id->idx_derived->raw_io_pipe_length);

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

      //PIDX_variable var = var_grp->variable[v];
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



PIDX_return_code PIDX_multi_patch_rst_buf_read_and_aggregate(PIDX_multi_patch_rst_id rst_id)
{
  int v;
  MPI_File fh;
  char *directory_path;
  char *data_set_path;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, rst_id->idx->current_time_step);

  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

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

      sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->grank, g);

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


PIDX_return_code PIDX_multi_patch_rst_buf_aggregated_write(PIDX_multi_patch_rst_id rst_id)
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

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->grank, g);
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


PIDX_return_code PIDX_multi_patch_rst_buf_aggregated_read(PIDX_multi_patch_rst_id rst_id)
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

    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->grank, g);
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
      ssize_t write_count = pread(fp, var_start->rst_patch_group[g]->reg_patch->buffer, buffer_size, data_offset);
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
