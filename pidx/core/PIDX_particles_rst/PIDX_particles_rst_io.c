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


// Writes out the restructured data (super patch)
PIDX_return_code PIDX_particles_rst_buf_aggregated_write(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // If the process does not hold a super patch
  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx_metadata->filename, strlen(rst_id->idx_metadata->filename) - 4);

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  sprintf(file_name, "%s/time%09d/%d_0", directory_path, rst_id->idx_metadata->current_time_step, rst_id->idx_c->simulation_rank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  for (int v = rst_id->first_index; v < rst_id->last_index + 1; v = v + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = rst_id->idx_metadata->variable[v];
    PIDX_patch out_patch = var_start->restructured_super_patch->restructured_patch;

    int bits = 0;
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    bits = (var->bpv/CHAR_BIT) * var->vps;

    int data_offset = 0;
    for (int v1 = 0; v1 < v; v1++)
      data_offset = data_offset + (out_patch->particle_count * (rst_id->idx_metadata->variable[v1]->vps * (rst_id->idx_metadata->variable[v1]->bpv/CHAR_BIT)));

    // TODO reshuffle here and benchmark  <<<<

    int buffer_size =  out_patch->particle_count * bits;
    uint64_t write_count = pwrite(fp, var_start->restructured_super_patch->restructured_patch->buffer, buffer_size, data_offset);
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


PIDX_return_code PIDX_particles_rst_buf_aggregated_read(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx_metadata->filename, strlen(rst_id->idx_metadata->filename) - 4);

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  sprintf(file_name, "%s/time%09d/%d_0", directory_path, rst_id->idx_metadata->current_time_step, rst_id->idx_c->simulation_rank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  for (int v = rst_id->first_index; v < rst_id->last_index + 1; v = v + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = rst_id->idx_metadata->variable[v];
    PIDX_patch out_patch = var_start->restructured_super_patch->restructured_patch;

    int bits = 0;
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    bits = (var->bpv/CHAR_BIT) * var->vps;

    int data_offset = 0;
    for (int v1 = 0; v1 < v; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx_metadata->variable[v1]->vps * (rst_id->idx_metadata->variable[v1]->bpv/CHAR_BIT)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
    uint64_t write_count = pread(fp, var_start->restructured_super_patch->restructured_patch->buffer, buffer_size, data_offset);
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
