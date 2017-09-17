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


// Writes out the restructured data (super patch)
PIDX_return_code PIDX_idx_rst_buf_aggregated_write(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

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

  sprintf(file_name, "%s/time%09d/%d_0", directory_path, rst_id->idx_metadata->current_time_step, rst_id->idx_comm_metadata->grank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  int v = 0;
  int svi = rst_id->first_index;
  int evi = rst_id->last_index + 1;
  for (v = svi; v < evi; v = v + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = var_grp->variable[v];
    PIDX_patch out_patch = var_start->restructured_super_patch->restructured_patch;

    int bits = 0;
    PIDX_variable var = var_grp->variable[v];
    bits = (var->bpv/8) * var->vps;

    int data_offset = 0, v1 = 0;
    for (v1 = 0; v1 < v; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
    ssize_t write_count = pwrite(fp, var_start->restructured_super_patch->restructured_patch->buffer, buffer_size, data_offset);
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


PIDX_return_code PIDX_idx_rst_buf_aggregated_read(PIDX_idx_rst_id rst_id)
{

  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx_metadata->filename, strlen(rst_id->idx_metadata->filename) - 4);


  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  sprintf(file_name, "%s/time%09d/%d_0", directory_path, rst_id->idx_metadata->current_time_step, rst_id->idx_comm_metadata->grank);
  int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

  int v = 0;
  int svi = rst_id->first_index;
  int evi = rst_id->last_index + 1;
  for (v = svi; v < evi; v = v + 1)
  {
    // copy the size and offset to output
    PIDX_variable var_start = var_grp->variable[v];
    PIDX_patch out_patch = var_start->restructured_super_patch->restructured_patch;

    int bits = 0;
    PIDX_variable var = var_grp->variable[v];
    bits = (var->bpv/8) * var->vps;

    int data_offset = 0, v1 = 0;
    for (v1 = 0; v1 < v; v1++)
      data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var_grp->variable[v1]->vps * (var_grp->variable[v1]->bpv/8)));

    int buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
    ssize_t write_count = pread(fp, var_start->restructured_super_patch->restructured_patch->buffer, buffer_size, data_offset);
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
