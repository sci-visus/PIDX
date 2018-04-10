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
 * \file PIDX.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX.h
 *
 */

#include "PIDX_file_handler.h"


PIDX_return_code PIDX_variable_create(char* variable_name, unsigned int bits_per_sample, PIDX_data_type type_name, PIDX_variable* variable)
{
  if (!variable_name)
    return PIDX_err_name;

  if (bits_per_sample <= 0)
    return PIDX_err_size;

  if (!type_name)
    return PIDX_err_type;

  *variable = malloc(sizeof *(*variable));
  memset(*variable, 0, sizeof *(*variable));

  int bits = 0;
  PIDX_default_bits_per_datatype(type_name, &bits);
  if (bits !=0 && bits != bits_per_sample)
    return PIDX_err_type;

  (*variable)->vps = 1;
  (*variable)->bpv = (bits_per_sample/1);

  /*
  if (strcmp(type_name, FLOAT64)  == 0)
  {
    (*variable)->vps = 1;
    (*variable)->bpv = (bits_per_sample/1);
  }
  else
  {
    (*variable)->vps = 3;
    (*variable)->bpv = (bits_per_sample/3);
  }
  */

  strcpy((*variable)->type_name, type_name);
  strcpy((*variable)->var_name, variable_name);

  (*variable)->sim_patch_count = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_variable_write_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;

  const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  temp_buffer = read_from_this_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = (unsigned char*)temp_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  return PIDX_success;
}




PIDX_return_code PIDX_variable_write_particle_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, int number_of_particles, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;

  const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  variable->sim_patch[variable->sim_patch_count]->particle_count = number_of_particles;

  temp_buffer = read_from_this_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = (unsigned char*)temp_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  variable->is_particle = 1;

  return PIDX_success;
}




PIDX_return_code PIDX_append_and_write_variable(PIDX_file file, PIDX_variable variable)
{

  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];

  if (var_grp->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_variable;

  var_grp->variable[var_grp->variable_index_tracker] = variable;

  var_grp->variable_index_tracker++;
  var_grp->local_variable_count++;

  return PIDX_success;
}



PIDX_return_code PIDX_get_next_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];
  *variable = var_grp->variable[var_grp->variable_index_tracker];

  return PIDX_success;
}



PIDX_return_code PIDX_variable_read_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, void* write_to_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;

  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  variable->sim_patch[variable->sim_patch_count]->buffer = write_to_this_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  return PIDX_success;
}



PIDX_return_code PIDX_read_next_variable(PIDX_file file, PIDX_variable variable)
{
  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];

  if (var_grp->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_variable;

  variable = var_grp->variable[var_grp->variable_index_tracker];

  var_grp->variable_index_tracker++;
  var_grp->local_variable_count++;

  return PIDX_success;
}



PIDX_return_code PIDX_reset_variable_counter(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx->variable_grp[0]->variable_index_tracker = 0;
  file->idx->variable_grp[0]->local_variable_count = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int variable_index)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_index < 0)
    return PIDX_err_size;
  
  PIDX_variable_group var_grp = file->idx->variable_grp[0];
  if(var_grp->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  var_grp->variable_index_tracker = variable_index;
  var_grp->local_variable_count = 1;
  var_grp->local_variable_index = variable_index;

  return PIDX_success;
}



PIDX_return_code PIDX_set_current_variable_by_name(PIDX_file file, const char* variable_name)
{
  if(!file)
    return PIDX_err_file;

  int variable_index = -1;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];
  if(var_grp->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;

  int i = 0;
  for (i = 0; i < file->idx->variable_count; i++)
  {
    //fprintf(stderr, "str %s %s\n", variable_name, file->idx->variable_grp[0]->variable[i]->var_name);
    if (strcmp(variable_name, file->idx->variable_grp[0]->variable[i]->var_name) == 0)
    {
      variable_index = i;
      break;
    }
  }

  var_grp->variable_index_tracker = variable_index;
  var_grp->local_variable_count = 1;
  var_grp->local_variable_index = variable_index;

  if (variable_index == -1)
    return PIDX_err_variable;
  else
    return PIDX_success;
}



PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable)
{
  if(!file)
    return PIDX_err_file;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];

  if(var_grp->variable_index_tracker >= var_grp->variable_count)
    return PIDX_err_count;

  var_grp->variable[var_grp->variable_index_tracker] = variable;

  return PIDX_success;
}



PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  PIDX_variable_group var_grp = file->idx->variable_grp[0];

  if(var_grp->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  (*variable) = var_grp->variable[var_grp->variable_index_tracker];
  
  return PIDX_success;
}




PIDX_return_code PIDX_read_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_write_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  temp_buffer = read_from_this_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = (unsigned char*)temp_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;
  //variable->io_state = 1;

  PIDX_variable_group var_grp = file->idx->variable_grp[0];

  var_grp->variable[var_grp->variable_index_tracker] = variable;

  var_grp->variable_index_tracker++;
  var_grp->local_variable_count++;

  return PIDX_success;
}
