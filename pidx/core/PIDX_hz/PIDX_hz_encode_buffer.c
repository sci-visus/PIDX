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


#include "../../PIDX_inc.h"


PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id)
{

  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  PIDX_variable var0 = var_grp->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int c = 0, v = 0, bytes_for_datatype = 0;

  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  // Allocate actual HZ buffer for the variables
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    var->hz_buffer->buffer = (unsigned char**)malloc( maxH * sizeof (unsigned char*));
    memset(var->hz_buffer->buffer, 0,  maxH * sizeof (unsigned char*));

    bytes_for_datatype = ((var->bpv / 8) * chunk_size * var->vps) / id->idx->compression_factor;
    for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
    {
      unsigned long long samples_per_level = (var->hz_buffer->end_hz_index[c] - var->hz_buffer->start_hz_index[c] + 1);

      //if (id->idx_c->grank == 1)
      //  fprintf(stderr, "[V %d] [C %d] %lld\n", v, c, samples_per_level);

      var->hz_buffer->buffer[c] = malloc(bytes_for_datatype * samples_per_level);
      memset(var->hz_buffer->buffer[c], 0, bytes_for_datatype * samples_per_level);
    }
  }

  return PIDX_success;
}


/// tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id)
{
  int itr = 0, v = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  PIDX_variable var0 = var_grp->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  if(var_grp->variable[id->first_index]->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->sim_patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_d->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->maxh (%d) not set.\n", __FILE__, __LINE__, id->idx_d->maxh);
    return 1;
  }

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (itr = id->resolution_from; itr < id->idx_d->maxh - id->resolution_to; itr++)
    {
      free(var->hz_buffer->buffer[itr]);
      var->hz_buffer->buffer[itr] = 0;
    }
  }

  return PIDX_success;
}
