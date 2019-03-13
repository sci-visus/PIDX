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


#include "../../PIDX_inc.h"


/// Allocate HZ buffer for every level
/// Buffer size at every hz level = end_hz_index[level] - start_hz_index[level] + 1
PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id)
{
  PIDX_variable var0 = id->idx->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  int c = 0, v = 0, bytes_for_datatype = 0;

  int maxH = id->idx->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  // Allocate actual HZ buffer for the variables
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    var->hz_buffer->buffer = (unsigned char**)malloc( maxH * sizeof (unsigned char*));
    memset(var->hz_buffer->buffer, 0,  maxH * sizeof (unsigned char*));

    bytes_for_datatype = ((var->bpv / 8) * chunk_size * var->vps) / id->idx->compression_factor;
    for (c = 0; c < maxH - id->resolution_to; c++)
    {
      uint64_t samples_per_level = (var->hz_buffer->end_hz_index[c] - var->hz_buffer->start_hz_index[c] + 1);

      var->hz_buffer->buffer[c] = malloc(bytes_for_datatype * samples_per_level);
      memset(var->hz_buffer->buffer[c], 0, bytes_for_datatype * samples_per_level);
    }
  }

  return PIDX_success;
}


/// Free the memory buffer for every hz level
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id)
{
  int itr = 0, v = 0;
  PIDX_variable var0 = id->idx->variable[id->first_index];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  if (id->idx->variable[id->first_index]->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->sim_patch_count not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  if (id->idx->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx->maxh (%d) not set.\n", __FILE__, __LINE__, id->idx->maxh);
    return PIDX_err_hz;
  }

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (itr = 0; itr < id->idx->maxh - id->resolution_to; itr++)
    {
      free(var->hz_buffer->buffer[itr]);
      var->hz_buffer->buffer[itr] = 0;
    }
  }

  return PIDX_success;
}
