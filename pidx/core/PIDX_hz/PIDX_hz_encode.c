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

#include "../../PIDX_inc.h"



PIDX_return_code PIDX_hz_encode_set_resolution(PIDX_hz_encode_id id, int resolution_from, int resolution_to)
{
  if (resolution_from < 0 || resolution_from < 0)
    return PIDX_err_hz;

  id->resolution_from = resolution_from;
  id->resolution_to = resolution_to;

  return PIDX_success;
}


PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_metadata_cache meta_data_cache, int fs_block_size, int first_index, int last_index)
{
  PIDX_hz_encode_id hz_id;
  hz_id = (PIDX_hz_encode_id)malloc(sizeof (*hz_id));
  memset(hz_id, 0, sizeof (*hz_id));

  hz_id->idx = idx_meta_data;
  hz_id->idx_c = idx_c;
  hz_id->idx_dbg = idx_dbg;
  hz_id->meta_data_cache = meta_data_cache;

  hz_id->first_index = first_index;
  hz_id->last_index = last_index;
  
  hz_id->fs_block_size = fs_block_size;

  hz_id->resolution_from = 0;
  hz_id->resolution_to = 0;

  return hz_id;
}



PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id)
{
  
  free(id);
  id = 0;
  
  return PIDX_success;
}

