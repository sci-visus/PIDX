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


PIDX_return_code create_agg_io_buffer(PIDX_io file)
{
  uint32_t lc = 0;
  assert (file->idx_b->file0_agg_group_from_index == 0);

  uint32_t vc =  file->idx->variable_count;
  if (vc <= 0)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  idx_dataset idx = file->idx;

  file->agg_id = malloc(sizeof(*(file->agg_id)) * vc);
  file->io_id = malloc(sizeof(*(file->io_id)) * vc);
  memset(file->agg_id, 0, sizeof(*(file->agg_id)) * vc);
  memset(file->io_id, 0, sizeof(*(file->io_id)) * vc);

  idx->agg_buffer = malloc(sizeof(*(idx->agg_buffer)) * vc);
  memset(idx->agg_buffer, 0, sizeof(*(idx->agg_buffer)) * vc);

  for (uint32_t v = 0; v < vc; v++)
  {
    lc = (file->idx_b->agg_level - file->idx_b->file0_agg_group_from_index);
    file->agg_id[v] = malloc(sizeof(*(file->agg_id[v])) * lc);
    file->io_id[v] = malloc(sizeof(*(file->io_id[v])) * lc);
    memset(file->agg_id[v], 0, sizeof(*(file->agg_id[v])) * lc);
    memset(file->io_id[v], 0, sizeof(*(file->io_id[v])) * lc);

    idx->agg_buffer[v] = malloc(sizeof(*(idx->agg_buffer[v])) * lc);
    memset(idx->agg_buffer[v], 0, sizeof(*(idx->agg_buffer[v])) * lc);
  }

  return PIDX_success;
}


PIDX_return_code destroy_agg_io_buffer(PIDX_io file, int svi, int evi)
{
  uint32_t vc =  file->idx->variable_count;
  idx_dataset idx = file->idx;

  for (uint32_t v = 0; v < vc; v++)
  {
    free(file->agg_id[v]);
    free(file->io_id[v]);
    free(idx->agg_buffer[v]);
  }

  free(file->agg_id);
  free(file->io_id);
  free(idx->agg_buffer);

  return PIDX_success;
}
