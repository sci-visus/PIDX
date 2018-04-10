/*
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


PIDX_return_code create_agg_io_buffer(PIDX_io file, int gi, int svi, int evi)
{
  int lc = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  assert (var_grp->shared_start_layout_index == 0);

  int vc =  file->idx->variable_count;// (evi - svi);
  if (vc <= 0)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  idx_dataset_derived_metadata idx = file->idx_d;

  file->agg_id = malloc(sizeof(*(file->agg_id)) * vc);
  file->io_id = malloc(sizeof(*(file->io_id)) * vc);
  memset(file->agg_id, 0, sizeof(*(file->agg_id)) * vc);
  memset(file->io_id, 0, sizeof(*(file->io_id)) * vc);

  idx->agg_buffer = malloc(sizeof(*(idx->agg_buffer)) * vc);
  memset(idx->agg_buffer, 0, sizeof(*(idx->agg_buffer)) * vc);

  int v = 0;
  for (v = 0; v < vc; v++)
  {
    lc = (var_grp->agg_level - var_grp->shared_start_layout_index);
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
  //int vc = evi - svi;
  int vc =  file->idx->variable_count;// (evi - svi);
  idx_dataset_derived_metadata idx = file->idx_d;

  int v = 0;
  for (v = 0; v < vc; v++)
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



PIDX_return_code finalize_aggregation(PIDX_io file, int gi, int start_index)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  int ret;
  int i = 0;
  int i_1 = 0;
  //start_index = start_index - local_var_index;

  int sli = var_grp->shared_start_layout_index;
  int agg_i = var_grp->agg_level;

  //fprintf(stderr, "[%d] sli and agg_i %d %d si %d\n", file->idx_c->grank, sli, agg_i, start_index);
  for (i = sli; i < agg_i; i++)
  {
    i_1 = i - sli;
    ret = PIDX_agg_buf_destroy(file->idx_d->agg_buffer[start_index][i_1]);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx_d->agg_buffer[start_index][i_1]);
    PIDX_agg_finalize(file->agg_id[start_index][i_1]);
    PIDX_file_io_finalize(file->io_id[start_index][i_1]);
  }

  return PIDX_success;
}
