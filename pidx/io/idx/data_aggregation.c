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



PIDX_return_code aggregation_setup(PIDX_io file, int svi, int evi)
{
  int ret = 0;
  idx_dataset idx = file->idx;
  PIDX_time time = file->time;

  assert(file->idx_b->file0_agg_group_from_index == 0);

  // aggregation happens in epochs of aggregation groups
  for (int j = file->idx_b->file0_agg_group_from_index; j < file->idx_b->agg_level; j++)
  {
    time->agg_init_start[svi][j] = PIDX_get_time();

    file->agg_id[svi][j] = PIDX_agg_init(file->idx, file->idx_c, file->idx_b, svi, evi);
    idx->agg_buffer[svi][j] = malloc(sizeof(*(idx->agg_buffer[svi][j])));
    memset(idx->agg_buffer[svi][j], 0, sizeof(*(idx->agg_buffer[svi][j])));

    idx->agg_buffer[svi][j]->file_number = -1;
    idx->agg_buffer[svi][j]->var_number = -1;

    time->agg_init_end[svi][j] = PIDX_get_time();

    time->agg_meta_start[svi][j] = PIDX_get_time();
    ret = PIDX_agg_meta_data_create(file->agg_id[svi][j], idx->agg_buffer[svi][j], file->idx_b->block_layout_by_agg_group[j]);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
    time->agg_meta_end[svi][j] = PIDX_get_time();

    time->agg_buf_start[svi][j] = PIDX_get_time();

    //ret = PIDX_agg_create_global_partition_localized_aggregation_buffer(file->agg_id[svi][j], idx->agg_buffer[svi][j], file->idx_b->block_layout_by_agg_group[j], j);
    //ret = PIDX_agg_create_local_partition_localized_aggregation_buffer(file->agg_id[svi][j], idx->agg_buffer[svi][j], file->idx_b->block_layout_by_agg_group[j], j);
    ret = PIDX_agg_buf_create_local_uniform_dist(file->agg_id[svi][j], idx->agg_buffer[svi][j], file->idx_b->block_layout_by_agg_group[j]);

    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
    time->agg_buf_end[svi][j] = PIDX_get_time();
  }

  return PIDX_success;
}



PIDX_return_code aggregation(PIDX_io file, int svi, int mode )
{
  PIDX_time time = file->time;
  assert(file->idx_b->file0_agg_group_from_index == 0);

  for (int j = file->idx_b->file0_agg_group_from_index; j < file->idx_b->agg_level; j++)
  {
    if (file->idx_dbg->debug_do_agg == 1)
    {
      time->agg_start[svi][j] = PIDX_get_time();
      if (PIDX_agg_global_and_local(file->agg_id[svi][j], file->idx->agg_buffer[svi][j], j, file->idx_b->block_layout_by_agg_group[j], mode) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
      time->agg_end[svi][j] = PIDX_get_time();
    }

    time->agg_meta_cleanup_start[svi][j] = PIDX_get_time();
    if (PIDX_agg_meta_data_destroy(file->agg_id[svi][j], file->idx_b->block_layout_by_agg_group[j]) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
    time->agg_meta_cleanup_end[svi][j] = PIDX_get_time();
  }

  return PIDX_success;
}



PIDX_return_code aggregation_cleanup(PIDX_io file, int start_index)
{
  for (uint32_t i = file->idx_b->file0_agg_group_from_index; i < file->idx_b->agg_level; i++)
  {
    uint32_t i_1 = i - file->idx_b->file0_agg_group_from_index;
    if (PIDX_agg_buf_destroy(file->idx->agg_buffer[start_index][i_1]) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx->agg_buffer[start_index][i_1]);
    PIDX_agg_finalize(file->agg_id[start_index][i_1]);
  }

  return PIDX_success;
}


