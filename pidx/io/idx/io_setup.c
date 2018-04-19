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


PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  file->hz_from_shared = 0;
  file->hz_to_shared = idx->total_partiton_level;

  file->hz_from_non_shared = idx->total_partiton_level;
  file->hz_to_non_shared =  idx->maxh;

  if (file->hz_from_shared == file->hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
  }

  if (file->hz_from_non_shared == file->hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
  }

  if (file->hz_from_shared == file->hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
    var_grp->shared_layout_count = 0;
  }
  else
  {
    var_grp->shared_start_layout_index = (file->hz_from_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->shared_start_layout_index <= 0)
      var_grp->shared_start_layout_index = 0;

    var_grp->shared_end_layout_index = (file->hz_to_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->shared_end_layout_index <= 0)
      var_grp->shared_end_layout_index = 1;

    var_grp->shared_layout_count = var_grp->shared_end_layout_index - var_grp->shared_start_layout_index;
  }

  if (file->hz_from_non_shared == file->hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
    var_grp->nshared_layout_count = 0;
  }
  else
  {
    var_grp->nshared_start_layout_index = (file->hz_from_non_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->nshared_start_layout_index <= 0)
      var_grp->nshared_start_layout_index = 0;

    var_grp->nshared_end_layout_index = (file->hz_to_non_shared - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (var_grp->nshared_end_layout_index <= 0)
      var_grp->nshared_end_layout_index = 1;

    var_grp->nshared_layout_count = var_grp->nshared_end_layout_index - var_grp->nshared_start_layout_index;
  }

  return PIDX_success;
}




PIDX_return_code find_agg_level(PIDX_io file, int gi, int svi, int evi)
{
  //fprintf(stderr, "svi and evi = %d and %d\n", svi, evi);
  int i = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  int total_aggregator = 0;
  int var_count = evi - svi;

#if 0
  if (file->idx->enable_agg == 0)
    var_grp->agg_level = var_grp->shared_start_layout_index;
  else
  {
    for (i = 0; i < var_grp->shared_layout_count + var_grp->nshared_layout_count ; i++)
    {
      no_of_aggregators = var_grp->block_layout_by_level[i]->efc;
      total_aggregator = total_aggregator + no_of_aggregators;
      if (no_of_aggregators <= file->idx_c->partition_nprocs)
        var_grp->agg_level = i + 1;
    }
  }

  if (total_aggregator > file->idx_c->partition_nprocs)
    var_grp->agg_level = var_grp->shared_start_layout_index;
#endif


  if (file->idx->enable_agg == 0)
  {
    var_grp->agg_level = var_grp->shared_start_layout_index;
    file->idx_d->variable_pipe_length = var_count - 1;
  }
  else
  {
    for (i = 0; i < var_grp->shared_layout_count + var_grp->nshared_layout_count ; i++)
      total_aggregator = total_aggregator + var_grp->block_layout_by_level[i]->efc;

    //fprintf(stderr, "npocs %d agg %d vc %d\n", file->idx_c->partition_nprocs, total_aggregator, var_count);
    if (file->idx_c->partition_nprocs >= total_aggregator * var_count)
    {
      var_grp->agg_level = var_grp->shared_layout_count + var_grp->nshared_layout_count;
      file->idx_d->variable_pipe_length = var_count - 1;

      //if (file->idx_c->partition_rank == 0)
      //  fprintf(stderr, "[A] agg level %d pipe length %d\n", var_grp->agg_level, file->idx_d->variable_pipe_length);
    }
    else
    {
      if (file->idx_c->partition_nprocs < total_aggregator)
      {
        var_grp->agg_level = var_grp->shared_start_layout_index;
        file->idx_d->variable_pipe_length = var_count - 1;

        //if (file->idx_c->partition_rank == 0)
        //  fprintf(stderr, "[B] agg level %d pipe length %d\n", var_grp->agg_level, file->idx_d->variable_pipe_length);
      }
      else
      {
        assert(var_count > 1);
        for (i = 0; i < var_count; i++)
        {
          if ((i + 1) * total_aggregator > file->idx_c->partition_nprocs)
            break;
        }
        file->idx_d->variable_pipe_length = i - 1;
        var_grp->agg_level = var_grp->shared_layout_count + var_grp->nshared_layout_count;
        //if (file->idx_c->partition_rank == 0)
        //  fprintf(stderr, "[C] agg level %d pipe length %d\n", var_grp->agg_level, file->idx_d->variable_pipe_length);
      }
    }
  }

  //if (file->idx_c->partition_rank == 0)
  //  fprintf(stderr, "agg level %d pipe length %d\n", var_grp->agg_level, file->idx_d->variable_pipe_length);

  return PIDX_success;
}
