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


PIDX_return_code select_io_mode(PIDX_io file)
{
  // HZ level range for file 0
  file->idx_b->hz_file0_from = 0;       // This is set to 0, but it can be non-zero x, if we want to skip writing the first x hz levels
  file->idx_b->hz_file0_to = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1;

  // aggregation group range for file 0 [0, 1)
  file->idx_b->file0_agg_group_from_index = (file->idx_b->hz_file0_from - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_b->file0_agg_group_from_index <= 0)
    file->idx_b->file0_agg_group_from_index = 0;

  file->idx_b->file0_agg_group_to_index = (file->idx_b->hz_file0_to - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_b->file0_agg_group_to_index <= 0)
    file->idx_b->file0_agg_group_to_index = 1;

  file->idx_b->file0_agg_group_count = file->idx_b->file0_agg_group_to_index - file->idx_b->file0_agg_group_from_index;


  // HZ level range for all other files
  file->idx_b->hz_n_file0_from = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1;
  file->idx_b->hz_n_file0_to =  file->idx->maxh;

  // aggregation group range for all other files [1, n)
  if (file->idx_b->hz_n_file0_from == file->idx_b->hz_n_file0_to)
  {
    file->idx_b->nfile0_agg_group_from_index = 0;
    file->idx_b->nfile0_agg_group_to_index = 0;
    file->idx_b->nfile0_agg_group_count = 0;
  }
  else
  {
    file->idx_b->nfile0_agg_group_from_index = (file->idx_b->hz_n_file0_from - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (file->idx_b->nfile0_agg_group_from_index <= 0)
      file->idx_b->nfile0_agg_group_from_index = 0;

    file->idx_b->nfile0_agg_group_to_index = (file->idx_b->hz_n_file0_to - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
    if (file->idx_b->nfile0_agg_group_to_index <= 0)
      file->idx_b->nfile0_agg_group_to_index = 1;

    file->idx_b->nfile0_agg_group_count = file->idx_b->nfile0_agg_group_to_index - file->idx_b->nfile0_agg_group_from_index;
  }

  return PIDX_success;
}



PIDX_return_code find_agg_level(PIDX_io file, int svi, int evi)
{
  int total_aggregator = 0;
  int var_count = evi - svi;

  if (file->idx_dbg->enable_agg == 0)
  {
    file->idx_b->agg_level = file->idx_b->file0_agg_group_from_index;
    file->idx->variable_pipe_length = var_count - 1;
  }
  else
  {
    for (uint32_t i = 0; i < file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count ; i++)
      total_aggregator = total_aggregator + file->idx_b->block_layout_by_agg_group[i]->efc;

    if (file->idx_c->partition_nprocs >= total_aggregator * var_count)
    {
      file->idx_b->agg_level = file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count;
      file->idx->variable_pipe_length = var_count - 1;

      //if (file->idx_c->partition_rank == 0)
      //  fprintf(stderr, "[A] agg level %d pipe length %d\n", file->idx_b->agg_level, file->idx->variable_pipe_length);
    }
    else
    {
      if (file->idx_c->partition_nprocs < total_aggregator)
      {
        file->idx_b->agg_level = file->idx_b->file0_agg_group_from_index;
        file->idx->variable_pipe_length = var_count - 1;

        //if (file->idx_c->partition_rank == 0)
        //  fprintf(stderr, "[B] agg level %d pipe length %d\n", file->idx_b->agg_level, file->idx->variable_pipe_length);
      }
      else
      {
        assert(var_count > 1);
        uint32_t i;
        for (i = 0; i < var_count; i++)
        {
          if ((i + 1) * total_aggregator > file->idx_c->partition_nprocs)
            break;
        }
        file->idx->variable_pipe_length = i - 1;
        file->idx_b->agg_level = file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count;
        //if (file->idx_c->partition_rank == 0)
        //  fprintf(stderr, "[C] agg level %d pipe length %d\n", file->idx_b->agg_level, file->idx->variable_pipe_length);
      }
    }
  }

  //if (file->idx_c->partition_rank == 0)
  //  fprintf(stderr, "agg level %d pipe length %d\n", file->idx_b->agg_level, file->idx->variable_pipe_length);

  return PIDX_success;
}
