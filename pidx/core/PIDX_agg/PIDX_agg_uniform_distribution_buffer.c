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

// distrubutes aggregators across rank space uniformly

PIDX_return_code PIDX_agg_buf_create_local_uniform_dist(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl)
{
  int rank_counter = 0;
  int aggregator_interval = id->idx_c->partition_nprocs / ((id->li - id->fi + 1) * lbl->efc);  // Distance between aggregators (in terms of number of processes)
  assert(aggregator_interval != 0);

  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  // loop through all the files in the particular aggregation group
  for (int k = 0; k < lbl->efc; k++)
  {
    // loop through all the variables (in the aggregation epoch)
    for (int i = id->fi; i <= id->li; i++)
    {
      id->agg_r[k][i - id->fi] = rank_counter;
      rank_counter = rank_counter + aggregator_interval;

      // if my rank is equal to the rank associated with file k and var number i, then I am the aggregator for file k variable i
      if (id->idx_c->partition_rank == id->agg_r[k][i - id->fi])
      {
        ab->file_number = lbl->existing_file_index[k];
        ab->var_number = i;

        int bpdt = (chunk_size * id->idx->variable[ab->var_number]->bpv/8) / (id->idx->compression_factor);
        uint64_t sample_count = lbl->bcpf[ab->file_number] * id->idx->samples_per_block;
        ab->buffer_size = sample_count * bpdt;

        ab->buffer = malloc(ab->buffer_size);
        memset(ab->buffer, 0, ab->buffer_size);
        if (ab->buffer == NULL)
        {
          fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) ab->buffer_size, __LINE__, __FILE__);
          return PIDX_err_agg;
        }

#if DEBUG_OUTPUT
        fprintf(stderr, "[File %d] [Variable %d] [Color %d] [n %d r %d p %d] Aggregator Partition rank %d\n", ab->file_number, ab->var_number, id->idx_c->color, id->idx_c->simulation_nprocs, id->idx_c->rnprocs, id->idx_c->partition_nprocs, id->idx_c->partition_rank);
#endif

      }
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_destroy(Agg_buffer ab)
{
  if (ab->buffer_size != 0)
  {
    free(ab->buffer);
    ab->buffer = 0;
  }

  return PIDX_success;
}
