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

#define PIDX_ACTIVE_TARGET 1

static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab);
static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode);
static int write_samples(PIDX_agg_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout, int mode);

PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl,  int MODE)
{
  if (create_window(id, ab) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

#ifdef PIDX_ACTIVE_TARGET
  if (MPI_Win_fence(0, id->win) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

  if (one_sided_data_com(id, ab, layout_id, lbl, MODE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

#ifdef PIDX_ACTIVE_TARGET
  if (MPI_Win_fence(0, id->win) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

  if (MPI_Win_free(&(id->win)) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

  return PIDX_success;
}


static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab)
{
  PIDX_variable var = id->idx->variable[ab->var_number];

  if (ab->buffer_size != 0)
  {
    int tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
    int bpdt = tcs * (var->bpv/8) / (id->idx->compression_factor);

    if (MPI_Win_create(ab->buffer, ab->buffer_size, bpdt, MPI_INFO_NULL, id->idx_c->partition_comm, &(id->win)) != MPI_SUCCESS)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }
  else
  {
    if (MPI_Win_create(0, 0, 1, MPI_INFO_NULL, id->idx_c->partition_comm, &(id->win)) != MPI_SUCCESS)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }

  return PIDX_success;
}



static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode)
{
  int i, v, ret = 0;
  uint64_t index = 0, count = 0;

  PIDX_variable var0 = id->idx->variable[id->fi];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (v = id->fi; v <= id->li; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    index = 0, count = 0;
    HZ_buffer hz_buf = var->hz_buffer;

    if (hz_buf->is_boundary_HZ_buffer == 1)
    {
      for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
      {
        if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
        {
          index = 0;
          count = hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;

          ret = write_samples(id, v, hz_buf->start_hz_index[i], count, hz_buf->buffer[i], 0, lbl, mode);
          if (ret != PIDX_success)
          {
            fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
    else if (hz_buf->is_boundary_HZ_buffer == 2)
    {
      for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
      {
        if (var0->hz_buffer->nsamples_per_level[i][0] * var0->hz_buffer->nsamples_per_level[i][1] * var0->hz_buffer->nsamples_per_level[i][2] != 0)
        {
          int start_block_index = hz_buf->start_hz_index[i] / id->idx->samples_per_block;
          int end_block_index = hz_buf->end_hz_index[i] / id->idx->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          // Block 0
          if (end_block_index == start_block_index)
          {
            count = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1);
            ret = write_samples(id, v, hz_buf->start_hz_index[i], count, hz_buf->buffer[i], 0, lbl, mode);
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
          // the remaining blocks
          else
          {
            int send_index = 0;
            for (int bl = start_block_index; bl <= end_block_index; bl++)
            {
              if (PIDX_blocks_is_block_present(bl, id->idx->bits_per_block, lbl))
              {
                if (bl == start_block_index)
                {
                  index = 0;
                  count = ((start_block_index + 1) * id->idx->samples_per_block) - hz_buf->start_hz_index[i];
                }
                else if (bl == end_block_index)
                {
                  index = (end_block_index * id->idx->samples_per_block - hz_buf->start_hz_index[i]);
                  count = hz_buf->end_hz_index[i] - ((end_block_index) * id->idx->samples_per_block) + 1;
                }
                else
                {
                  index = (bl * id->idx->samples_per_block - hz_buf->start_hz_index[i]);
                  count = id->idx->samples_per_block;
                }

                ret = write_samples(id, v, index + hz_buf->start_hz_index[i], count, hz_buf->buffer[i], send_index, lbl, mode);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_agg;
                }
                send_index = send_index + count;
              }
              else // if a block is skipped
                send_index = send_index + id->idx->samples_per_block;
            }
          }
        }
      }
    }
  }

  return PIDX_success;
}



static int write_samples(PIDX_agg_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout, int mode)
{
  int block_number, file_index, file_count, block_negative_offset = 0, file_number;
  uint64_t samples_per_file = id->idx->samples_per_block * id->idx->blocks_per_file;
  uint64_t data_offset = 0;

  PIDX_variable var = id->idx->variable[variable_index];

  int bytes_per_datatype = (var->bpv / 8) * var->vps * (id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2]) / id->idx->compression_factor;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype;

  while (hz_count)
  {
    block_number = hz_start_index / id->idx->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;

    assert(file_count > hz_count);

    if ((uint64_t)file_count > hz_count)
      file_count = hz_count;

    data_offset = 0;
    data_offset = file_index * bytes_per_datatype;

    // Adjusting for missing blocks
    block_negative_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, id->idx->bits_per_block, block_number, layout);
    data_offset -= block_negative_offset * id->idx->samples_per_block * bytes_per_datatype;

    uint64_t samples_per_file = (uint64_t) id->idx->samples_per_block * id->idx->blocks_per_file;
    int file_no = hz_start_index / samples_per_file;
    int target_rank = id->agg_r[layout->inverse_existing_file_index[file_no]][variable_index - id->fi];

    if (mode == PIDX_WRITE)
    {
      if (MPI_Put(hz_buffer, file_count * bytes_per_datatype, MPI_BYTE, target_rank, data_offset / bytes_per_datatype, file_count * bytes_per_datatype, MPI_BYTE, id->win) != MPI_SUCCESS)
      {
        fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_agg;
      }
    }
    else
    {
      if (MPI_Get(hz_buffer, file_count * bytes_per_datatype, MPI_BYTE, target_rank, data_offset / bytes_per_datatype, file_count * bytes_per_datatype, MPI_BYTE, id->win) != MPI_SUCCESS)
      {
        fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_agg;
      }
    }

    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * bytes_per_datatype;
  }

  return PIDX_success;
}
