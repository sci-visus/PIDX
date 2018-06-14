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


// There are two modes of MPI One-sided RMA communication
// 1 - Active target communication (using fences)
// 2 - Passive target communication (using locks)

// My observation is that RMA works best with active target communication except for some machines
// like ash, this might be a hit and trial thing.
// Active target communication hangs on some machines (example ash), so, on those machines it is important to
// uncomment the following line
// #undef PIDX_ACTIVE_TARGET

#define PIDX_ACTIVE_TARGET

static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab);
static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode);
static int write_samples(PIDX_agg_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout, int mode);


// Perform aggregation
PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl,  int MODE)
{
  // Steps for aggregation
  // Step 1: Create one sided Windows
  // Step 2: RMA fence for synchronization - begin data transfer
  // Step 3: Transfer data (one sided)
  // Step 4: RMA fence for synchronization - end data transfer
  // Step 5: Free the MPI windows

  // Step 1
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

  // If I am an aggregator then create a window of size buffer size
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
  else // If I am not an aggregator then create an empty window
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
  int ret = 0;
  uint64_t index = 0, count = 0;

  PIDX_variable var0 = id->idx->variable[id->fi];

  if (var0->restructured_super_patch_count == 0)
    return PIDX_success;

  for (int v = id->fi; v <= id->li; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    index = 0, count = 0;
    HZ_buffer hz_buf = var->hz_buffer;

    // if a block is power in two dimension then iterate through all samples in a level at once
    if (hz_buf->is_boundary_HZ_buffer == power_two_block)
    {
      for (int i = lbl->resolution_from; i < lbl->resolution_to; i++)
      {
        if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
        {
          index = 0;
          count = hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;  // all samples in the hz level local to the process

          ret = write_samples(id, v, hz_buf->start_hz_index[i], count, hz_buf->buffer[i], 0, lbl, mode);
          if (ret != PIDX_success)
          {
            fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
    // if a block is non-power in two then iterte through samples in intervals of blocks
    else if (hz_buf->is_boundary_HZ_buffer == non_power_two_block)
    {
      for (int i = lbl->resolution_from; i < lbl->resolution_to; i++)
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
  int block_number, file_index, file_count, block_negative_offset = 0;
  uint64_t samples_per_file = id->idx->samples_per_block * id->idx->blocks_per_file;
  uint64_t data_offset = 0;

  PIDX_variable var = id->idx->variable[variable_index];

  int bytes_per_datatype = (var->bpv / 8) * var->vps * (id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2]) / id->idx->compression_factor;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype;

  // This while loop is redundant, it will only be executed once
  // this code is an exact replica of file per process io look at function write_samples PIDX_hz_encode_io.c
  // this is because file-per-process io a process is allowed to write to multiple files but here we restrict
  // one aggregator per file (the setup is already done in such a manner), look at function find_agg_level in
  // io_setup.c

  while (hz_count)
  {
    block_number = hz_start_index / id->idx->samples_per_block;
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


#ifndef PIDX_ACTIVE_TARGET
    MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , id->win);
#endif

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

#ifndef PIDX_ACTIVE_TARGET
    MPI_Win_unlock(target_rank, id->win);
#endif


    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * bytes_per_datatype;
  }

  return PIDX_success;
}
