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

//#define PIDX_ACTIVE_TARGET 1

static void log_status(char* log_message, int step, int line_number, MPI_Comm comm);
static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab);
static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode);
static PIDX_return_code aggregate(PIDX_agg_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer ab, PIDX_block_layout lbl, int MODE, int layout_id);

PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl,  int MODE)
{
  if (create_window(id, ab) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
  char log_message[1024];
  sprintf(log_message, "[AGG Phase] Creat window for layout %d\n", layout_id);
  log_status(log_message, id->idx_c->color,  __LINE__, id->idx_c->partition_comm);

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
  sprintf(log_message, "[AGG Phase] One sided call for for layout %d\n", layout_id);
  log_status(log_message, id->idx_c->color,  __LINE__, id->idx_c->partition_comm);

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
  sprintf(log_message, "[AGG Phase] Comm free for layout %d\n", layout_id);
  log_status(log_message, id->idx_c->color,  __LINE__, id->idx_c->partition_comm);

  return PIDX_success;
}


static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab)
{
  int ret = 0;
  PIDX_variable var = id->idx->variable[ab->var_number];

  if (ab->buffer_size != 0)
  {
    int tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
    int bpdt = tcs * (var->bpv/8) / (id->idx->compression_factor);

    ret = MPI_Win_create(ab->buffer, ab->buffer_size, bpdt, MPI_INFO_NULL, id->idx_c->partition_comm, &(id->win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }
  else
  {
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, id->idx_c->partition_comm, &(id->win));
    if (ret != MPI_SUCCESS)
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
#ifdef PIDX_DUMP_AGG
      if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
      {
        fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
        fflush(agg_dump_fp);
      }
#endif

      for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
      {
        if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
        {
          index = 0;
          count = hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;

#ifdef PIDX_DUMP_AGG
          if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[%d]: ", i);
            fflush(agg_dump_fp);
          }
#endif
          ret = aggregate(id, v, hz_buf->start_hz_index[i], count, hz_buf->buffer[i], 0, ab, lbl, mode, layout_id);
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
#ifdef PIDX_DUMP_AGG
      if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
      {
        fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
        fflush(agg_dump_fp);
      }
#endif

      for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
      {
        if (var0->hz_buffer->nsamples_per_level[i][0] * var0->hz_buffer->nsamples_per_level[i][1] * var0->hz_buffer->nsamples_per_level[i][2] != 0)
        {
#ifdef PIDX_DUMP_AGG
          if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[%d]: ", i);
            fflush(agg_dump_fp);
          }
#endif
          int start_block_index = hz_buf->start_hz_index[i] / id->idx->samples_per_block;
          int end_block_index = hz_buf->end_hz_index[i] / id->idx->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          if (end_block_index == start_block_index)
          {
            count = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1);
            ret = aggregate(id, v, var0->hz_buffer->start_hz_index[i], count, hz_buf->buffer[i], 0, ab, lbl, mode, layout_id);
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
          else
          {
            int send_index = 0;
            int bl;
            for (bl = start_block_index; bl <= end_block_index; bl++)
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

                ret = aggregate(id, v, index + hz_buf->start_hz_index[i], count, hz_buf->buffer[i], send_index, ab, lbl, mode, layout_id);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_agg;
                }
                send_index = send_index + count;
              }
              else
                send_index = send_index + id->idx->samples_per_block;
            }
          }
        }
      }
    }
  }

  return PIDX_success;
}

static PIDX_return_code aggregate(PIDX_agg_id id, int variable_index, uint64_t hz_start, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer ab, PIDX_block_layout lbl, int MODE, int layout_id)
{
  int ret;
  int itr;
  int bpdt;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, vps;
  int target_rank = 0;
  uint64_t start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, samples_in_file = 0;
  uint64_t samples_per_file = (uint64_t) id->idx->samples_per_block * id->idx->blocks_per_file;
  uint64_t tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

  PIDX_variable var = id->idx->variable[variable_index];

  vps = var->vps; //number of samples for variable j

  // hz_start is the starting HZ index for the data buffer at level "level" and for regular patch number "patch"
  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / id->idx->samples_per_block;


  //number of empty blocks befor block "block_no" in the file "file_no"
  //negative_block_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, block_no, id->idx->variable[id->ini]->global_block_layout);
  negative_block_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, id->idx->bits_per_block, block_no, lbl);
  if (negative_block_offset < 0)
    return PIDX_err_agg;

  //number of samples in file "file_no"
  //samples_in_file = id->idx->variable[id->ini]->bcpf[file_no] * id->idx->samples_per_block;
  samples_in_file = lbl->bcpf[file_no] * id->idx->samples_per_block;
  if (samples_in_file > samples_per_file)
    return PIDX_err_agg;

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * id->idx->samples_per_block))) * vps)
      %
      (samples_in_file * vps);

  sample_index = target_disp / (samples_in_file / ab->agg_f);
  if (sample_index >= var->vps * ab->agg_f)
    return PIDX_err_agg;

  target_disp = target_disp % (samples_in_file / ab->agg_f);

  target_rank = id->agg_r[lbl->inverse_existing_file_index[file_no]][variable_index - id->fi][sample_index];

  //fprintf(stderr, "[F %d] [V %d] target rank %d\n", file_no, (variable_index - id->fi), target_rank);
  target_count = hz_count * vps;
  bpdt = ((var->bpv / 8) * tcs) / (id->idx->compression_factor);
  hz_buffer = hz_buffer + buffer_offset * bpdt * vps;

  start_agg_index = target_disp / (uint64_t) (samples_in_file / ab->agg_f);
  end_agg_index = ((target_disp + target_count - 1) / (uint64_t) (samples_in_file / ab->agg_f));

  if (start_agg_index != end_agg_index)
  {
    if (target_rank != id->idx_c->partition_rank)
    {
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , id->win);
#endif
      //target_disp_address = target_disp;
      if (MODE == PIDX_WRITE)
      {
        ret = MPI_Put(hz_buffer, ((samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        fprintf(stderr, "[1] TR %d\n", target_rank);
        ret = MPI_Get(hz_buffer, ((samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }

#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, id->win);
#endif
    }
    else
    {
      if (MODE == PIDX_WRITE)
      {
        memcpy( ab->buffer + target_disp * bpdt, hz_buffer, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt);
      }
      else
      {
        memcpy( hz_buffer, ab->buffer + target_disp * bpdt, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt);
      }
    }

    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++)
    {
      if (target_rank != id->idx_c->partition_rank)
      {
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank + ab->aggregator_interval, 0, id->win);
#endif
        if (MODE == PIDX_WRITE)
        {
          ret = MPI_Put(hz_buffer + (( (samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
        else
        {
          fprintf(stderr, "[2] TR %d\n", target_rank + ab->aggregator_interval);
          ret = MPI_Get(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(target_rank + ab->aggregator_interval, id->win);
#endif
      }
      else
      {
        if (MODE == PIDX_WRITE)
        {
          memcpy( ab->buffer, hz_buffer + (( (samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt);
        }
        else
        {
          memcpy( hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, ab->buffer, (samples_in_file / ab->agg_f) * bpdt);
        }
      }
    }

    if (target_rank + ab->aggregator_interval != id->idx_c->partition_rank)
    {
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank + ab->aggregator_interval, 0, id->win);
#endif
      if (MODE == PIDX_WRITE)
      {
        ret = MPI_Put(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / ab->agg_f))) + (((samples_in_file / ab->agg_f)) - target_disp))) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        fprintf(stderr, "[3] TR %d (%d %d)\n", target_rank + ab->aggregator_interval, target_rank, ab->aggregator_interval);
        ret = MPI_Get(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / ab->agg_f))) + (((samples_in_file / ab->agg_f)) - target_disp))) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank + ab->aggregator_interval, id->win);
#endif
    }
    else
    {
      if (MODE == PIDX_WRITE)
      {
        memcpy( ab->buffer, hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt);
      }
      else
      {
        memcpy( hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, ab->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt);
      }
    }
  }
  else
  {
    if (target_rank != id->idx_c->partition_rank)
    {
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , id->win);
#endif
      //target_disp_address = target_disp;
      if (MODE == PIDX_WRITE)
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        ret = MPI_Put(hz_buffer, hz_count * vps * bpdt, MPI_BYTE, target_rank, target_disp, hz_count * vps * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Get(hz_buffer, hz_count * vps * bpdt, MPI_BYTE, target_rank, target_disp, hz_count * vps * bpdt, MPI_BYTE, id->win);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, id->win);
#endif
    }
    else
    {
      if (MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( ab->buffer + target_disp * bpdt, hz_buffer, hz_count * vps * bpdt);
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( hz_buffer, ab->buffer + target_disp * bpdt, hz_count * vps * bpdt);
      }
    }
  }

  return PIDX_success;
}


static void log_status(char* log_message, int step, int line_number, MPI_Comm comm)
{
  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Barrier(comm);

  if (rank == 0)
    fprintf(stderr, "[nprocs %d] R%d Color %d [%d] Log message: %s", size, rank, step, line_number, log_message);

  return;
}
