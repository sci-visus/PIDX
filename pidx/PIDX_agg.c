/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

#include "PIDX_inc.h"
#include <sys/time.h>

#define PIDX_ACTIVE_TARGET

#define PIDX_DUMP_AGG

#ifdef PIDX_DUMP_AGG
static FILE* agg_dump_fp;
#endif

struct PIDX_agg_struct
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Comm global_comm;
  MPI_Win win;
#endif

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;
  

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int init_index;
  int first_index;
  int last_index;

  int aggregator_interval;
};

enum IO_MODE {PIDX_READ, PIDX_WRITE};

PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  agg_id->idx = idx_meta_data;
  agg_id->idx_d = idx_d;

  agg_id->init_index = init_index;
  agg_id->first_index = first_index;
  agg_id->last_index = last_index;

  return agg_id;
}

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_agg_set_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->comm = comm;

  return PIDX_success;
}

PIDX_return_code PIDX_agg_set_global_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->global_comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, int MODE)
{
  int rank = 0, itr;// nrank = 0;
  int bytes_per_datatype;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, values_per_sample;
  int target_rank = 0;
  int64_t start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, hz_start = 0, samples_in_file = 0;
  int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  //MPI_Aint target_disp_address;

  int64_t total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

#if PIDX_HAVE_MPI
  int ret;
  MPI_Comm_rank(agg_id->comm, &rank);
  //MPI_Comm_rank(agg_id->global_comm, &nrank);
#endif

  values_per_sample = agg_id->idx->variable[variable_index]->values_per_sample; //number of samples for variable j

  //starting HZ index for the data buffer at level "level" and for regular patch number "patch"
  hz_start = hz_start_index;

  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / agg_id->idx_d->samples_per_block;

  //number of empty blocks befor block "block_no" in the file "file_no"
  negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, agg_id->idx->variable[agg_id->init_index]->global_block_layout);
  assert(negative_block_offset >= 0);

  //number of samples in file "file_no"
  samples_in_file = agg_id->idx->variable[agg_id->init_index]->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;
  assert(samples_in_file <= samples_per_file);

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  assert(hz_start >= (samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block));
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);
  assert(target_disp >= 0);

  sample_index = target_disp / (samples_in_file / agg_id->idx_d->aggregation_factor);
  assert(sample_index < agg_id->idx->variable[variable_index]->values_per_sample * agg_id->idx_d->aggregation_factor);

  target_disp = target_disp % (samples_in_file / agg_id->idx_d->aggregation_factor);


  target_rank = agg_id->idx_d->agg_buffer->rank_holder[file_no][variable_index - agg_id->first_index][sample_index];
  target_count = hz_count * values_per_sample;

  bytes_per_datatype = ((agg_id->idx->variable[variable_index]->bits_per_value / 8) * total_chunk_size) / (64 / agg_id->idx->compression_bit_rate);

  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;

  start_agg_index = target_disp / (int64_t) (samples_in_file / agg_id->idx_d->aggregation_factor);
  end_agg_index = ((target_disp + target_count - 1) / (int64_t) (samples_in_file / agg_id->idx_d->aggregation_factor));
  assert(start_agg_index >= 0 && end_agg_index >= 0 && end_agg_index >= start_agg_index);

  if (start_agg_index != end_agg_index)
  {
    if (target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      //target_disp_address = target_disp;
      if (MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[A] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank, (long long)((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp), 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Put(hz_buffer, ((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, ((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }

#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    }
    else
      if (MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MA] Count %lld Local Disp %d Target Disp %lld\n", (long long)((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp), 0, (long long) target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( agg_id->idx_d->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, ( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer, agg_id->idx_d->agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) * bytes_per_datatype);

    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++)
    {
      if (target_rank != rank)
      {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
        if (MODE == PIDX_WRITE)
        {

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[B] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_id->aggregator_interval), (long long)(samples_in_file / agg_id->idx_d->aggregation_factor), (long long)(( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
          ret = MPI_Put(hz_buffer + (( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }

        }
        else
        {
          ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
      }
      else
      {
        if (MODE == PIDX_WRITE)
        {

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[MB] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_id->idx_d->aggregation_factor), (long long)(((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
          memcpy( agg_id->idx_d->agg_buffer->buffer, hz_buffer + (( (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, ( samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype);
        }
        else
          memcpy( hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, agg_id->idx_d->agg_buffer->buffer, (samples_in_file / agg_id->idx_d->aggregation_factor) * bytes_per_datatype);
      }
    }

    if (target_rank + agg_id->aggregator_interval != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
      if (MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[C] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_id->aggregator_interval), (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_d->aggregation_factor))) + (((samples_in_file / agg_id->idx_d->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif


        ret = MPI_Put(hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_d->aggregation_factor))) + (((samples_in_file / agg_id->idx_d->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp)) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_d->aggregation_factor))) + (((samples_in_file / agg_id->idx_d->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp)) * bytes_per_datatype,
                      MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
    }
    else
      if(MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MC] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_d->aggregation_factor))) + (((samples_in_file / agg_id->idx_d->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

        memcpy( agg_id->idx_d->agg_buffer->buffer, hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp)) * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer + (((samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_d->aggregation_factor))) * bytes_per_datatype, agg_id->idx_d->agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_d->aggregation_factor) - target_disp)) * bytes_per_datatype);
  }
  else
  {
    if(target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      //target_disp_address = target_disp;
      if(MODE == PIDX_WRITE)
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    }
    else
    {
      if(MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif        
        memcpy( agg_id->idx_d->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
      }
      else
      {
        //double x;
        //memcpy(&x, agg_id->idx_d->agg_buffer->buffer + target_disp * bytes_per_datatype, bytes_per_datatype);
        //printf("count %d = %f\n", hz_count, x);
        memcpy( hz_buffer, agg_id->idx_d->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }
  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_create(PIDX_agg_id agg_id)
{
  int v, p;
  if (agg_id->idx->enable_agg == 0)
  {
    PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
    for (v = agg_id->first_index; v <= agg_id->last_index; v++)
    {
      PIDX_variable var = agg_id->idx->variable[v];
      for (p = 0; p < var0->patch_group_count; p++)
      {
        var->hz_buffer[p]->HZ_io_from = 0;
        var->hz_buffer[p]->HZ_io_to = agg_id->idx_d->maxh;

        var->hz_buffer[p]->HZ_agg_from = 0;
        var->hz_buffer[p]->HZ_agg_to = 0;
      }
    }
    return PIDX_success;
  }

  int no_of_aggregators = 0;
  int per_file_aggregator = 0;
  int rank= 0, nprocs = 1;
  int i, j, k;
  int rank_counter = 0;
  int level;

#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  for (v = agg_id->first_index; v <= agg_id->last_index; v++)
    no_of_aggregators = no_of_aggregators + agg_id->idx->variable[v]->values_per_sample * agg_id->idx->variable[agg_id->init_index]->existing_file_count;

  for (v = agg_id->first_index; v <= agg_id->last_index; v++)
    per_file_aggregator = per_file_aggregator + agg_id->idx->variable[v]->values_per_sample;

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  if (nprocs < no_of_aggregators * agg_id->idx_d->aggregation_factor)
  {
    if ((agg_id->idx_d->max_file_count == 1 || nprocs == 1))
    {
      for (v = agg_id->first_index; v <= agg_id->last_index; v++)
      {
        PIDX_variable var = agg_id->idx->variable[v];
        for (p = 0; p < var0->patch_group_count; p++)
        {
          var->hz_buffer[p]->HZ_io_from = 0;
          var->hz_buffer[p]->HZ_io_to = agg_id->idx_d->maxh;

          var->hz_buffer[p]->HZ_agg_from = 0;
          var->hz_buffer[p]->HZ_agg_to = 0;
        }
      }
      agg_id->idx->enable_agg = 0;
      return PIDX_success;
    }
    else
    {
      while (no_of_aggregators > nprocs)
      {
        agg_id->idx_d->max_file_count = agg_id->idx_d->max_file_count / 2;
        agg_id->idx->variable[agg_id->init_index]->existing_file_count = 0;
        for (i = 0; i < agg_id->idx_d->max_file_count; i++)
          if (agg_id->idx->variable[agg_id->init_index]->file_index[i] == 1)
            agg_id->idx->variable[agg_id->init_index]->existing_file_count++;
        no_of_aggregators = per_file_aggregator * agg_id->idx->variable[agg_id->init_index]->existing_file_count * agg_id->idx_d->aggregation_factor;
      }
      level = getLeveL((agg_id->idx->variable[agg_id->init_index]->existing_file_count * agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file) - 1);


      for (v = agg_id->first_index; v <= agg_id->last_index; v++)
      {
        PIDX_variable var = agg_id->idx->variable[v];
        for (p = 0; p < var0->patch_group_count; p++)
        {
          var->hz_buffer[p]->HZ_agg_from = 0;
          var->hz_buffer[p]->HZ_agg_to = level;

          var->hz_buffer[p]->HZ_io_from = level + 1;
          var->hz_buffer[p]->HZ_io_to = agg_id->idx_d->maxh;
        }
      }
      agg_id->idx->enable_agg = 1;
    }
  }
  else
  {
    agg_id->idx->enable_agg = 2;
    for (v = agg_id->first_index; v <= agg_id->last_index; v++)
    {
      PIDX_variable var = agg_id->idx->variable[v];
      for (p = 0; p < var0->patch_group_count; p++)
      {
        var->hz_buffer[p]->HZ_io_from = 0;
        var->hz_buffer[p]->HZ_io_to = 0;

        var->hz_buffer[p]->HZ_agg_from = 0;
        var->hz_buffer[p]->HZ_agg_to = agg_id->idx_d->maxh;
      }
    }
  }

  agg_id->aggregator_interval = nprocs / (no_of_aggregators * agg_id->idx_d->aggregation_factor);
  assert(agg_id->aggregator_interval != 0);

  agg_id->idx_d->agg_buffer = malloc(sizeof(*agg_id->idx_d->agg_buffer));
  memset(agg_id->idx_d->agg_buffer, 0, sizeof(*agg_id->idx_d->agg_buffer));

  Agg_buffer agg_buffer = agg_id->idx_d->agg_buffer;

  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_buffer->rank_holder = malloc(agg_id->idx_d->max_file_count * sizeof (int**));
  for (i = 0; i < agg_id->idx_d->max_file_count; i++)
  {
    agg_buffer->rank_holder[i] = malloc((agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      agg_buffer->rank_holder[i][j - agg_id->first_index] = malloc(agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_d->aggregation_factor);
      memset(agg_buffer->rank_holder[i][j - agg_id->first_index], 0, agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_d->aggregation_factor);
    }
  }

  rank_counter = 0;
  for (k = 0; k < agg_id->idx->variable[agg_id->init_index]->existing_file_count; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_id->idx_d->aggregation_factor; j++)
      {
        agg_buffer->rank_holder[agg_id->idx->variable[agg_id->init_index]->existing_file_index[k]][i - agg_id->first_index][j] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;

        if(rank == agg_buffer->rank_holder[agg_id->idx->variable[agg_id->init_index]->existing_file_index[k]][i - agg_id->first_index][j])
        {
          agg_buffer->file_number = agg_id->idx->variable[agg_id->init_index]->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;

          uint64_t sample_count = agg_id->idx->variable[agg_id->init_index]->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_id->idx_d->aggregation_factor;

          int total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

          int bytes_per_datatype = (total_chunk_size * agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (64 / agg_id->idx->compression_bit_rate);

          agg_buffer->buffer_size = sample_count * bytes_per_datatype;
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          if (agg_buffer->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) agg_buffer->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_agg_write(PIDX_agg_id agg_id)
{
  if (agg_id->idx->enable_agg == 0)
    return PIDX_success;

  int i, p, e1, v, ret = 0;
  int send_index = 0;
  int64_t index = 0, count = 0, hz_index = 0;
  int rank = 0;

#if PIDX_HAVE_MPI
  Agg_buffer agg_buffer = agg_id->idx_d->agg_buffer;
  MPI_Comm_rank(agg_id->comm, &rank);

  if (agg_id->idx_d->agg_buffer->buffer_size != 0)
  {
    int total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
    int bytes_per_datatype = total_chunk_size * (agg_id->idx->variable[agg_id->idx_d->agg_buffer->var_number]->bits_per_value/8) / (64/agg_id->idx->compression_bit_rate);

    ret = MPI_Win_create(agg_buffer->buffer, agg_buffer->buffer_size, bytes_per_datatype, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }
  else
  {
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }

#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, agg_id->win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif

#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    char agg_file_name[1024];
    ret = mkdir(agg_id->idx_d->agg_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, agg_id->idx_d->agg_dump_dir_name);
      return PIDX_err_agg;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(agg_id->comm);
#endif

    sprintf(agg_file_name, "%s/rank_%d", agg_id->idx_d->agg_dump_dir_name, rank);
    agg_dump_fp = fopen(agg_file_name, "a+");
    if (!agg_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] agg_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, agg_file_name);
      return PIDX_err_agg;
    }
  }
#endif

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  for (p = 0; p < var0->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    HZ_buffer hz_buf = var0->hz_buffer[p];
    if(var0->hz_buffer[p]->type == 0)
    {  
      for (i = 0; i < hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i++)
        hz_index = hz_index + hz_buf->samples_per_level[i];

      for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
      {
        if (hz_buf->samples_per_level[i] != 0)
        {
          for(e1 = 0; e1 < hz_buf->samples_per_level[i] ; e1++)
          {
            if(e1 == 0)
            {
              index = hz_buf->buffer_index[hz_index];
              send_index = e1;
              count = 1;

              if(hz_buf->samples_per_level[i] == 1)
              {
                for(v = agg_id->first_index; v <= agg_id->last_index; v++)
                {
                  ret = aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
                  }
                }
              }
            }
            else
            {
              if(hz_buf->buffer_index[hz_index] - hz_buf->buffer_index[hz_index - 1] == 1)
              {
                count++;
                if(e1 == hz_buf->samples_per_level[i] - 1)
                {
                  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
                  {
                    aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
                    }
                  }
                }
              }
              else
              {
                for(v = agg_id->first_index; v <= agg_id->last_index; v++)
                {
                  aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
                  }
                }

                if(e1 == hz_buf->samples_per_level[i] - 1)
                {
                  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
                  {
                    aggregate_write_read(agg_id, v, hz_buf->buffer_index[hz_index], 1, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, PIDX_WRITE);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
                    }
                  }
                }
                index = hz_buf->buffer_index[hz_index];
                count = 1;
                send_index = e1;
              }
            }
            hz_index++;
          }
        }
      }
    }
    else
    {
      for(v = agg_id->first_index; v <= agg_id->last_index; v++)
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Variable %d\n", v);
          fflush(agg_dump_fp);
        }
#endif
        for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        {
          if (hz_buf->samples_per_level[i] != 0)
          {
            index = 0;
            count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1 - (var0->hz_buffer[p]->missing_block_count_per_level[i] * agg_id->idx_d->samples_per_block);


#ifdef PIDX_DUMP_AGG
            if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif
            ret = aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, PIDX_WRITE);
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
        }
      }
    }
  }

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, agg_id->win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#else
  //MPI_Win_create has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
  ret = MPI_Win_free(&(agg_id->win));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_agg_read(PIDX_agg_id agg_id)
{
  int i, p, var, ret = 0;
  int64_t count = 0;

  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    char agg_file_name[1024];

    ret = mkdir(agg_id->idx_d->agg_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, agg_id->idx_d->agg_dump_dir_name);
      return -1;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(agg_id->comm);
#endif

    sprintf(agg_file_name, "%s/rank_%d", agg_id->idx_d->agg_dump_dir_name, rank);
    agg_dump_fp = fopen(agg_file_name, "a+");
    if (!agg_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] agg_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, agg_file_name);
      return PIDX_err_agg;
    }
  }
#endif

#if PIDX_HAVE_MPI
  //agg_id->idx_d->win_time_start = MPI_Wtime();
  if (agg_id->idx_d->agg_buffer->buffer_size != 0)
    MPI_Win_create(agg_id->idx_d->agg_buffer->buffer, agg_id->idx_d->agg_buffer->buffer_size, agg_id->idx->variable[agg_id->idx_d->agg_buffer->var_number]->bits_per_value/8, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  else
    MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  //agg_id->idx_d->win_time_end = MPI_Wtime();
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
#endif

  for (p = 0; p < agg_id->idx->variable[agg_id->first_index]->patch_group_count; p++)
  {
    count = 0;
    for (var = agg_id->first_index; var <= agg_id->last_index; var++)
    {
#ifdef PIDX_DUMP_AGG
      if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
      {
        fprintf(agg_dump_fp, "Variable %d\n", var);
        fflush(agg_dump_fp);
      }
#endif
      for (i = agg_id->idx->variable[agg_id->first_index]->hz_buffer[p]->HZ_agg_from + agg_id->idx_d->res_from; i < agg_id->idx->variable[agg_id->first_index]->hz_buffer[p]->HZ_agg_to - agg_id->idx_d->res_to; i++)
      {
        count =  agg_id->idx->variable[var]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i] + 1 - (agg_id->idx->variable[var]->hz_buffer[p]->missing_block_count_per_level[i] * agg_id->idx_d->samples_per_block);

        //if (agg_id->idx->variable[agg_id->first_index]->hz_buffer[p]->samples_per_level[i] != 0)
        if (count != 0)
        {
#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[%d]: ", i);
            fflush(agg_dump_fp);
          }
#endif
          //agg_id->idx_d->agg_level_start[p][var][i] = MPI_Wtime();
          ret = aggregate_write_read(agg_id, var, agg_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[var]->hz_buffer[p]->buffer[i], 0, PIDX_READ);
          if (ret != PIDX_success)
          {
            fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
          //agg_id->idx_d->agg_level_end[p][var][i] = MPI_Wtime();
        }
      }
    }
  }

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_create has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
  //agg_id->idx_d->win_free_time_start = MPI_Wtime();
  MPI_Win_free(&(agg_id->win));
  //agg_id->idx_d->win_free_time_end = MPI_Wtime();
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_destroy(PIDX_agg_id agg_id)
{
  int i = 0, j = 0, v = 0, p = 0, itr = 0;
  if (agg_id->idx->enable_agg == 1 || agg_id->idx->enable_agg == 2)
  {
    if (agg_id->idx_d->agg_buffer->buffer_size != 0)
    {
      free(agg_id->idx_d->agg_buffer->buffer);
      agg_id->idx_d->agg_buffer->buffer = 0;
    }

    for (i = 0; i < agg_id->idx_d->max_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        free(agg_id->idx_d->agg_buffer->rank_holder[i][j - agg_id->first_index]);
        agg_id->idx_d->agg_buffer->rank_holder[i][j - agg_id->first_index] = 0;
      }
      free(agg_id->idx_d->agg_buffer->rank_holder[i]);
    }

    free(agg_id->idx_d->agg_buffer->rank_holder);
    agg_id->idx_d->agg_buffer->rank_holder = 0;

    free(agg_id->idx_d->agg_buffer);
  }

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  for (p = 0; p < var0->patch_group_count; p++)
  {
    if (agg_id->idx->enable_agg == 0  || agg_id->idx->enable_agg == 1)
    {
      free(var0->hz_buffer[p]->samples_per_level);
      var0->hz_buffer[p]->samples_per_level = 0;

      free(var0->hz_buffer[p]->start_hz_index);
      free(var0->hz_buffer[p]->end_hz_index);

      free(var0->hz_buffer[p]->missing_block_count_per_level);

      if (var0->hz_buffer[p]->type == 0)
        free(var0->hz_buffer[p]->buffer_index);

      for (itr = 0; itr < agg_id->idx_d->maxh; itr++)
        free(var0->hz_buffer[p]->nsamples_per_level[itr]);
      free(var0->hz_buffer[p]->nsamples_per_level);
    }
  }

  for (v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    for (p = 0; p < var0->patch_group_count; p++)
    {
      if (agg_id->idx->enable_agg == 0 || agg_id->idx->enable_agg == 1)
      {
        for (itr = var->hz_buffer[p]->HZ_io_from ; itr < var->hz_buffer[p]->HZ_io_to ; itr++)
        {
          free(var->hz_buffer[p]->buffer[itr]);
          var->hz_buffer[p]->buffer[itr] = 0;
        }
      }

      free(var->hz_buffer[p]->buffer);
      var->hz_buffer[p]->buffer = 0;

      free(var->hz_buffer[p]);
      var->hz_buffer[p] = 0;

    }
  }


  return PIDX_success;
}


PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id)
{

  free(agg_id);
  agg_id = 0;

  return PIDX_success;
}
