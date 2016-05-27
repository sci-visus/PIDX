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

#include "../PIDX_agg.h"

#define PIDX_ACTIVE_TARGET
#define PIDX_DUMP_AGG

struct PIDX_local_agg_struct
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Comm global_comm;
  MPI_Comm local_comm;
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

  int ***rank_holder2;

};

#ifdef PIDX_DUMP_AGG
static FILE* agg_dump_fp;
#endif

//static PIDX_return_code create_window(PIDX_local_agg_id agg_id);
//static PIDX_return_code one_sided_data_com(PIDX_local_agg_id agg_id, int mode);

#if PIDX_HAVE_MPI
static PIDX_return_code local_aggregate_write_read(PIDX_local_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE, MPI_Comm comm, int layout_id);

static PIDX_return_code create_open_log_file (PIDX_local_agg_id agg_id);

#endif


#if PIDX_HAVE_MPI
static PIDX_return_code create_window(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, MPI_Comm comm)
{
  int rank = 0, ret = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(comm, &rank);

  if (agg_buffer->buffer_size != 0)
  {
    int total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
    int bytes_per_datatype = total_chunk_size * (agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (agg_id->idx->compression_factor);

#if !SIMULATE_IO
    ret = MPI_Win_create(agg_buffer->buffer, agg_buffer->buffer_size, bytes_per_datatype, MPI_INFO_NULL, comm, &(agg_id->win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
#endif
  }
  else
  {
#if !SIMULATE_IO
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, comm, &(agg_id->win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
#endif
  }
#endif
  return PIDX_success;
}
#endif


#if PIDX_HAVE_MPI
static PIDX_return_code layout_zero(PIDX_local_agg_id agg_id, PIDX_block_layout block_layout)
{

  int rank = 0;
#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif
  int v = 0, p = 0;
  int ret = 0;

  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    int values_per_sample = var->values_per_sample;
    int64_t total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
    int bytes_per_datatype = ((var->bits_per_value / 8) * total_chunk_size) / (agg_id->idx->compression_factor);
    int target_rank = agg_id->rank_holder2[block_layout->inverse_existing_file_index[0]][v - agg_id->first_index][0];

    for (p = 0; p < var->patch_group_count; p++)
    {
      HZ_buffer hz_buf = agg_id->idx->variable[v]->hz_buffer[p];
#if !SIMULATE_IO
      //printf("buffer size = %d (%d %d %d) TR %d TD %d\n", hz_buf->lower_hz_buffer_size * values_per_sample * bytes_per_datatype, hz_buf->lower_hz_buffer_size, values_per_sample, bytes_per_datatype, target_rank, hz_buf->lower_hz_taret_disp * values_per_sample * bytes_per_datatype);

      ret = MPI_Put(hz_buf->lower_hz_buffer, hz_buf->lower_hz_buffer_size * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, 0, 1, hz_buf->lower_level_datatype, agg_id->win);
      if(ret != MPI_SUCCESS)
      {
        fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_agg;
      }
#endif
      if (ret != PIDX_success)
      {
        fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_agg;
      }
    }
  }

  return PIDX_success;
}
#endif



#if PIDX_HAVE_MPI
static PIDX_return_code local_one_sided_data_com(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout block_layout, int mode)
{
  int i, p, v, ret = 0;
  int64_t index = 0, count = 0;
  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      index = 0, count = 0;
      HZ_buffer hz_buf = var->hz_buffer[p];


      if (hz_buf->type == 1)
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
          fflush(agg_dump_fp);
        }
#endif

        //printf("From to To :: %d %d\n", block_layout->resolution_from, block_layout->resolution_to);
        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            index = 0;
            count =  hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;

#ifdef PIDX_DUMP_AGG
            if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif

#if !SIMULATE_IO

            ret = local_aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
#else
            ret = local_aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
#endif
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
        }
      }

      else if (hz_buf->type == 2)
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
          fflush(agg_dump_fp);
        }
#endif

        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
#ifdef PIDX_DUMP_AGG
            if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif
            int start_block_index = agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] / agg_id->idx_d->samples_per_block;
            int end_block_index = agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] / agg_id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              count = (agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);

#if !SIMULATE_IO
              ret = local_aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
#else
              ret = local_aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
#endif
            }
            else
            {
              int send_index = 0;
              int bl;
              for (bl = start_block_index; bl <= end_block_index; bl++)
              {
                if (PIDX_blocks_is_block_present(bl, block_layout))
                {
                  if (bl == start_block_index)
                  {
                    index = 0;
                    count = ((start_block_index + 1) * agg_id->idx_d->samples_per_block) - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i];
                  }
                  else if (bl == end_block_index)
                  {
                    index = (end_block_index * agg_id->idx_d->samples_per_block - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - ((end_block_index) * agg_id->idx_d->samples_per_block) + 1;
                  }
                  else
                  {
                    index = (bl * agg_id->idx_d->samples_per_block - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = agg_id->idx_d->samples_per_block;
                  }

#if !SIMULATE_IO
                  ret = local_aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
#else
                  ret = local_aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, NULL, send_index, agg_buffer, block_layout, mode, agg_id->comm, layout_id);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
#endif
                  send_index = send_index + count;
                }
                else
                  send_index = send_index + agg_id->idx_d->samples_per_block;
              }
            }
          }
        }
      }
    }
  }

  return PIDX_success;
}
#endif


#if PIDX_HAVE_MPI
static PIDX_return_code local_aggregate_write_read(PIDX_local_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE, MPI_Comm comm, int layout_id)
{
  int rank = 0, nprocs = 1, itr;// nrank = 0;nprocs
  int bytes_per_datatype;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, values_per_sample;
  int target_rank = 0;
  int64_t start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, hz_start = 0, samples_in_file = 0;
  int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  //MPI_Aint target_disp_address;

  int64_t total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

#if PIDX_HAVE_MPI
  int ret;

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
  //negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, agg_id->idx->variable[agg_id->init_index]->global_block_layout);
  negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
  if (negative_block_offset < 0)
    return PIDX_err_agg;

  //number of samples in file "file_no"
  //samples_in_file = agg_id->idx->variable[agg_id->init_index]->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;
  samples_in_file = block_layout->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;
  if (samples_in_file > samples_per_file)
    return PIDX_err_agg;

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  //assert(hz_start >= (samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block));


  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);
  if (target_disp < 0)
    return PIDX_err_agg;


  sample_index = target_disp / (samples_in_file / agg_buffer->aggregation_factor);
  if (sample_index >= agg_id->idx->variable[variable_index]->values_per_sample * agg_buffer->aggregation_factor)
    return PIDX_err_agg;

  target_disp = target_disp % (samples_in_file / agg_buffer->aggregation_factor);

  target_rank = agg_id->rank_holder2[block_layout->inverse_existing_file_index[file_no]][variable_index - agg_id->first_index][sample_index];

  MPI_Comm_rank(agg_id->comm, &rank);
  MPI_Comm_size(agg_id->comm, &nprocs);
  //printf("COMM size = %d\n", nprocs);

  //
  if (layout_id != 0 && agg_id->idx->current_time_step == 0)
  {
    MPI_Comm agg_comm;
    int max_rank = 0;
    int min_rank = 0;

    MPI_Comm_split(agg_id->comm, target_rank, rank, &agg_comm);
    int nrank = 0;
    MPI_Comm_rank(agg_comm, &nrank);

    MPI_Allreduce(&rank, &max_rank, 1, MPI_INT, MPI_MAX, agg_comm);
    MPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, agg_comm);

    //if (nrank == 0 && variable_index == 0)
    //  printf("[%d] [%d %d]: %d [%d - %d]\n", layout_id, block_layout->inverse_existing_file_index[file_no], sample_index, target_rank, min_rank, max_rank);
    MPI_Comm_free(&agg_comm);

    if (target_rank < min_rank || target_rank > max_rank || rank < min_rank || rank > max_rank)
    {
      printf("[TR %d] [%d] V %d P %d A %d FN %d [%d %d]\n", target_rank, rank, variable_index, layout_id, block_layout->inverse_existing_file_index[file_no], file_no, min_rank, max_rank);

      //printf("[TR %d] [%d] V %d P %d A %d = %d %d\n", target_rank, rank, variable_index, layout_id, block_layout->inverse_existing_file_index[file_no], agg_id->idx_d->layout_agg_range[variable_index][layout_id][block_layout->inverse_existing_file_index[file_no]][0], agg_id->idx_d->layout_agg_range[variable_index][layout_id][block_layout->inverse_existing_file_index[file_no]][1]);
    }

    assert(target_rank >= min_rank);
    assert(target_rank <= max_rank);
    assert(rank >= min_rank);
    assert(rank <= max_rank);
  }
  //

  target_count = hz_count * values_per_sample;

  bytes_per_datatype = ((agg_id->idx->variable[variable_index]->bits_per_value / 8) * total_chunk_size) / (agg_id->idx->compression_factor);

#if !SIMULATE_IO
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;
#endif

  start_agg_index = target_disp / (int64_t) (samples_in_file / agg_buffer->aggregation_factor);
  end_agg_index = ((target_disp + target_count - 1) / (int64_t) (samples_in_file / agg_buffer->aggregation_factor));
  //assert(start_agg_index >= 0 && end_agg_index >= 0 && end_agg_index >= start_agg_index);


  /*
  if (start_agg_index != end_agg_index)
  {

#if PIDX_HAVE_MPI
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[A] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank, (long long)((samples_in_file / agg_buffer->aggregation_factor) - target_disp), 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        ret = MPI_Put(hz_buffer, ((samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif

#endif


    if (rank == 31)
      printf("ABCD [%d] %d %d\n", end_agg_index - start_agg_index - 1, end_agg_index, start_agg_index);

    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++)
    {

#if PIDX_HAVE_MPI
#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[B] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), (long long)(samples_in_file / agg_buffer->aggregation_factor), (long long)(( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
#if !SIMULATE_IO
          ret = MPI_Put(hz_buffer + (( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
#endif
#endif
    }

#if PIDX_HAVE_MPI

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[C] Target Rank %d [%d + %d] Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), target_rank, agg_buffer->aggregator_interval, (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        ret = MPI_Put(hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif

#endif


  }
  else
  {
#if PIDX_HAVE_MPI

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif


#endif
  }
  */

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
          fprintf(agg_dump_fp, "[A] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank, (long long)((samples_in_file / agg_buffer->aggregation_factor) - target_disp), 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        ret = MPI_Put(hz_buffer, ((samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif
      }
      else
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[A] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank, (long long)((samples_in_file / agg_buffer->aggregation_factor) - target_disp), 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Get(hz_buffer, ((samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
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
          fprintf(agg_dump_fp, "[MA] Count %lld Local Disp %d Target Disp %lld\n", (long long)((samples_in_file / agg_buffer->aggregation_factor) - target_disp), 0, (long long) target_disp);
          fflush(agg_dump_fp);
        }
#endif
#if !SIMULATE_IO
        memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype);
#endif
      }
      else
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MA] Count %lld Local Disp %d Target Disp %lld\n", (long long)((samples_in_file / agg_buffer->aggregation_factor) - target_disp), 0, (long long) target_disp);
          fflush(agg_dump_fp);
        }
#endif

        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype);
      }

    //if (rank == 31)
    //  printf("ABCD [%d] %d %d\n", end_agg_index - start_agg_index - 1, end_agg_index, start_agg_index);

    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++)
    {
      if (target_rank != rank)
      {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_buffer->aggregator_interval, 0, agg_id->win);
#endif
        if (MODE == PIDX_WRITE)
        {

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[B] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), (long long)(samples_in_file / agg_buffer->aggregation_factor), (long long)(( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
#if !SIMULATE_IO
          ret = MPI_Put(hz_buffer + (( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
#endif
        }
        else
        {

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[B] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), (long long)(samples_in_file / agg_buffer->aggregation_factor), (long long)(( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif

          ret = MPI_Get(hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(target_rank + agg_buffer->aggregator_interval, agg_id->win);
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
            fprintf(agg_dump_fp, "[MB] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_buffer->aggregation_factor), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif

#if !SIMULATE_IO
          memcpy( agg_buffer->buffer, hz_buffer + (( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, ( samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype);
#endif
        }
        else
        {
#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[MB] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_buffer->aggregation_factor), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
          memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype);
        }
      }
    }

    if (target_rank + agg_buffer->aggregator_interval != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_buffer->aggregator_interval, 0, agg_id->win);
#endif
      if (MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[C] Target Rank %d [%d + %d] Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), target_rank, agg_buffer->aggregator_interval, (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        ret = MPI_Put(hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif
      }
      else
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[C] Target Rank %d [%d + %d] Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_buffer->aggregator_interval), target_rank, agg_buffer->aggregator_interval, (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Get(hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_buffer->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype,
                      MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank + agg_buffer->aggregator_interval, agg_id->win);
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
          fprintf(agg_dump_fp, "[MC] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

#if !SIMULATE_IO
        memcpy( agg_buffer->buffer, hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype);
#endif
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MC] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_buffer->aggregation_factor))) + (((samples_in_file / agg_buffer->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

        memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype);
      }
    }
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

#if !SIMULATE_IO
        //printf("[%d %d] XX Count %d Offset %d\n", rank, target_rank, hz_count * values_per_sample * bytes_per_datatype, target_disp);
        ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
#endif
      }
      else
      {

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

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
#if !SIMULATE_IO
        memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
#endif
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }

  return PIDX_success;
}
#endif


PIDX_local_agg_id PIDX_local_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_local_agg_id agg_id;

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
PIDX_return_code PIDX_local_agg_set_communicator(PIDX_local_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->comm = comm;

  return PIDX_success;
}

PIDX_return_code PIDX_local_agg_set_global_communicator(PIDX_local_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->global_comm = comm;


  return PIDX_success;
}
#endif


PIDX_return_code PIDX_local_agg_meta_data_create(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout)
{
#if PIDX_HAVE_MPI
  int i = 0, j = 0, d = 0;
  int rank = 0, nprocs = 1;
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);

  //printf("local_block_layout->existing_file_count * agg_buffer->aggregation_factor = %d x %d\n", local_block_layout->existing_file_count, agg_buffer->aggregation_factor);
  agg_buffer->aggregator_interval = nprocs / ((agg_id->last_index - agg_id->first_index + 1) * local_block_layout->existing_file_count * agg_buffer->aggregation_factor);
  //assert(agg_buffer->aggregator_interval != 0);

  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_id->rank_holder2 = malloc(local_block_layout->existing_file_count * sizeof (int**));
  memset(agg_id->rank_holder2, 0, local_block_layout->existing_file_count * sizeof (int**));
  //if (rank == 0)
  //  printf("FC = %d\n", local_block_layout->existing_file_count);
  for (i = 0; i < local_block_layout->existing_file_count; i++)
  {
    agg_id->rank_holder2[i] = malloc((agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    memset(agg_id->rank_holder2[i], 0, (agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    //if (rank == 0)
    //  printf("VC = %d %d\n", agg_id->first_index, agg_id->last_index);
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      agg_id->rank_holder2[i][j - agg_id->first_index] = malloc(agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      memset(agg_id->rank_holder2[i][j - agg_id->first_index], 0, agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      for (d = 0 ; d < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; d++)
        agg_id->rank_holder2[i][j - agg_id->first_index][d] = -1;
    }
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_meta_data_destroy(PIDX_local_agg_id agg_id, PIDX_block_layout local_block_layout)
{
  int i = 0, j = 0;
  for (i = 0; i < local_block_layout->existing_file_count; i++)
  {
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      free(agg_id->rank_holder2[i][j - agg_id->first_index]);
      agg_id->rank_holder2[i][j - agg_id->first_index] = 0;
    }
    free(agg_id->rank_holder2[i]);
  }
  free(agg_id->rank_holder2);
  agg_id->rank_holder2 = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_local_agg_buf_create_multiple_level(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int agg_offset, int var_offset, int file_status)
{
#if PIDX_HAVE_MPI
  int rank = 0, nprocs = 1, grank = 0, gnprocs = 1;
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
  MPI_Comm_rank(agg_id->global_comm, &grank);
  MPI_Comm_size(agg_id->global_comm, &gnprocs);
  //MPI_Allreduce(&rank, global_rank, 1, MPI_INT, )



  /*
  PIDX_variable var = agg_id->idx->variable[agg_id->first_index];
  HZ_buffer hz_buf = var->hz_buffer[0];
  int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  int values_per_sample = agg_id->idx->variable[0]->values_per_sample;
  int hz_start = hz_buf->start_hz_index[local_block_layout->resolution_from];
  int file_no = hz_start / samples_per_file;
  int block_no = hz_start / agg_id->idx_d->samples_per_block;
  int negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, local_block_layout);
  int samples_in_file = local_block_layout->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;
  int target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);

  int sample_index = target_disp / (samples_in_file / agg_buffer->aggregation_factor);
  int target_rank = sample_index;//agg_id->rank_holder2[local_block_layout->inverse_existing_file_index[file_no]][0][sample_index];
  //printf("[%d] [%d %d] Target Rank = %d\n", rank, (local_block_layout->inverse_existing_file_index[file_no]), sample_index, target_rank);
  MPI_Comm_split(agg_id->comm, target_rank, rank, &(agg_id->local_comm));

  //int tsize;
  //MPI_Comm_size(agg_id->local_comm, &tsize);
  //printf("SSSSSSSSSSS222 = %d\n", tsize);


  MPI_Comm_rank(agg_id->local_comm, &rank);
  */

  int i = 0, j = 0, k = 0;
  for (k = 0; k < local_block_layout->existing_file_count; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_buffer->aggregation_factor; j++)
      {
        uint64_t file_index = 0;

        int *first;
        first = malloc(sizeof(*first) * PIDX_MAX_DIMENSIONS);
        memset(first, 0, sizeof(*first) * PIDX_MAX_DIMENSIONS);

        int file_count = 1;

        int negative_file_index = 0;
        if (agg_offset == 0 || agg_offset == 1)
        {
          file_count = 1;
          file_index = 0 * agg_buffer->aggregation_factor + j;
        }
        else
        {
          file_count = (int)pow(2, agg_offset - 1);
          negative_file_index = (int)pow(2, agg_offset - 1);
          file_index = (local_block_layout->existing_file_index[k] - negative_file_index) * agg_buffer->aggregation_factor + j;
        }

        int old_rank = 0;
        if (agg_offset == 0)
          old_rank = file_index * agg_buffer->aggregator_interval;
        else if (agg_offset == 1)
          old_rank = (file_index + 1) * agg_buffer->aggregator_interval - 1;
        else
          old_rank = (file_index * agg_buffer->aggregator_interval) + (agg_buffer->aggregator_interval / 2);

#if 1
        int bits = (agg_id->idx_d->maxh - 1 - (int)log2(file_count * agg_buffer->aggregation_factor));

        //if (rank == 0)
        //  printf("[B %d] [%d] [%d - %d (%d) %d %d] [%s] -- %d %d\n", grank, agg_offset, k, local_block_layout->existing_file_index[k], local_block_layout->existing_file_count, i, j, agg_id->idx->idx_cl_bitSequence, agg_id->idx_d->maxh, file_index);

        file_index = file_index << bits;

        //if (rank == 0)
        //  printf("[A %d] [%d] [%d - %d (%d) %d %d] [%s] -- %d %d PB %d\n", grank, agg_offset, k, local_block_layout->existing_file_index[k], local_block_layout->existing_file_count, i, j, agg_id->idx->bitSequence, agg_id->idx_d->maxh, file_index, partition_bits);

        Deinterleave(agg_id->idx->bitPattern, (agg_id->idx_d->maxh - 1), file_index, first);

        int calculated_rank = 0;
        int rank_x = first[0] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[0]);
        int rank_y = first[1] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[1]);
        int rank_z = first[2] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[2]);

        int nrank_x = (agg_id->idx->bounds[0] / agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[0]);
        int nrank_y = (agg_id->idx->bounds[1] / agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[1]);

        calculated_rank = rank_x + (rank_y * nrank_x) + (rank_z * nrank_x * nrank_y);

        int trank = 0;
        if (file_status == 0)
        {
          if (agg_offset == 0)
            trank = agg_id->idx_d->rank_buffer[calculated_rank];
          else
            trank = agg_id->idx_d->rank_buffer[calculated_rank + (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2))];
          agg_id->rank_holder2[k][i - agg_id->first_index][j] = trank;//((trank / 2) + 1) * 2;
        }
        else
        {
          if (agg_offset == 0)
            agg_id->rank_holder2[k][i - agg_id->first_index][j] = agg_id->idx_d->rank_buffer[calculated_rank];
          else
            agg_id->rank_holder2[k][i - agg_id->first_index][j] = agg_id->idx_d->rank_buffer[calculated_rank + (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2))];

          trank = agg_id->rank_holder2[k][i - agg_id->first_index][j];
        }
        //
        //printf("%d -> %d : %d -> %d + %d\n", grank, trank, calculated_rank + (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2)), calculated_rank, (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2)) );

        free(first);
#endif
        //if (grank == 0)
        //  printf("%d -> %d : %d -> %d + %d\n", grank, trank, calculated_rank + (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2)), calculated_rank, (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2)) );
        if(rank == agg_id->rank_holder2[k][i - agg_id->first_index][j])
        {
          agg_buffer->file_number = local_block_layout->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;

          uint64_t sample_count = local_block_layout->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_buffer->aggregation_factor;

          int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

          int bytes_per_datatype = 0;
          bytes_per_datatype = (chunk_size * agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (agg_id->idx->compression_factor);

          agg_buffer->buffer_size = sample_count * bytes_per_datatype;

          if (i == 0)
            printf("[G %d] [%d] [%d] [F [%d %d] V %d  S %d] -> UR %d TR %d CR %d OR %d Size %d [%d (%d x %d / %d) x %d]\n", grank, agg_offset, agg_offset, k, local_block_layout->existing_file_index[k], i, j, rank, trank, calculated_rank, old_rank, (int)agg_buffer->buffer_size, agg_buffer->file_number, local_block_layout->block_count_per_file[agg_buffer->file_number], agg_id->idx_d->samples_per_block, agg_buffer->aggregation_factor, bytes_per_datatype);//, first[0], first[1], first[2], rank_x, rank_y, rank_z);

#if !SIMULATE_IO
          //double bs_time = MPI_Wtime();
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
          if (agg_buffer->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) agg_buffer->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
#endif
        }
      }
    }
  }

#if 0
  if (rank == 0)
  {
    for (i = 0; i < local_block_layout->existing_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        for (k = 0; k < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; k++)
        {
          printf("[%d] [%d %d %d] -> %d\n", agg_offset, k, i, j, agg_id->rank_holder2[i][j-agg_id->first_index][k]);
        }
      }
    }
  }
#endif
#endif

  return PIDX_success;
}




PIDX_return_code PIDX_local_agg_buf_create(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int agg_offset, int var_offset)
{
#if PIDX_HAVE_MPI
  int rank = 0, nprocs = 1, grank = 0, gnprocs = 1;
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
  MPI_Comm_rank(agg_id->global_comm, &grank);
  MPI_Comm_size(agg_id->global_comm, &gnprocs);
  //MPI_Allreduce(&rank, global_rank, 1, MPI_INT, )



  /*
  PIDX_variable var = agg_id->idx->variable[agg_id->first_index];
  HZ_buffer hz_buf = var->hz_buffer[0];
  int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  int values_per_sample = agg_id->idx->variable[0]->values_per_sample;
  int hz_start = hz_buf->start_hz_index[local_block_layout->resolution_from];
  int file_no = hz_start / samples_per_file;
  int block_no = hz_start / agg_id->idx_d->samples_per_block;
  int negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, local_block_layout);
  int samples_in_file = local_block_layout->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;
  int target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_d->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);

  int sample_index = target_disp / (samples_in_file / agg_buffer->aggregation_factor);
  int target_rank = sample_index;//agg_id->rank_holder2[local_block_layout->inverse_existing_file_index[file_no]][0][sample_index];
  //printf("[%d] [%d %d] Target Rank = %d\n", rank, (local_block_layout->inverse_existing_file_index[file_no]), sample_index, target_rank);
  MPI_Comm_split(agg_id->comm, target_rank, rank, &(agg_id->local_comm));

  //int tsize;
  //MPI_Comm_size(agg_id->local_comm, &tsize);
  //printf("SSSSSSSSSSS222 = %d\n", tsize);


  MPI_Comm_rank(agg_id->local_comm, &rank);
  */

  int i = 0, j = 0, k = 0;
  for (k = 0; k < local_block_layout->existing_file_count; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_buffer->aggregation_factor; j++)
      {
        uint64_t file_index = 0;

        int *first;
        first = malloc(sizeof(*first) * PIDX_MAX_DIMENSIONS);
        memset(first, 0, sizeof(*first) * PIDX_MAX_DIMENSIONS);

        int file_count = 1;

        int negative_file_index = 0;
        if (agg_offset == 0 || agg_offset == 1)
        {
          file_count = 1;
          file_index = 0 * agg_buffer->aggregation_factor + j;
        }
        else
        {
          file_count = (int)pow(2, agg_offset - 1);
          negative_file_index = (int)pow(2, agg_offset - 1);
          file_index = (local_block_layout->existing_file_index[k] - negative_file_index) * agg_buffer->aggregation_factor + j;
        }

        int old_rank = 0;
        if (agg_offset == 0)
          old_rank = file_index * agg_buffer->aggregator_interval;
        else if (agg_offset == 1)
          old_rank = (file_index + 1) * agg_buffer->aggregator_interval - 1;
        else
          old_rank = (file_index * agg_buffer->aggregator_interval) + (agg_buffer->aggregator_interval / 2);

#if 1
        int bits = (agg_id->idx_d->maxh - 1 - (int)log2(file_count * agg_buffer->aggregation_factor));

        //if (rank == 0)
        //  printf("[B %d] [%d] [%d - %d (%d) %d %d] [%s] -- %d %d\n", grank, agg_offset, k, local_block_layout->existing_file_index[k], local_block_layout->existing_file_count, i, j, agg_id->idx->idx_cl_bitSequence, agg_id->idx_d->maxh, file_index);

        file_index = file_index << bits;

        //if (rank == 0)
        //  printf("[A %d] [%d] [%d - %d (%d) %d %d] [%s] -- %d %d PB %d\n", grank, agg_offset, k, local_block_layout->existing_file_index[k], local_block_layout->existing_file_count, i, j, agg_id->idx->bitSequence, agg_id->idx_d->maxh, file_index, partition_bits);

        Deinterleave(agg_id->idx->bitPattern, (agg_id->idx_d->maxh - 1), file_index, first);

        int calculated_rank = 0;
        int rank_x = first[0] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[0]);
        int rank_y = first[1] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[1]);
        int rank_z = first[2] / (agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[2]);

        int nrank_x = (agg_id->idx->bounds[0] / agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[0]);
        int nrank_y = (agg_id->idx->bounds[1] / agg_id->idx->variable[agg_id->first_index]->sim_patch[0]->size[1]);

        calculated_rank = rank_x + (rank_y * nrank_x) + (rank_z * nrank_x * nrank_y);

        if (agg_offset == 0)
          agg_id->rank_holder2[k][i - agg_id->first_index][j] = agg_id->idx_d->rank_buffer[calculated_rank];
        else
          agg_id->rank_holder2[k][i - agg_id->first_index][j] = agg_id->idx_d->rank_buffer[calculated_rank + (nprocs/ (local_block_layout->existing_file_count * agg_buffer->aggregation_factor * 2))];

        free(first);
#endif
        if(rank == agg_id->rank_holder2[k][i - agg_id->first_index][j])
        {
          agg_buffer->file_number = local_block_layout->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;

          uint64_t sample_count = local_block_layout->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_buffer->aggregation_factor;

          int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

          int bytes_per_datatype = 0;
          bytes_per_datatype = (chunk_size * agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (agg_id->idx->compression_factor);

          agg_buffer->buffer_size = sample_count * bytes_per_datatype;

          if (i == 0)
            printf("[G %d] [%d] [%d] [F [%d %d] V %d  S %d] -> UR %d CR %d OR %d Size %d [%d (%d x %d / %d) x %d]\n", grank, agg_offset, agg_offset, k, local_block_layout->existing_file_index[k], i, j, rank, calculated_rank, old_rank, (int)agg_buffer->buffer_size, agg_buffer->file_number, local_block_layout->block_count_per_file[agg_buffer->file_number], agg_id->idx_d->samples_per_block, agg_buffer->aggregation_factor, bytes_per_datatype);//, first[0], first[1], first[2], rank_x, rank_y, rank_z);

#if !SIMULATE_IO
          //double bs_time = MPI_Wtime();
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
          if (agg_buffer->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) agg_buffer->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
          //double be_time = MPI_Wtime();
          //if (rank == 0)
          //  printf("[XX] %d time at %d = %f\n", var_offset, agg_offset - agg_id->idx_d->start_layout_index, be_time - bs_time);
#endif
        }
      }
    }
  }


  //double e_time = MPI_Wtime();
  //if (rank == 0)
  //  printf("[YY] %d time at %d = %f\n", var_offset, agg_offset - agg_id->idx_d->start_layout_index, e_time - s_time);

#if 0
  if (rank == 0)
  {
    for (i = 0; i < local_block_layout->existing_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        for (k = 0; k < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; k++)
        {
          printf("[%d %d %d] -> %d\n", k, i, j, agg_id->rank_holder2[i][j-agg_id->first_index][k]);
        }
      }
    }
  }
#endif
#endif

  return PIDX_success;
}



PIDX_return_code PIDX_local_agg_buf_destroy(Agg_buffer agg_buffer)
{
  //MPI_Comm_free(&(agg_id->local_comm));
#if !SIMULATE_IO
  if (agg_buffer->buffer_size != 0)
  {
    free(agg_buffer->buffer);
    agg_buffer->buffer = 0;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code create_open_log_file (PIDX_local_agg_id agg_id)
{
  int rank = 0;
  int ret = 0;
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

  return PIDX_success;
}

static PIDX_return_code close_log_file (PIDX_local_agg_id agg_id)
{
#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}

#if PIDX_HAVE_MPI
static PIDX_return_code create_file_zero_buffer(PIDX_local_agg_id agg_id, PIDX_block_layout block_layout, int start_index, int end_index)
{
  int rank = 0;
#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  int64_t send_index = 0;
  int p = 0, i = 0, bytes_for_datatype = 0, v = 0;
  //int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  int count = 0;
  int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
  int buffer_count = 0;
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    bytes_for_datatype = ((var->bits_per_value / 8) * chunk_size * var->values_per_sample) / (agg_id->idx->compression_factor);

    for (p = 0; p < var->patch_group_count; p++)
    {
      HZ_buffer hz_buf = var->hz_buffer[p];

      if (hz_buf->type == 2)
      {
        buffer_count = 0;
        hz_buf->lower_hz_buffer_size = 0;

        for (i = /*block_layout->resolution_from*/start_index; i < /*block_layout->resolution_to*/end_index; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            int start_block_index = hz_buf->start_hz_index[i] / agg_id->idx_d->samples_per_block;
            int end_block_index = hz_buf->end_hz_index[i] / agg_id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              count = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1);
              hz_buf->lower_hz_buffer_size = hz_buf->lower_hz_buffer_size + count;
              buffer_count++;
            }
            else
            {
              send_index = 0;
              int bl;
              for (bl = start_block_index; bl <= end_block_index; bl++)
              {
                if (PIDX_blocks_is_block_present(bl, block_layout))
                {
                  if (bl == start_block_index)
                  {
                    count = ((start_block_index + 1) * agg_id->idx_d->samples_per_block) - hz_buf->start_hz_index[i];
                    hz_buf->lower_hz_buffer_size = hz_buf->lower_hz_buffer_size + count;
                    buffer_count++;
                  }
                  else if (bl == end_block_index)
                  {
                    count = hz_buf->end_hz_index[i] - ((end_block_index) * agg_id->idx_d->samples_per_block) + 1;
                    hz_buf->lower_hz_buffer_size = hz_buf->lower_hz_buffer_size + count;
                    buffer_count++;
                  }
                  else
                  {
                    count = agg_id->idx_d->samples_per_block;
                    hz_buf->lower_hz_buffer_size = hz_buf->lower_hz_buffer_size + count;
                    buffer_count++;
                  }
                  send_index = send_index + count;
                }
                else
                  send_index = send_index + agg_id->idx_d->samples_per_block;
              }
            }
          }
        }

        hz_buf->lower_hz_buffer = malloc(bytes_for_datatype * hz_buf->lower_hz_buffer_size);
        memset(hz_buf->lower_hz_buffer, 0, bytes_for_datatype * hz_buf->lower_hz_buffer_size);

        hz_buf->lower_hz_disp = malloc(sizeof(int) * buffer_count);
        memset(hz_buf->lower_hz_disp, 0, sizeof(int) * buffer_count);

        hz_buf->lower_hz_count = malloc(sizeof(int) * buffer_count);
        memset(hz_buf->lower_hz_count, 0, sizeof(int) * buffer_count);

        int block_no = 0;
        int negative_block_offset = 0;
        buffer_count = 0;
        int target_offset = 0;
        send_index = 0;
        for (i = /*block_layout->resolution_from*/start_index; i < /*block_layout->resolution_to*/end_index; i++)
        //for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            int start_block_index = hz_buf->start_hz_index[i] / agg_id->idx_d->samples_per_block;
            int end_block_index = hz_buf->end_hz_index[i] / agg_id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              hz_buf->lower_hz_disp[buffer_count] = hz_buf->start_hz_index[i] * bytes_for_datatype;
              hz_buf->lower_hz_count[buffer_count] = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1)  * bytes_for_datatype;

              block_no = hz_buf->start_hz_index[i] / (agg_id->idx_d->samples_per_block);
              negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
              hz_buf->lower_hz_disp[buffer_count] = hz_buf->lower_hz_disp[buffer_count] - (negative_block_offset * bytes_for_datatype * agg_id->idx_d->samples_per_block);

              //if (rank == 4)
              //  printf("[%d] O %d C %d SI %d\n", buffer_count, hz_buf->lower_hz_disp[buffer_count] / bytes_for_datatype, hz_buf->lower_hz_count[buffer_count] / bytes_for_datatype, send_index/bytes_for_datatype);

              memcpy(hz_buf->lower_hz_buffer + target_offset, hz_buf->buffer[i],  hz_buf->lower_hz_count[buffer_count]);
              target_offset = target_offset + hz_buf->lower_hz_count[buffer_count];
              buffer_count++;
            }
            else
            {
              send_index = 0;
              int bl;
              count = 0;
              for (bl = start_block_index; bl <= end_block_index; bl++)
              {
                if (PIDX_blocks_is_block_present(bl, block_layout))
                {
                  if (bl == start_block_index)
                  {
                    hz_buf->lower_hz_disp[buffer_count] = hz_buf->start_hz_index[i] * bytes_for_datatype;

                    block_no = hz_buf->start_hz_index[i] / agg_id->idx_d->samples_per_block;
                    negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
                    hz_buf->lower_hz_disp[buffer_count] = hz_buf->lower_hz_disp[buffer_count] - (negative_block_offset * bytes_for_datatype * agg_id->idx_d->samples_per_block);

                    hz_buf->lower_hz_count[buffer_count] = (((start_block_index + 1) * agg_id->idx_d->samples_per_block) - hz_buf->start_hz_index[i]) * bytes_for_datatype;
                    count = hz_buf->lower_hz_count[buffer_count];

                    //if (rank == 4)
                    //  printf("[%d] O %d C %d SI %d NI %d\n", buffer_count, hz_buf->lower_hz_disp[buffer_count] / bytes_for_datatype, hz_buf->lower_hz_count[buffer_count] / bytes_for_datatype, send_index/bytes_for_datatype, negative_block_offset);

                    memcpy(hz_buf->lower_hz_buffer + target_offset, hz_buf->buffer[i] + (send_index),  hz_buf->lower_hz_count[buffer_count]);
                    target_offset = target_offset + hz_buf->lower_hz_count[buffer_count];
                    buffer_count++;
                  }
                  else if (bl == end_block_index)
                  {
                    hz_buf->lower_hz_disp[buffer_count] = (( (end_block_index * agg_id->idx_d->samples_per_block - hz_buf->start_hz_index[i]) * bytes_for_datatype ) + hz_buf->start_hz_index[i] * bytes_for_datatype);

                    block_no = hz_buf->lower_hz_disp[buffer_count] / (agg_id->idx_d->samples_per_block * bytes_for_datatype);
                    negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
                    hz_buf->lower_hz_disp[buffer_count] = hz_buf->lower_hz_disp[buffer_count] - (negative_block_offset * bytes_for_datatype * agg_id->idx_d->samples_per_block);

                    hz_buf->lower_hz_count[buffer_count] = (hz_buf->end_hz_index[i] - ((end_block_index) * agg_id->idx_d->samples_per_block) + 1) * bytes_for_datatype;
                    count = hz_buf->lower_hz_count[buffer_count];

                    //if (rank == 4)
                    //  printf("[%d] O %d C %d SI %d NI %d\n", buffer_count, hz_buf->lower_hz_disp[buffer_count] / bytes_for_datatype, hz_buf->lower_hz_count[buffer_count] / bytes_for_datatype, send_index/bytes_for_datatype, negative_block_offset);

                    memcpy(hz_buf->lower_hz_buffer + target_offset, hz_buf->buffer[i] + (send_index),  hz_buf->lower_hz_count[buffer_count]);
                    target_offset = target_offset + hz_buf->lower_hz_count[buffer_count];
                    buffer_count++;
                  }
                  else
                  {
                    hz_buf->lower_hz_disp[buffer_count] = (bl * agg_id->idx_d->samples_per_block * bytes_for_datatype);

                    block_no = hz_buf->lower_hz_disp[buffer_count] / (agg_id->idx_d->samples_per_block * bytes_for_datatype);
                    negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
                    hz_buf->lower_hz_disp[buffer_count] = hz_buf->lower_hz_disp[buffer_count] - (negative_block_offset * bytes_for_datatype * agg_id->idx_d->samples_per_block);

                    hz_buf->lower_hz_count[buffer_count] = agg_id->idx_d->samples_per_block  * bytes_for_datatype;
                    count = hz_buf->lower_hz_count[buffer_count];

                    //if (rank == 4)
                    //  printf("[%d] O %d C %d SI %d NI %d\n", buffer_count, hz_buf->lower_hz_disp[buffer_count] / bytes_for_datatype, hz_buf->lower_hz_count[buffer_count] / bytes_for_datatype, send_index/bytes_for_datatype, negative_block_offset);

                    memcpy(hz_buf->lower_hz_buffer + target_offset, hz_buf->buffer[i] + (send_index),  hz_buf->lower_hz_count[buffer_count]);
                    target_offset = target_offset + hz_buf->lower_hz_count[buffer_count];
                    buffer_count++;
                  }
                  send_index = send_index + count;
                }
                else
                {
                  send_index = send_index + agg_id->idx_d->samples_per_block * bytes_for_datatype;
                }
              }
            }
          }
        }
      }

      //
      if (hz_buf->type == 1)
      {
        buffer_count = 0;
        hz_buf->lower_hz_buffer_size = 0;

        for (i = /*block_layout->resolution_from*/start_index; i < /*block_layout->resolution_to*/end_index; i++)
        //for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            hz_buf->lower_hz_buffer_size = hz_buf->lower_hz_buffer_size + (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1);
            buffer_count = buffer_count + 1;
          }
        }

        hz_buf->lower_hz_disp = malloc(sizeof(int) * buffer_count);
        memset(hz_buf->lower_hz_disp, 0, sizeof(int) * buffer_count);

        hz_buf->lower_hz_count = malloc(sizeof(int) * buffer_count);
        memset(hz_buf->lower_hz_count, 0, sizeof(int) * buffer_count);

        hz_buf->lower_hz_buffer = malloc(bytes_for_datatype * hz_buf->lower_hz_buffer_size);
        memset(hz_buf->lower_hz_buffer, 0, bytes_for_datatype * hz_buf->lower_hz_buffer_size);

        buffer_count = 0;
        int prev_count = 0;
        int block_no = 0;
        int negative_block_offset = 0;
        for (i = /*block_layout->resolution_from*/start_index; i < /*block_layout->resolution_to*/end_index; i++)
        //for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            hz_buf->lower_hz_disp[buffer_count] = hz_buf->start_hz_index[i] * bytes_for_datatype;
            hz_buf->lower_hz_count[buffer_count] = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1) * bytes_for_datatype;

            block_no = hz_buf->start_hz_index[i] / (agg_id->idx_d->samples_per_block);
            negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx->blocks_per_file, block_no, block_layout);
            hz_buf->lower_hz_disp[buffer_count] = hz_buf->lower_hz_disp[buffer_count] - (negative_block_offset * bytes_for_datatype * agg_id->idx_d->samples_per_block);

            //if (rank == 0)
            //  printf("PC %d [%d %d] [%d %d] O %d C %d\n", prev_count, i, buffer_count, block_layout->resolution_from, block_layout->resolution_to, hz_buf->lower_hz_disp[buffer_count], hz_buf->lower_hz_count[buffer_count]);

            memcpy(hz_buf->lower_hz_buffer + prev_count, hz_buf->buffer[i],  (hz_buf->lower_hz_count[buffer_count]));
            prev_count = prev_count + hz_buf->lower_hz_count[buffer_count];

            buffer_count++;
          }
        }
      }


      //int b = 0;
      //if (rank == 0)
      //  for (b = 0; b < buffer_count; b++)
      //    printf("[%d] [%d] O %d C %d\n", rank, b, hz_buf->lower_hz_disp[b] / bytes_for_datatype, hz_buf->lower_hz_count[b] / bytes_for_datatype);

      MPI_Type_indexed( buffer_count, hz_buf->lower_hz_count, hz_buf->lower_hz_disp, MPI_BYTE, &(hz_buf->lower_level_datatype));
      MPI_Type_commit(&(hz_buf->lower_level_datatype));
    }
  }

  return PIDX_success;
}
#endif


static PIDX_return_code destroy_file_zero_buffer(PIDX_local_agg_id agg_id)
{
  int p = 0, v = 0;
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      HZ_buffer hz_buf = var->hz_buffer[p];

      free(hz_buf->lower_hz_disp);
      hz_buf->lower_hz_disp = 0;

      free(hz_buf->lower_hz_count);
      hz_buf->lower_hz_count = 0;

      free(hz_buf->lower_hz_buffer);
      hz_buf->lower_hz_buffer = 0;
#if PIDX_HAVE_MPI
      MPI_Type_free(&hz_buf->lower_level_datatype);
#endif
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg(PIDX_local_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout block_layout, int MODE)
{
#if PIDX_HAVE_MPI
  //int file_zero = 1;
  int ret;

  ret = create_window(agg_id, agg_buffer, agg_id->comm);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

  ret = create_open_log_file(agg_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] PIDX error (create_open_log_file).\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

  if (agg_id->idx_d->file_zero_optimization == 1)
  {
    if (layout_id == 0)
    {
      ret = create_file_zero_buffer(agg_id, block_layout, block_layout->resolution_from, block_layout->resolution_to);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] PIDX error (create_file_zero_buffer).\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
  }

#if !SIMULATE_IO
#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, agg_id->win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif
#endif

  if (MODE == PIDX_WRITE)
  {
    if (agg_id->idx_d->file_zero_optimization == 1)
    {
      if (layout_id == 0)
      {
        ret = layout_zero(agg_id, block_layout);
        if (ret != PIDX_success)
        {
          fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        ret = local_one_sided_data_com(agg_id, agg_buffer, layout_id, block_layout, PIDX_WRITE);
        if (ret != PIDX_success)
        {
          fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
          return PIDX_err_agg;
        }
      }
    }
    else
    {
      ret = local_one_sided_data_com(agg_id, agg_buffer, layout_id, block_layout, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
  }
  else
  {
    ret = local_one_sided_data_com(agg_id, agg_buffer, layout_id, block_layout, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }

#if !SIMULATE_IO
#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, agg_id->win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif
#endif
#endif

#if !SIMULATE_IO
#if PIDX_HAVE_MPI
  ret = MPI_Win_free(&(agg_id->win));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " [%s] [%d] Window create error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif
#endif

  if (agg_id->idx_d->file_zero_optimization == 1)
  {
    if (layout_id == 0)
    {
      ret = destroy_file_zero_buffer(agg_id);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] PIDX error (create_file_zero_buffer).\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
  }

  ret = close_log_file(agg_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] PIDX error (create_open_log_file).\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_finalize(PIDX_local_agg_id agg_id)
{
  free(agg_id);
  return PIDX_success;
}
