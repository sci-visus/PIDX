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

//static PIDX_return_code create_window(PIDX_agg_id agg_id);
//static PIDX_return_code one_sided_data_com(PIDX_agg_id agg_id, int mode);

static PIDX_return_code aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout,  int MODE);

static PIDX_return_code local_aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE, MPI_Comm comm);

static PIDX_return_code local_aggregate_local_comm_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout,  int MODE, MPI_Comm comm);

static PIDX_return_code local_aggregate_write_read_layout_zero(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, MPI_Datatype, int shuffle, int MODE);

struct PIDX_agg_struct
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

  //int agg_existing_file_count;

  int ***rank_holder;
  int ***rank_holder2;
  int ***rank_holder3;
};

static PIDX_return_code create_window(PIDX_agg_id agg_id, Agg_buffer agg_buffer, MPI_Comm comm)
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


static PIDX_return_code one_sided_data_com(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int mode)
{
  int i, p, e1, v, ret = 0;
  int send_index = 0;
  int64_t index = 0, count = 0, hz_index = 0;
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
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      hz_index = 0, index = 0, count = 0, send_index = 0;
      HZ_buffer hz_buf = var->hz_buffer[p];
      //printf("hz_buf->type = %d\n", hz_buf->type);
      if (hz_buf->type == 0)
      {
        for (i = 0; i < hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i++)
          hz_index = hz_index + hz_buf->samples_per_level[i];

        for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        //for (i = /*hz_buf->HZ_agg_from + agg_id->idx_d->res_from*/(agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)); i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
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
                  ret = aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
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
                    aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
                    }

                  }
                }
                else
                {
                  aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
                  }

                  if(e1 == hz_buf->samples_per_level[i] - 1)
                  {
                    aggregate_write_read(agg_id, v, hz_buf->buffer_index[hz_index], 1, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, agg_buffer, block_layout, PIDX_WRITE);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
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

      else if (hz_buf->type == 1)
      {
#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
          fflush(agg_dump_fp);
        }
#endif
        //for (i = 0; i < agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file) + 1; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 1; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 1; i < (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 2; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 2; i < (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 3; i++)
        //for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        //printf("[%d] from %d to %d\n", rank, block_layout->resolution_from, block_layout->resolution_to);
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

            //double x;
            //memcpy(&x, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], sizeof(double));
            //printf("offset %d count %d val %f\n", hz_buf->start_hz_index[i], count, x);
            ret = aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, agg_buffer, block_layout, mode);
#else
            ret = aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode);
#endif
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
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
        //for (i = 0; i < agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file) + 1; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 1; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 1; i < (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 2; i++)
        //for (i = (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 2; i < (agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)) + 3; i++)
        //for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
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
              index = 0;
              count = (agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);

#if !SIMULATE_IO
              ret = aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, agg_buffer, block_layout, mode);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
#else
              ret = aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
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
                //if (PIDX_blocks_is_block_present(bl, agg_id->idx->variable[agg_id->init_index]->global_block_layout))
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
                  ret = aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, mode);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
#else
                  ret = aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, NULL, send_index, agg_buffer, block_layout, mode);
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

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}


static PIDX_return_code layout_zero(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int shuffle, MPI_Datatype chunk_data_type, int mode)
{

  int rank = 0;
#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif
  int v = 0, p = 0, i = 0;
  int total_count = 0, count = 0, bytes_for_datatype = 0, index = 0;
  int ret = 0;

  int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    total_count = 0;
    for (p = 0; p < var->patch_group_count; p++)
    {
      bytes_for_datatype = ((var->bits_per_value / 8) * chunk_size * var->values_per_sample) / (agg_id->idx->compression_factor);
      HZ_buffer hz_buf = var->hz_buffer[p];

      index = 0;
      count = 0;
      total_count = 0;
      for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
      {
        if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
        {
          count =  hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;
          memcpy(hz_buf->lower_hz_buffer + index * bytes_for_datatype, hz_buf->buffer[i],  count * bytes_for_datatype);
          total_count = total_count + count;
          index = total_count;
        }
      }

#if !SIMULATE_IO
      ret = local_aggregate_write_read_layout_zero(agg_id, v, /*hz_buf->start_hz_index[i]*/0, total_count, agg_id->idx->variable[v]->hz_buffer[p]->lower_hz_buffer, 0, agg_buffer, block_layout, chunk_data_type, shuffle, mode);
#else
      ret = local_aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode);
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


static PIDX_return_code local_one_sided_data_com(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int mode)
{
  int i, p, e1, v, ret = 0;
  int send_index = 0;
  int64_t index = 0, count = 0, hz_index = 0;
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
      hz_index = 0, index = 0, count = 0, send_index = 0;
      HZ_buffer hz_buf = var->hz_buffer[p];
      //printf("hz_buf->type = %d\n", hz_buf->type);
      if (hz_buf->type == 0)
      {
        for (i = 0; i < hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i++)
          hz_index = hz_index + hz_buf->samples_per_level[i];

        for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
        //for (i = /*hz_buf->HZ_agg_from + agg_id->idx_d->res_from*/(agg_id->idx->bits_per_block + log2(agg_id->idx->blocks_per_file)); i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
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
                  ret = local_aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE, agg_id->comm);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
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
                    local_aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE, agg_id->comm);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
                    }

                  }
                }
                else
                {
                  local_aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_WRITE, agg_id->comm);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
                  }

                  if(e1 == hz_buf->samples_per_level[i] - 1)
                  {
                    local_aggregate_write_read(agg_id, v, hz_buf->buffer_index[hz_index], 1, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, agg_buffer, block_layout, PIDX_WRITE, agg_id->comm);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
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

      else if (hz_buf->type == 1)
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

            ret = local_aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->comm);
#else
            ret = local_aggregate_write_read(agg_id, v, hz_buf->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode);
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
              index = 0;
              count = (agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);

#if !SIMULATE_IO
              ret = local_aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->comm);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
#else
              ret = local_aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode, agg_id->comm);
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
                //if (PIDX_blocks_is_block_present(bl, agg_id->idx->variable[agg_id->init_index]->global_block_layout))
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
                  ret = local_aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, mode, agg_id->comm);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
#else
                  ret = local_aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, NULL, send_index, agg_buffer, block_layout, mode, agg_id->comm);
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



static PIDX_return_code local_one_sided_data_com_local_comm(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int mode)
{
  int i, p, e1, v, ret = 0;
  int send_index = 0;
  int64_t index = 0, count = 0, hz_index = 0;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->local_comm, &rank);
  MPI_Comm_size(agg_id->local_comm, &nprocs);
  //printf("[R %d] Comm Size %d\n", rank, nprocs);
#endif

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      hz_index = 0, index = 0, count = 0, send_index = 0;
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

            ret = local_aggregate_local_comm_write_read(agg_id, v, hz_buf->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->local_comm);
#else
            ret = local_aggregate_local_comm_write_read(agg_id, v, hz_buf->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode);
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
              index = 0;
              count = (agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);

#if !SIMULATE_IO
              ret = local_aggregate_local_comm_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, agg_buffer, block_layout, mode, agg_id->local_comm);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in local_aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
#else
              ret = local_aggregate_local_comm_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, mode, agg_id->local_comm);
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
                //if (PIDX_blocks_is_block_present(bl, agg_id->idx->variable[agg_id->init_index]->global_block_layout))
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
                  ret = local_aggregate_local_comm_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, mode, agg_id->local_comm);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
#else
                  ret = local_aggregate_local_comm_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, NULL, send_index, agg_buffer, block_layout, mode, agg_id->local_comm);
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



static PIDX_return_code aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE)
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

  target_rank = agg_id->rank_holder[file_no][variable_index - agg_id->first_index][sample_index];

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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype);

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
          memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype);
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
        memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype);

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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }

  return PIDX_success;
}



static PIDX_return_code local_aggregate_write_read_layout_zero(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, MPI_Datatype chunk_type, int shuffle, int MODE)
{
  int rank = 0;// nrank = 0;
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

  target_count = hz_count * values_per_sample;

  bytes_per_datatype = ((agg_id->idx->variable[variable_index]->bits_per_value / 8) * total_chunk_size) / (agg_id->idx->compression_factor);

#if !SIMULATE_IO
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;
#endif

  start_agg_index = target_disp / (int64_t) (samples_in_file / agg_buffer->aggregation_factor);
  end_agg_index = ((target_disp + target_count - 1) / (int64_t) (samples_in_file / agg_buffer->aggregation_factor));
  //assert(start_agg_index >= 0 && end_agg_index >= 0 && end_agg_index >= start_agg_index);

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
    //printf("[%d] Count %d Offset %d\n", rank, hz_count * values_per_sample * bytes_per_datatype, rank * hz_count * values_per_sample, target_rank);

    if (shuffle == 1)
      ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, 1, chunk_type, agg_id->win);
    else
      ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, rank * hz_count * values_per_sample, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
    if(ret != MPI_SUCCESS)
    {
      fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_agg;
    }
#endif
  }
  else
  {
    ret = MPI_Get(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, 1, chunk_type, agg_id->win);
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


  return PIDX_success;
}


static PIDX_return_code local_aggregate_local_comm_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE, MPI_Comm comm)
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
  MPI_Comm_rank(comm, &rank);
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

  //printf("file index = %d\n", block_layout->inverse_existing_file_index[file_no]);
  //target_rank = agg_id->rank_holder3[block_layout->inverse_existing_file_index[file_no]][variable_index - agg_id->first_index][sample_index];
  target_rank = agg_id->rank_holder3[0][variable_index - agg_id->first_index][sample_index];

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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype);

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
          memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype);
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
        memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype);

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
        //printf("[%d] XX Count %d Offset %d\n", rank, hz_count * values_per_sample * bytes_per_datatype, target_disp);
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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }

  return PIDX_success;
}



static PIDX_return_code local_aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer agg_buffer, PIDX_block_layout block_layout, int MODE, MPI_Comm comm)
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
  MPI_Comm_rank(comm, &rank);
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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_buffer->aggregation_factor) - target_disp) * bytes_per_datatype);

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
          memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (samples_in_file / agg_buffer->aggregation_factor) * bytes_per_datatype);
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
        memcpy( hz_buffer + (((samples_in_file / agg_buffer->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_buffer->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_buffer->aggregation_factor) - target_disp)) * bytes_per_datatype);

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
        //printf("[%d] XX Count %d Offset %d\n", rank, hz_count * values_per_sample * bytes_per_datatype, target_disp);
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
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }

  return PIDX_success;
}


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


PIDX_return_code PIDX_create_local_aggregation_comm(PIDX_agg_id agg_id)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  int color;
  int rank;
  int nprocs;

  MPI_Comm temp_comm = agg_id->comm;
  MPI_Comm_rank(temp_comm, &rank);
  MPI_Comm_size(agg_id->comm, &nprocs);
  //printf("[Before] Agg size = %d\n", nprocs);

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  if (var0->patch_group_count == 0)
    color = 0;
  else
    color = 1;

  MPI_Comm_split(temp_comm, color, rank, &(agg_id->comm));
  MPI_Comm_size(agg_id->comm, &nprocs);

  return PIDX_success;
}

PIDX_return_code PIDX_destroy_local_aggregation_comm(PIDX_agg_id agg_id)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  MPI_Comm_free(&(agg_id->comm));

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


PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout)
{
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  int i, j, d = 0;

  agg_buffer->aggregator_interval = nprocs / ((agg_id->last_index - agg_id->first_index + 1) * local_block_layout->existing_file_count * agg_buffer->aggregation_factor);

  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_id->rank_holder = malloc(agg_id->idx_d->max_file_count * sizeof (int**));
  memset(agg_id->rank_holder, 0, agg_id->idx_d->max_file_count * sizeof (int**));
  for (i = 0; i < agg_id->idx_d->max_file_count; i++)
  {
    agg_id->rank_holder[i] = malloc((agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    memset(agg_id->rank_holder[i], 0, (agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      agg_id->rank_holder[i][j - agg_id->first_index] = malloc(agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      memset(agg_id->rank_holder[i][j - agg_id->first_index], 0, agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      for (d = 0 ; d < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; d++)
        agg_id->rank_holder[i][j - agg_id->first_index][d] = -1;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout)
{
  int i = 0, j = 0, d = 0;
  int rank = 0, nprocs = 1;
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  agg_buffer->aggregator_interval = nprocs / ((agg_id->last_index - agg_id->first_index + 1) * local_block_layout->existing_file_count * agg_buffer->aggregation_factor);
  //assert(agg_buffer->aggregator_interval != 0);

  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_id->rank_holder2 = malloc(local_block_layout->existing_file_count * sizeof (int**));
  memset(agg_id->rank_holder2, 0, local_block_layout->existing_file_count * sizeof (int**));
  for (i = 0; i < local_block_layout->existing_file_count; i++)
  {
    agg_id->rank_holder2[i] = malloc((agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    memset(agg_id->rank_holder2[i], 0, (agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      agg_id->rank_holder2[i][j - agg_id->first_index] = malloc(agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      memset(agg_id->rank_holder2[i][j - agg_id->first_index], 0, agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      for (d = 0 ; d < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; d++)
        agg_id->rank_holder2[i][j - agg_id->first_index][d] = -1;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_local_comm_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout)
{
  int i = 0, j = 0, d = 0;
  int rank = 0, nprocs = 1;
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  agg_buffer->aggregator_interval = nprocs / ((agg_id->last_index - agg_id->first_index + 1) * 1 * agg_buffer->aggregation_factor);
  //assert(agg_buffer->aggregator_interval != 0);

  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_id->rank_holder3 = malloc(1 * sizeof (int**));
  memset(agg_id->rank_holder3, 0, 1 * sizeof (int**));
  for (i = 0; i < 1; i++)
  {
    agg_id->rank_holder3[i] = malloc((agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    memset(agg_id->rank_holder3[i], 0, (agg_id->last_index - agg_id->first_index + 1) * sizeof (int*));
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      agg_id->rank_holder3[i][j - agg_id->first_index] = malloc(agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      memset(agg_id->rank_holder3[i][j - agg_id->first_index], 0, agg_id->idx->variable[j]->values_per_sample * sizeof (int) * agg_buffer->aggregation_factor);
      for (d = 0 ; d < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; d++)
        agg_id->rank_holder3[i][j - agg_id->first_index][d] = -1;
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id agg_id)
{
  int i = 0, j = 0;
  for (i = 0; i < agg_id->idx_d->max_file_count; i++)
  {
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      free(agg_id->rank_holder[i][j - agg_id->first_index]);
      agg_id->rank_holder[i][j - agg_id->first_index] = 0;
    }
    free(agg_id->rank_holder[i]);
  }
  free(agg_id->rank_holder);
  agg_id->rank_holder = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_local_agg_meta_data_destroy(PIDX_agg_id agg_id, PIDX_block_layout local_block_layout)
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


PIDX_return_code PIDX_local_agg_local_comm_meta_data_destroy(PIDX_agg_id agg_id, PIDX_block_layout local_block_layout)
{
  int i = 0, j = 0;
  for (i = 0; i < 1; i++)
  {
    for (j = agg_id->first_index; j <= agg_id->last_index; j++)
    {
      free(agg_id->rank_holder3[i][j - agg_id->first_index]);
      agg_id->rank_holder3[i][j - agg_id->first_index] = 0;
    }
    free(agg_id->rank_holder3[i]);
  }
  free(agg_id->rank_holder3);
  agg_id->rank_holder3 = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_agg_buf_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout, int i1, int j1)
{

  int rank = 0, nprocs = 1;
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  int rank_counter = 0, i = 0, j = 0, k = 0;
  for (k = 0; k < local_block_layout->existing_file_count; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_buffer->aggregation_factor; j++)
      {
        //agg_id->rank_holder[agg_id->idx->variable[agg_id->init_index]->existing_file_index[k]][i - agg_id->first_index][j] = rank_counter;
        agg_id->rank_holder[global_block_layout->existing_file_index[k]][i - agg_id->first_index][j] = rank_counter;
        rank_counter = rank_counter + agg_buffer->aggregator_interval;

        //if(rank == agg_id->rank_holder[agg_id->idx->variable[agg_id->init_index]->existing_file_index[k]][i - agg_id->first_index][j])
        if(rank == agg_id->rank_holder[global_block_layout->existing_file_index[k]][i - agg_id->first_index][j])
        {
          //agg_buffer->file_number = agg_id->idx->variable[agg_id->init_index]->existing_file_index[k];
          agg_buffer->file_number = global_block_layout->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;

          //uint64_t sample_count = agg_id->idx->variable[agg_id->init_index]->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_buffer->aggregation_factor;
          uint64_t sample_count = local_block_layout->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_buffer->aggregation_factor;

          int total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

          int bytes_per_datatype = (total_chunk_size * agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / ( agg_id->idx->compression_factor);

          agg_buffer->buffer_size = sample_count * bytes_per_datatype;
          printf("O [%d] [%d %d %d] buffer_size = %d (%d %d)\n", rank, agg_buffer->file_number, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->buffer_size, i1, j1);

#if !SIMULATE_IO
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
    for (i = 0; i < agg_id->idx_d->agg_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        for (k = 0; k < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; k++)
        {
          printf("[%d %d %d] -> %d\n", k, i, j, agg_id->rank_holder[i][j-agg_id->first_index][k]);
        }
      }
    }
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_buf_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int agg_offset)
{
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  int rank_counter = 0, i = 0, j = 0, k = 0;
  rank_counter = agg_offset;
  for (k = 0; k < local_block_layout->existing_file_count; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_buffer->aggregation_factor; j++)
      {
        agg_id->rank_holder2[k][i - agg_id->first_index][j] = rank_counter;
        rank_counter = rank_counter + agg_buffer->aggregator_interval;

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

#if !SIMULATE_IO
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
    for (i = 0; i < agg_id->idx_d->agg_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        for (k = 0; k < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; k++)
        {
          printf("[%d %d %d] -> %d\n", k, i, j, agg_id->rank_holder[i][j-agg_id->first_index][k]);
        }
      }
    }
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg_local_comm_buf_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, int agg_offset)
{
  int rank = 0, nprocs = 1;

  int rank_counter = 0, i = 0, j = 0, k = 0, p = 0, v = 0;
  rank_counter = agg_offset;

  int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  int file_no = 0, prev_file_no = 0;

  if (var0->patch_group_count > 1)
  {
    printf("Not Implemented");
    return PIDX_err_agg;
  }

  if (var0->patch_group_count != 0)
  {
    HZ_buffer hz_buf0 = var0->hz_buffer[0];
    prev_file_no = hz_buf0->start_hz_index[local_block_layout->resolution_from] / samples_per_file;
    file_no = hz_buf0->start_hz_index[local_block_layout->resolution_from] / samples_per_file;

    for (p = 1; p < var0->patch_group_count; p++)
    {
      HZ_buffer hz_buf = var0->hz_buffer[p];
      file_no = hz_buf->start_hz_index[local_block_layout->resolution_from] / samples_per_file;

      if (prev_file_no != file_no)
      {
        printf("Not Implemented");
        return PIDX_err_agg;
      }
    }
  }

  //printf("[%d] File Number %d\n", rank, file_no);
  MPI_Comm_split(agg_id->comm, file_no, rank, &(agg_id->local_comm));
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->local_comm, &nprocs);
  MPI_Comm_rank(agg_id->local_comm, &rank);
#endif

  for (k = 0; k < 1; k++)
  {
    for (i = agg_id->first_index; i <= agg_id->last_index; i++)
    {
      for (j = 0; j < agg_id->idx->variable[i]->values_per_sample * agg_buffer->aggregation_factor; j++)
      {
        agg_id->rank_holder3[k][i - agg_id->first_index][j] = rank_counter;
        rank_counter = rank_counter + agg_buffer->aggregator_interval;

        if(rank == agg_id->rank_holder3[k][i - agg_id->first_index][j])
        {
          agg_buffer->file_number = file_no;
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;

          uint64_t sample_count = local_block_layout->block_count_per_file[agg_buffer->file_number] * agg_id->idx_d->samples_per_block / agg_buffer->aggregation_factor;

          int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];

          int bytes_per_datatype = 0;
          bytes_per_datatype = (chunk_size * agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (agg_id->idx->compression_factor);

          agg_buffer->buffer_size = sample_count * bytes_per_datatype;

#if !SIMULATE_IO
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

  //int target_rank = agg_id->rank_holder3[local_block_layout->inverse_existing_file_index[file_no]][v - agg_id->first_index][0];
  //MPI_Comm_split(agg_id->comm, file_no, rank, &(agg_id->local_comm));

#if 0
  if (rank == 0)
  {
    for (i = 0; i < agg_id->idx_d->agg_file_count; i++)
    {
      for (j = agg_id->first_index; j <= agg_id->last_index; j++)
      {
        for (k = 0; k < agg_id->idx->variable[j]->values_per_sample * agg_buffer->aggregation_factor; k++)
        {
          printf("[%d %d %d] -> %d\n", k, i, j, agg_id->rank_holder[i][j-agg_id->first_index][k]);
        }
      }
    }
  }
#endif

  return PIDX_success;
}



PIDX_return_code PIDX_agg_buf_destroy(PIDX_agg_id agg_id, Agg_buffer agg_buffer)
{

#if !SIMULATE_IO
  if (agg_buffer->buffer_size != 0)
  {
    free(agg_buffer->buffer);
    agg_buffer->buffer = 0;
  }
#endif

  MPI_Comm_free(&(agg_id->local_comm));

  return PIDX_success;
}


PIDX_return_code PIDX_agg_write(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout)
{
  int ret;
  ret = create_window(agg_id, agg_buffer, agg_id->comm);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
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

  //printf("Final Block Bitmap\n");
  //PIDX_blocks_print_layout(block_layout);

  ret = one_sided_data_com(agg_id, agg_buffer, block_layout, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
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

  return PIDX_success;
}


PIDX_return_code PIDX_local_agg(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout block_layout, int MODE)
{
  int file_zero = 0;
  int shuffle = 1;
  int ret;
  MPI_Datatype chunk_data_type;
  ret = create_window(agg_id, agg_buffer, agg_id->comm);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

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

  if (file_zero == 1)
  {
    if (layout_id == 0)
    {
      int p = 0, i = 0, samples_in_file = 0, bytes_for_datatype = 0, v = 0;
      int64_t samples_per_file = (int64_t) agg_id->idx_d->samples_per_block * agg_id->idx->blocks_per_file;
      int count = 0, total_count = 0, file_no = 0;
      int chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
      for(v = agg_id->first_index; v <= agg_id->last_index; v++)
      {
        PIDX_variable var = agg_id->idx->variable[v];
        total_count = 0;
        for (p = 0; p < var->patch_group_count; p++)
        {
          bytes_for_datatype = ((var->bits_per_value / 8) * chunk_size * var->values_per_sample) / (agg_id->idx->compression_factor);
          HZ_buffer hz_buf = var->hz_buffer[p];

          if (shuffle == 1)
          {
            hz_buf->lower_hz_disp = malloc(sizeof(int) * (block_layout->resolution_to - block_layout->resolution_from));
            memset(hz_buf->lower_hz_disp, 0, sizeof(int) * (block_layout->resolution_to - block_layout->resolution_from));

            hz_buf->lower_hz_count = malloc(sizeof(int) * (block_layout->resolution_to - block_layout->resolution_from));
            memset(hz_buf->lower_hz_count, 0, sizeof(int) * (block_layout->resolution_to - block_layout->resolution_from));
          }

          for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
          {
            if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
            {
              count =  hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;
              total_count = total_count + count;

              file_no = hz_buf->start_hz_index[i] / samples_per_file;
              samples_in_file = block_layout->block_count_per_file[file_no] * agg_id->idx_d->samples_per_block;

              if (shuffle == 1)
              {
                hz_buf->lower_hz_disp[i] = (hz_buf->start_hz_index[i] - samples_per_file * file_no) % samples_in_file;
                hz_buf->lower_hz_count[i] = count;
              }
            }
          }

          hz_buf->lower_hz_buffer = malloc(bytes_for_datatype * total_count);
          memset(hz_buf->lower_hz_buffer, 0, bytes_for_datatype * total_count);
        }
      }

      if (shuffle == 1)
      {
        PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
        HZ_buffer hz_buf0 = var0->hz_buffer[0];
        MPI_Type_indexed( (block_layout->resolution_to - block_layout->resolution_from), hz_buf0->lower_hz_count, hz_buf0->lower_hz_disp, MPI_DOUBLE, &chunk_data_type);
        MPI_Type_commit(&chunk_data_type);
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

  //printf("Final Block Bitmap\n");
  //PIDX_blocks_print_layout(block_layout);

  if (file_zero == 1)
  {
    if (layout_id == 0)
    {
      ret = layout_zero(agg_id, agg_buffer, block_layout, shuffle, chunk_data_type, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
    else
    {
      ret = local_one_sided_data_com(agg_id, agg_buffer, block_layout, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
  }
  else
  {
    if (MODE == PIDX_WRITE)
    {
      ret = local_one_sided_data_com(agg_id, agg_buffer, block_layout, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
    }
    else
    {
      ret = local_one_sided_data_com(agg_id, agg_buffer, block_layout, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }
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

  if (file_zero == 1)
  {
    if (layout_id == 0)
    {
      int p = 0, v = 0;
      for(v = agg_id->first_index; v <= agg_id->last_index; v++)
      {
        PIDX_variable var = agg_id->idx->variable[v];
        for (p = 0; p < var->patch_group_count; p++)
        {
          HZ_buffer hz_buf = var->hz_buffer[p];

          if (shuffle == 1)
          {
            free(hz_buf->lower_hz_disp);
            hz_buf->lower_hz_disp = 0;

            free(hz_buf->lower_hz_count);
            hz_buf->lower_hz_count = 0;
          }

          free(hz_buf->lower_hz_buffer);
          hz_buf->lower_hz_buffer = 0;
        }
      }

      if (shuffle == 1)
        MPI_Type_free(&chunk_data_type);
    }
  }

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}



PIDX_return_code PIDX_local_agg_local_comm(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout block_layout, int MODE)
{
  int ret;
  ret = create_window(agg_id, agg_buffer, agg_id->local_comm);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

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

  //printf("Final Block Bitmap\n");
  //PIDX_blocks_print_layout(block_layout);

  if (MODE == PIDX_WRITE)
  {
    ret = local_one_sided_data_com_local_comm(agg_id, agg_buffer, block_layout, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }
  else
  {
    ret = local_one_sided_data_com_local_comm(agg_id, agg_buffer, block_layout, PIDX_READ);
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

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}




PIDX_return_code PIDX_agg_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout)
{

  int ret;
  ret = create_window(agg_id, agg_buffer, agg_id->comm);
  if (ret != PIDX_success)
  {
    fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }


  agg_id->idx_d->staged_aggregation = 0;
  if (agg_id->idx_d->staged_aggregation == 0)
  {

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

    ret = one_sided_data_com(agg_id, agg_buffer, block_layout, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr, " [%s] [%d] Fence error.\n", __FILE__, __LINE__);
      return PIDX_err_agg;
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

  }

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

  return PIDX_success;
}


PIDX_return_code PIDX_agg_read1(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout block_layout)
{
  int i, p, v, ret = 0;
#if 0
  int send_index = 0;
  int e1;
  int64_t hz_index = 0, index = 0;
#endif
  int64_t count = 0;
  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);

  if (agg_buffer->buffer_size != 0)
  {
    int total_chunk_size = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2] * agg_id->idx->chunk_size[3] * agg_id->idx->chunk_size[4];
    int bytes_per_datatype = total_chunk_size * (agg_id->idx->variable[agg_buffer->var_number]->bits_per_value/8) / (agg_id->idx->compression_factor);

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
  for(v = agg_id->first_index; v <= agg_id->last_index; v++)
  {
    PIDX_variable var = agg_id->idx->variable[agg_id->first_index];
    for (p = 0; p < var->patch_group_count; p++)
    {
      count = 0;
#if 0
      send_index = 0;
      hz_index = 0, index = 0;
#endif
      HZ_buffer hz_buf = var->hz_buffer[p];
      if(var->hz_buffer[p]->type == 0)
      {
#if 0
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
                  ret = aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, PIDX_READ);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
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
                    aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, PIDX_READ);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
                    }
                  }
                }
                else
                {
                  aggregate_write_read(agg_id, v, index, count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, PIDX_READ);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return PIDX_err_agg;
                  }

                  if(e1 == hz_buf->samples_per_level[i] - 1)
                  {
                    aggregate_write_read(agg_id, v, hz_buf->buffer_index[hz_index], 1, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, agg_buffer,  PIDX_READ);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return PIDX_err_agg;
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
#endif
      }
      else if (hz_buf->type == 1)
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
          if ((hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2]) != 0)
          {
            count =  var->hz_buffer[p]->end_hz_index[i] - var->hz_buffer[p]->start_hz_index[i] + 1;
#ifdef PIDX_DUMP_AGG
            if (agg_id->idx_d->dump_agg_info == 1 && agg_id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif
            //printf("offset %d count %d\n", var->hz_buffer[p]->start_hz_index[i], count);
#if !SIMULATE_IO
            ret = aggregate_write_read(agg_id, v, var->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], 0, agg_buffer, block_layout, PIDX_READ);
#else
            ret = aggregate_write_read(agg_id, v, var->hz_buffer[p]->start_hz_index[i], count, NULL, 0, agg_buffer, block_layout, PIDX_READ);
#endif
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
        }
      }
      else if (hz_buf->type == 2)
      {
        uint64_t index = 0, count = 0;
        for(v = agg_id->first_index; v <= agg_id->last_index; v++)
        {
          HZ_buffer hz_buf = agg_id->idx->variable[v]->hz_buffer[p];
          for (i = hz_buf->HZ_agg_from + agg_id->idx_d->res_from; i < hz_buf->HZ_agg_to - agg_id->idx_d->res_to; i++)
          {
            if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
            {
              int start_block_index = agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] / agg_id->idx_d->samples_per_block;
              int end_block_index = agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] / agg_id->idx_d->samples_per_block;
              assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

              if (end_block_index == start_block_index)
              {
                index = 0;
                count = (agg_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);
                //printf("A [%d] offset 0 Count %lld\n", i, (unsigned long long)count);
                ret = aggregate_write_read(agg_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, agg_buffer, block_layout, PIDX_READ);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                  return PIDX_err_io;
                }
              }
              else
              {
                int send_index = 0;
                int bl;
                for (bl = start_block_index; bl <= end_block_index; bl++)
                {
                  //if (PIDX_blocks_is_block_present(bl, agg_id->idx->variable[agg_id->init_index]->global_block_layout))
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

                    //printf("B [%d] offset %lld send offset %lld Count %lld\n", i, (unsigned long long)index, (unsigned long long)send_index, (unsigned long long)count);
                    ret = aggregate_write_read(agg_id, v, index + agg_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, agg_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, agg_buffer, block_layout, PIDX_READ);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_io;
                    }
                    send_index = send_index + count;
                  }
                }
              }
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

  //double x;
  //memcpy(&x, agg_buf->buffer, sizeof(double));
  //printf("first elemet = %f\n", x);

  return PIDX_success;
}


PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id)
{

  free(agg_id);
  agg_id = 0;

  return PIDX_success;
}
