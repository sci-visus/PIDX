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

#undef PIDX_RECORD_TIME
//#define RANK_ORDER 1
static uint32_t *cached_header_copy;
static int enable_caching = 0;

#define PIDX_DUMP_IO

#ifdef PIDX_DUMP_IO
static FILE* io_dump_fp;
#endif

enum IO_MODE { PIDX_READ, PIDX_WRITE};

struct PIDX_io_struct
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Win win;
#endif

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int init_index;
  int first_index;
  int last_index;
};

static int write_read_samples(PIDX_io_id io_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int64_t buffer_offset, int MODE);


PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_io_id)malloc(sizeof (*io_id));
  memset(io_id, 0, sizeof (*io_id));

  io_id->idx = idx_meta_data;
  io_id->idx_d = idx_d;

  io_id->init_index = init_index;
  io_id->first_index = first_index;
  io_id->last_index = last_index;

  return io_id;
}

#if PIDX_HAVE_MPI
int PIDX_io_set_communicator(PIDX_io_id io_id, MPI_Comm comm)
{
  if (io_id == NULL)
    return PIDX_err_id;

  io_id->comm = comm;

  return PIDX_success;
}
#endif

int PIDX_io_cached_data(uint32_t* cached_header)
{
  cached_header_copy = cached_header;
  enable_caching = 1;

  return PIDX_success;
}


int PIDX_io_aggregated_write(PIDX_io_id io_id)
{
  int64_t data_offset = 0;  
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  uint32_t *headers;
  int total_header_size = 0;
  double t1, t2, t3, t4, t5;
  int rank = 0;

#if PIDX_HAVE_MPI
  int mpi_ret;
  MPI_File fh;
  MPI_Status status;
  MPI_Comm_rank(io_id->comm, &rank);
#else
  int fh;
#endif

  int total_chunk_size = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]);

  if (io_id->idx->enable_agg != 0)
  {
    Agg_buffer agg_buf = io_id->idx_d->agg_buffer;

    if (enable_caching == 1 && agg_buf->var_number == io_id->init_index && agg_buf->sample_number == 0)
    {
      t1 = PIDX_get_time();

      generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) agg_buf->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed filename %s.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
#else
      fh = open(file_name, O_WRONLY);
#endif

      t2 = PIDX_get_time();

      data_offset = 0;
      total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
      headers = (uint32_t*)malloc(total_header_size);
      memset(headers, 0, total_header_size);


      if (enable_caching == 1)
        memcpy (headers, cached_header_copy, total_header_size);
      else
      {
        //TODO
      }

      t3 = PIDX_get_time();

      uint64_t header_size = (io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size);

      unsigned char* temp_buffer = (unsigned char*)realloc(agg_buf->buffer, agg_buf->buffer_size  + header_size);

      if (temp_buffer == NULL)
      {
        fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      else
      {
        agg_buf->buffer = temp_buffer;
        memmove(agg_buf->buffer + header_size, agg_buf->buffer, agg_buf->buffer_size);
        memcpy(agg_buf->buffer, headers, total_header_size);
        memset(agg_buf->buffer + total_header_size, 0, (header_size - total_header_size));
      }
      free(headers);
#if PIDX_HAVE_MPI

      mpi_ret = MPI_File_write_at(fh, 0, agg_buf->buffer, agg_buf->buffer_size + header_size, MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }

      int write_count;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != agg_buf->buffer_size + header_size)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

#else
      ssize_t write_count = pwrite(fh, agg_buf->buffer, agg_buf->buffer_size + header_size, 0);
      if (write_count != agg_buf->buffer_size + header_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#endif

      t4 = PIDX_get_time();

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#else
      close(fh);
#endif

      t5 = PIDX_get_time();

#ifdef PIDX_RECORD_TIME
      printf("V0. [R %d] [O 0 C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long)agg_buf->buffer_size + header_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t3-t2), (t4-t3), (t5-t4));
#else
#endif
    }
    else if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
    {
      t1 = PIDX_get_time();

      generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) agg_buf->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
#else
      fh = open(file_name, O_WRONLY);
#endif

      t2 = PIDX_get_time();

      data_offset = 0;
      data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;

      //TODO
      for (k = 0; k < agg_buf->var_number; k++)
      {
        PIDX_variable vark = io_id->idx->variable[k];
        int bytes_per_datatype =  ((vark->bits_per_value/8) * total_chunk_size) / (64/io_id->idx->compression_bit_rate);
        int64_t prev_var_sample = (int64_t) io_id->idx->variable[io_id->init_index]->block_count_per_file[agg_buf->file_number] * io_id->idx_d->samples_per_block * bytes_per_datatype * io_id->idx->variable[k]->values_per_sample;

        data_offset = (int64_t) data_offset + prev_var_sample;
      }

      for (i = 0; i < agg_buf->sample_number; i++)
        data_offset = (int64_t) data_offset + agg_buf->buffer_size;

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_write_at(fh, data_offset, agg_buf->buffer, agg_buf->buffer_size , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }

      int write_count = 0;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != agg_buf->buffer_size)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

#else
      ssize_t write_count = pwrite(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
      if (write_count != agg_buf->buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#endif

      t3 = PIDX_get_time();

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#else
      close(fh);
#endif

      t4 = PIDX_get_time();

#ifdef PIDX_RECORD_TIME
      printf("V. [R %d] [O %lld C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, (long long)agg_buf->buffer_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
    }
  }

  int ret;
  ret = PIDX_io_per_process_write(io_id);
  if (ret != PIDX_success)
    return PIDX_err_io;

  return PIDX_success;
}


int PIDX_io_aggregated_read(PIDX_io_id io_id)
{
  int64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0;
  uint32_t *headers;
  int total_header_size = 0;
  double t1, t2, t3, t4;
  int rank = 0;

#if PIDX_HAVE_MPI
  int mpi_ret;
  MPI_File fh;
  MPI_Status status;
  MPI_Comm_rank(io_id->comm, &rank);
#else
  int fh;
#endif

  int total_chunk_size = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]);

  if (io_id->idx->enable_agg != 0)
  {
    Agg_buffer agg_buf = io_id->idx_d->agg_buffer;

    if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
    {
      t1 = PIDX_get_time();

      generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) agg_buf->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
#else
      fh = open(file_name, O_WRONLY);
#endif

      t2 = PIDX_get_time();

      data_offset = 0;
      //TODO
      headers = malloc(total_header_size);
#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_read_at(fh, 0, headers, total_header_size , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
#endif
      free(headers);

      int var_num = io_id->idx_d->agg_buffer->var_number;
      data_offset = htonl(headers[12 + ((0 + (io_id->idx->blocks_per_file * var_num))*10 )]);
      for (i = 1; i < io_id->idx->blocks_per_file; i++)
      {
        if (htonl(headers[12 + ((i + (io_id->idx->blocks_per_file * var_num))*10 )]) < data_offset)
          data_offset = htonl(headers[12 + ((i + (io_id->idx->blocks_per_file * var_num))*10 )]);
      }

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_read_at(fh, data_offset, agg_buf->buffer, agg_buf->buffer_size , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }

      int write_count = 0;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != agg_buf->buffer_size)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

#else
      ssize_t write_count = pwrite(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
      if (write_count != agg_buf->buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#endif

      t3 = PIDX_get_time();

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#else
      close(fh);
#endif

      t4 = PIDX_get_time();

#ifdef PIDX_RECORD_TIME
      printf("V. [R %d] [O %lld C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, (long long)agg_buf->buffer_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
    }
  }

  int ret;
  ret = PIDX_io_per_process_read(io_id);
  if (ret != PIDX_success)
    return PIDX_err_io;

   return PIDX_success;
}


int PIDX_io_per_process_write(PIDX_io_id io_id)
{
  int e1 = 0, i = 0, p = 0, v = 0, ret;
  int send_index = 0;
  int64_t hz_index = 0;
  int64_t index = 0, count = 0;
  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(io_id->comm, &rank);
#endif

#ifdef PIDX_DUMP_IO
  if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
  {
    char io_file_name[1024];
    ret = mkdir(io_id->idx_d->io_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, io_id->idx_d->io_dump_dir_name);
      return PIDX_err_io;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(io_id->comm);
#endif

    sprintf(io_file_name, "%s/rank_%d", io_id->idx_d->io_dump_dir_name, rank);
    io_dump_fp = fopen(io_file_name, "a+");
    if (!io_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, io_file_name);
      return PIDX_err_io;
    }
  }
#endif

  PIDX_variable var0 = io_id->idx->variable[io_id->first_index];
  for (p = 0; p < var0->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(var0->hz_buffer[p]->type == 0)
    {
      for (i = 0; i < var0->hz_buffer[p]->HZ_io_from; i++)
        hz_index = hz_index + var0->hz_buffer[p]->samples_per_level[i];

      for (i = var0->hz_buffer[p]->HZ_io_from; i < var0->hz_buffer[p]->HZ_io_to; i++)
      {
        for(e1 = 0; e1 < var0->hz_buffer[p]->samples_per_level[i] ; e1++)
        {
          if(e1 == 0)
          {
            index = var0->hz_buffer[p]->buffer_index[hz_index];
            send_index = e1;
            count = 1;

            if(var0->hz_buffer[p]->samples_per_level[i] == 1)
            {
              for(v = io_id->first_index; v <= io_id->last_index; v++)
              {
                ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], (send_index), PIDX_WRITE);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_io;
                }
              }
            }
          }
          else
          {
            if(var0->hz_buffer[p]->buffer_index[hz_index] - var0->hz_buffer[p]->buffer_index[hz_index - 1] == 1)
            {
              count++;
              if (e1 == var0->hz_buffer[p]->samples_per_level[i] - 1)
              {
                for(v = io_id->first_index; v <= io_id->last_index; v++)
                {
                  ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                }
              }
            }
            else
            {
              for(v = io_id->first_index; v <= io_id->last_index; v++)
              {
                ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_io;
                }
              }

              if(e1 == var0->hz_buffer[p]->samples_per_level[i] - 1)
              {
                for(v = io_id->first_index; v <= io_id->last_index; v++)
                {
                  ret = write_read_samples(io_id, v, var0->hz_buffer[p]->buffer_index[hz_index], 1, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, PIDX_WRITE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                }
              }
              index = var0->hz_buffer[p]->buffer_index[hz_index];
              count = 1;
              send_index = e1;
            }
          }
          hz_index++;
        }
      }
    }

    else if (var0->hz_buffer[p]->type == 1)
    {
      for (i = var0->hz_buffer[p]->HZ_io_from; i < var0->hz_buffer[p]->HZ_io_to; i++)
      {
        if (var0->hz_buffer[p]->samples_per_level[i] != 0)
        {
          for(v = io_id->first_index; v <= io_id->last_index; v++)
          {
#ifdef PIDX_DUMP_IO
            //if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
            //{
            //  fprintf(io_dump_fp, "Variable %d\n", v);
            //  fflush(io_dump_fp);
            //}
#endif
            HZ_buffer hz_buf = io_id->idx->variable[v]->hz_buffer[p];
            index = 0;
            count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1;

#ifdef PIDX_DUMP_IO
            if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
            {
              fprintf(io_dump_fp, "[%d]: ", i);
              fflush(io_dump_fp);
            }
#endif
            ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, PIDX_WRITE);
            if (ret != PIDX_success)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return PIDX_err_io;
            }
          }
        }
      }
    }

    else if (var0->hz_buffer[p]->type == 2)
    {
      for(v = io_id->first_index; v <= io_id->last_index; v++)
      {
#ifdef PIDX_DUMP_IO
        if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
        {
          fprintf(io_dump_fp, "Variable %d\n", v);
          fflush(io_dump_fp);
        }
#endif
        HZ_buffer hz_buf = io_id->idx->variable[v]->hz_buffer[p];
        for (i = hz_buf->HZ_io_from + io_id->idx_d->res_from; i < hz_buf->HZ_io_to - io_id->idx_d->res_to; i++)
        {
          if (var0->hz_buffer[p]->samples_per_level[i] != 0)
          {
            index = 0;
            count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1 - (var0->hz_buffer[p]->missing_block_count_per_level[i] * io_id->idx_d->samples_per_block);
            //printf("%d [%d] count = %d %d\n", rank, i, var0->hz_buffer[p]->missing_block_count_per_level[i], count);

#ifdef PIDX_DUMP_IO
            if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
            {
              fprintf(io_dump_fp, "[%d]: ", i);
              fflush(io_dump_fp);
            }
#endif

            ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, PIDX_WRITE);
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_io;
            }
          }
        }
      }

#if 0
      int start_block_index, end_block_index, bl;
      for (i = var0->hz_buffer[p]->HZ_io_from; i < var0->hz_buffer[p]->HZ_io_to; i++)
      {
        for (var = io_id->first_index; var <= io_id->last_index; var++)
        {
          start_block_index = io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i] / io_id->idx_d->samples_per_block;
          end_block_index = io_id->idx->variable[var]->hz_buffer[p]->end_hz_index[i] / io_id->idx_d->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          send_index = 0;
          for (bl = start_block_index; bl <= end_block_index; bl++)
          {
            if (end_block_index == start_block_index)
            {
              index = 0;
              count = (io_id->idx->variable[var]->hz_buffer[p]->end_hz_index[i] - io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i] + 1);
            }
            else
            {
              if (bl == start_block_index)
              {
                index = 0;
                count = ((start_block_index + 1) * io_id->idx_d->samples_per_block) - io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i];
              }
              else if (bl == end_block_index)
              {
                index = (end_block_index * io_id->idx_d->samples_per_block - io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i]);
                count = io_id->idx->variable[var]->hz_buffer[p]->end_hz_index[i] - ((end_block_index) * io_id->idx_d->samples_per_block) + 1;
              }
              else
              {
                index = (bl * io_id->idx_d->samples_per_block - io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i]);
                count = io_id->idx_d->samples_per_block;
              }
            }
            ret = write_read_samples(io_id, var, index + io_id->idx->variable[var]->hz_buffer[p]->start_hz_index[i], count, io_id->idx->variable[var]->hz_buffer[p]->buffer[i], send_index, PIDX_WRITE);
            if (ret != PIDX_success)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return PIDX_err_io;
            }

            send_index = send_index + count;
          }
        }
      }
#endif
    }
  }
#ifdef PIDX_DUMP_IO
  if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
  {
    fprintf(io_dump_fp, "\n");
    fclose(io_dump_fp);
  }
#endif
  return PIDX_success;
}


int PIDX_io_per_process_read(PIDX_io_id io_id)
{
    int e1 = 0, i = 0, p = 0, v = 0, ret;
    int send_index = 0;
    int64_t hz_index = 0;
    int64_t index = 0, count = 0;
    int rank = 0;

  #if PIDX_HAVE_MPI
    MPI_Comm_rank(io_id->comm, &rank);
  #endif

  #ifdef PIDX_DUMP_IO
    if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
    {
      char io_file_name[1024];
      ret = mkdir(io_id->idx_d->io_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
      if (ret != 0 && errno != EEXIST)
      {
        perror("mkdir");
        fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, io_id->idx_d->io_dump_dir_name);
        return PIDX_err_io;
      }

  #if PIDX_HAVE_MPI
      MPI_Barrier(io_id->comm);
  #endif

      sprintf(io_file_name, "%s/rank_%d", io_id->idx_d->io_dump_dir_name, rank);
      io_dump_fp = fopen(io_file_name, "a+");
      if (!io_dump_fp)
      {
        fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, io_file_name);
        return PIDX_err_io;
      }
    }
  #endif

    PIDX_variable var0 = io_id->idx->variable[io_id->first_index];
    for (p = 0; p < var0->patch_group_count; p++)
    {
      hz_index = 0, index = 0, count = 0, send_index = 0;
      if(var0->hz_buffer[p]->type == 0)
      {
        for (i = 0; i < var0->hz_buffer[p]->HZ_io_from; i++)
          hz_index = hz_index + var0->hz_buffer[p]->samples_per_level[i];

        for (i = var0->hz_buffer[p]->HZ_io_from; i < var0->hz_buffer[p]->HZ_io_to; i++)
        {
          for(e1 = 0; e1 < var0->hz_buffer[p]->samples_per_level[i] ; e1++)
          {
            if(e1 == 0)
            {
              index = var0->hz_buffer[p]->buffer_index[hz_index];
              send_index = e1;
              count = 1;

              if(var0->hz_buffer[p]->samples_per_level[i] == 1)
              {
                for(v = io_id->first_index; v <= io_id->last_index; v++)
                {
                  ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], (send_index), PIDX_READ);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                }
              }
            }
            else
            {
              if(var0->hz_buffer[p]->buffer_index[hz_index] - var0->hz_buffer[p]->buffer_index[hz_index - 1] == 1)
              {
                count++;
                if (e1 == var0->hz_buffer[p]->samples_per_level[i] - 1)
                {
                  for(v = io_id->first_index; v <= io_id->last_index; v++)
                  {
                    ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_READ);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_io;
                    }
                  }
                }
              }
              else
              {
                for(v = io_id->first_index; v <= io_id->last_index; v++)
                {
                  ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, PIDX_READ);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                }

                if(e1 == var0->hz_buffer[p]->samples_per_level[i] - 1)
                {
                  for(v = io_id->first_index; v <= io_id->last_index; v++)
                  {
                    ret = write_read_samples(io_id, v, var0->hz_buffer[p]->buffer_index[hz_index], 1, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1, PIDX_READ);
                    if (ret != PIDX_success)
                    {
                      fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_io;
                    }
                  }
                }
                index = var0->hz_buffer[p]->buffer_index[hz_index];
                count = 1;
                send_index = e1;
              }
            }
            hz_index++;
          }
        }
      }

      else if (var0->hz_buffer[p]->type == 1)
      {
        for (i = var0->hz_buffer[p]->HZ_io_from; i < var0->hz_buffer[p]->HZ_io_to; i++)
        {
          if (var0->hz_buffer[p]->samples_per_level[i] != 0)
          {
            for(v = io_id->first_index; v <= io_id->last_index; v++)
            {
  #ifdef PIDX_DUMP_IO
              //if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
              //{
              //  fprintf(io_dump_fp, "Variable %d\n", v);
              //  fflush(io_dump_fp);
              //}
  #endif
              HZ_buffer hz_buf = io_id->idx->variable[v]->hz_buffer[p];
              index = 0;
              count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1;

  #ifdef PIDX_DUMP_IO
              if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
              {
                fprintf(io_dump_fp, "[%d]: ", i);
                fflush(io_dump_fp);
              }
  #endif
              ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, PIDX_READ);
              if (ret != PIDX_success)
              {
                fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                return PIDX_err_io;
              }
            }
          }
        }
      }

      else if (var0->hz_buffer[p]->type == 2)
      {
        for(v = io_id->first_index; v <= io_id->last_index; v++)
        {
  #ifdef PIDX_DUMP_IO
          if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
          {
            fprintf(io_dump_fp, "Variable %d\n", v);
            fflush(io_dump_fp);
          }
  #endif
          HZ_buffer hz_buf = io_id->idx->variable[v]->hz_buffer[p];
          for (i = hz_buf->HZ_io_from + io_id->idx_d->res_from; i < hz_buf->HZ_io_to - io_id->idx_d->res_to; i++)
          {
            if (var0->hz_buffer[p]->samples_per_level[i] != 0)
            {
              index = 0;
              count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1 - (var0->hz_buffer[p]->missing_block_count_per_level[i] * io_id->idx_d->samples_per_block);
              //printf("%d [%d] count = %d %d\n", rank, i, var0->hz_buffer[p]->missing_block_count_per_level[i], count);

  #ifdef PIDX_DUMP_IO
              if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
              {
                fprintf(io_dump_fp, "[%d]: ", i);
                fflush(io_dump_fp);
              }
  #endif

              ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, PIDX_READ);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_io;
              }
            }
          }
        }
      }
    }
  #ifdef PIDX_DUMP_IO
    if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
    {
      fprintf(io_dump_fp, "\n");
      fclose(io_dump_fp);
    }
  #endif
    return PIDX_success;
}


int PIDX_io_finalize(PIDX_io_id io_id) 
{

  free(io_id);
  io_id = 0;

  return PIDX_success;
}


static int write_read_samples(PIDX_io_id io_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int64_t buffer_offset, int MODE)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, block_negative_offset = 0, file_number;
  int mpi_ret;
  int bytes_per_sample, bytes_per_datatype;
  int i = 0, l = 0;
  char file_name[PATH_MAX];
  off_t data_offset = 0;

#if PIDX_HAVE_MPI
  MPI_File fh;
  MPI_Status status;
#else
  int fh;
#endif
    
  samples_per_file = io_id->idx_d->samples_per_block * io_id->idx->blocks_per_file;

  bytes_per_datatype = (io_id->idx->variable[variable_index]->bits_per_value / 8) * (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]) / (64/io_id->idx->compression_bit_rate);
  
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * io_id->idx->variable[variable_index]->values_per_sample;
  
  while (hz_count) 
  {
    block_number = hz_start_index / io_id->idx_d->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;
    
    if ((int64_t)file_count > hz_count)
      file_count = hz_count;

    // build file name
    ret = generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

#if PIDX_HAVE_MPI
    if(MODE == PIDX_WRITE)
    {
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS) 
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed. (%s)\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
    }
    else
    {
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS) 
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
#else
    fh = open(file_name, O_WRONLY);
#endif
    
    data_offset = 0;
    bytes_per_sample = io_id->idx->variable[variable_index]->bits_per_value / 8;
    data_offset = file_index * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample;
    data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
    

    block_negative_offset = PIDX_blocks_find_negative_offset(io_id->idx->blocks_per_file, block_number, io_id->idx->variable[io_id->init_index]->global_block_layout);
      
    data_offset -= block_negative_offset * io_id->idx_d->samples_per_block * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample;
      
    for (l = 0; l < variable_index; l++) 
    {
      bytes_per_sample = io_id->idx->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (io_id->idx->blocks_per_file * file_number)), io_id->idx->variable[io_id->init_index]->global_block_layout))
          data_offset = data_offset + (io_id->idx->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_d->samples_per_block);
    }
    
    if(MODE == PIDX_WRITE)
    {
#if PIDX_HAVE_MPI
      
      //printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->idx_d->start_fs_block, io_id->idx_d->fs_block_size, data_offset, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), file_number, hz_start_index, samples_per_file);
      
#ifdef PIDX_DUMP_IO
      if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
      {
        fprintf(io_dump_fp, "[A] Count %lld Target Disp %d (%d %d)\n", (long long)file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), (file_index * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample - block_negative_offset * io_id->idx_d->samples_per_block * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample)/8, (int)io_id->idx_d->start_fs_block, (int)io_id->idx_d->fs_block_size);
        fflush(io_dump_fp);
      }
#endif

      mpi_ret = MPI_File_write_at(fh, data_offset, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS) 
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      int write_count;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#else
      ssize_t write_count = pwrite(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
      if (write_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#endif
    }
    if(MODE == PIDX_READ)
    {
#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_read_at(fh, data_offset, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS) 
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#else
      ssize_t read_count = pread(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
      if (read_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
#endif
    }
        
    //file_count = file_count / io_id->idx->variable[variable_index]->values_per_sample;
    
    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * io_id->idx->variable[variable_index]->values_per_sample * bytes_per_datatype;

#if PIDX_HAVE_MPI
    MPI_File_close(&fh);
#else
    close(fh);
#endif

  }
  return PIDX_success;
}
