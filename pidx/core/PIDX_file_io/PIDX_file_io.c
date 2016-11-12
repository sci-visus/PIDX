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

#include "../../PIDX_inc.h"
#include "async_io.h"
#include "compressed_io.h"
#include "PIDX_file_io.h"
#include "sync_io.h"


static uint32_t *cached_header_copy;
static int enable_caching = 0;

//#define PIDX_DUMP_IO 1

#ifdef PIDX_DUMP_IO
static FILE* io_dump_fp;
#endif



PIDX_file_io_id PIDX_file_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int first_index, int last_index)
{
  PIDX_file_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_file_io_id)malloc(sizeof (*io_id));
  memset(io_id, 0, sizeof (*io_id));

  io_id->idx = idx_meta_data;
  io_id->idx_d = idx_d;
  io_id->idx_c = idx_c;

  io_id->group_index = 0;
  io_id->first_index = first_index;
  io_id->last_index = last_index;

  return io_id;
}


int PIDX_file_io_cached_data(uint32_t* cached_header)
{
  cached_header_copy = cached_header;
  enable_caching = 1;

  return PIDX_success;
}


int PIDX_aggregated_io(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, int MODE)
{
  unsigned long long data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  uint32_t *headers;
  int total_header_size = 0;
  int ret;
  MPI_File fh;
  MPI_Status status;
  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);

  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
    /*
    int adjusted_file_index = 0;
    int l = pow(2, ((int)log2((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file)));
    adjusted_file_index = (l * (io_id->idx_d->partition_count[0] * io_id->idx_d->partition_count[1] * io_id->idx_d->partition_count[2]) + (((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file) - l) + (io_id->idx_d->color * l)) / io_id->idx->blocks_per_file;
    */

    generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) /*adjusted_file_index*/agg_buf->file_number, file_name, PATH_MAX);

    int use_compression = 0;
    if (use_compression == 0)
    {
      if (MODE == PIDX_WRITE)
      {
        ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }
      else
      {
        ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }

      PIDX_variable_group var_grp = io_id->idx->variable_grp[io_id->group_index];
      data_offset = 0;
      data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
      if (MODE == PIDX_WRITE)
      {
        for (k = 0; k < agg_buf->var_number; k++)
        {
          PIDX_variable vark = var_grp->variable[k];
          int bytes_per_datatype =  ((vark->bpv/8) * tck) / (io_id->idx->compression_factor);
          unsigned long long prev_var_sample = (unsigned long long) block_layout->bcpf[agg_buf->file_number] * io_id->idx_d->samples_per_block * bytes_per_datatype * var_grp->variable[k]->vps;

          data_offset = (unsigned long long) data_offset + prev_var_sample;
        }

        for (i = 0; i < agg_buf->sample_number; i++)
          data_offset = (unsigned long long) data_offset + agg_buf->buffer_size;
      }
      else
      {
        int total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
        headers = malloc(total_header_size);
        memset(headers, 0, total_header_size);

        ret = MPI_File_read_at(fh, 0, headers, total_header_size , MPI_BYTE, &status);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }

      if (MODE == PIDX_WRITE)
      {
        /*
        int rank;
        MPI_Comm_rank(io_id->comm, &rank);
        double x1, x2;
        memcpy(&x1, agg_buf->buffer, sizeof(double));
        memcpy(&x2, agg_buf->buffer + sizeof(double), sizeof(double));
        printf("W [%d] [%d %d %d] size = %d and offset = %d [%f %f]\n", rank, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, agg_buf->buffer_size, data_offset, x1, x2);
        */
        ret = MPI_File_write_at(fh, data_offset, agg_buf->buffer, agg_buf->buffer_size , MPI_BYTE, &status);
        if (ret != MPI_SUCCESS)
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
      }
      else
      {
        int data_size = 0;
        int block_count = 0;
        for (i = 0; i < io_id->idx->blocks_per_file; i++)
        {
          if (PIDX_blocks_is_block_present(agg_buf->file_number * io_id->idx->blocks_per_file + i, block_layout))
          {
            data_offset = htonl(headers[12 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);
            data_size = htonl(headers[14 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);

            ret = MPI_File_read_at(fh, data_offset, agg_buf->buffer + (block_count * io_id->idx_d->samples_per_block * (var_grp->variable[agg_buf->var_number]->bpv/8) * var_grp->variable[agg_buf->var_number]->vps * io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]) / io_id->idx->compression_factor, /*agg_buf->buffer_size*/data_size , MPI_BYTE, &status);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
              return PIDX_err_io;
            }

            int read_count = 0;
            MPI_Get_count(&status, MPI_BYTE, &read_count);
            if (read_count != /*agg_buf->buffer_size*/data_size)
            {
              fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %lldd\n", __FILE__, __LINE__, read_count, (long long)agg_buf->buffer_size);
              return PIDX_err_io;
            }
            block_count++;
          }
        }
        free(headers);
      }

      ret = MPI_File_close(&fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
       return PIDX_err_io;
      }
    }
    else
    {
      MPI_Comm new_comm;
      int new_rank = 0, new_count = 0;
      MPI_Comm_split(io_id->idx_c->comm, agg_buf->file_number, io_id->idx_c->rank, &new_comm);

      MPI_Comm_rank(new_comm, &new_rank);
      MPI_Comm_size(new_comm, &new_count);

      unsigned long long *aggregator_process_offset;
      int *aggregator_process_var_number;
      int *aggregator_process_sample_number;

      aggregator_process_offset = malloc(sizeof(*aggregator_process_offset) * new_count);
      memset(aggregator_process_offset, 0, sizeof(*aggregator_process_offset) * new_count);
      aggregator_process_var_number = malloc(sizeof(*aggregator_process_var_number) * new_count);
      memset(aggregator_process_var_number, 0, sizeof(*aggregator_process_var_number) * new_count);
      aggregator_process_sample_number = malloc(sizeof(*aggregator_process_sample_number) * new_count);
      memset(aggregator_process_sample_number, 0, sizeof(*aggregator_process_sample_number) * new_count);

      MPI_Allgather(&agg_buf->compressed_buffer_size, 1, MPI_UNSIGNED_LONG_LONG, aggregator_process_offset, 1, MPI_UNSIGNED_LONG_LONG, new_comm);
      MPI_Allgather(&agg_buf->var_number, 1, MPI_INT, aggregator_process_var_number, 1, MPI_INT, new_comm);
      MPI_Allgather(&agg_buf->sample_number, 1, MPI_INT, aggregator_process_sample_number, 1, MPI_INT, new_comm);

      //MPI_Gather(agg_buf->compressed_block_size, agg_buf->num_idx_blocks, MPI_UNSIGNED_LONG, void *recvbuf, agg_buf->num_idx_blocks, MPI_UNSIGNED_LONG, 0, new_comm);
      //global_block_layout->bcpf[agg_buffer->file_number];

      if (MODE == PIDX_WRITE)
      {
        ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }

      total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
      data_offset = total_header_size;
      for (i = 0; i < new_rank; i++)
        data_offset = data_offset + aggregator_process_offset[i];

      printf("data offset = %lld and length = %lld\n", (long long)data_offset, (long long)agg_buf->compressed_buffer_size);
      ret = MPI_File_write_at(fh, data_offset, agg_buf->buffer, agg_buf->compressed_buffer_size , MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
        return PIDX_err_io;
      }

      ret = MPI_File_close(&fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      free(aggregator_process_offset);
      free(aggregator_process_var_number);
      free(aggregator_process_sample_number);
      MPI_Comm_free(&new_comm);
    }
  }

  return PIDX_success;
}





PIDX_return_code PIDX_async_aggregated_io(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template, int mode)
{
  int ret;
  if (mode == PIDX_WRITE)
    ret = PIDX_async_aggregated_write(io_id, agg_buf, block_layout, request, fh, filename_template);
  else if (mode == PIDX_READ)
    ret = PIDX_async_aggregated_read(io_id, agg_buf, block_layout, request, fh, filename_template);

  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error in file: [%s] line: [%d]\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  return PIDX_success;
}



int PIDX_file_io_finalize(PIDX_file_io_id io_id)
{

  free(io_id);
  io_id = 0;

  return PIDX_success;
}
