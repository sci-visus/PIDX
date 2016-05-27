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

#undef PIDX_RECORD_TIME
//#define RANK_ORDER 1
static uint32_t *cached_header_copy;
static int enable_caching = 0;

#define PIDX_DUMP_IO 1

#ifdef PIDX_DUMP_IO
static FILE* io_dump_fp;
#endif

struct PIDX_file_io_struct
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

static int write_read_samples(PIDX_file_io_id io_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int64_t buffer_offset, PIDX_block_layout layout, int MODE);


PIDX_file_io_id PIDX_file_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_file_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_file_io_id)malloc(sizeof (*io_id));
  memset(io_id, 0, sizeof (*io_id));

  io_id->idx = idx_meta_data;
  io_id->idx_d = idx_d;

  io_id->init_index = init_index;
  io_id->first_index = first_index;
  io_id->last_index = last_index;

  return io_id;
}

#if PIDX_HAVE_MPI
int PIDX_file_io_set_communicator(PIDX_file_io_id io_id, MPI_Comm comm)
{
  if (io_id == NULL)
    return PIDX_err_id;

  io_id->comm = comm;

  return PIDX_success;
}
#endif

int PIDX_file_io_cached_data(uint32_t* cached_header)
{
  cached_header_copy = cached_header;
  enable_caching = 1;

  return PIDX_success;
}


int PIDX_aggregated_io(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, int MODE)
{
  int64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  uint32_t *headers;
  int total_header_size = 0;
#ifdef PIDX_RECORD_TIME
  double t1, t2, t3, t4, t5;
#endif

#if PIDX_HAVE_MPI
  int mpi_ret;
  MPI_File fh;
  MPI_Status status;
#else
  int fh;
#endif

  int total_chunk_size = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]);

  if (enable_caching == 1 && agg_buf->var_number == io_id->init_index && agg_buf->sample_number == 0)
  {
#ifdef PIDX_RECORD_TIME
    t1 = PIDX_get_time();
#endif

    /*
    int adjusted_file_index = 0;
    int l = pow(2, ((int)log2((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file)));
    adjusted_file_index = (l * (io_id->idx_d->idx_count[0] * io_id->idx_d->idx_count[1] * io_id->idx_d->idx_count[2]) + (((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file) - l) + (io_id->idx_d->color * l)) / io_id->idx->blocks_per_file;
    */

    generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) agg_buf->file_number /*adjusted_file_index*/, file_name, PATH_MAX);

#if !SIMULATE_IO
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
#endif

#ifdef PIDX_RECORD_TIME
    t2 = PIDX_get_time();
#endif

    data_offset = 0;
    total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
    headers = (uint32_t*)malloc(total_header_size);
    memset(headers, 0, total_header_size);

#if !SIMULATE_IO
    if (enable_caching == 1)
      memcpy (headers, cached_header_copy, total_header_size);
    else
    {
      //TODO
    }
#endif

#ifdef PIDX_RECORD_TIME
    t3 = PIDX_get_time();
#endif

    uint64_t header_size = (io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size);

#if !SIMULATE_IO
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
#endif
    free(headers);

#if !SIMULATE_IO
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
#endif

#ifdef PIDX_RECORD_TIME
    t4 = PIDX_get_time();
#endif

#if !SIMULATE_IO
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
#endif

#ifdef PIDX_RECORD_TIME
    t5 = PIDX_get_time();
#endif

#ifdef PIDX_RECORD_TIME
    printf("V0. [R %d] [O 0 C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long)agg_buf->buffer_size + header_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t3-t2), (t4-t3), (t5-t4));
#else
#endif
  }
  else if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
#ifdef PIDX_RECORD_TIME
    t1 = PIDX_get_time();
#endif

    /*
    int adjusted_file_index = 0;
    int l = pow(2, ((int)log2((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file)));
    adjusted_file_index = (l * (io_id->idx_d->idx_count[0] * io_id->idx_d->idx_count[1] * io_id->idx_d->idx_count[2]) + (((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file) - l) + (io_id->idx_d->color * l)) / io_id->idx->blocks_per_file;
    */

    generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, (unsigned int) /*adjusted_file_index*/agg_buf->file_number, file_name, PATH_MAX);

    int use_compression = 0;
    if (use_compression == 0)
    {
#if !SIMULATE_IO
#if PIDX_HAVE_MPI
      if (MODE == PIDX_WRITE)
      {
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }
      else
      {
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }
#else
      if (MODE == PIDX_WRITE)
        fh = open(file_name, O_WRONLY);
      else
        fh = open(file_name, O_RDONLY);
#endif
#endif

 #ifdef PIDX_RECORD_TIME
      t2 = PIDX_get_time();
#endif

      data_offset = 0;
      data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
      if (MODE == PIDX_WRITE)
      {
        for (k = 0; k < agg_buf->var_number; k++)
        {
          PIDX_variable vark = io_id->idx->variable[k];
          int bytes_per_datatype =  ((vark->bits_per_value/8) * total_chunk_size) / (io_id->idx->compression_factor);
          int64_t prev_var_sample = (int64_t) block_layout->block_count_per_file[agg_buf->file_number] * io_id->idx_d->samples_per_block * bytes_per_datatype * io_id->idx->variable[k]->values_per_sample;

          data_offset = (int64_t) data_offset + prev_var_sample;
        }

        for (i = 0; i < agg_buf->sample_number; i++)
          data_offset = (int64_t) data_offset + agg_buf->buffer_size;
      }
      else
      {
        int total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
        headers = malloc(total_header_size);
        memset(headers, 0, total_header_size);

  #if PIDX_HAVE_MPI
        mpi_ret = MPI_File_read_at(fh, 0, headers, total_header_size , MPI_BYTE, &status);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
  #endif
      }

#if !SIMULATE_IO
#if PIDX_HAVE_MPI
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

            mpi_ret = MPI_File_read_at(fh, data_offset, agg_buf->buffer + (block_count * io_id->idx_d->samples_per_block * (io_id->idx->variable[agg_buf->var_number]->bits_per_value/8) * io_id->idx->variable[agg_buf->var_number]->values_per_sample * io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]) / io_id->idx->compression_factor, /*agg_buf->buffer_size*/data_size , MPI_BYTE, &status);
            if (mpi_ret != MPI_SUCCESS)
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

#else
      if (MODE == PIDX_WRITE)
      {
        ssize_t write_count = pwrite(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
        if (write_count != agg_buf->buffer_size)
        {
          fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
      else
      {
        ssize_t read_count = pread(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
        if (read_count != agg_buf->buffer_size)
        {
          fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
#endif
#endif

#ifdef PIDX_RECORD_TIME
      t3 = PIDX_get_time();
#endif

#if !SIMULATE_IO
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
#endif

#ifdef PIDX_RECORD_TIME
      t4 = PIDX_get_time();
#endif

#ifdef PIDX_RECORD_TIME
      printf("V. [R %d] [O %lld C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, (long long)agg_buf->buffer_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
    }
    else
    {
      MPI_Comm new_comm;
      int rank = 0, new_rank = 0, new_count = 0;
      MPI_Comm_rank(io_id->comm, &rank);
      MPI_Comm_split(io_id->comm, agg_buf->file_number, rank, &new_comm);

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

      MPI_Allgather(&agg_buf->compressed_buffer_size, 1, MPI_LONG_LONG, aggregator_process_offset, 1, MPI_LONG_LONG, new_comm);
      MPI_Allgather(&agg_buf->var_number, 1, MPI_INT, aggregator_process_var_number, 1, MPI_INT, new_comm);
      MPI_Allgather(&agg_buf->sample_number, 1, MPI_INT, aggregator_process_sample_number, 1, MPI_INT, new_comm);

      //MPI_Gather(agg_buf->compressed_block_size, agg_buf->num_idx_blocks, MPI_UNSIGNED_LONG, void *recvbuf, agg_buf->num_idx_blocks, MPI_UNSIGNED_LONG, 0, new_comm);
      //global_block_layout->block_count_per_file[agg_buffer->file_number];

      if (MODE == PIDX_WRITE)
      {
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (mpi_ret != MPI_SUCCESS)
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
      mpi_ret = MPI_File_write_at(fh, data_offset, agg_buf->buffer, agg_buf->compressed_buffer_size , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
        return PIDX_err_io;
      }

      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
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



PIDX_return_code PIDX_async_aggregated_io(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, int MODE, MPI_Request* request, MPI_File* fh, char* filename_template)
{
  int64_t data_offset = 0;  
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  uint32_t *headers;
#ifdef PIDX_RECORD_TIME
  double t1, t2, t3, t4, t5;
#endif

#if PIDX_HAVE_MPI
  int mpi_ret;
  //MPI_File fh;
  MPI_Status status;
#else
  int fh;
#endif

  int total_chunk_size = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]);

  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
#ifdef PIDX_RECORD_TIME
    t1 = PIDX_get_time();
#endif

    /*
    int adjusted_file_index = 0;
    int l = pow(2, ((int)log2((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file)));
    adjusted_file_index = (l * (io_id->idx_d->idx_count[0] * io_id->idx_d->idx_count[1] * io_id->idx_d->idx_count[2]) + (((unsigned int) agg_buf->file_number * io_id->idx->blocks_per_file) - l) + (io_id->idx_d->color * l)) / io_id->idx->blocks_per_file;
    */

    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int) /*adjusted_file_index*/agg_buf->file_number, file_name, PATH_MAX);

    int use_compression = 0;
    if (use_compression == 0)
    {
#if !SIMULATE_IO
#if PIDX_HAVE_MPI
      if (MODE == PIDX_WRITE)
      {
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, (fh));
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }
      else
      {
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }
      }
#else
      if (MODE == PIDX_WRITE)
        fh = open(file_name, O_WRONLY);
      else
        fh = open(file_name, O_RDONLY);
#endif
#endif

 #ifdef PIDX_RECORD_TIME
      t2 = PIDX_get_time();
#endif

      data_offset = 0;
      data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
      if (MODE == PIDX_WRITE)
      {
        for (k = 0; k < agg_buf->var_number; k++)
        {
          PIDX_variable vark = io_id->idx->variable[k];
          int bytes_per_datatype =  ((vark->bits_per_value/8) * total_chunk_size) / (io_id->idx->compression_factor);
          int64_t prev_var_sample = (int64_t) block_layout->block_count_per_file[agg_buf->file_number] * io_id->idx_d->samples_per_block * bytes_per_datatype * io_id->idx->variable[k]->values_per_sample;

          data_offset = (int64_t) data_offset + prev_var_sample;
        }

        for (i = 0; i < agg_buf->sample_number; i++)
          data_offset = (int64_t) data_offset + agg_buf->buffer_size;
      }
      else
      {
        int total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
        headers = malloc(total_header_size);
        memset(headers, 0, total_header_size);

  #if PIDX_HAVE_MPI
        mpi_ret = MPI_File_read_at(*fh, 0, headers, total_header_size , MPI_BYTE, &status);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          //return PIDX_err_io;
        }
  #endif
      }

#if !SIMULATE_IO
#if PIDX_HAVE_MPI
      if (MODE == PIDX_WRITE)
      {

        //MPI_Status status;

        /*
        int rank;
        MPI_Comm_rank(io_id->comm, &rank);
        double x1, x2;
        memcpy(&x1, agg_buf->buffer, sizeof(double));
        memcpy(&x2, agg_buf->buffer + sizeof(double), sizeof(double));
        printf("W [%d] [%d %d %d] size = %d and offset = %d [%f %f]\n", rank, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, agg_buf->buffer_size, data_offset, x1, x2);
        */

        mpi_ret = MPI_File_iwrite_at(*fh, data_offset, agg_buf->buffer, agg_buf->buffer_size , MPI_BYTE, request);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          //return PIDX_err_io;
        }

        //printf("[I] %d %p\n", request, request);
        //MPI_Wait((request), &status);
        //MPI_Wait(&(io_id->idx_d->request), &status);
        //if (ret != MPI_SUCCESS)
        //{
        //  fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        //  return PIDX_err_file;
        //}

        /*
        int write_count = 0;
        MPI_Get_count(&status, MPI_BYTE, &write_count);
        if (write_count != agg_buf->buffer_size)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
        */
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

            mpi_ret = MPI_File_read_at(*fh, data_offset, agg_buf->buffer + (block_count * io_id->idx_d->samples_per_block * (io_id->idx->variable[agg_buf->var_number]->bits_per_value/8) * io_id->idx->variable[agg_buf->var_number]->values_per_sample * io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]) / io_id->idx->compression_factor, /*agg_buf->buffer_size*/data_size , MPI_BYTE, &status);
            if (mpi_ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
              //return PIDX_err_io;
            }

            int read_count = 0;
            MPI_Get_count(&status, MPI_BYTE, &read_count);
            if (read_count != /*agg_buf->buffer_size*/data_size)
            {
              fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %lldd\n", __FILE__, __LINE__, read_count, (long long)agg_buf->buffer_size);
              //return PIDX_err_io;
            }
            block_count++;
          }
        }
        free(headers);
      }

#else
      if (MODE == PIDX_WRITE)
      {
        ssize_t write_count = pwrite(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
        if (write_count != agg_buf->buffer_size)
        {
          fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
      else
      {
        ssize_t read_count = pread(fh, agg_buf->buffer, agg_buf->buffer_size, data_offset);
        if (read_count != agg_buf->buffer_size)
        {
          fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
#endif
#endif

#ifdef PIDX_RECORD_TIME
      t3 = PIDX_get_time();
#endif

#if !SIMULATE_IO
#if PIDX_HAVE_MPI
      if (MODE != PIDX_WRITE)
      {
        mpi_ret = MPI_File_close(fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
#else
      close(fh);
#endif
#endif

#ifdef PIDX_RECORD_TIME
      t4 = PIDX_get_time();
#endif

#ifdef PIDX_RECORD_TIME
      printf("V. [R %d] [O %lld C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, (long long)agg_buf->buffer_size, agg_buf->file_number, agg_buf->var_number, agg_buf->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
    }
  }

  return PIDX_success;
}



int PIDX_file_io_per_process(PIDX_file_io_id io_id, PIDX_block_layout block_layout, int MODE)
{
  int i = 0, p = 0, v = 0, ret, e1 = 0;
  int send_index = 0;
  int hz_index = 0;
  int64_t index = 0, count = 0;
  int rank = 0;

#if PIDX_HAVE_MPI
  if (io_id->idx_d->parallel_mode == 1)
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
    if (io_id->idx_d->parallel_mode == 1)
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
      for(v = io_id->first_index; v <= io_id->last_index; v++)
      {
        hz_index = 0, index = 0, count = 0, send_index = 0;
        for (i = 0; i < block_layout->resolution_from; i++)
          hz_index = hz_index + var0->hz_buffer[p]->samples_per_level[i];

        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
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
                ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, block_layout, MODE);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_io;
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
                  ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, block_layout, MODE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                }
              }
              else
              {
                ret = write_read_samples(io_id, v, index, count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, block_layout, MODE);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_io;
                }

                if(e1 == var0->hz_buffer[p]->samples_per_level[i] - 1)
                {
                  ret = write_read_samples(io_id, v, var0->hz_buffer[p]->buffer_index[hz_index], 1, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], e1,  block_layout, MODE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
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
    }

    index = 0, count = 0, send_index = 0;
    if (var0->hz_buffer[p]->type == 1)
    {
      for(v = io_id->first_index; v <= io_id->last_index; v++)
      {
        HZ_buffer hz_buf = io_id->idx->variable[v]->hz_buffer[p];
        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
#ifdef PIDX_DUMP_IO
            //if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
            //{
            //  fprintf(io_dump_fp, "Variable %d\n", v);
            //  fflush(io_dump_fp);
            //}
#endif
            
            index = 0;
            count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1;

#ifdef PIDX_DUMP_IO
            if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
            {
              fprintf(io_dump_fp, "[%d]: ", i);
              fflush(io_dump_fp);
            }
#endif
#if !SIMULATE_IO
            ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout, MODE);
#else
            ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, block_layout, MODE);
#endif
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
        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
            int start_block_index = io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] / io_id->idx_d->samples_per_block;
            int end_block_index = io_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] / io_id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              index = 0;
              count = (io_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);
              //printf("A [%d] offset 0 Count %lld\n", i, (unsigned long long)count);
#if !SIMULATE_IO
              ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout, MODE);
#else
              ret = write_read_samples(io_id, v, var0->hz_buffer[p]->start_hz_index[i], count, NULL, 0, block_layout, MODE);
#endif
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_io;
              }
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
                    index = 0;
                    count = ((start_block_index + 1) * io_id->idx_d->samples_per_block) - io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i];
                  }
                  else if (bl == end_block_index)
                  {
                    index = (end_block_index * io_id->idx_d->samples_per_block - io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = io_id->idx->variable[v]->hz_buffer[p]->end_hz_index[i] - ((end_block_index) * io_id->idx_d->samples_per_block) + 1;
                  }
                  else
                  {
                    index = (bl * io_id->idx_d->samples_per_block - io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = io_id->idx_d->samples_per_block;
                  }

                  //printf("B [%d] offset %lld send offset %lld Count %lld\n", i, (unsigned long long)index, (unsigned long long)send_index, (unsigned long long)count);
#if !SIMULATE_IO
                  ret = write_read_samples(io_id, v, index + io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, io_id->idx->variable[v]->hz_buffer[p]->buffer[i], send_index, block_layout, MODE);
#else
                  ret = write_read_samples(io_id, v, index + io_id->idx->variable[v]->hz_buffer[p]->start_hz_index[i], count, NULL, send_index, block_layout, MODE);
#endif
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                  send_index = send_index + count;
                }
                else
                  send_index = send_index + io_id->idx_d->samples_per_block;
              }
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



int PIDX_file_io_finalize(PIDX_file_io_id io_id)
{

  free(io_id);
  io_id = 0;

  return PIDX_success;
}


static int write_read_samples(PIDX_file_io_id io_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int64_t buffer_offset, PIDX_block_layout layout, int MODE)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, block_negative_offset = 0, file_number;
  int bytes_per_sample, bytes_per_datatype;
  int i = 0;
  char file_name[PATH_MAX];
  off_t data_offset = 0;

  samples_per_file = io_id->idx_d->samples_per_block * io_id->idx->blocks_per_file;

  bytes_per_datatype = (io_id->idx->variable[variable_index]->bits_per_value / 8) * (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2] * io_id->idx->chunk_size[3] * io_id->idx->chunk_size[4]) / (io_id->idx->compression_factor);
  
#if !SIMULATE_IO
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * io_id->idx->variable[variable_index]->values_per_sample;
#endif
  
  while (hz_count) 
  {
    block_number = hz_start_index / io_id->idx_d->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;
    
    if ((int64_t)file_count > hz_count)
      file_count = hz_count;

    // build file name
    int adjusted_file_index = 0;
    int l = pow(2, ((int)log2((unsigned int) file_number * io_id->idx->blocks_per_file)));
    adjusted_file_index = (l * (io_id->idx_d->idx_count[0] * io_id->idx_d->idx_count[1] * io_id->idx_d->idx_count[2]) + (((unsigned int) file_number * io_id->idx->blocks_per_file) - l) + (io_id->idx_d->color * l)) / io_id->idx->blocks_per_file;

    ret = generate_file_name(io_id->idx->blocks_per_file, io_id->idx->filename_template, /*file_number*/adjusted_file_index, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

    data_offset = 0;
    bytes_per_sample = io_id->idx->variable[variable_index]->bits_per_value / 8;
    data_offset = file_index * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample;
    data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;

    block_negative_offset = PIDX_blocks_find_negative_offset(io_id->idx->blocks_per_file, block_number, layout);
      
    data_offset -= block_negative_offset * io_id->idx_d->samples_per_block * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample;
      
    for (l = 0; l < variable_index; l++) 
    {
      bytes_per_sample = io_id->idx->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (io_id->idx->blocks_per_file * file_number)), layout))
          data_offset = data_offset + (io_id->idx->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_d->samples_per_block);
    }
    
    if(MODE == PIDX_WRITE)
    {
#if !SIMULATE_IO
#if PIDX_HAVE_MPI

#ifdef PIDX_DUMP_IO
      if (io_id->idx_d->dump_io_info == 1 && io_id->idx->current_time_step == 0)
      {
        fprintf(io_dump_fp, "[A] Count %lld Target Disp %d (%d %d)\n", (long long)file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), (file_index * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample - block_negative_offset * io_id->idx_d->samples_per_block * bytes_per_sample * io_id->idx->variable[variable_index]->values_per_sample)/8, (int)io_id->idx_d->start_fs_block, (int)io_id->idx_d->fs_block_size);
        fflush(io_dump_fp);
      }
#endif

      if (io_id->idx_d->parallel_mode == 1)
      {
        int rank = 0;
        MPI_Comm_rank(io_id->comm, &rank);
        MPI_File fh;
        MPI_Status status;
        int mpi_ret;
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() failed. (%s) [%d]\n", __FILE__, __LINE__, file_name, file_number);
          return PIDX_err_io;
        }

        /*
        printf("[%d] Data Offset %d Count %d\n", rank, data_offset, (file_count));
        int x = 0;
        for (x = 0; x < file_count; x++)
        {
          double x1;
          memcpy(&x1, hz_buffer + x * sizeof(double), sizeof(double));
          printf("Values %d %f\n", x, x1);
        }
        */

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
        MPI_File_close(&fh);
      }
      else
      {
        int fh;
        fh = open(file_name, O_WRONLY);
        ssize_t write_count = pwrite(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
        if (write_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
        {
          fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
        close(fh);
      }
#else
      int fh;
      fh = open(file_name, O_WRONLY);
      /*
      double x1, x2, x3, x4;
      memcpy(&x1, hz_buffer, sizeof(double));
      memcpy(&x2, hz_buffer + sizeof(double), sizeof(double));
      memcpy(&x3, hz_buffer + 2*sizeof(double), sizeof(double));
      memcpy(&x4, hz_buffer + 3*sizeof(double), sizeof(double));
      printf("[%d] [%d %d] Values %f %f %f %f\n", variable_index, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset, x1, x2, x3, x4);
      */
      ssize_t write_count = pwrite(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
      if (write_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      close(fh);
#endif
#endif
    }
    if(MODE == PIDX_READ)
    {
#if PIDX_HAVE_MPI
      if (io_id->idx_d->parallel_mode == 1)
      {
        MPI_File fh;
        MPI_Status status;
        int mpi_ret;
        mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }

        mpi_ret = MPI_File_read_at(fh, data_offset, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }

        MPI_File_close(&fh);
      }
      else
      {
        int fh;
        fh = open(file_name, O_RDONLY);
        ssize_t read_count = pread(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
        if (read_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
        {
          fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
        close(fh);
      }
#else
      int fh;
      fh = open(file_name, O_RDONLY);
      ssize_t read_count = pread(fh, hz_buffer, file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8), data_offset);
      if (read_count != file_count * io_id->idx->variable[variable_index]->values_per_sample * (io_id->idx->variable[variable_index]->bits_per_value/8))
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      close(fh);
#endif
    }

    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * io_id->idx->variable[variable_index]->values_per_sample * bytes_per_datatype;

#if PIDX_HAVE_MPI

#else

#endif

  }
  return PIDX_success;
}
