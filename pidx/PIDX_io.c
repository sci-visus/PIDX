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

enum IO_MODE { PIDX_READ, PIDX_WRITE};

struct PIDX_io_struct
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Win win;
#endif

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;

  int start_var_index;
  int end_var_index;
};

static int write_read_samples(PIDX_io_id io_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int64_t buffer_offset, int MODE);

//static int generate_file_name(PIDX_io_id io_id, int file_number, char* filename, int maxlen);

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

  samples_per_file = io_id->idx_derived_ptr->samples_per_block * io_id->idx_ptr->blocks_per_file;

  bytes_per_datatype = (io_id->idx_ptr->variable[variable_index]->bits_per_value / 8) * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4]);

  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * io_id->idx_ptr->variable[variable_index]->values_per_sample;

  while (hz_count)
  {
    block_number = hz_start_index / io_id->idx_derived_ptr->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;

    if ((int64_t)file_count > hz_count)
      file_count = hz_count;

    // build file name
    ret = generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return -1;
    }
    //printf("File name = %s\n", file_name);
#if PIDX_HAVE_MPI
    if(MODE == PIDX_WRITE)
    {
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed. (%s)\n", __FILE__, __LINE__, file_name);
	return -1;
      }
    }
    else
    {
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
	return -1;
      }
    }
#else
    fh = open(file_name, O_WRONLY);
#endif

    data_offset = 0;
    bytes_per_sample = io_id->idx_ptr->variable[variable_index]->bits_per_value / 8;
    data_offset = file_index * bytes_per_sample * io_id->idx_ptr->variable[variable_index]->values_per_sample;
    data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;

#ifdef PIDX_VAR_SLOW_LOOP
    block_negative_offset = PIDX_blocks_find_negative_offset(io_id->idx_ptr->blocks_per_file, block_number, io_id->idx_ptr->variable[variable_index]->VAR_global_block_layout);
#else
    block_negative_offset = PIDX_blocks_find_negative_offset(io_id->idx_ptr->blocks_per_file, block_number, io_id->idx_derived_ptr->global_block_layout);
#endif

    data_offset -= block_negative_offset * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[variable_index]->values_per_sample;

#ifdef PIDX_VAR_SLOW_LOOP
    for (l = 0; l < variable_index; l++)
    {
      bytes_per_sample = io_id->idx_ptr->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_ptr->variable[l]->VAR_global_block_layout))
          data_offset = data_offset + (io_id->idx_ptr->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_derived_ptr->samples_per_block);
    }
#else
    for (l = 0; l < variable_index; l++)
    {
      bytes_per_sample = io_id->idx_ptr->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_derived_ptr->global_block_layout))
          data_offset = data_offset + (io_id->idx_ptr->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_derived_ptr->samples_per_block);
    }
#endif

#if 0
    int u = 0;
    printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->idx_derived_ptr->start_fs_block, io_id->idx_derived_ptr->fs_block_size, data_offset, file_count, file_number, hz_start_index, samples_per_file);
    if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(double) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "float64") == 0)
    {
      double value;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&value, hz_buffer + (u * sizeof(double)), sizeof(double));
	printf("[DOUBLE] Value at %d = %f\n", u,  value);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(int) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "float32") == 0)
    {
      float fvalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&fvalue, hz_buffer + (u * sizeof(float)), sizeof(float));
	printf("[FLOAT] Value at %d = %f\n", u,  fvalue);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(float) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "int32") == 0)
    {
      int ivalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&ivalue, hz_buffer + (u * sizeof(int)), sizeof(int));
	printf("[INT] Value at %d = %f\n", u,  ivalue);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(unsigned char) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "uint8") == 0)
    {
      unsigned char uivalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&uivalue, hz_buffer + (u * sizeof(unsigned char)), sizeof(unsigned char));
	printf("[U INT] Value at %d = %f\n", u,  uivalue);
      }
    }
#endif

    if(MODE == PIDX_WRITE)
    {
#if PIDX_HAVE_MPI

      //printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->idx_derived_ptr->start_fs_block, io_id->idx_derived_ptr->fs_block_size, data_offset, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), file_number, hz_start_index, samples_per_file);

      mpi_ret = MPI_File_write_at(fh, data_offset, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }
#else
      pwrite(fh, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), data_offset);
#endif
    }
    if(MODE == PIDX_READ)
    {
#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_read_at(fh, data_offset, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
	return -1;
      }
#else
      pread(fh, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), data_offset);
#endif
    }

    //file_count = file_count / io_id->idx_ptr->variable[variable_index]->values_per_sample;

    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * bytes_per_datatype;

#if PIDX_HAVE_MPI
    MPI_File_close(&fh);
#else
    close(fh);
#endif

  }
  return (0);
}

PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  PIDX_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_io_id)malloc(sizeof (*io_id));
  memset(io_id, 0, sizeof (*io_id));

  io_id->idx_ptr = idx_meta_data;
  io_id->idx_derived_ptr = idx_derived_ptr;

  /*
  io_id->idx_ptr = (idx_dataset)malloc(sizeof(*(io_id->idx_ptr)));
  memcpy(io_id->idx_ptr, idx_meta_data, sizeof(*(io_id->idx_ptr)));

  io_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(io_id->idx_derived_ptr)));
  memcpy(io_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(io_id->idx_derived_ptr)));
  */

  io_id->start_var_index = start_var_index;
  io_id->end_var_index = end_var_index;

  return io_id;
}

#if PIDX_HAVE_MPI
int PIDX_io_set_communicator(PIDX_io_id io_id, MPI_Comm comm)
{
  io_id->comm = comm;
  //MPI_Comm_dup(comm, &io_id->comm);
  return 0;
}
#endif

void print_buffer(unsigned char *buffer, int offset, int count, int sample_size, int datatype_size, char* name)
{
  int i;
  if(datatype_size/8 == sizeof(int) && strcmp(name, "int32") == 0)
  {
    int ival = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&ival, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("Int value at %d = %d\n", i, ival);
    }
  }
  else if(datatype_size/8 == sizeof(double) && strcmp(name, "float64") == 0)
  {
    double dval = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&dval, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Double value at %d = %f\n", count, i, dval);
    }
  }
  else if (datatype_size/8 == sizeof(float) && strcmp(name, "float32") == 0)
  {
    float fval = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&fval, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Float value at %d = %f\n", count, i, fval);
    }
  }
  else if (datatype_size/8 == sizeof(unsigned char) && strcmp(name, "uint8") == 0)
  {
    unsigned int uival = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&uival, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Unsigned int value at %d = %d\n", count, i, uival);
    }
  }
}

int PIDX_io_cached_data(uint32_t* cached_header)
{
  cached_header_copy = cached_header;
  enable_caching = 1;

  return 0;
}

#if 1
int compress_aggregation_buffer(PIDX_io_id io_id, unsigned char* agg_buffer, int buffer_size)
{
  int i = 0, j = 0, double_elements_per_block = 0, memmove_count = 0;
  double value = 0;
  int nan_counter = 0, total_nan_count = 0, block_nan_count = 0, existing_block_count = 0;
  int bytes_per_datatype =  (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) /* (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4]) */;

  for (i = 0; i < io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]; i++)
  {
    if (PIDX_blocks_is_block_present((i + (io_id->idx_ptr->blocks_per_file * io_id->idx_derived_ptr->agg_buffer->file_number)), io_id->idx_derived_ptr->global_block_layout))
    {
      block_nan_count = 0;

      double_elements_per_block = io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4]);
      printf("samples per block = %d double elements per block = %d block count %d\n", io_id->idx_derived_ptr->samples_per_block, double_elements_per_block, existing_block_count);
      for (j = 0; j < double_elements_per_block; j++)
      {
        memcpy(&value, io_id->idx_derived_ptr->agg_buffer->buffer + ((existing_block_count * double_elements_per_block) + (j /*- memmove_count*/)) * bytes_per_datatype, bytes_per_datatype );

        /// if value is nan
        //if (isnan(value) != 0)
        if (value == 9999999999)
        {
          nan_counter++;
          total_nan_count++;
          block_nan_count++;
          //printf("YYYYYYYYY %d %d\n", ((existing_block_count * double_elements_per_block) + (j /*- memmove_count*/)), total_nan_count);
        }
        /// when performing memmove
        else
        {
          //printf("XXXXXXXXX %d %d\n", j, total_nan_count);
          if( nan_counter != 0)
          {
            printf("[%d] XXXXXXXXX %d %d\n", i, j, nan_counter);
            //printf("Destination %d (%d - %d) Source %d (%d - %d + %d) Count %d (%d - %d)\n", (int)(((existing_block_count * double_elements_per_block) + j) - total_nan_count), (int)((existing_block_count * double_elements_per_block) + j), (int)total_nan_count, (int)(((existing_block_count * double_elements_per_block) + j) - total_nan_count + nan_counter), (int)((existing_block_count * double_elements_per_block) + j), (int)total_nan_count, (int)nan_counter, (int)io_id->idx_derived_ptr->agg_buffer->buffer_size - ((existing_block_count * double_elements_per_block) + j) * bytes_per_datatype, (int)io_id->idx_derived_ptr->agg_buffer->buffer_size, (int)((existing_block_count * double_elements_per_block) + j) * bytes_per_datatype);

            //memmove(io_id->idx_derived_ptr->agg_buffer->buffer + (((existing_block_count * double_elements_per_block) + j) - total_nan_count) * bytes_per_datatype, io_id->idx_derived_ptr->agg_buffer->buffer + (((existing_block_count * double_elements_per_block) + j) - total_nan_count + nan_counter) * bytes_per_datatype, io_id->idx_derived_ptr->agg_buffer->buffer_size - ((existing_block_count * double_elements_per_block) + j) * bytes_per_datatype);
            memmove_count = memmove_count + nan_counter;
            nan_counter = 0;
          }
        }
        value = 0;
      }
      io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i] = (io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])) - (block_nan_count * (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));
      printf("Compressed block size %d = %lld\n", i, (long long)io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i]/(io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));
      existing_block_count++;
    }
    else
      io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i] = 0;
  }

  unsigned char* temp_buffer = (unsigned char*)realloc(io_id->idx_derived_ptr->agg_buffer->buffer, io_id->idx_derived_ptr->agg_buffer->buffer_size - (total_nan_count * bytes_per_datatype));
  if (temp_buffer == NULL)
  {
    fprintf(stderr, "[%s] [%d] realloc() of aggregation buffer for compression failed.\n", __FILE__, __LINE__);
    return -1;
  }
  else
    io_id->idx_derived_ptr->agg_buffer->buffer = temp_buffer;

  io_id->idx_derived_ptr->agg_buffer->compressed_buffer_size = ((existing_block_count * io_id->idx_derived_ptr->samples_per_block) * (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])) - (total_nan_count * (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));

#if 1
  printf("Total NAN count = %d %lld %d %d %d %d %d\n",
         (int)total_nan_count,
         (long long)io_id->idx_derived_ptr->agg_buffer->compressed_buffer_size,
         (int)existing_block_count,
         (int)io_id->idx_derived_ptr->samples_per_block,
         (int)total_nan_count,
         (int)(io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8),
         (int)(io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4]));
#endif

  return 0;
}
#endif

int PIDX_io_aggregated_write(PIDX_io_id io_id)
{
  int64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0, rank, mpi_ret, new_rank, new_count;
  uint32_t *headers;
  uint64_t *aggregator_process_offset;
  int *aggregator_process_var_number;
  int *aggregator_process_sample_number;
  int total_header_size = 0;
  int write_count;
  int bytes_per_datatype;
  int aggregated_header_write = 0;

#ifdef PIDX_RECORD_TIME
  double t1, t2, t3, t4, t5;
#endif

#if PIDX_HAVE_MPI
  MPI_Comm new_comm;
  MPI_File fh;
  MPI_Status status;
  MPI_Comm_rank(io_id->comm, &rank);
#else
  int fh;
#endif

  if (io_id->idx_ptr->enable_compression == 1 && io_id->idx_ptr->compression_type == 0)
  {
    uint32_t *header;
    if (io_id->idx_derived_ptr->agg_buffer->file_number == -1)
      io_id->idx_derived_ptr->agg_buffer->file_number = io_id->idx_derived_ptr->existing_file_count + 1;

    int var_used_in_binary_file = 0;
    var_used_in_binary_file = (io_id->idx_ptr->variable_count < 0) ? 64 : io_id->idx_ptr->variable_count;
    total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * var_used_in_binary_file;

    MPI_Comm_split(io_id->comm, io_id->idx_derived_ptr->agg_buffer->file_number, rank, &new_comm);

    if (io_id->idx_derived_ptr->agg_buffer->var_number != -1 && io_id->idx_derived_ptr->agg_buffer->sample_number != -1 && io_id->idx_derived_ptr->agg_buffer->file_number != -1)
    {
      compress_aggregation_buffer (io_id, io_id->idx_derived_ptr->agg_buffer->buffer, io_id->idx_derived_ptr->agg_buffer->buffer_size);

      MPI_Comm_rank(new_comm, &new_rank);
      MPI_Comm_size(new_comm, &new_count);
      //printf("Rank %d and new rank %d\n", rank, new_rank);

      aggregator_process_offset = malloc(sizeof(*aggregator_process_offset) * new_count);
      memset(aggregator_process_offset, 0, sizeof(*aggregator_process_offset) * new_count);
      aggregator_process_var_number = malloc(sizeof(*aggregator_process_var_number) * new_count);
      memset(aggregator_process_var_number, 0, sizeof(*aggregator_process_var_number) * new_count);
      aggregator_process_sample_number = malloc(sizeof(*aggregator_process_sample_number) * new_count);
      memset(aggregator_process_sample_number, 0, sizeof(*aggregator_process_sample_number) * new_count);

      MPI_Allgather(&io_id->idx_derived_ptr->agg_buffer->compressed_buffer_size, 1, MPI_LONG_LONG, aggregator_process_offset, 1, MPI_LONG_LONG, new_comm);
      MPI_Allgather(&io_id->idx_derived_ptr->agg_buffer->var_number, 1, MPI_INT, aggregator_process_var_number, 1, MPI_INT, new_comm);
      MPI_Allgather(&io_id->idx_derived_ptr->agg_buffer->sample_number, 1, MPI_INT, aggregator_process_sample_number, 1, MPI_INT, new_comm);

      //TODO calculate data_offset
      generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);

      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return -1;
      }

      data_offset = total_header_size;
      for (i = 0; i < new_rank; i++)
        data_offset = data_offset + aggregator_process_offset[i];

      printf("data offset = %lld and length = %lld\n", (long long)data_offset, (long long)io_id->idx_derived_ptr->agg_buffer->compressed_buffer_size);
      mpi_ret = MPI_File_write_at(fh, data_offset, io_id->idx_derived_ptr->agg_buffer->buffer, io_id->idx_derived_ptr->agg_buffer->compressed_buffer_size , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
        return -1;
      }


      if (aggregated_header_write == 1)
      {
#if 0
        uint32_t *aggregated_header = 0;
        int var_used_in_binary_file = 0, total_header_size = 0;
        if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
        {
          var_used_in_binary_file = (io_id->idx_ptr->variable_count < 0) ? 64 : io_id->idx_ptr->variable_count;
          total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * var_used_in_binary_file;
          io_id->idx_derived_ptr->start_fs_block = total_header_size / io_id->idx_derived_ptr->fs_block_size;
          if (total_header_size % io_id->idx_derived_ptr->fs_block_size)
            io_id->idx_derived_ptr->start_fs_block++;

          aggregated_header = (uint32_t*)malloc(total_header_size);
          memset(aggregated_header, 0, total_header_size);

#if PIDX_HAVE_MPI
          mpi_ret = MPI_Win_create(aggregated_header, total_header_size, sizeof (uint32_t), MPI_INFO_NULL, new_comm, &(io_id->win));
          if (mpi_ret != MPI_SUCCESS)
          {
            fprintf(stderr, "[%s] [%d] MPI_Win_create() failed.\n", __FILE__, __LINE__);
            return -1;
          }
#endif
        }
        else
        {
#if PIDX_HAVE_MPI
          mpi_ret =  MPI_Win_create(0, 0, 1, MPI_INFO_NULL, new_comm, &(io_id->win));
          if (mpi_ret != MPI_SUCCESS)
          {
            fprintf(stderr, "[%s] [%d] MPI_Win_create() failed.\n", __FILE__, __LINE__);
            return -1;
          }
        }
#endif

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
        mpi_ret = MPI_Win_fence(0, io_id->win);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_Win_fence() failed.\n", __FILE__, __LINE__);
          return -1;
        }
#else
        //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif

#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, io_id->idx_derived_ptr->agg_buffer->rank_holder[io_id->start_var_index][0][io_id->idx_derived_ptr->agg_buffer->file_number], 0 , io_id->win);
#endif

#if RANK_ORDER
        mpi_ret = MPI_Put(const void *origin_addr, (10 * io_id->idx_ptr->blocks_per_file), MPI_INT, io_id->idx_derived_ptr->agg_buffer->rank_holder[io_id->start_var_index][0][io_id->idx_derived_ptr->agg_buffer->file_number], MPI_Aint target_disp, (10 * io_id->idx_ptr->blocks_per_file), MPI_INT, io_id->win);
#else
        mpi_ret = MPI_Put(const void *origin_addr, (10 * io_id->idx_ptr->blocks_per_file), MPI_INT, io_id->idx_derived_ptr->agg_buffer->rank_holder[io_id->idx_derived_ptr->agg_buffer->file_number][0][io_id->start_var_index], MPI_Aint target_disp, (10 * io_id->idx_ptr->blocks_per_file), MPI_INT, io_id->win);
#endif
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_Put() failed.\n", __FILE__, __LINE__);
          return -1;
        }

#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(io_id->idx_derived_ptr->agg_buffer->rank_holder[io_id->start_var_index][0][io_id->idx_derived_ptr->agg_buffer->file_number], io_id->win);
#endif

#ifdef PIDX_ACTIVE_TARGET
        mpi_ret = MPI_Win_fence(0, io_id->win);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_Win_fence() failed.\n", __FILE__, __LINE__);
          return -1;
        }
#else
        //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif

        mpi_ret = MPI_File_write_at(fh, 0, aggregated_header, total_header_size , MPI_BYTE, &status);
        if (mpi_ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
          return -1;
        }

#endif

        if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
        {
          free(aggregated_header);
          aggregated_header = 0;
        }
#endif
      }
      else
      {
        data_offset = total_header_size;
        for (i = 0; i < new_rank; i++)
          data_offset = data_offset + aggregator_process_offset[i];
        printf("init offset = %ld\n", data_offset);
        if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
        {
          int header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t);
          header = (uint32_t*)malloc(header_size);
          memset(header, 0, header_size);

          for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
          {
            data_offset = data_offset + io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i];
            //printf("Block %d : Offset %ld Count %lld block size %lld\n", i, data_offset, (long long)io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i], io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i]);

            header[12 + (i*10)] = htonl(data_offset);
            header[14 + (i*10)] = htonl(io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i]);
          }

          mpi_ret = MPI_File_write_at(fh, 0, header, (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) , MPI_BYTE, &status);
          if (mpi_ret != MPI_SUCCESS)
          {
            fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
            return -1;
          }
        }
        else
        {
          int header_size = (10 * io_id->idx_ptr->blocks_per_file) * sizeof (uint32_t);
          header = (uint32_t*)malloc(header_size);
          memset(header, 0, header_size);

          for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
          {
            data_offset = data_offset + io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i];
            header[12 + (i*10)] = htonl(data_offset);
            header[14 + (i*10)] = htonl(io_id->idx_derived_ptr->agg_buffer->compressed_block_size[i]);
          }

          mpi_ret = MPI_File_write_at(fh, (10 + (10 * io_id->idx_ptr->blocks_per_file * new_rank)) * sizeof(uint32_t), header, (10 * io_id->idx_ptr->blocks_per_file) * sizeof (uint32_t) , MPI_BYTE, &status);
          if (mpi_ret != MPI_SUCCESS)
          {
            fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
            return -1;
          }
        }
        free(header);
      }

      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }

      free(aggregator_process_offset);
      free(aggregator_process_var_number);
      free(aggregator_process_sample_number);
    }
    MPI_Comm_free(&new_comm);
  }
  else
  {
    if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
    {
      bytes_per_datatype =  ((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8)  * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])) / (64/io_id->idx_ptr->compression_bit_rate);

#ifdef PIDX_RECORD_TIME
      t1 = MPI_Wtime();
#endif
      generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed filename %s.\n", __FILE__, __LINE__, file_name);
        return -1;
      }
#else
      fh = open(file_name, O_WRONLY);
#endif

#ifdef PIDX_RECORD_TIME
      t2 = MPI_Wtime();
#endif

      data_offset = 0;
      total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * io_id->idx_ptr->variable_count;

      headers = (uint32_t*)malloc(total_header_size);
      memset(headers, 0, total_header_size);

      if (enable_caching == 1)
        memcpy (headers, cached_header_copy, total_header_size);

#ifdef PIDX_RECORD_TIME
      t3 = MPI_Wtime();
#endif

#ifdef PIDX_VAR_SLOW_LOOP
      unsigned char* temp_buffer = realloc(io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size));

      if (temp_buffer == NULL)
      {
        ;
      }
      else
      {
        io_id->idx_derived_ptr->agg_buffer->buffer = temp_buffer;
        memmove(io_id->idx_derived_ptr->agg_buffer->buffer + io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) );
        memcpy(io_id->idx_derived_ptr->agg_buffer->buffer, headers, total_header_size);
        memset(io_id->idx_derived_ptr->agg_buffer->buffer + total_header_size, 0, (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size - total_header_size));
      }
      free(headers);


#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_write_at(fh, 0, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }

      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != (((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block  / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return -1;
      }

#else
      pwrite(fh, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), 0);
#endif

#else
      unsigned char* temp_buffer = (unsigned char*)realloc(io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size));

      if (temp_buffer == NULL)
      {
        fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
        return -1;
      }
      else
      {
        io_id->idx_derived_ptr->agg_buffer->buffer = temp_buffer;
        memmove(io_id->idx_derived_ptr->agg_buffer->buffer + io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) );
        memcpy(io_id->idx_derived_ptr->agg_buffer->buffer, headers, total_header_size);

        memset(io_id->idx_derived_ptr->agg_buffer->buffer + total_header_size, 0, (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size - total_header_size));
      }
      free(headers);

#if PIDX_HAVE_MPI
      /*
      double val[24];
      int t = 0;
      memcpy(val, io_id->idx_derived_ptr->agg_buffer->buffer + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), 24 * sizeof(double));
      for (t = 0; t < 24; t++)
        printf("Agg value at %d = %f\n", t, val[t]);
      */
      mpi_ret = MPI_File_write_at(fh, 0, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
        return -1;
      }

      double x, y, z;
      memcpy(&x, io_id->idx_derived_ptr->agg_buffer->buffer, sizeof(double));
      memcpy(&y, io_id->idx_derived_ptr->agg_buffer->buffer + sizeof(double), sizeof(double));
      memcpy(&z, io_id->idx_derived_ptr->agg_buffer->buffer + sizeof(double) + sizeof(double), sizeof(double));

      MPI_Get_count(&status, MPI_BYTE, &write_count);
      //printf("[A] Elemets to write %d\n", write_count);
      if (write_count != (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block  / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return -1;
      }

#else
      pwrite(fh, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), 0);
#endif

#endif

#ifdef PIDX_RECORD_TIME
      t4 = MPI_Wtime();
#endif

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }
#else
        close(fh);
#endif

#ifdef PIDX_RECORD_TIME
      t5 = MPI_Wtime();
#endif

#ifdef PIDX_RECORD_TIME
      printf("A. [R %d] [O 0 C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank,
             (long long)(((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), io_id->idx_derived_ptr->agg_buffer->file_number, io_id->idx_derived_ptr->agg_buffer->var_number, io_id->idx_derived_ptr->agg_buffer->sample_number, (t2-t1), (t3-t2), (t4-t3), (t5-t4));
#endif
    }
    else if (io_id->idx_derived_ptr->agg_buffer->var_number != -1 && io_id->idx_derived_ptr->agg_buffer->sample_number != -1 && io_id->idx_derived_ptr->agg_buffer->file_number != -1)
    {
#ifdef PIDX_RECORD_TIME
      t1 = MPI_Wtime();
#endif
      bytes_per_datatype =  (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])  / (64/io_id->idx_ptr->compression_bit_rate);

      generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return -1;
      }
#else
      fh = open(file_name, O_WRONLY);
#endif

#ifdef PIDX_RECORD_TIME
      t2 = MPI_Wtime();
#endif

      data_offset = 0;
      data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;

#ifdef PIDX_VAR_SLOW_LOOP
      for (k = 0; k < io_id->idx_derived_ptr->agg_buffer->var_number; k++)
        for (i = 0; i < io_id->idx_ptr->variable[k]->values_per_sample; i++)
          data_offset = data_offset + io_id->idx_ptr->variable[k]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[k]->bits_per_value/8) /* io_id->idx_derived_ptr->aggregation_factor*/;

      for (i = 0; i < io_id->idx_derived_ptr->agg_buffer->sample_number; i++)
        data_offset = data_offset + io_id->idx_ptr->variable[k]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_write_at(fh, data_offset, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)) , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }

      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != ((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)) )
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return -1;
      }

#else
      pwrite(fh, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->blocks_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)), data_offset);
#endif

#else
      for (k = 0; k < io_id->idx_derived_ptr->agg_buffer->var_number; k++)
        for (i = 0; i < io_id->idx_ptr->variable[k]->values_per_sample; i++)
          data_offset = (int64_t) data_offset + (int64_t) io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (bytes_per_datatype) /* io_id->idx_derived_ptr->aggregation_factor */;

      for (i = 0; i < io_id->idx_derived_ptr->agg_buffer->sample_number; i++)
        data_offset = (int64_t) data_offset + (int64_t) io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype);

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_write_at(fh, data_offset, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor)  * (bytes_per_datatype)) , MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long) data_offset, __FILE__, __LINE__, file_name);
        return -1;
      }

      MPI_Get_count(&status, MPI_BYTE, &write_count);
      //printf("[B] Elemets to write %d\n", write_count);
      if (write_count != ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return -1;
      }

#else
      pwrite(fh, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)), data_offset);
#endif

#endif

#ifdef PIDX_RECORD_TIME
      t3 = MPI_Wtime();
#endif

#if PIDX_HAVE_MPI
      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }
#else
      close(fh);
#endif

#ifdef PIDX_RECORD_TIME
      t4 = MPI_Wtime();
#endif

#ifdef PIDX_RECORD_TIME
      printf("B. [R %d] [O %lld C %lld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, (long long)((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor)  * (bytes_per_datatype)), io_id->idx_derived_ptr->agg_buffer->file_number, io_id->idx_derived_ptr->agg_buffer->var_number, io_id->idx_derived_ptr->agg_buffer->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
    }
  }

  if (io_id->idx_ptr->enable_compression == 1 && io_id->idx_ptr->compression_type == 2) // zfp done without block restructuring
  {
    if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
    {
      bytes_per_datatype =  ((io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8)  * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])) / (64/io_id->idx_ptr->compression_bit_rate);
      generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);
      mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed filename %s.\n", __FILE__, __LINE__, file_name);
        return -1;
      }

      data_offset = 0;
      total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * io_id->idx_ptr->variable_count;

      headers = (uint32_t*)malloc(total_header_size);
      memset(headers, 0, total_header_size);

      if (enable_caching == 1)
        memcpy (headers, cached_header_copy, total_header_size);

      uint64_t buffer_size = (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)));
      uint64_t header_size = (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size);
      //unsigned char* temp_buffer = (unsigned char*)realloc(io_id->idx_derived_ptr->agg_buffer->buffer, bufer_size+header_size);

      zfp_params zfp;
      zfp.type = ZFP_TYPE_DOUBLE;
      zfp.nx = 4;
      zfp.ny = 4;
      zfp.nz = 4;
      uint64_t input_block_size = zfp.nx * zfp.ny * zfp.nz * sizeof(double);
      uint64_t compressed_block_size = zfp_estimate_compressed_size(&zfp);
      assert(compressed_block_size > 0);
      uint64_t input_offset = 0;
      uint64_t output_offset = 0;
      uint64_t numBlocks = buffer_size / input_block_size;
      uint64_t compressed_buffer_size = compressed_block_size * numBlocks;
      unsigned char* compressed_buffer = (unsigned char*)malloc(compressed_buffer_size);
      uint64_t num_blocks_encoded = 0;
      while (input_offset < buffer_size)
      {
        assert(output_offset < compressed_buffer_size);
        uint64_t bytes_written = zfp_compress(&zfp, io_id->idx_derived_ptr->agg_buffer->buffer + input_offset,
                                          compressed_buffer + output_offset, compressed_block_size);
        assert(bytes_written > 0);
        input_offset += input_block_size;
        output_offset += bytes_written;
        ++num_blocks_encoded;
      }
      assert(num_blocks_encoded == numBlocks);
      assert(output_offset == compressed_buffer_size);

      unsigned char* temp_buffer = (unsigned char*)realloc(compressed_buffer, compressed_buffer_size+header_size);

      if (temp_buffer == NULL)
      {
        fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
        return -1;
      }
      else
      {
        io_id->idx_derived_ptr->agg_buffer->buffer = temp_buffer;
        memmove(io_id->idx_derived_ptr->agg_buffer->buffer + io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) );
        memcpy(io_id->idx_derived_ptr->agg_buffer->buffer, headers, total_header_size);
        memset(io_id->idx_derived_ptr->agg_buffer->buffer + total_header_size, 0, (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size - total_header_size));
      }
      free(headers);

      mpi_ret = MPI_File_write_at(fh, 0, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), MPI_BYTE, &status);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
        return -1;
      }

      MPI_Get_count(&status, MPI_BYTE, &write_count);
      //printf("[A] Elemets to write %d\n", write_count);
      if (write_count != (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block  / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return -1;
      }

      mpi_ret = MPI_File_close(&fh);
      if (mpi_ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return -1;
      }
    }
  }

  return 0;
}

int PIDX_io_aggregated_read(PIDX_io_id io_id)
{
  int64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0, rank, mpi_ret;
  //uint32_t *headers;
  //int total_header_size;
  int write_count;
  int bytes_per_datatype;

#ifdef PIDX_RECORD_TIME
  double t1, t2, t3, t4, t5;
#endif

#if PIDX_HAVE_MPI
  MPI_File fh;
  MPI_Status status;
  MPI_Comm_rank(io_id->comm, &rank);
#else
  int fh;
#endif

  if (io_id->idx_derived_ptr->agg_buffer->var_number == io_id->start_var_index && io_id->idx_derived_ptr->agg_buffer->sample_number == 0)
  {
    bytes_per_datatype =  (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8)  * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4])  / (64/io_id->idx_ptr->compression_bit_rate);

#ifdef PIDX_RECORD_TIME
    t1 = MPI_Wtime();
#endif
    generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
    mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed filename %s.\n", __FILE__, __LINE__, file_name);
      return -1;
    }
#else
    fh = open(file_name, O_RDONLY);
#endif

#ifdef PIDX_RECORD_TIME
    t2 = MPI_Wtime();
#endif


#ifdef PIDX_RECORD_TIME
    t3 = MPI_Wtime();
#endif


#if PIDX_HAVE_MPI
    mpi_ret = MPI_File_read_at(fh, (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))), MPI_BYTE, &status);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
      return -1;
    }


    MPI_Get_count(&status, MPI_BYTE, &write_count);
    //printf("[A] Elemets to write %d\n", write_count);
    if (write_count != (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * ( io_id->idx_derived_ptr->samples_per_block  / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))))
    {
      fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
      return -1;
    }

#else
    pread(fh, io_id->idx_derived_ptr->agg_buffer->buffer, (((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))), (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size));

    //double x[8];
    //memcpy(x, io_id->idx_derived_ptr->agg_buffer->buffer, sizeof(double) * 8);
    //printf("count = %f %f %f %f %f %f %f %f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);

#endif


#ifdef PIDX_RECORD_TIME
    t4 = MPI_Wtime();
#endif

#if PIDX_HAVE_MPI
    mpi_ret = MPI_File_close(&fh);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
      return -1;
    }
#else
      close(fh);
#endif

#ifdef PIDX_RECORD_TIME
    t5 = MPI_Wtime();
#endif

#ifdef PIDX_RECORD_TIME
    printf("A. [R %d] [OS 0 %ld %ld (%d %d %d) %ld] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank,
            (long int)(((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype))) + (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size),
            (long int)io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype,
            (int)io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number], (int)(io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor), (int)bytes_per_datatype,
            (io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size), io_id->idx_derived_ptr->agg_buffer->file_number, io_id->idx_derived_ptr->agg_buffer->var_number, io_id->idx_derived_ptr->agg_buffer->sample_number, (t2-t1), (t3-t2), (t4-t3), (t5-t4));
#endif
  }
  else if (io_id->idx_derived_ptr->agg_buffer->var_number != -1 && io_id->idx_derived_ptr->agg_buffer->sample_number != -1 && io_id->idx_derived_ptr->agg_buffer->file_number != -1)
  {
#ifdef PIDX_RECORD_TIME
    t1 = MPI_Wtime();
#endif
    bytes_per_datatype =  (io_id->idx_ptr->variable[io_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) * (io_id->idx_ptr->compression_block_size[0] * io_id->idx_ptr->compression_block_size[1] * io_id->idx_ptr->compression_block_size[2] * io_id->idx_ptr->compression_block_size[3] * io_id->idx_ptr->compression_block_size[4]);

    generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) io_id->idx_derived_ptr->agg_buffer->file_number, file_name, PATH_MAX);


#if PIDX_HAVE_MPI
    mpi_ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return -1;
    }
#else
    fh = open(file_name, O_WRONLY);
#endif

#ifdef PIDX_RECORD_TIME
    t2 = MPI_Wtime();
#endif

    data_offset = 0;
    data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;


    for (k = 0; k < io_id->idx_derived_ptr->agg_buffer->var_number; k++)
      for (i = 0; i < io_id->idx_ptr->variable[k]->values_per_sample; i++)
        data_offset = (int64_t) data_offset + (int64_t) io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[k]->bits_per_value/8) /* io_id->idx_derived_ptr->aggregation_factor */;

    for (i = 0; i < io_id->idx_derived_ptr->agg_buffer->sample_number; i++)
      data_offset = (int64_t) data_offset + (int64_t) io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number] * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype);

#if PIDX_HAVE_MPI

    mpi_ret = MPI_File_read_at(fh, data_offset, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor)  * (bytes_per_datatype)) , MPI_BYTE, &status);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_open() failed.\n", (long long) data_offset, __FILE__, __LINE__);
      return -1;
    }

    MPI_Get_count(&status, MPI_BYTE, &write_count);
    //printf("[B] Elemets to write %d\n", write_count);
    if (write_count != ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)))
    {
      fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
      return -1;
    }

#else
    pread(fh, io_id->idx_derived_ptr->agg_buffer->buffer, ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block / io_id->idx_derived_ptr->aggregation_factor) * (bytes_per_datatype)), data_offset);
#endif


#ifdef PIDX_RECORD_TIME
    t3 = MPI_Wtime();
#endif

#if PIDX_HAVE_MPI
    mpi_ret = MPI_File_close(&fh);
    if (mpi_ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
      return -1;
    }
#else
    close(fh);
#endif

#ifdef PIDX_RECORD_TIME
    t4 = MPI_Wtime();
#endif

#ifdef PIDX_RECORD_TIME
    //printf("B. [R %d] [OS %lld %d] [FVS %d %d %d] Time: O %f H %f W %f C %f\n", rank, (long long) data_offset, ((io_id->idx_derived_ptr->existing_blocks_index_per_file[io_id->idx_derived_ptr->agg_buffer->file_number]) * (io_id->idx_derived_ptr->samples_per_block/io_id->idx_derived_ptr->aggregation_factor)  * (bytes_per_datatype)), io_id->idx_derived_ptr->agg_buffer->file_number, io_id->idx_derived_ptr->agg_buffer->var_number, io_id->idx_derived_ptr->agg_buffer->sample_number, (t2-t1), (t2-t2), (t3-t2), (t4-t3));
#endif
  }
  return 0;
}


int PIDX_io_per_process_write(PIDX_io_id io_id)
{
  int e1 = 0, i = 0, p = 0, var = 0, ret;
  int send_index = 0;
  int64_t hz_index = 0;
  int64_t index = 0, count = 0;

  for (p = 0; p < io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 0)
    {
      for (i = 0; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i++)
        hz_index = hz_index + io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i];

      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        for(e1 = 0; e1 < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
        {
          if(e1 == 0)
          {
            index = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
            send_index = e1;
            count = 1;

            if(io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
            {
              for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
              {
                ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], (send_index), PIDX_WRITE);
                if (ret == -1)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return -1;
                }
                //printf("X");
                //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
              }
            }
          }
          else
          {
            if(io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
            {
              count++;
              if(e1 == io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
              {
                for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
                {
                  ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret == -1)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return -1;
                  }
                  //printf("Y");
                  //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
                }
              }
            }
            else
            {
              for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
              {
                ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
                if (ret == -1)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return -1;
                }
                //printf("Z");
                //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
              }

              if(e1 == io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
              {
                for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
                {
                  ret = write_read_samples(io_id, var, io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, PIDX_WRITE);
                  if (ret == -1)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return -1;
                  }
                  //printf("M");
                  //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, 1, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
                }
              }
              index = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
              count = 1;
              send_index = e1;
            }
          }
          hz_index++;
        }
      }
    }

    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 1)
    {
      for (i = /*io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from*/io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to - 1; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        if (io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
        {
          for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
          {
            index = 0;
            count =  io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1;
            ret = write_read_samples(io_id, var, io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_WRITE);
            if (ret == -1)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return -1;
            }
          }
        }
      }
    }

    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 2)
    {
      int start_block_index, end_block_index, bl;
      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        for (var = io_id->start_var_index; var <= io_id->end_var_index; var++)
        {
          start_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
          end_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          send_index = 0;
          for (bl = start_block_index; bl <= end_block_index; bl++)
          {
            if (end_block_index == start_block_index)
            {
              index = 0;
              count = (io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1);
            }
            else
            {
              if (bl == start_block_index)
              {
                index = 0;
                count = ((start_block_index + 1) * io_id->idx_derived_ptr->samples_per_block) - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i];
              }
              else if (bl == end_block_index)
              {
                index = (end_block_index * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
                count = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - ((end_block_index) * io_id->idx_derived_ptr->samples_per_block) + 1;
              }
              else
              {
                index = (bl * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
                count = io_id->idx_derived_ptr->samples_per_block;
              }
            }


            ret = write_read_samples(io_id, var, index + io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
            if (ret == -1)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return -1;
            }

            send_index = send_index + count;
          }
        }
      }
    }
  }
  return 0;
}

int PIDX_io_per_process_read(PIDX_io_id io_id)
{
  int e1 = 0, i = 0, p = 0, var = 0, ret;
  int send_index = 0;
  int64_t hz_index = 0;
  int64_t index = 0, count = 0;

  for (p = 0; p < io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 0)
    {
      for (i = 0; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i++)
        hz_index = hz_index + io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i];

      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        for(e1 = 0; e1 < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
        {
          if(e1 == 0)
          {
            index = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
            send_index = e1;
            count = 1;

            if(io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
            {
              for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
              {
                ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], (send_index), PIDX_READ);
                if (ret == -1)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return -1;
                }
                //printf("X");
                //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
              }
            }
          }
          else
          {
            if(io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
            {
              count++;
              if(e1 == io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
              {
                for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
                {
                  ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_READ);
                  if (ret == -1)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return -1;
                  }
                  //printf("Y");
                  //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
                }
              }
            }
            else
            {
              for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
              {
                ret = write_read_samples(io_id, var, index, count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_READ);
                if (ret == -1)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return -1;
                }
                //printf("Z");
                //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, count, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
              }

              if(e1 == io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
              {
                for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
                {
                  ret = write_read_samples(io_id, var, io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, PIDX_READ);
                  if (ret == -1)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return -1;
                  }
                  //printf("M");
                  //print_buffer(io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, 1, io_id->idx_ptr->variable[var]->values_per_sample, io_id->idx_ptr->variable[var]->bits_per_value, io_id->idx_ptr->variable[var]->type_name);
                }
              }
              index = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
              count = 1;
              send_index = e1;
            }
          }
          hz_index++;
        }
      }
    }

    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 1)
    {
      for (i = /*io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from*/io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to - 1; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        if (io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
        {
          for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
          {
            index = 0;
            count =  io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1;
            ret = write_read_samples(io_id, var, io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_READ);
            if (ret == -1)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return -1;
            }
          }
        }
      }
    }

    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->box_group_type == 2)
    {
      int start_block_index, end_block_index, bl;
      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        for (var = io_id->start_var_index; var <= io_id->end_var_index; var++)
        {
          start_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
          end_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          send_index = 0;
          for (bl = start_block_index; bl <= end_block_index; bl++)
          {
            if (end_block_index == start_block_index)
            {
              index = 0;
              count = (io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1);
            }
            else
            {
              if (bl == start_block_index)
              {
                index = 0;
                count = ((start_block_index + 1) * io_id->idx_derived_ptr->samples_per_block) - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i];
              }
              else if (bl == end_block_index)
              {
                index = (end_block_index * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
                count = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - ((end_block_index) * io_id->idx_derived_ptr->samples_per_block) + 1;
              }
              else
              {
                index = (bl * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
                count = io_id->idx_derived_ptr->samples_per_block;
              }
            }
            //printf("[%d] Offset %d Count %d\n", i, index + io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count);
            //printf(" %d [B%d] INDEX and COUNT and SEND INDEX %d %d %d\n", bl, i, index, count, send_index);
            //if(bl == 14)
            //{
              //double check;
              //memcpy(&check, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i] + send_index * (io_id->idx_ptr->variable[var]->bits_per_value / 8) * io_id->idx_ptr->variable[var]->values_per_sample, sizeof(double));
              //printf("%d %d %d (%d) VALUE = %f\n", var, p, i, (send_index  * io_id->idx_ptr->variable[var]->values_per_sample), check);
            //}

            ret = write_read_samples(io_id, var, index + io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_READ);
            if (ret == -1)
            {
              fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
              return -1;
            }

            send_index = send_index + count;
          }
        }
      }
    }
  }
  return 0;
}


int PIDX_io_finalize(PIDX_io_id io_id)
{

  free(io_id);
  io_id = 0;

  return 0;
}
