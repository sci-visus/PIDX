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
#undef PIDX_PRINT_AGG
//#define RANK_ORDER 1

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
  /// Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;

  int start_var_index;
  int end_var_index;

  int aggregator_interval;
};

int decompress_buffer(PIDX_agg_id agg_id, void* buffer, int length);

enum IO_MODE { PIDX_READ, PIDX_WRITE};

PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  agg_id->idx_ptr = idx_meta_data;
  agg_id->idx_derived_ptr = idx_derived_ptr;
  agg_id->start_var_index = start_var_index;
  agg_id->end_var_index = end_var_index;

  return agg_id;
}

#if PIDX_HAVE_MPI
int PIDX_agg_set_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  agg_id->comm = comm;
  //MPI_Comm_dup(comm, &agg_id->comm);
  return 0;
}

int PIDX_agg_set_global_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  agg_id->global_comm = comm;
  //MPI_Comm_dup(comm, &agg_id->comm);
  return 0;
}
#endif

/*
static double now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}
*/

#ifdef PIDX_HAVE_LOSSY_ZFP
// compress 'inbytes' bytes from 'data' to stream 'FPZ'
static int compress(FPZ* fpz, const void* data, size_t inbytes, const char* medium)
{
  //fprintf(stderr, "compressing to %s\n", medium);
  //double t = now();
#ifndef WITHOUT_HEADER
  // write header (optional)
  if (!fpzip_write_header(fpz))
  {
    fprintf(stderr, "cannot write header: %s\n", fpzip_errstr[fpzip_errno]);
    return 0;
  }
#endif
  // perform actual compression
  size_t outbytes = fpzip_write(fpz, data);
  if (!outbytes)
  {
    fprintf(stderr, "compression failed: %s\n", fpzip_errstr[fpzip_errno]);
    return 0;
  }
  //t = now() - t;
  //fprintf(stderr, "in=%zu out=%zu ratio=%.2f seconds=%.3f MB/s=%.3f\n", inbytes, outbytes, (double)inbytes / outbytes, t, (double)inbytes / (1024 * 1024 * t));
  return (int)outbytes;
}

int compress_buffer(PIDX_agg_id agg_id, void* buffer, int length)
{
  int i = 0, j = 0;
  int prec = 0;
  int nf = 1;
  int type = FPZIP_TYPE_DOUBLE;

  int64_t total_compression_block_size = agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4];

  //LOSSLESS Compression
  if (agg_id->idx_ptr->compression_type == 0)
  {
    size_t size = (type == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
    if (prec == 0)
      prec = CHAR_BIT * size;

    size_t total_size = 0;

    for ( i = 0; i < length; i = i + total_compression_block_size * size)
    {
      void* decompressed_buffer = malloc(total_compression_block_size * size + 1024);
      //printf("inint decomp = %d\n", (int)(total_compression_block_size * size + 1024));
      FPZ* fpz = fpzip_write_to_buffer(decompressed_buffer, (total_compression_block_size * size + 1024));
      fpz->type = type;
      fpz->prec = prec;
      fpz->nx = agg_id->idx_ptr->compression_block_size[0];
      fpz->ny = agg_id->idx_ptr->compression_block_size[1];
      fpz->nz = agg_id->idx_ptr->compression_block_size[2];
      fpz->nf = nf;

      //printf("offset = %d input size = %d\n", i, (int)total_compression_block_size * size);
      size_t x = compress(fpz, buffer + i, total_compression_block_size * size, "memory");

      if (x / sizeof(double) != 0)
        x = ((x / sizeof(double)) + 1) * sizeof(double);
      //printf("[%d] x = %ld\n", i, x/sizeof(double));
      //printf("decompressed buffe size = %ld %ld\n", x);
      //if (!compress(fpz, data, inbytes, "memory"))
      //  return EXIT_FAILURE;
      fpzip_write_close(fpz);

      //TODO: big problem
      double value = 0;
      for (j = 0; j < x ; j = j + 8)
      {
        memcpy(&value, decompressed_buffer + j, 8);
        if (isnan(value) != 0)
          printf("NAN\n");
      }
      //printf("reallocating to %d\n", (int)x);
      void* temp_buffer = (double*)realloc(decompressed_buffer, x);
      if (temp_buffer == NULL)
        printf("realloc error\n!!!");
      else
        decompressed_buffer = temp_buffer;
      //printf("copy location %d\n", (int)total_size);
      memcpy(buffer + total_size, decompressed_buffer, x);

      total_size = total_size + x;
      free(decompressed_buffer);
    }
    return (int)total_size;
  }
  //LOSSY Compression
  else if (agg_id->idx_ptr->compression_type == 1)
  {
    size_t typesize, outsize;
    typesize = sizeof(double);
    zfp_params params;
    params.type = ZFP_TYPE_DOUBLE;
    params.nx = agg_id->idx_ptr->compression_block_size[0];
    params.ny = agg_id->idx_ptr->compression_block_size[1];
    params.nz = agg_id->idx_ptr->compression_block_size[2];
    int total_size = 0;
    zfp_set_rate(&params, agg_id->idx_ptr->compression_bit_rate);
    outsize = zfp_estimate_compressed_size(&params);
    unsigned char* zip = malloc(outsize);
    for ( i = 0; i < length; i = i + total_compression_block_size * typesize)
    {
      memset(zip, 0, outsize);
      outsize = zfp_compress(&params, buffer + i , zip, outsize);
      if (outsize == 0)
      {
        fprintf(stderr, "compression failed\n");
        return -1;
      }
      memcpy(buffer + total_size, zip, outsize);
      total_size = total_size + outsize;
    }
    free(zip);
    return total_size;
  }

  return -1;
}

int decompress_buffer(PIDX_agg_id agg_id, void* buffer, int length)
{
  int64_t compression_block_num_elems = agg_id->idx_ptr->compression_block_size[0] *
                                        agg_id->idx_ptr->compression_block_size[1] *
                                        agg_id->idx_ptr->compression_block_size[2];
  unsigned int type_size = sizeof(double); // DUONG_HARDCODE
  if (agg_id->idx_ptr->compression_type == 1) // zfp lossy compression
  {
    zfp_params zfp;
    zfp.type = ZFP_TYPE_DOUBLE; // DUONG_HARDCODE
    zfp.nx = agg_id->idx_ptr->compression_block_size[0];
    zfp.ny = agg_id->idx_ptr->compression_block_size[1];
    zfp.nz = agg_id->idx_ptr->compression_block_size[2];
    unsigned int bit_rate = (unsigned int) agg_id->idx_ptr->compression_bit_rate;
    zfp_set_rate(&zfp, bit_rate);
    unsigned int input_offset = 0;
    unsigned int output_offset = 0;
    unsigned int input_buffer_size = (unsigned int) length;
    unsigned int output_buffer_size = (length / (bit_rate / 8.0)) * type_size; // DUONG_HARDCODE
    unsigned char* srcPtr = buffer;
    unsigned char* dstPtr = (unsigned char*) malloc(output_buffer_size);
    unsigned int compressed_block_size = compression_block_num_elems * (bit_rate / 8.0);
    unsigned int uncompressed_block_size = compression_block_num_elems * type_size;
    while (input_offset < input_buffer_size) // decompress one compression block at a time
    {
      int success = zfp_decompress(&zfp, dstPtr + output_offset, srcPtr + input_offset, compressed_block_size);
      if (success == 0) // failure
      {
        free(dstPtr);
        return -1;
      }
      input_offset += compressed_block_size;
      output_offset += uncompressed_block_size;
      assert(output_offset <= output_buffer_size);
    }
    memcpy(buffer, dstPtr, output_buffer_size);
    free(dstPtr);
    return 0;
  }
  return -1;
}


#endif

int aggregate_write_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, int MODE)
{
  int ret;
  int rank = 0, itr;// nrank = 0;
  int bytes_per_datatype;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, values_per_sample;
  int target_rank = 0;
  int64_t start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, hz_start = 0, samples_in_file = 0;
  int64_t samples_per_file = (int64_t) agg_id->idx_derived_ptr->samples_per_block * agg_id->idx_ptr->blocks_per_file;
  //MPI_Aint target_disp_address;

  int64_t total_compression_block_size = agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4];

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
  //MPI_Comm_rank(agg_id->global_comm, &nrank);
#endif

  values_per_sample = agg_id->idx_ptr->variable[variable_index]->values_per_sample; //number of samples for variable j

  //starting HZ index for the data buffer at level "level" and for regular box number "box"
  hz_start = hz_start_index;

  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / agg_id->idx_derived_ptr->samples_per_block;

  //number of empty blocks befor block "block_no" in the file "file_no"
#ifdef PIDX_VAR_SLOW_LOOP
  negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx_ptr->blocks_per_file, block_no, agg_id->idx_ptr->variable[variable_index]->VAR_global_block_layout);
  assert(negative_block_offset >= 0);

  //number of samples in file "file_no"
  samples_in_file = agg_id->idx_ptr->variable[variable_index]->VAR_blocks_per_file[file_no] * agg_id->idx_derived_ptr->samples_per_block;
  assert(samples_in_file <= samples_per_file);
#else
  negative_block_offset = PIDX_blocks_find_negative_offset(agg_id->idx_ptr->blocks_per_file, block_no, agg_id->idx_derived_ptr->global_block_layout);
  assert(negative_block_offset >= 0);

  //number of samples in file "file_no"
  samples_in_file = agg_id->idx_derived_ptr->existing_blocks_index_per_file[file_no] * agg_id->idx_derived_ptr->samples_per_block;
  assert(samples_in_file <= samples_per_file);
#endif

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  assert(hz_start >= (samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block));
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);
  assert(target_disp >= 0);

  sample_index = target_disp / (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);
  assert(sample_index < agg_id->idx_ptr->variable[variable_index]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor);

  target_disp = target_disp % (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);

#if RANK_ORDER
  target_rank = agg_id->idx_derived_ptr->agg_buffer->rank_holder[variable_index - agg_id->start_var_index][sample_index][file_no];
#else
  target_rank = agg_id->idx_derived_ptr->agg_buffer->rank_holder[file_no][variable_index - agg_id->start_var_index][sample_index];
#endif
  target_count = hz_count * values_per_sample;

  bytes_per_datatype = ((agg_id->idx_ptr->variable[variable_index]->bits_per_value / 8) * total_compression_block_size) / (64 / agg_id->idx_ptr->compression_bit_rate);

  //printf("CCCCCCCCCCCc : %d = %d * %d / %d\n", bytes_per_datatype, (agg_id->idx_ptr->variable[variable_index]->bits_per_value / 8), total_compression_block_size, (64 / agg_id->idx_ptr->compression_bit_rate));
  
  
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;

  start_agg_index = target_disp / (int64_t) (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);
  end_agg_index = ((target_disp + target_count - 1) / (int64_t) (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor));
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[A] Count %lld Local Disp %d Target Disp %lld\n", (long long)((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp), 0, (long long)target_disp);
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[A] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank, (long long)((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp), 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        if (agg_id->idx_ptr->enable_compression == 1)
        {
#ifdef PIDX_HAVE_LOSSY_ZFP
          int length = 0;
          length = compress_buffer(agg_id, hz_buffer, ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate));

          ret = MPI_Put(hz_buffer, length, MPI_BYTE, target_rank, target_disp, length, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
#else
          if (rank == 0)
            printf("Compression Library not available\n");
#endif
        }
        else
        {
          ret = MPI_Put(hz_buffer, ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MA] Count %lld Local Disp %d Target Disp %lld\n", (long long)target_disp, 0, (long long)((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp));
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MA] Count %lld Local Disp %d Target Disp %lld\n", (long long)((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp), 0, (long long) target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype);

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
#ifdef PIDX_PRINT_AGG
          if (rank == 0)
            printf("[B] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[B] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_id->aggregator_interval), (long long)(samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (long long)(( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
          if (agg_id->idx_ptr->enable_compression == 1)
          {
#ifdef PIDX_HAVE_LOSSY_ZFP
            int length = 0;
            unsigned char* offseted_buffer = hz_buffer + (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate);
            length = compress_buffer(agg_id, offseted_buffer, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate));

            ret = MPI_Put(offseted_buffer, length, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, length, MPI_BYTE, agg_id->win);
            if(ret != MPI_SUCCESS)
            {
              fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
              return (-1);
            }
#else
            if (rank == 0)
              printf("Compression Library not available\n");
#endif
          }
          else
          {
            ret = MPI_Put(hz_buffer + (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
              return (-1);
            }
          }
        }
        else
        {
          ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
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
#ifdef PIDX_PRINT_AGG
          if (rank == 0)
            printf("[MB] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[MB] Count %lld Local Disp %lld Target Disp %d\n", (long long)(samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
            fflush(agg_dump_fp);
          }
#endif
          memcpy( agg_id->idx_derived_ptr->agg_buffer->buffer, hz_buffer + (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, ( samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype);
        }
        else
          memcpy( hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, agg_id->idx_derived_ptr->agg_buffer->buffer, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype);
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[C] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[C] Target Rank %d Count %lld Local Disp %lld Target Disp %d\n", (target_rank + agg_id->aggregator_interval), (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif
        if (agg_id->idx_ptr->enable_compression == 1)
        {
#ifdef PIDX_HAVE_LOSSY_ZFP
          int length = 0;
          unsigned char* offseted_buffer = hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate);
          length = compress_buffer(agg_id, offseted_buffer, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))) * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate));

          ret = MPI_Put(offseted_buffer, length, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, length, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
#else
          if (rank == 0)
            printf("Compression Library not available\n");
#endif
        }
        else
        {
          ret = MPI_Put(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype,
                      MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MC] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MC] Count %lld Local Disp %lld Target Disp %d\n", (long long)(target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (long long)(((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
          fflush(agg_dump_fp);
        }
#endif

        memcpy( agg_id->idx_derived_ptr->agg_buffer->buffer, hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, agg_id->idx_derived_ptr->agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype);
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[D] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        if (agg_id->idx_ptr->enable_compression == 1)
        {
#ifdef PIDX_HAVE_LOSSY_ZFP
          int length = 0;
          unsigned char* offseted_buffer = hz_buffer;
          length = compress_buffer(agg_id, offseted_buffer, hz_count * values_per_sample * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate));

          ret = MPI_Put(offseted_buffer, length, MPI_BYTE, target_rank, target_disp, length, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
#else
          if (rank == 0)
            printf("Compression Library not available\n");
#endif
        }
        else
        {
          ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
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
#ifdef PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
#endif

#ifdef PIDX_DUMP_AGG
        if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        if (agg_id->idx_ptr->enable_compression == 1)
        {
#ifdef PIDX_HAVE_LOSSY_ZFP
          int length = 0;
          unsigned char* offseted_buffer = hz_buffer;
          /*
          double test1;
          double x, y, z;
          memcpy(&x, offseted_buffer, sizeof(double));
          memcpy(&y, offseted_buffer + sizeof(double), sizeof(double));
          memcpy(&z, offseted_buffer + sizeof(double) + sizeof(double), sizeof(double));
          printf("MMMMMMMMMMM: %f %f %f [%d]\n", x, y, z, hz_count * values_per_sample * bytes_per_datatype);
          */
          length = compress_buffer(agg_id, offseted_buffer, hz_count * values_per_sample * bytes_per_datatype * (64 / agg_id->idx_ptr->compression_bit_rate));

          printf("Original : Compressed :: %d (%d %d %d) : %d\n", (int)hz_count * values_per_sample * bytes_per_datatype, (int)hz_count, (int)values_per_sample, (int)bytes_per_datatype, (int)length);
          memcpy( agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, offseted_buffer, length);
          /*
          memcpy(&x, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, sizeof(double));
          memcpy(&y, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype + sizeof(double), sizeof(double));
          memcpy(&z, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype + sizeof(double) + sizeof(double), sizeof(double));
          printf("UUUUUUUUUUU: %f %f %f\n", x, y, z);
          */
#else
          if (rank == 0)
            printf("Compression Library not available\n");
#endif
        }
        else
          memcpy( agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
      }
      else
      {
        //double x;
        //memcpy(&x, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, bytes_per_datatype);
        //printf("count %d = %f\n", hz_count, x);
        memcpy( hz_buffer, agg_id->idx_derived_ptr->agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
      }
    }
  }
  return PIDX_success;
}

int decompess_read(PIDX_agg_id agg_id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, int buffer_offset, int MODE)
{
  int bytes_per_datatype;
  int values_per_sample;
  
  int64_t total_compression_block_size = agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4];

  bytes_per_datatype = ((agg_id->idx_ptr->variable[variable_index]->bits_per_value / 8) * total_compression_block_size) / (64 / agg_id->idx_ptr->compression_bit_rate);
  
  values_per_sample = agg_id->idx_ptr->variable[variable_index]->values_per_sample; //number of samples for variable j
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;
  
  //DUONG_HARDCODE
  
#ifdef PIDX_HAVE_LOSSY_ZFP
  if (agg_id->idx_ptr->enable_compression == 1)
  {
    decompress_buffer(agg_id, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
  }
#else
  printf("Compression build required !");
#endif

  return PIDX_success;
}

int PIDX_agg_buf_create(PIDX_agg_id agg_id)
{
  int i, j, k, var;
  int rank_counter = 0, no_of_aggregators = 0, nprocs = 1, rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  agg_id->idx_derived_ptr->agg_buffer->buffer_size = 0;
  agg_id->idx_derived_ptr->agg_buffer->sample_number = -1;
  agg_id->idx_derived_ptr->agg_buffer->var_number = -1;
  agg_id->idx_derived_ptr->agg_buffer->file_number = -1;

  agg_id->idx_derived_ptr->agg_buffer->compressed_block_size = malloc( sizeof(*(agg_id->idx_derived_ptr->agg_buffer->compressed_block_size)) * agg_id->idx_ptr->blocks_per_file);
  memset(agg_id->idx_derived_ptr->agg_buffer->compressed_block_size, 0, sizeof(*(agg_id->idx_derived_ptr->agg_buffer->compressed_block_size)) * agg_id->idx_ptr->blocks_per_file);


#ifdef PIDX_VAR_SLOW_LOOP
  for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    no_of_aggregators = no_of_aggregators + agg_id->idx_ptr->variable[var]->values_per_sample * agg_id->idx_ptr->variable[var]->VAR_existing_file_count;
#else
  for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    no_of_aggregators = no_of_aggregators + agg_id->idx_ptr->variable[var]->values_per_sample * agg_id->idx_derived_ptr->existing_file_count;
#endif
  agg_id->aggregator_interval = nprocs/ (no_of_aggregators * agg_id->idx_derived_ptr->aggregation_factor);

  assert(agg_id->aggregator_interval != 0);
#if RANK_ORDER
  agg_id->idx_derived_ptr->agg_buffer->rank_holder = malloc((agg_id->end_var_index - agg_id->start_var_index + 1) * sizeof (int**));
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index] = malloc( agg_id->idx_ptr->variable[i]->values_per_sample  * sizeof (int*) * agg_id->idx_derived_ptr->aggregation_factor);
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
    {
      agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j] = malloc(/*agg_id->idx_ptr->variable[i]->existing_file_count*/ agg_id->idx_derived_ptr->max_file_count * sizeof (int));
      memset(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j], 0, agg_id->idx_derived_ptr->max_file_count * sizeof (int));
    }
  }
#else
  agg_id->idx_derived_ptr->agg_buffer->rank_holder = malloc(agg_id->idx_derived_ptr->max_file_count * sizeof (int**));
  for (i = 0; i < agg_id->idx_derived_ptr->max_file_count; i++)
  {
    agg_id->idx_derived_ptr->agg_buffer->rank_holder[i] = malloc( (agg_id->end_var_index - agg_id->start_var_index + 1)  * sizeof (int*));
    for (j = agg_id->start_var_index; j <= agg_id->end_var_index; j++)
    {
      agg_id->idx_derived_ptr->agg_buffer->rank_holder[i][j - agg_id->start_var_index] = malloc( agg_id->idx_ptr->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_derived_ptr->aggregation_factor);
      memset(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i][j - agg_id->start_var_index], 0, agg_id->idx_ptr->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_derived_ptr->aggregation_factor);
    }
  }
#endif

  rank_counter = 0;
#if RANK_ORDER

#ifdef PIDX_VAR_SLOW_LOOP
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
    {
      for (k = 0; k < agg_id->idx_ptr->variable[i]->VAR_existing_file_count; k++)
      {
        agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;

        if(rank == agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]])
        {
          agg_id->idx_derived_ptr->agg_buffer->file_number = agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k];
          agg_id->idx_derived_ptr->agg_buffer->var_number = i;
          agg_id->idx_derived_ptr->agg_buffer->sample_number = j;

          agg_id->idx_derived_ptr->agg_buffer->buffer_size = agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number] * (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor) * (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8);
          agg_id->idx_derived_ptr->agg_buffer->buffer = malloc(agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          memset(agg_id->idx_derived_ptr->agg_buffer->buffer, 0, agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          //printf("Aggregator Rank %d Buffer Size %d (Var no: %d) (Sample no: %d) (File no: %d) (%d x %d x %d)\n", rank, agg_id->idx_derived_ptr->agg_buffer->buffer_size, agg_id->idx_derived_ptr->agg_buffer->var_number, agg_id->idx_derived_ptr->agg_buffer->sample_number, agg_id->idx_derived_ptr->agg_buffer->file_number, agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->blocks_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number], agg_id->idx_derived_ptr->samples_per_block, (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));
        }
      }
    }
  }
#else
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
    {
      for (k = 0; k < agg_id->idx_derived_ptr->existing_file_count; k++)
      {
        agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_derived_ptr->existing_file_index[k]] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;

        if(rank == agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_derived_ptr->existing_file_index[k]])
        {
          agg_id->idx_derived_ptr->agg_buffer->file_number = agg_id->idx_derived_ptr->existing_file_index[k];
          agg_id->idx_derived_ptr->agg_buffer->var_number = i;
          agg_id->idx_derived_ptr->agg_buffer->sample_number = j;

          agg_id->idx_derived_ptr->agg_buffer->buffer_size = agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number] * (agg_id->idx_derived_ptr->samples_per_block/ agg_id->idx_derived_ptr->aggregation_factor) * (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8);
          agg_id->idx_derived_ptr->agg_buffer->buffer = malloc(agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          memset(agg_id->idx_derived_ptr->agg_buffer->buffer, 0, agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          //printf("Aggregator Rank %d Buffer Size %d (Var no: %d) (Sample no: %d) (File no: %d) (%d x %d x %d)\n", rank, agg_id->idx_derived_ptr->agg_buffer->buffer_size, agg_id->idx_derived_ptr->agg_buffer->var_number, agg_id->idx_derived_ptr->agg_buffer->sample_number, agg_id->idx_derived_ptr->agg_buffer->file_number, agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number], agg_id->idx_derived_ptr->samples_per_block, (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));
        }
      }
    }
  }
#endif

#else

#ifdef PIDX_VAR_SLOW_LOOP
  for (k = 0; k < agg_id->idx_ptr->variable[i]->VAR_existing_file_count; k++)
  {
    for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
    {
      for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
      {
        agg_id->idx_derived_ptr->agg_buffer->rank_holder[agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]][i - agg_id->start_var_index][j] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;

        if(rank == agg_id->idx_derived_ptr->agg_buffer->rank_holder[agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]][i - agg_id->start_var_index][j])
        {
          agg_id->idx_derived_ptr->agg_buffer->file_number = agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k];
          agg_id->idx_derived_ptr->agg_buffer->var_number = i;
          agg_id->idx_derived_ptr->agg_buffer->sample_number = j;

          agg_id->idx_derived_ptr->agg_buffer->buffer_size = agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->VAR_blocks_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number] * (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor) * (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8);
          agg_id->idx_derived_ptr->agg_buffer->buffer = malloc(agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          memset(agg_id->idx_derived_ptr->agg_buffer->buffer, 0, agg_id->idx_derived_ptr->agg_buffer->buffer_size);
        }
      }
    }
  }
#else
  for (k = 0; k < agg_id->idx_derived_ptr->existing_file_count; k++)
  {
    for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
    {
      for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
      {
        agg_id->idx_derived_ptr->agg_buffer->rank_holder[agg_id->idx_derived_ptr->existing_file_index[k]][i - agg_id->start_var_index][j] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;

        if(rank == agg_id->idx_derived_ptr->agg_buffer->rank_holder[agg_id->idx_derived_ptr->existing_file_index[k]][i - agg_id->start_var_index][j])
        {
          agg_id->idx_derived_ptr->agg_buffer->file_number = agg_id->idx_derived_ptr->existing_file_index[k];
          agg_id->idx_derived_ptr->agg_buffer->var_number = i;
          agg_id->idx_derived_ptr->agg_buffer->sample_number = j;

          agg_id->idx_derived_ptr->agg_buffer->buffer_size = ((agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number] * (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor) * (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8)) / (64 / agg_id->idx_ptr->compression_bit_rate)) * agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4];
          
          if (agg_id->idx_ptr->current_time_step == 0)
            printf("[Aggregator] [%d] [%d %d %d] : %lld (%d %d %d %d %d)\n", rank, i, j, k, 
                 (long long) agg_id->idx_derived_ptr->agg_buffer->buffer_size, 
                 (int)agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number], 
                 (int)(agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor), 
                 (int)(agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8),
                 (int)(64 / agg_id->idx_ptr->compression_bit_rate), (int)(agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4])
                );

          agg_id->idx_derived_ptr->agg_buffer->buffer = malloc(agg_id->idx_derived_ptr->agg_buffer->buffer_size);
          if (agg_id->idx_derived_ptr->agg_buffer->buffer == NULL)
          {
            printf("[Aggregation] [%d] [%d %d %d] : %lld (%d %d (%d/%d) %d)\n", rank, i, j, k, (long long) agg_id->idx_derived_ptr->agg_buffer->buffer_size, agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_id->idx_derived_ptr->agg_buffer->file_number], (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor), agg_id->idx_derived_ptr->samples_per_block, agg_id->idx_derived_ptr->aggregation_factor, (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));

            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) agg_id->idx_derived_ptr->agg_buffer->buffer_size, __LINE__, __FILE__);
            return (-1);
          }
#if 0
          if (agg_id->idx_ptr->enable_compression == 1)
          {
            double nan_d = NAN;
            //printf("Aggregatio buffer size = %d\n", (int)agg_id->idx_derived_ptr->agg_buffer->buffer_size);
            for (l = 0; l < agg_id->idx_derived_ptr->agg_buffer->buffer_size; l = l + (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8))
              memcpy (agg_id->idx_derived_ptr->agg_buffer->buffer + l, &nan_d, (agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8));
              //agg_id->idx_derived_ptr->agg_buffer->buffer[l] = NAN;
          }
          else
            memset(agg_id->idx_derived_ptr->agg_buffer->buffer, 0, agg_id->idx_derived_ptr->agg_buffer->buffer_size);
#endif
        }
      }
    }
  }
#endif

#endif

  return PIDX_success;
}

int PIDX_agg_write(PIDX_agg_id agg_id)
{
  int i, p, e1, var, ret = 0;
  int send_index = 0;
  int64_t index = 0, count = 0, hz_index = 0;
  int variable_order = 1;
  //int element_count = 0;

  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
  {
    char agg_file_name[1024];

    ret = mkdir(agg_id->idx_derived_ptr->agg_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, agg_id->idx_derived_ptr->agg_dump_dir_name);
      return -1;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(agg_id->comm);
#endif

    sprintf(agg_file_name, "%s/rank_%d", agg_id->idx_derived_ptr->agg_dump_dir_name, rank);
    agg_dump_fp = fopen(agg_file_name, "a+");
    if (!agg_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] agg_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, agg_file_name);
      return (-1);
    }
  }
#endif

#if PIDX_HAVE_MPI
  //agg_id->idx_derived_ptr->win_time_start = MPI_Wtime();
  if (agg_id->idx_derived_ptr->agg_buffer->buffer_size != 0)
    MPI_Win_create(agg_id->idx_derived_ptr->agg_buffer->buffer, agg_id->idx_derived_ptr->agg_buffer->buffer_size, ((agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8) * agg_id->idx_ptr->compression_block_size[0] * agg_id->idx_ptr->compression_block_size[1] * agg_id->idx_ptr->compression_block_size[2] * agg_id->idx_ptr->compression_block_size[3] * agg_id->idx_ptr->compression_block_size[4]) / (64/agg_id->idx_ptr->compression_bit_rate), MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  else
    MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  //agg_id->idx_derived_ptr->win_time_end = MPI_Wtime();
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
#endif


  for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;

    if(agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_ptr[p]->box_group_type == 0)
    {
      for (i = 0; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i++)
        hz_index = hz_index + agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i];

      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to - agg_id->idx_derived_ptr->resolution_to; i++)
      {
        if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
        {
          for(e1 = 0; e1 < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
          {
            if(e1 == 0)
            {
              index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
              send_index = e1;
              count = 1;

              if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
              {
                for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                {
                  //printf("[A] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
                  ret = aggregate_write_read(agg_id, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret == -1)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return (-1);
                  }
                }
              }
            }
            else
            {
              if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
              {
                count++;
                if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
                {
                  for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                  {
                    //printf("[B] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
                    aggregate_write_read(agg_id, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
                    if (ret == -1)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return (-1);
                    }
                  }
                }
              }
              else
              {
                for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                {
                  //printf("[C] Size %lld Offset %lld\n", count, index);
                  aggregate_write_read(agg_id, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, PIDX_WRITE);
                  if (ret == -1)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return (-1);
                  }
                }

                if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
                {
                  for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                  {
                    //printf("[D] Size %lld Offset %lld\n", count, index);
                    aggregate_write_read(agg_id, var, agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, PIDX_WRITE);
                    if (ret == -1)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return (-1);
                    }
                  }
                }
                index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
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
      if (variable_order == 0)
      {
        for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to  - agg_id->idx_derived_ptr->resolution_to; i++)
        {
          if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
          {
            for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
            {
              index = 0;
              count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);

              ret = aggregate_write_read(agg_id, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_WRITE);
              if (ret == -1)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return (-1);
              }
            }
          }
        }
      }
      else
      {
        for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
        {
#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "Variable %d\n", var);
            fflush(agg_dump_fp);
          }
#endif
          for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to - agg_id->idx_derived_ptr->resolution_to; i++)
          {
            if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
            {
              index = 0;
              count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);

#ifdef PIDX_PRINT_AGG
              if (rank == 0)
                printf("[AGG] [Color %d] [VAR %d] [HZ %d] Size %lld Send Offset %lld\n", agg_id->idx_derived_ptr->color, var, i, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
#endif

#ifdef PIDX_DUMP_AGG
              if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
              {
                fprintf(agg_dump_fp, "[%d]: ", i);
                fflush(agg_dump_fp);
              }
#endif
              //agg_id->idx_derived_ptr->agg_level_start[p][var][i] = MPI_Wtime();
              ret = aggregate_write_read(agg_id, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_WRITE);
              if (ret == -1)
              {
                fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                return (-1);
              }
              //agg_id->idx_derived_ptr->agg_level_end[p][var][i] = MPI_Wtime();
            }
          }
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
  //agg_id->idx_derived_ptr->win_free_time_start = MPI_Wtime();
  MPI_Win_free(&(agg_id->win));
  //agg_id->idx_derived_ptr->win_free_time_end = MPI_Wtime();
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}

int PIDX_agg_read(PIDX_agg_id agg_id)
{
  int i, p, var, ret = 0;
  int64_t count = 0;

  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
  {
    char agg_file_name[1024];

    ret = mkdir(agg_id->idx_derived_ptr->agg_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, agg_id->idx_derived_ptr->agg_dump_dir_name);
      return -1;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(agg_id->comm);
#endif

    sprintf(agg_file_name, "%s/rank_%d", agg_id->idx_derived_ptr->agg_dump_dir_name, rank);
    agg_dump_fp = fopen(agg_file_name, "a+");
    if (!agg_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] agg_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, agg_file_name);
      return (-1);
    }
  }
#endif

#if PIDX_HAVE_MPI
  //agg_id->idx_derived_ptr->win_time_start = MPI_Wtime();
  if (agg_id->idx_derived_ptr->agg_buffer->buffer_size != 0)
    MPI_Win_create(agg_id->idx_derived_ptr->agg_buffer->buffer, agg_id->idx_derived_ptr->agg_buffer->buffer_size, agg_id->idx_ptr->variable[agg_id->idx_derived_ptr->agg_buffer->var_number]->bits_per_value/8, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  else
    MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  //agg_id->idx_derived_ptr->win_time_end = MPI_Wtime();
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
#endif


  for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
  {
    count = 0;

    for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    {
#ifdef PIDX_DUMP_AGG
      if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
      {
        fprintf(agg_dump_fp, "Variable %d\n", var);
        fflush(agg_dump_fp);
      }
#endif
      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to - agg_id->idx_derived_ptr->resolution_to; i++)
      {
        count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);

        //if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
        if (count != 0)
        {
#ifdef PIDX_PRINT_AGG
          if (rank == 0)
            printf("[AGG] [Color %d] [VAR %d] [HZ %d] Size %lld Send Offset %lld\n", agg_id->idx_derived_ptr->color, var, i, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
#endif

#ifdef PIDX_DUMP_AGG
          if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
          {
            fprintf(agg_dump_fp, "[%d]: ", i);
            fflush(agg_dump_fp);
          }
#endif
          //agg_id->idx_derived_ptr->agg_level_start[p][var][i] = MPI_Wtime();
          ret = aggregate_write_read(agg_id, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_READ);
          if (ret == -1)
          {
            fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
          //agg_id->idx_derived_ptr->agg_level_end[p][var][i] = MPI_Wtime();

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
  //agg_id->idx_derived_ptr->win_free_time_start = MPI_Wtime();
  MPI_Win_free(&(agg_id->win));
  //agg_id->idx_derived_ptr->win_free_time_end = MPI_Wtime();
#endif

#ifdef PIDX_DUMP_AGG
  if (agg_id->idx_derived_ptr->dump_agg_info == 1 && agg_id->idx_ptr->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  if (agg_id->idx_ptr->enable_compression == 1)
  {
    //printf("XXXXXXXX %d: %d %d\n", agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count, agg_id->start_var_index, agg_id->end_var_index);
    for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
    {
      for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
      {
        for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from + agg_id->idx_derived_ptr->resolution_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to - agg_id->idx_derived_ptr->resolution_to; i++)
        {
          //printf("YYYYYY %d: %d\n", i, agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i]);
          //if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
          {
            count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);

            double test;
            //memcpy(&test, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], sizeof(double));
            //printf("[%d] Before decompression = %f\n", i, test);
            
            ret = decompess_read(agg_id, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, PIDX_WRITE);
            if (ret == -1)
            {
              fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
              return (-1);
            }
            
            memcpy(&test, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], sizeof(double));
            printf("[%d] After decompression = %f\n", i, test);
          }
        }
      }
    }
  }
  
  return PIDX_success;
}

int PIDX_agg_buf_destroy(PIDX_agg_id agg_id)
{
  if (agg_id->idx_derived_ptr->agg_buffer->buffer_size != 0)
  {
    free(agg_id->idx_derived_ptr->agg_buffer->buffer);
    agg_id->idx_derived_ptr->agg_buffer->buffer = 0;
  }
  free(agg_id->idx_derived_ptr->agg_buffer->compressed_block_size);
  agg_id->idx_derived_ptr->agg_buffer->compressed_block_size = 0;
  int i = 0, j = 0;
#if RANK_ORDER
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      free(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j]);
      agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index][j] = 0;
    }
    free(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index]);
    agg_id->idx_derived_ptr->agg_buffer->rank_holder[i - agg_id->start_var_index] = 0;
  }
#else
  for (i = 0; i < agg_id->idx_derived_ptr->max_file_count; i++)
  {
    for (j = agg_id->start_var_index; j <= agg_id->end_var_index; j++)
    {
      free(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i][j - agg_id->start_var_index]);
      agg_id->idx_derived_ptr->agg_buffer->rank_holder[i][j - agg_id->start_var_index] = 0;
    }
    free(agg_id->idx_derived_ptr->agg_buffer->rank_holder[i]);
  }
#endif

  free(agg_id->idx_derived_ptr->agg_buffer->rank_holder);
  agg_id->idx_derived_ptr->agg_buffer->rank_holder = 0;

  return PIDX_success;
}

int PIDX_agg_finalize(PIDX_agg_id agg_id)
{

  free(agg_id);
  agg_id = 0;

  return 0;
}
