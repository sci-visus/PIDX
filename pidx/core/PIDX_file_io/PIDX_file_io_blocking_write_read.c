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

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);



PIDX_return_code PIDX_file_io_blocking_read(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, char* filename_template)
{
  uint64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0;
  MPI_File fp;
  uint32_t *headers;
  int ret;
  MPI_Status status;

  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);
  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int) agg_buf->file_number, file_name, PATH_MAX);

    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    int total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
    headers = malloc(total_header_size);
    memset(headers, 0, total_header_size);

    ret = MPI_File_read_at(fp, 0, headers, total_header_size , MPI_BYTE, &status);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Data offset = [%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }
    int read_count = 0;
    MPI_Get_count(&status, MPI_BYTE, &read_count);
    if (read_count != total_header_size)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %dd\n", __FILE__, __LINE__, read_count, total_header_size);
      return PIDX_err_io;
    }

    int data_size = 0;
    int block_count = 0;
    for (i = 0; i < io_id->idx->blocks_per_file; i++)
    {
      if (PIDX_blocks_is_block_present(agg_buf->file_number * io_id->idx->blocks_per_file + i, io_id->idx->bits_per_block, block_layout))
      {
        data_offset = htonl(headers[12 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);
        data_size = htonl(headers[14 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);

        int buffer_index = (block_count * io_id->idx->samples_per_block * (io_id->idx->variable[agg_buf->var_number]->bpv/8) * io_id->idx->variable[agg_buf->var_number]->vps * tck) / io_id->idx->compression_factor;

        //fprintf(stderr, "DO and DS %d %d\n", data_offset, data_size);
        ret = MPI_File_read_at(fp, data_offset, agg_buf->buffer + buffer_index, data_size, MPI_BYTE, &status);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }

#if 0
        if (io_id->idx_d->color == 1 && agg_buf->buffer_size != 0)
        {
          double x1, x2, x3, x4, x5, x6, x7, x8;
          memcpy(&x1, agg_buf->buffer, sizeof(double));
          memcpy(&x2, agg_buf->buffer + sizeof(double), sizeof(double));
          memcpy(&x3, agg_buf->buffer + sizeof(double) * 2, sizeof(double));
          memcpy(&x4, agg_buf->buffer + sizeof(double) * 3, sizeof(double));
          memcpy(&x5, agg_buf->buffer + sizeof(double) * 4, sizeof(double));
          memcpy(&x6, agg_buf->buffer + sizeof(double) * 5, sizeof(double));
          memcpy(&x7, agg_buf->buffer + sizeof(double) * 6, sizeof(double));
          memcpy(&x8, agg_buf->buffer + sizeof(double) * 7, sizeof(double));
          fprintf(stderr, "x %f %f %f %f %f %f %f %f\n", x1, x2, x3, x4, x5, x6, x7, x8);
        }
#endif

#if 0
        int read_count = 0;
        MPI_Get_count(&status, MPI_BYTE, &read_count);
        if (read_count != data_size)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %d\n", __FILE__, __LINE__, read_count, total_header_size);
          return PIDX_err_io;
        }
#endif

        if (io_id->idx->flip_endian == 1)
        {
          PIDX_variable curr_var = io_id->idx->variable[agg_buf->var_number];
          if (curr_var->bpv/8 == 4 || curr_var->bpv/8 == 12)
          {
            int y = 0;
            float temp;
            float temp2;

            for (y = 0; y < data_size / sizeof(float); y++)
            {
              memcpy(&temp, agg_buf->buffer + buffer_index + (y * sizeof(float)), sizeof(float));
              bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
              memcpy(agg_buf->buffer + buffer_index + (y * sizeof(float)), &temp2, sizeof(float));
            }
          }
          else if (curr_var->bpv/8 == 8 || curr_var->bpv/8 == 24)
          {
            int y = 0;
            double temp;
            double temp2;

            for (y = 0; y < data_size / sizeof(double); y++)
            {
              memcpy(&temp, agg_buf->buffer + buffer_index + (y * sizeof(double)), sizeof(double));
              bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
              memcpy(agg_buf->buffer + buffer_index + (y * sizeof(double)), &temp2, sizeof(double));
            }
          }
        }

        block_count++;
      }
    }

    MPI_File_close(&fp);
    free(headers);
  }

  return PIDX_success;
}



PIDX_return_code PIDX_file_io_blocking_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, char* filename_template)
{
  uint64_t data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  int ret;
  MPI_File fh;
  MPI_Status status;
  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);

  uint64_t total_header_size;
  total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
  uint64_t start_fs_block = total_header_size / io_id->fs_block_size;
  if (total_header_size % io_id->fs_block_size)
    start_fs_block++;

  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int) agg_buf->file_number, file_name, PATH_MAX);
    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    data_offset = 0;
    data_offset += start_fs_block * io_id->fs_block_size;

    for (k = 0; k < agg_buf->var_number; k++)
    {
      PIDX_variable vark = io_id->idx->variable[k];
      int bytes_per_datatype =  ((vark->bpv/8) * tck) / (io_id->idx->compression_factor);
      uint64_t prev_var_sample = (uint64_t) block_layout->bcpf[agg_buf->file_number] * io_id->idx->samples_per_block * bytes_per_datatype * io_id->idx->variable[k]->vps;

      data_offset = (uint64_t) data_offset + prev_var_sample;
    }

    for (i = 0; i < agg_buf->sample_number; i++)
      data_offset = (uint64_t) data_offset + agg_buf->buffer_size;

    //fprintf(stderr, "DO %d DS %d\n", data_offset, agg_buf->buffer_size);
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

    ret = MPI_File_close(&fh);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}



static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
  unsigned char *data = ((unsigned char *)val) + 3;
  unsigned char *out = (unsigned char *)outbuf;

  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out = *data;

  return;
}



static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
  unsigned char *data = ((unsigned char *)val) + 7;
  unsigned char *out = (unsigned char *)outbuf;

  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out = *data;

  return;
}
