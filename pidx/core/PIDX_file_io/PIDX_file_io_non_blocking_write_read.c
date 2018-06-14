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

#ifdef PIDX_DUMP_IO
static FILE* io_dump_fp;
#endif

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);


PIDX_return_code PIDX_file_io_async_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template)
{
  uint64_t data_offset = 0;
  char file_name[PATH_MAX];
  int ret;

  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);

  if (agg_buf->var_number != -1 && agg_buf->file_number != -1)
  {
    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int)agg_buf->file_number, file_name, PATH_MAX);

    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, (fh));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    uint64_t total_header_size;
    total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
    uint64_t start_fs_block = total_header_size / io_id->fs_block_size;
    if (total_header_size % io_id->fs_block_size)
      start_fs_block++;

    data_offset += start_fs_block * io_id->fs_block_size;

    for (int k = 0; k < agg_buf->var_number; k++)
    {
      PIDX_variable vark = io_id->idx->variable[k];
      int bytes_per_datatype =  ((vark->bpv/8) * tck) / (io_id->idx->compression_factor);
      uint64_t prev_var_sample = (uint64_t) block_layout->bcpf[agg_buf->file_number] * io_id->idx->samples_per_block * bytes_per_datatype * io_id->idx->variable[k]->vps;

      data_offset = (uint64_t) data_offset + prev_var_sample;
    }

    //for (i = 0; i < agg_buf->sample_number; i++)
    //  data_offset = (uint64_t) data_offset + agg_buf->buffer_size;

    if (io_id->idx->flip_endian == 1)
    {
      PIDX_variable curr_var = io_id->idx->variable[agg_buf->var_number];
      if (curr_var->bpv/8 == 4 || curr_var->bpv/8 == 12)
      {
        int y = 0;
        float temp;
        float temp2;

        for (y = 0; y < agg_buf->buffer_size / sizeof(float); y++)
        {
          memcpy(&temp, agg_buf->buffer + (y * sizeof(float)), sizeof(float));
          bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(agg_buf->buffer + (y * sizeof(float)), &temp2, sizeof(float));
        }
      }
      else if (curr_var->bpv/8 == 8 || curr_var->bpv/8 == 24)
      {
        int y = 0;
        double temp;
        double temp2;

        for (y = 0; y < agg_buf->buffer_size / sizeof(double); y++)
        {
          memcpy(&temp, agg_buf->buffer + (y * sizeof(double)), sizeof(double));
          bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(agg_buf->buffer + (y * sizeof(double)), &temp2, sizeof(double));
        }
      }
    }

    ret = MPI_File_iwrite_at(*fh, data_offset, agg_buf->buffer, agg_buf->buffer_size, MPI_BYTE, request);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
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
