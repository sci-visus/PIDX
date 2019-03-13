/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

static int write_samples(PIDX_hz_encode_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout);
static int read_samples(PIDX_hz_encode_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout);

static int opened_file_number = -1;
static uint32_t *headers;
static int total_header_size = 0;
static MPI_File fp = 0;
static MPI_Status status;


int PIDX_file_io_per_process(PIDX_hz_encode_id id, PIDX_block_layout block_layout, int MODE)
{
  int i = 0, v = 0, ret;
  int send_index = 0;
  uint64_t index = 0, count = 0;

  PIDX_variable var0 = id->idx->variable[id->first_index];

  total_header_size = (10 + (10 * id->idx->blocks_per_file)) * sizeof (uint32_t) * id->idx->variable_count;
  headers = malloc(total_header_size);
  memset(headers, 0, total_header_size);

  index = 0, count = 0, send_index = 0;

  if (var0->hz_buffer->is_boundary_HZ_buffer == 1)
  {
    for (v = id->first_index; v <= id->last_index; v++)
    {
      if (id->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || id->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        fprintf(id->idx_dbg->debug_file_output_fp, "Variable %d\n", v);
        fflush(id->idx_dbg->debug_file_output_fp);
      }

      HZ_buffer hz_buf = id->idx->variable[v]->hz_buffer;
      for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
      {
        if (var0->hz_buffer->nsamples_per_level[i][0] * var0->hz_buffer->nsamples_per_level[i][1] * var0->hz_buffer->nsamples_per_level[i][2] != 0)
        {
          index = 0;
          count =  var0->hz_buffer->end_hz_index[i] - var0->hz_buffer->start_hz_index[i] + 1;

          if (id->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || id->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
          {
            fprintf(id->idx_dbg->debug_file_output_fp, "[%d]: ", i);
            fflush(id->idx_dbg->debug_file_output_fp);
          }

          if (MODE == PIDX_WRITE)
            ret = write_samples(id, v, var0->hz_buffer->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout);
          else
            ret = read_samples(id, v, var0->hz_buffer->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout);
          if (ret != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
        }
      }
    }
  }
  else if (var0->hz_buffer->is_boundary_HZ_buffer == 2)
  {
    for (v = id->first_index; v <= id->last_index; v++)
    {
      if (id->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || id->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        fprintf(id->idx_dbg->debug_file_output_fp, "Variable %d\n", v);
        fflush(id->idx_dbg->debug_file_output_fp);
      }

      HZ_buffer hz_buf = id->idx->variable[v]->hz_buffer;
      for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
      {
        if (var0->hz_buffer->nsamples_per_level[i][0] * var0->hz_buffer->nsamples_per_level[i][1] * var0->hz_buffer->nsamples_per_level[i][2] != 0)
        {
          int start_block_index = id->idx->variable[v]->hz_buffer->start_hz_index[i] / id->idx->samples_per_block;
          int end_block_index = id->idx->variable[v]->hz_buffer->end_hz_index[i] / id->idx->samples_per_block;
          assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

          if (end_block_index == start_block_index)
          {
            index = 0;
            count = (id->idx->variable[v]->hz_buffer->end_hz_index[i] - id->idx->variable[v]->hz_buffer->start_hz_index[i] + 1);

            if (MODE == PIDX_WRITE)
              ret = write_samples(id, v, var0->hz_buffer->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout);
            else
              ret = read_samples(id, v, var0->hz_buffer->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout);
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
              if (PIDX_blocks_is_block_present(bl, id->idx->bits_per_block, block_layout))
              {
                if (bl == start_block_index)
                {
                  index = 0;
                  count = ((start_block_index + 1) * id->idx->samples_per_block) - id->idx->variable[v]->hz_buffer->start_hz_index[i];
                }
                else if (bl == end_block_index)
                {
                  index = (end_block_index * id->idx->samples_per_block - id->idx->variable[v]->hz_buffer->start_hz_index[i]);
                  count = id->idx->variable[v]->hz_buffer->end_hz_index[i] - ((end_block_index) * id->idx->samples_per_block) + 1;
                }
                else
                {
                  index = (bl * id->idx->samples_per_block - id->idx->variable[v]->hz_buffer->start_hz_index[i]);
                  count = id->idx->samples_per_block;
                }

                if (MODE == PIDX_WRITE)
                  ret = write_samples(id, v, index + id->idx->variable[v]->hz_buffer->start_hz_index[i], count, id->idx->variable[v]->hz_buffer->buffer[i], send_index, block_layout);
                else
                  ret = read_samples(id, v, index + id->idx->variable[v]->hz_buffer->start_hz_index[i], count, id->idx->variable[v]->hz_buffer->buffer[i], send_index, block_layout);
                if (ret != PIDX_success)
                {
                  fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_io;
                }
                send_index = send_index + count;
              }
              else
                send_index = send_index + id->idx->samples_per_block;
            }
          }
        }
      }
    }
  }

  free(headers);

  if (fp != 0)
    MPI_File_close(&fp);

  opened_file_number = -1;

  return PIDX_success;
}


static int write_samples(PIDX_hz_encode_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, block_negative_offset = 0, file_number;
  int bytes_per_datatype;
  int i = 0;
  char file_name[PATH_MAX];
  uint64_t data_offset = 0;

  PIDX_variable curr_var = id->idx->variable[variable_index];

  samples_per_file = id->idx->samples_per_block * id->idx->blocks_per_file;

  bytes_per_datatype = (curr_var->bpv / 8) * curr_var->vps * (id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2]) / id->idx->compression_factor;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype;

  while (hz_count)
  {
    block_number = hz_start_index / id->idx->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;

    ret = generate_file_name(id->idx->blocks_per_file, id->idx->filename_template_partition, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

    if (opened_file_number != file_number)
    {
      opened_file_number = file_number;

      if (fp != 0)
        MPI_File_close(&fp);

      //fprintf(stderr, "Opening file %s\n", file_name);
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }
    }

    if ((uint64_t)file_count > hz_count)
      file_count = hz_count;

    uint64_t total_header_size;
    total_header_size = (10 + (10 * id->idx->blocks_per_file)) * sizeof (uint32_t) * id->idx->variable_count;
    uint64_t start_fs_block = total_header_size / id->fs_block_size;
    if (total_header_size % id->fs_block_size)
      start_fs_block++;

    data_offset = 0;
    data_offset = file_index * bytes_per_datatype;
    data_offset += start_fs_block * id->fs_block_size;

    // Adjusting for missing blocks
    block_negative_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, id->idx->bits_per_block, block_number, layout);

    data_offset -= block_negative_offset * id->idx->samples_per_block * bytes_per_datatype;

    int l = 0;
    for (l = 0; l < variable_index; l++)
    {
      for (i = 0; i < id->idx->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (id->idx->blocks_per_file * file_number)), id->idx->bits_per_block, layout))
          data_offset = data_offset + (id->idx->variable[l]->vps * (id->idx->variable[l]->bpv / 8) * ((id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2]) / id->idx->compression_factor) * id->idx->samples_per_block);
    }


    if (id->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || id->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
    {
      fprintf(id->idx_dbg->debug_file_output_fp, "[A] Count %lld Target Disp %lld (%d %d)\n", (unsigned long long)file_count * bytes_per_datatype, (unsigned long long)(file_index * bytes_per_datatype - block_negative_offset * id->idx->samples_per_block * bytes_per_datatype), (int)start_fs_block, (int)id->fs_block_size);
      fflush(id->idx_dbg->debug_file_output_fp);
    }

    int ret;
    if (id->idx->flip_endian == 1)
    {
      if (curr_var->bpv/8 == 4 || curr_var->bpv/8 == 12)
      {
        int y = 0;
        float temp, temp2;

        for (y = 0; y < (file_count * curr_var->vps * (curr_var->bpv/8)) / sizeof(float); y++)
        {
          memcpy(&temp, hz_buffer + (y * sizeof(float)), sizeof(float));
          bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(hz_buffer + (y * sizeof(float)), &temp2, sizeof(float));
        }
      }
      else if (curr_var->bpv/8 == 8 || curr_var->bpv/8 == 24)
      {
        int y = 0;
        double temp, temp2;

        for (y = 0; y < (file_count * curr_var->vps * (curr_var->bpv/8)) / sizeof(double); y++)
        {
          memcpy(&temp, hz_buffer + (y * sizeof(double)), sizeof(double));
          bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
          memcpy(hz_buffer + (y * sizeof(double)), &temp2, sizeof(double));
        }
      }
    }

    //if (id->idx_c->partition_rank == 0)
    //  fprintf(stderr, "[%d] Data offset %d data size %lld\n", variable_index, data_offset, file_count * bytes_per_datatype);
    ret = MPI_File_write_at(fp, data_offset, hz_buffer, file_count * bytes_per_datatype, MPI_BYTE, &status);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed. wc %d %s\n", __FILE__, __LINE__, file_count, file_name);
      return PIDX_err_io;
    }

    int write_count;
    MPI_Get_count(&status, MPI_BYTE, &write_count);
    if (write_count != file_count * bytes_per_datatype)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * bytes_per_datatype;
  }
  return PIDX_success;
}



static int read_samples(PIDX_hz_encode_id id, int variable_index, uint64_t hz_start_index, uint64_t hz_count, unsigned char* hz_buffer, uint64_t buffer_offset, PIDX_block_layout layout)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, file_number, bytes_per_datatype;
  char file_name[PATH_MAX];
  uint64_t data_offset = 0;
  uint64_t data_size = 0;

  PIDX_variable curr_var = id->idx->variable[variable_index];

  samples_per_file = id->idx->samples_per_block * id->idx->blocks_per_file;

  bytes_per_datatype = (curr_var->bpv / 8) * curr_var->vps * (id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2]) / (id->idx->compression_factor);
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype;

  while (hz_count)
  {
    block_number = hz_start_index / id->idx->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;

    ret = generate_file_name(id->idx->blocks_per_file, id->idx->filename_template_partition, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

    if (opened_file_number != file_number)
    {
      opened_file_number = file_number;

      if (fp != 0)
        MPI_File_close(&fp);

      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }

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

    }

    if ((uint64_t)file_count > hz_count)
      file_count = hz_count;

    int block_size_bytes = id->idx->samples_per_block * bytes_per_datatype;
    unsigned char *temp_buffer = malloc(block_size_bytes);
    int bl = 0;
    int blocks_to_read = 0;
    if (file_count % id->idx->samples_per_block == 0)
      blocks_to_read = file_count / id->idx->samples_per_block;
    else
      blocks_to_read = (file_count / id->idx->samples_per_block) + 1;

    for (bl = 0; bl < blocks_to_read; bl++)
    {
      memset(temp_buffer, 0, block_size_bytes);
      data_offset = ntohl(headers[12 + ((((block_number % id->idx->blocks_per_file) + bl) + (id->idx->blocks_per_file * variable_index))*10 )]);
      data_size = ntohl(headers[14 + ((((block_number % id->idx->blocks_per_file) + bl) + (id->idx->blocks_per_file * variable_index))*10 )]);

      if (data_size == 0)
        continue;

      int ret;
      ret = MPI_File_read_at(fp, data_offset, temp_buffer, block_size_bytes, MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      if (bl == blocks_to_read - 1)
      {
        if (file_count % id->idx->samples_per_block == 0)
          memcpy (hz_buffer + bl * block_size_bytes, temp_buffer + (hz_start_index % id->idx->samples_per_block) * bytes_per_datatype, block_size_bytes);
        else
          memcpy (hz_buffer + bl * block_size_bytes, temp_buffer + (hz_start_index % id->idx->samples_per_block) * bytes_per_datatype, (file_count % id->idx->samples_per_block) * bytes_per_datatype);
      }
      else
        memcpy (hz_buffer + bl * block_size_bytes, temp_buffer + (hz_start_index % id->idx->samples_per_block) * bytes_per_datatype, block_size_bytes);
    }
    free(temp_buffer);

    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * bytes_per_datatype;
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
