#include "../../PIDX_inc.h"

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);

PIDX_return_code PIDX_async_aggregated_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template)
{
  unsigned long long data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0, k = 0;
  int ret;

  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);

  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int)agg_buf->file_number, file_name, PATH_MAX);

    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, (fh));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    data_offset = 0;
    data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
    PIDX_variable_group var_grp = io_id->idx->variable_grp[io_id->group_index];

    for (k = 0; k < agg_buf->var_number; k++)
    {
      PIDX_variable vark = var_grp->variable[k];
      int bytes_per_datatype =  ((vark->bpv/8) * tck) / (io_id->idx->compression_factor);
      unsigned long long prev_var_sample = (unsigned long long) block_layout->bcpf[agg_buf->file_number] * io_id->idx_d->samples_per_block * bytes_per_datatype * var_grp->variable[k]->vps;

      data_offset = (unsigned long long) data_offset + prev_var_sample;
    }

    for (i = 0; i < agg_buf->sample_number; i++)
      data_offset = (unsigned long long) data_offset + agg_buf->buffer_size;

    if (io_id->idx->flip_endian == 1)
    {
      PIDX_variable curr_var = var_grp->variable[agg_buf->var_number];
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

    ret = MPI_File_iwrite_at(*fh, data_offset, agg_buf->buffer, agg_buf->buffer_size , MPI_BYTE, request);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_async_aggregated_read(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template)
{
  unsigned long long data_offset = 0;
  char file_name[PATH_MAX];
  int i = 0;
  MPI_File fp;
  uint32_t *headers;
  int ret;
  MPI_Status status;

#if 1
  int tck = (io_id->idx->chunk_size[0] * io_id->idx->chunk_size[1] * io_id->idx->chunk_size[2]);
  if (agg_buf->var_number != -1 && agg_buf->sample_number != -1 && agg_buf->file_number != -1)
  {
    generate_file_name(io_id->idx->blocks_per_file, filename_template, (unsigned int) /*adjusted_file_index*/agg_buf->file_number, file_name, PATH_MAX);

    ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
      return PIDX_err_io;
    }

    data_offset = 0;
    data_offset += io_id->idx_d->start_fs_block * io_id->idx_d->fs_block_size;
    PIDX_variable_group var_grp = io_id->idx->variable_grp[io_id->group_index];
    int total_header_size = (10 + (10 * io_id->idx->blocks_per_file)) * sizeof (uint32_t) * io_id->idx->variable_count;
    headers = malloc(total_header_size);
    memset(headers, 0, total_header_size);

    ret = MPI_File_read_at(fp, 0, headers, total_header_size , MPI_BYTE, &status);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
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
      if (PIDX_blocks_is_block_present(agg_buf->file_number * io_id->idx->blocks_per_file + i, block_layout))
      {
        data_offset = htonl(headers[12 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);
        data_size = htonl(headers[14 + ((i + (io_id->idx->blocks_per_file * agg_buf->var_number))*10 )]);

        int buffer_index = (block_count * io_id->idx_d->samples_per_block * (var_grp->variable[agg_buf->var_number]->bpv/8) * var_grp->variable[agg_buf->var_number]->vps * tck) / io_id->idx->compression_factor;
        ret = MPI_File_read_at(fp, data_offset, agg_buf->buffer + buffer_index, data_size, MPI_BYTE, &status);
        if (ret != MPI_SUCCESS)
        {
          fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }

        if (io_id->idx->flip_endian == 1)
        {
          PIDX_variable curr_var = var_grp->variable[agg_buf->var_number];
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
#endif

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
