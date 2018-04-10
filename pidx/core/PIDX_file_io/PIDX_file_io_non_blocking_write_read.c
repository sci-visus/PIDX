#include "../../PIDX_inc.h"

#ifdef PIDX_DUMP_IO
static FILE* io_dump_fp;
#endif

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);


PIDX_return_code PIDX_file_io_async_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template)
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
