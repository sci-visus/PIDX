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

static PIDX_return_code dump_debug_data_init(PIDX_hz_encode_id hz_id);
static PIDX_return_code dump_debug_data_finalie (PIDX_hz_encode_id id);
static int write_read_samples(PIDX_hz_encode_id hz_id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, unsigned long long buffer_offset, PIDX_block_layout layout, int MODE);


int PIDX_file_io_per_process(PIDX_hz_encode_id hz_id, PIDX_block_layout block_layout, int MODE)
{
  int i = 0, p = 0, v = 0, ret;
  int send_index = 0;
  unsigned long long index = 0, count = 0;

  dump_debug_data_init(hz_id);

  PIDX_variable_group var_grp = hz_id->idx->variable_grp[hz_id->group_index];
  PIDX_variable var0 = var_grp->variable[hz_id->first_index];
  for (p = 0; p < var0->patch_group_count; p++)
  {
    index = 0, count = 0, send_index = 0;
    index = 0, count = 0, send_index = 0;
    if (var0->hz_buffer[p]->type == 1)
    {
      for(v = hz_id->first_index; v <= hz_id->last_index; v++)
      {

        if (hz_id->idx_d->dump_io_info == 1 && hz_id->idx->current_time_step == 0)
        {
          fprintf(hz_id->idx_d->io_dump_fp, "Variable %d\n", v);
          fflush(hz_id->idx_d->io_dump_fp);
        }

        HZ_buffer hz_buf = var_grp->variable[v]->hz_buffer[p];
        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
            index = 0;
            count =  var0->hz_buffer[p]->end_hz_index[i] - var0->hz_buffer[p]->start_hz_index[i] + 1;

            if (hz_id->idx_d->dump_io_info == 1 && hz_id->idx->current_time_step == 0)
            {
              fprintf(hz_id->idx_d->io_dump_fp, "[%d]: ", i);
              fflush(hz_id->idx_d->io_dump_fp);
            }

            ret = write_read_samples(hz_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout, MODE);
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
      for(v = hz_id->first_index; v <= hz_id->last_index; v++)
      {
        if (hz_id->idx_d->dump_io_info == 1 && hz_id->idx->current_time_step == 0)
        {
          fprintf(hz_id->idx_d->io_dump_fp, "Variable %d\n", v);
          fflush(hz_id->idx_d->io_dump_fp);
        }

        HZ_buffer hz_buf = var_grp->variable[v]->hz_buffer[p];
        for (i = block_layout->resolution_from; i < block_layout->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
            int start_block_index = var_grp->variable[v]->hz_buffer[p]->start_hz_index[i] / hz_id->idx_d->samples_per_block;
            int end_block_index = var_grp->variable[v]->hz_buffer[p]->end_hz_index[i] / hz_id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              index = 0;
              count = (var_grp->variable[v]->hz_buffer[p]->end_hz_index[i] - var_grp->variable[v]->hz_buffer[p]->start_hz_index[i] + 1);

              ret = write_read_samples(hz_id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, block_layout, MODE);
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
                    count = ((start_block_index + 1) * hz_id->idx_d->samples_per_block) - var_grp->variable[v]->hz_buffer[p]->start_hz_index[i];
                  }
                  else if (bl == end_block_index)
                  {
                    index = (end_block_index * hz_id->idx_d->samples_per_block - var_grp->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = var_grp->variable[v]->hz_buffer[p]->end_hz_index[i] - ((end_block_index) * hz_id->idx_d->samples_per_block) + 1;
                  }
                  else
                  {
                    index = (bl * hz_id->idx_d->samples_per_block - var_grp->variable[v]->hz_buffer[p]->start_hz_index[i]);
                    count = hz_id->idx_d->samples_per_block;
                  }

                  ret = write_read_samples(hz_id, v, index + var_grp->variable[v]->hz_buffer[p]->start_hz_index[i], count, var_grp->variable[v]->hz_buffer[p]->buffer[i], send_index, block_layout, MODE);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_io;
                  }
                  send_index = send_index + count;
                }
                else
                  send_index = send_index + hz_id->idx_d->samples_per_block;
              }
            }
          }
        }
      }
    }
  }

  dump_debug_data_finalie(hz_id);

  return PIDX_success;
}


static int write_read_samples(PIDX_hz_encode_id hz_id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, unsigned long long buffer_offset, PIDX_block_layout layout, int MODE)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, block_negative_offset = 0, file_number;
  int bytes_per_sample, bytes_per_datatype;
  int i = 0;
  char file_name[PATH_MAX];
  off_t data_offset = 0;

  PIDX_variable_group var_grp = hz_id->idx->variable_grp[hz_id->group_index];
  samples_per_file = hz_id->idx_d->samples_per_block * hz_id->idx->blocks_per_file;

  bytes_per_datatype = (var_grp->variable[variable_index]->bpv / 8) * (hz_id->idx->chunk_size[0] * hz_id->idx->chunk_size[1] * hz_id->idx->chunk_size[2]) / (hz_id->idx->compression_factor);
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * var_grp->variable[variable_index]->vps;


  while (hz_count)
  {
    block_number = hz_start_index / hz_id->idx_d->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;

    if ((unsigned long long)file_count > hz_count)
      file_count = hz_count;

    ret = generate_file_name(hz_id->idx->blocks_per_file, hz_id->idx->filename_template, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }

    data_offset = 0;
    bytes_per_sample = var_grp->variable[variable_index]->bpv / 8;
    data_offset = file_index * bytes_per_sample * var_grp->variable[variable_index]->vps;
    data_offset += hz_id->idx_d->start_fs_block * hz_id->idx_d->fs_block_size;

    block_negative_offset = PIDX_blocks_find_negative_offset(hz_id->idx->blocks_per_file, block_number, layout);

    data_offset -= block_negative_offset * hz_id->idx_d->samples_per_block * bytes_per_sample * var_grp->variable[variable_index]->vps;

    int l = 0;
    for (l = 0; l < variable_index; l++)
    {
      bytes_per_sample = var_grp->variable[l]->bpv / 8;
      for (i = 0; i < hz_id->idx->blocks_per_file; i++)
        if (PIDX_blocks_is_block_present((i + (hz_id->idx->blocks_per_file * file_number)), layout))
          data_offset = data_offset + (var_grp->variable[l]->vps * bytes_per_sample * hz_id->idx_d->samples_per_block);
    }

    if(MODE == PIDX_WRITE)
    {
      if (hz_id->idx_d->dump_io_info == 1 && hz_id->idx->current_time_step == 0)
      {
        fprintf(hz_id->idx_d->io_dump_fp, "[A] Count %lld Target Disp %d (%d %d)\n", (long long)file_count * var_grp->variable[variable_index]->vps * (var_grp->variable[variable_index]->bpv/8), (file_index * bytes_per_sample * var_grp->variable[variable_index]->vps - block_negative_offset * hz_id->idx_d->samples_per_block * bytes_per_sample * var_grp->variable[variable_index]->vps)/8, (int)hz_id->idx_d->start_fs_block, (int)hz_id->idx_d->fs_block_size);
        fflush(hz_id->idx_d->io_dump_fp);
      }

      MPI_File fh;
      MPI_Status status;
      int ret;
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed. (%s) [%d]\n", __FILE__, __LINE__, file_name, file_number);
        return PIDX_err_io;
      }

      ret = MPI_File_write_at(fh, data_offset, hz_buffer, file_count * var_grp->variable[variable_index]->vps * (var_grp->variable[variable_index]->bpv/8), MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      int write_count;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      if (write_count != file_count * var_grp->variable[variable_index]->vps * (var_grp->variable[variable_index]->bpv/8))
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      MPI_File_close(&fh);

    }
    if(MODE == PIDX_READ)
    {

      MPI_File fh;
      MPI_Status status;
      int ret;
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      ret = MPI_File_read_at(fh, data_offset, hz_buffer, file_count * var_grp->variable[variable_index]->vps * (var_grp->variable[variable_index]->bpv/8), MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      MPI_File_close(&fh);
    }
    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * var_grp->variable[variable_index]->vps * bytes_per_datatype;
  }
  return PIDX_success;
}


static PIDX_return_code dump_debug_data_init(PIDX_hz_encode_id hz_id)
{
  if (hz_id->idx_d->dump_io_info == 1 && hz_id->idx->current_time_step == 0)
  {
    int ret;
    char io_file_name[1024];
    ret = mkdir(hz_id->idx_d->io_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, hz_id->idx_d->io_dump_dir_name);
      return PIDX_err_io;
    }

    MPI_Barrier(hz_id->idx_c->local_comm);

    sprintf(io_file_name, "%s/rank_%d", hz_id->idx_d->io_dump_dir_name, hz_id->idx_c->lrank);
    hz_id->idx_d->io_dump_fp = fopen(io_file_name, "a+");
    if (!hz_id->idx_d->io_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, io_file_name);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}


static PIDX_return_code dump_debug_data_finalie (PIDX_hz_encode_id id)
{

  if (id->idx_d->dump_io_info == 1 && id->idx->current_time_step == 0)
  {
    fprintf(id->idx_d->rst_dump_fp, "\n");
    fclose(id->idx_d->rst_dump_fp);
  }

  return PIDX_success;
}
