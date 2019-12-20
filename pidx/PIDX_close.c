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

#include "PIDX_file_handler.h"

static void PIDX_debug_output(PIDX_file file, int svi, int evi, int io_type);
static PIDX_return_code PIDX_dump_state_finalize (PIDX_file file);
static int approx_maxh(PIDX_file file);

int pidx_global_variable = 0;

PIDX_return_code PIDX_flush(PIDX_file file)
{
  PIDX_time time = file->time;

  // making sure that variables are added to the dataset
  if (file->idx->variable_count <= 0)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_variable;
  }

  file->io = PIDX_io_init(file->idx, file->idx_c, file->idx_dbg, file->meta_data_cache, file->idx_b, file->restructured_grid, file->time, file->fs_block_size, file->variable_index_tracker);
  if (file->io == NULL)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_flush;
  }

  // this is an approximate calculation of total number of aggregation
  // groups so that timming buffers can be populated properly
  int agg_group_count = approx_maxh(file);

  // Populate all timming buffers
  // they need to be the first ones to be populated
  PIDX_init_timming_buffers1(time, file->idx->variable_count, agg_group_count);


  // index range of variables within a flush
  int lvi = file->local_variable_index;
  int lvc = file->local_variable_count;

  // currently only two modes are supported, one for write and other for read
  if (file->flags == MPI_MODE_CREATE)
  {
    if (PIDX_write(file->io, lvi, (lvi + lvc), file->idx->io_type) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  else if (file->flags == PIDX_MODE_RDONLY)
  {
    if (PIDX_read(file->io, lvi, (lvi + lvc), file->idx->io_type) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  // Output to stderr the timmings of all the io phases
  PIDX_debug_output(file, lvi, (lvi + lvc), file->idx->io_type);

  // delete timming buffers
  PIDX_delete_timming_buffers1(time, file->idx->variable_count);


  // finalize step
  PIDX_io_finalize(file->io);


  // freeing buffers
  for (int j = file->local_variable_index; j < file->local_variable_index + file->local_variable_count; j++)
  {
    for (int p = 0; p < file->idx->variable[j]->sim_patch_count; p++)
    {
      free(file->idx->variable[j]->sim_patch[p]);
      file->idx->variable[j]->sim_patch[p] = 0;
    }
  }

  // getting ready for next phase of flush
  file->local_variable_index = file->variable_index_tracker;
  file->local_variable_count = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_close(PIDX_file file)
{
  if (PIDX_flush(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  PIDX_time time = file->time;
  time->sim_end = PIDX_get_time();

  for (uint32_t j = 0; j < file->idx->variable_count; j++)
  {
    for (uint32_t k = 0; k < file->idx->variable[j]->sim_patch_count; k++)
    {
      free(file->idx->variable[j]->sim_patch[k]);
      file->idx->variable[j]->sim_patch[k] = 0;
    }
    free(file->idx->variable[j]);
    file->idx->variable[j] = 0;
  }

  file->idx->variable_count = 0;

  PIDX_dump_state_finalize(file);

  free(file->idx);
  free(file->restructured_grid);
  free(file->time);
  free(file->idx_dbg);
  free(file->idx_c);
  free(file->idx_b);

  free(file);

  return PIDX_success;
}


static int approx_maxh(PIDX_file file)
{
  int maxh = log2(getPowerOf2(file->idx->bounds[0])) + log2(getPowerOf2(file->idx->bounds[1])) + log2(getPowerOf2(file->idx->bounds[2])) + 1;
  int lc = maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file));

  if (lc < 1)
    return 1;
  else
    return maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file));
}


static void PIDX_debug_output(PIDX_file file, int svi, int evi, int io_type)
{
  pidx_global_variable++;

#if DETAIL_OUTPUT
  if (file->idx_c->simulation_rank == 0)
  {
#if 0
    if (file->idx->io_type == PIDX_io_type::PIDX_RAW_IO)
      fprintf(stderr, "PIDX_RAW_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_IDX_IO)
      fprintf(stderr, "PIDX_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_LOCAL_PARTITION_IDX_IO)
      fprintf(stderr, "PIDX_LOCAL_PARTITION_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_PARTICLE_IO)
      fprintf(stderr, "PIDX_PARTICLE_IO %s\n", file->idx->filename);
#endif

    if (file->idx->io_type != PIDX_RAW_IO)
    {
      fprintf(stderr, "[%d : %d %d] [%d %d %d : %d]\n", file->idx->current_time_step, file->idx_c->simulation_rank, file->idx_c->simulation_nprocs, (int) file->idx->bounds[0], (int) file->idx->bounds[1], (int) file->idx->bounds[2], file->idx->variable_count);
      fprintf(stderr, "Box set by user (PIDX_USER_RST_BOX) %d %d %d\n", (int)file->restructured_grid->patch_size[0], (int)file->restructured_grid->patch_size[1], (int)file->restructured_grid->patch_size[2]);

      fprintf(stderr, "Compression Bit rate set to %f\n", file->idx->compression_bit_rate);

      if (file->idx->endian == PIDX_LITTLE_ENDIAN)
        fprintf(stderr, "Little Endian | ");
      else if (file->idx->endian == PIDX_BIG_ENDIAN)
        fprintf(stderr, "Big Endian | ");

      if (file->idx->flip_endian == 1)
        fprintf(stderr, "Endian Flipping Done\n");
      if (file->idx->flip_endian == 0)
        fprintf(stderr, "Endian Flipping Not Done\n");

      fprintf(stderr, "Partition count %d = %d x %d x %d Partitio size = %d x %d x %d\n", file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2], file->idx->partition_count[0], file->idx->partition_count[1], file->idx->partition_count[2], file->idx->partition_size[0], file->idx->partition_size[1], file->idx->partition_size[2]);
      fprintf(stderr, "Comp = %d\n", file->idx->compression_type);
      fprintf(stderr, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx->max_file_count);
      fprintf(stderr, "maxh = %d\n", file->idx->maxh);
    }
  }
#endif

  PIDX_time time = file->time;
  time->sim_end = MPI_Wtime();

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->idx_c->simulation_comm);

  if (io_type == PIDX_IDX_IO || io_type == PIDX_LOCAL_PARTITION_IDX_IO)
  {
    if (max_time == total_time)
    {
      double group_init = time->init_end - time->init_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;

      double pre_group_total = group_init + group_reg_box;
      double post_group_total = group_bitstring + group_block_layout + header_io;

#if DETAIL_OUTPUT
      fprintf(stderr, "\n[T %d R %d N %d G %d]\n", file->idx->current_time_step, file->idx_c->simulation_rank, file->idx_c->simulation_nprocs, pidx_global_variable);
      fprintf(stderr, "PRE INIT          :[%.4f + %.4f] = %.4f\n\n", group_init, group_reg_box, pre_group_total);
      fprintf(stderr, "POST INIT         :[%.4f + %.4f + %.4f] = %.4f\n\n", group_bitstring, group_block_layout, header_io, post_group_total);
#endif

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0;
      double hz_init = 0, hz_meta_data = 0, hz_buffer = 0, hz = 0, hz_compress = 0, hz_buffer_free = 0, hz_cleanup = 0, hz_total = 0, hz_all = 0;
      double chunk_init = 0, chunk_meta = 0, chunk_buffer = 0, chunk = 0, chunk_buffer_free = 0, chunk_cleanup = 0, chunk_total = 0, chunk_all = 0;
      double compression_init = 0, compression = 0, compression_total = 0, compression_all = 0;
      double io = 0, io_all = 0;

      for (int si = svi; si < evi; si++)
      {
        rst_init = time->rst_init_end[si] - time->rst_init_start[si];
        rst_meta_data_create = time->rst_meta_data_create_end[si] - time->rst_meta_data_create_start[si];
        rst_meta_data_io = time->rst_meta_data_io_end[si] - time->rst_meta_data_io_start[si];
        rst_buffer = time->rst_buffer_end[si] - time->rst_buffer_start[si];
        rst_write_read = time->rst_write_read_end[si] - time->rst_write_read_start[si];
        rst_buff_agg = time->rst_buff_agg_end[si] - time->rst_buff_agg_start[si];
        rst_buff_agg_free = time->rst_buff_agg_free_end[si] - time->rst_buff_agg_free_start[si];
        rst_buff_agg_io = time->rst_buff_agg_io_end[si] - time->rst_buff_agg_io_start[si];
        rst_cleanup = time->rst_cleanup_end[si] - time->rst_cleanup_start[si];
        rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
        rst_all = rst_all + rst_total;

#if DETAIL_OUTPUT
        fprintf(stderr, "RST                         :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
#endif
      }

      double partition_time = time->partition_end - time->partition_start;
#if DETAIL_OUTPUT
      fprintf(stderr, "Partition time %f\n", partition_time);
#endif

      double grp_rst_hz_chunk_agg_io = pre_group_total + post_group_total + partition_time;
      double agg_all = 0;
      double hz_io_all = 0;

      for (int si = svi; si < evi; si++)
      {
        if (file->idx->variable_tracker[si] == 1)
        {
#if 0
          rst_init = time->rst_init_end[si] - time->rst_init_start[si];
          rst_meta_data_create = time->rst_meta_data_create_end[si] - time->rst_meta_data_create_start[si];
          rst_meta_data_io = time->rst_meta_data_io_end[si] - time->rst_meta_data_io_start[si];
          rst_buffer = time->rst_buffer_end[si] - time->rst_buffer_start[si];
          rst_write_read = time->rst_write_read_end[si] - time->rst_write_read_start[si];
          rst_buff_agg = time->rst_buff_agg_end[si] - time->rst_buff_agg_start[si];
          rst_buff_agg_free = time->rst_buff_agg_free_end[si] - time->rst_buff_agg_free_start[si];
          rst_buff_agg_io = time->rst_buff_agg_io_end[si] - time->rst_buff_agg_io_start[si];
          rst_cleanup = time->rst_cleanup_end[si] - time->rst_cleanup_start[si];
          rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
          rst_all = rst_all + rst_total;
#endif


          hz_init = time->hz_init_end[si] - time->hz_init_start[si];
          hz_meta_data = time->hz_meta_end[si] - time->hz_meta_start[si];
          hz_buffer = time->hz_buffer_end[si] - time->hz_buffer_start[si];
          hz = time->hz_end[si] - time->hz_start[si];
          hz_compress = time->hz_compress_end[si] - time->hz_compress_start[si];
          hz_buffer_free = time->hz_buffer_free_end[si] - time->hz_buffer_free_start[si];
          hz_cleanup = time->hz_cleanup_end[si] - time->hz_cleanup_start[si];
          hz_total = hz_init + hz_meta_data + hz_buffer + hz + hz_compress + hz_buffer_free + hz_cleanup;
          hz_all = hz_all + hz_total;


          chunk_init = time->chunk_init_end[si] - time->chunk_init_start[si];
          chunk_meta = time->chunk_meta_end[si] - time->chunk_meta_start[si];
          chunk_buffer = time->chunk_buffer_end[si] - time->chunk_buffer_start[si];
          chunk = time->chunk_end[si] - time->chunk_start[si];
          chunk_buffer_free = time->chunk_buffer_free_end[si] - time->chunk_buffer_free_start[si];
          chunk_cleanup = time->chunk_cleanup_end[si] - time->chunk_cleanup_start[si];
          chunk_total = chunk_init + chunk_meta + chunk_buffer + chunk + chunk_buffer_free + chunk_cleanup;
          chunk_all = chunk_all + chunk_total;


          compression_init = time->compression_init_end[si] - time->compression_init_start[si];
          compression = time->compression_end[si] - time->compression_start[si];
          compression_total = compression_init + compression;
          compression_all = compression_all + compression_total;

#if DETAIL_OUTPUT
          //fprintf(stderr, "RST           :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
          fprintf(stderr, "HZ            :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, hz_init, hz_meta_data, hz_buffer, hz_compress, hz, hz_buffer_free, hz_cleanup, hz_total);
          fprintf(stderr, "CHUNK         :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, chunk_init, chunk_meta, chunk_buffer, chunk, chunk_buffer_free, chunk_cleanup, chunk_total);
          fprintf(stderr, "CMP           :[%d] [%.4f + %.4f] = %.4f\n", si, compression_init, compression, compression_total);
#endif

          double hz_io = 0;
          for (int i = file->idx_b->agg_level; i < file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count ; i++)
          {
            hz_io = time->hz_io_end[si][i] - time->hz_io_start[si][i];
            hz_io_all = hz_io_all + hz_io;

#if DETAIL_OUTPUT
            fprintf(stderr, "[HZ I/O S %d %d]  :[%d] [%d] %.4f [%.4f] \n", file->idx_b->agg_level, file->idx_b->file0_agg_group_to_index, si, i, hz_io, hz_io_all);
#endif
          }
        }

        double agg_init = 0, agg_meta = 0, agg_buf = 0, agg = 0, agg_meta_cleanup = 0, agg_total = 0, agg_cmp = 0;
        for (int i = file->idx_b->file0_agg_group_from_index; i < file->idx_b->agg_level ; i++)
        {
          agg_init = time->agg_init_end[si][i] - time->agg_init_start[si][i];
          agg_meta = time->agg_meta_end[si][i] - time->agg_meta_start[si][i];
          agg_buf = time->agg_buf_end[si][i] - time->agg_buf_start[si][i];
          agg = time->agg_end[si][i] - time->agg_start[si][i];
          agg_cmp = time->agg_compress_end[si][i] - time->agg_compress_start[si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[si][i] - time->agg_meta_cleanup_start[si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG S %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx_b->file0_agg_group_from_index, file->idx_b->file0_agg_group_to_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }

#if 0
        for (int i = file->idx_b->file0_agg_group_from_index; i < file->idx_b->agg_l_shared ; i++)
        {
          agg_init = time->agg_init_end[si][i] - time->agg_init_start[si][i];
          agg_meta = time->agg_meta_end[si][i] - time->agg_meta_start[si][i];
          agg_buf = time->agg_buf_end[si][i] - time->agg_buf_start[si][i];
          agg = time->agg_end[si][i] - time->agg_start[si][i];
          agg_cmp = time->agg_compress_end[si][i] - time->agg_compress_start[si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[si][i] - time->agg_meta_cleanup_start[si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG S %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx_b->file0_agg_group_from_index, file->idx_b->file0_agg_group_to_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }

        for (int i = file->idx_b->nfile0_agg_group_from_index; i < file->idx_b->agg_l_nshared ; i++)
        {
          agg_init = time->agg_init_end[si][i] - time->agg_init_start[si][i];
          agg_meta = time->agg_meta_end[si][i] - time->agg_meta_start[si][i];
          agg_buf = time->agg_buf_end[si][i] - time->agg_buf_start[si][i];
          agg = time->agg_end[si][i] - time->agg_start[si][i];
          agg_cmp = time->agg_compress_end[si][i] - time->agg_compress_start[si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[si][i] - time->agg_meta_cleanup_start[si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);
          //fprintf(stderr, "[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta + agg_buf, agg,  agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG N %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx_b->nfile0_agg_group_from_index, file->idx_b->nfile0_agg_group_to_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }
#endif

        io = time->io_end[si] - time->io_start[si];
        io_all = io_all + io;
#if DETAIL_OUTPUT
        fprintf(stderr, "IO [%d]        :[%d] %.4f\n", file->idx->variable_count, si, io);
#endif
      }

      grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all + hz_all + hz_io_all + chunk_all + compression_all + agg_all + io_all;

#if DETAIL_OUTPUT
//      fprintf(stderr, "XIRPICCHHAI      :[%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f = %.4f] + %.4f [%.4f %.4f]\n", pre_group_total, rst_all, partition_time, post_group_total, chunk_all, compression_all, hz_all, hz_io_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);
fprintf(stderr, "[%s %d %d (%d %d %d)] IRPICCHHAI      :[%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f = %.4f] + %.4f [%.4f %.4f]\n", file->idx->filename, file->idx->current_time_step, (evi - svi), (int)file->idx->bounds[0], (int)file->idx->bounds[1], (int)file->idx->bounds[2], pre_group_total, rst_all, partition_time, post_group_total, chunk_all, compression_all, hz_all, hz_io_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);
#endif
    }
  }
  else if (io_type == PIDX_PARTICLE_IO)
  {
    if (max_time == total_time)
    {
      double init_time = time->SX - time->sim_start;
      double idx_time = time->header_io_end - time->header_io_start;
      double meta_time = time->particle_meta_data_io_end - time->particle_meta_data_io_start;
      double data_time = time->particle_data_io_end - time->particle_data_io_start;
#if DETAIL_OUTPUT
      fprintf(stderr, "Timestep %d Variables %d Time [Init + idx file + meta + data] [%f + %f + %f + %f] = %f [%f]\n", file->idx->current_time_step, (evi - svi), init_time, idx_time, meta_time, data_time, (idx_time + meta_time + data_time + init_time), max_time );
#endif
    }
  }
  else if (io_type == PIDX_RST_PARTICLE_IO)
  {
    if (max_time == total_time)
    {
      double init_time = time->SX - time->sim_start;
      double set_box_size = time->set_reg_box_end - time->set_reg_box_start;
      double idx_time = time->header_io_end - time->header_io_start;

      double rst_total = 0;
      for (int si = 0; si < 1; si++)
      {
        double rst_init = time->rst_init_end[si] - time->rst_init_start[si];
        double rst_meta_data_create = time->rst_meta_data_create_end[si] - time->rst_meta_data_create_start[si];
        double rst_meta_data_io = time->rst_meta_data_io_end[si] - time->rst_meta_data_io_start[si];
        double rst_buffer = time->rst_buffer_end[si] - time->rst_buffer_start[si];
        double rst_write_read = time->rst_write_read_end[si] - time->rst_write_read_start[si];
        double rst_buff_agg = time->rst_buff_agg_end[si] - time->rst_buff_agg_start[si];
        double rst_buff_agg_free = time->rst_buff_agg_free_end[si] - time->rst_buff_agg_free_start[si];
        double rst_buff_agg_io = time->rst_buff_agg_io_end[si] - time->rst_buff_agg_io_start[si];
        double rst_cleanup = time->rst_cleanup_end[si] - time->rst_cleanup_start[si];
        rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
#if DETAIL_OUTPUT
        fprintf(stderr, "RST: [%d] [%f + %f + %f + %f + %f + %f + %f + %f + %f] = %f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
#endif
      }
#if DETAIL_OUTPUT
      fprintf(stderr, "Timestep %d Variables %d Time [Init + set box size + idx file + rst_io] [%f + %f + %f + %f] = %f [%f]\n", file->idx->current_time_step, (evi - svi), init_time, set_box_size,  idx_time, rst_total, (init_time + set_box_size + idx_time + rst_total), max_time);
#endif
    }
  }
  else
  {
    if (max_time == total_time)
    {
      double group_init = time->init_end - time->init_start;
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;
      double group_total = group_init + group_bitstring + group_reg_box + group_block_layout + header_io;

      //fprintf(stderr, "[T %d R %d N %d G %d] INIT          :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", file->idx->current_time_step, file->idx_c->simulation_rank, file->idx_c->simulation_nprocs, pidx_global_variable, group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      //fprintf(stderr, "INIT      :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0, r1 = 0, r2 = 0, r3 = 0;
      //double grp_rst_hz_chunk_agg_io = group_total;

      for (int si = svi; si < evi; si++)
      {
        if (file->idx->variable_tracker[si] == 1)
        {
          rst_init = time->rst_init_end[si] - time->rst_init_start[si];
          rst_meta_data_create = time->rst_meta_data_create_end[si] - time->rst_meta_data_create_start[si];
          rst_meta_data_io = time->rst_meta_data_io_end[si] - time->rst_meta_data_io_start[si];
          rst_buffer = time->rst_buffer_end[si] - time->rst_buffer_start[si];
          r1 = r1 + rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer;

          rst_write_read = time->rst_write_read_end[si] - time->rst_write_read_start[si];
          rst_buff_agg = time->rst_buff_agg_end[si] - time->rst_buff_agg_start[si];
          rst_buff_agg_free = time->rst_buff_agg_free_end[si] - time->rst_buff_agg_free_start[si];
          r2 = r2 + rst_write_read + rst_buff_agg + rst_buff_agg_free;

          rst_buff_agg_io = time->rst_buff_agg_io_end[si] - time->rst_buff_agg_io_start[si];
          rst_cleanup = time->rst_cleanup_end[si] - time->rst_cleanup_start[si];
          r3 = r3 + rst_buff_agg_io + rst_cleanup;

          rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
          rst_all = rst_all + rst_total;

          //fprintf(stderr, "RST       :[%d] [%.4f + [MD] %.4f + %.4f + [B] %.4f + [R] %.4f + [MC] %.4f + %.4f + [IO] %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
        }
      }

      //grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all;
#if DETAIL_OUTPUT
      fprintf(stderr, "[RAW] [%s] [%d %d %d : %d %d %d] [%d] [T %d R %d N %d V %d] : [%.4f + %.4f (%.4f = %.4f + %.4f + %.4f) = %.4f] + %.4f [%.4f %.4f]\n", file->idx->filename, (int)file->idx->bounds[0], (int)file->idx->bounds[1], (int)file->idx->bounds[2], (int)file->restructured_grid->patch_size[0], (int)file->restructured_grid->patch_size[1], (int)file->restructured_grid->patch_size[2],  pidx_global_variable, file->idx->current_time_step, file->idx_c->simulation_rank, file->idx_c->simulation_nprocs, (evi - svi),
             group_total,
             rst_all,
             r1 + r2 + r3,
             r1, r2, r3,
             (group_total + rst_all),
             (time->SX - time->sim_start),
             (group_total + rst_all) + (time->SX - time->sim_start), max_time);
#endif
    }
  }
}


static PIDX_return_code PIDX_dump_state_finalize (PIDX_file file)
{
  if (file->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
    fclose(file->idx_dbg->debug_file_output_fp);

  return PIDX_success;
}
