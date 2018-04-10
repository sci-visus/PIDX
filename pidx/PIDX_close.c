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

#include "PIDX_file_handler.h"

static void PIDX_debug_output(PIDX_file file, int gi, int svi, int evi, int io_type);
static PIDX_return_code PIDX_dump_state_finalize (PIDX_file file);
static int approx_maxh(PIDX_file file);

int pidx_global_variable = 0;

PIDX_return_code PIDX_flush(PIDX_file file)
{
  int i;
  int ret = PIDX_success;
  int vgc = file->idx->variable_group_count;
  PIDX_time time = file->idx_d->time;

  if (file->idx->variable_count <= 0)
    return PIDX_err_variable;

  file->io = PIDX_io_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, file->idx_cache);
  if (file->io == NULL)
    return PIDX_err_flush;

  int block_layout_count = approx_maxh(file);

  PIDX_init_timming_buffers1(time, vgc, file->idx->variable_count, block_layout_count, file->idx_d->wavelet_levels);
  for (i = file->local_group_index; i < file->local_group_index + file->local_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    int lvi = var_grp->local_variable_index;
    int lvc = var_grp->local_variable_count;

    if (file->flags == MPI_MODE_CREATE)
      ret = PIDX_write(file->io, i, lvi, (lvi + lvc), file->idx->io_type);

    else if (file->flags == PIDX_MODE_RDONLY)
      ret = PIDX_read(file->io, i, lvi, (lvi + lvc), file->idx->io_type);

    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  for (i = file->local_group_index; i < file->local_group_index + file->local_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];

    int lvi = var_grp->local_variable_index;
    int lvc = var_grp->local_variable_count;
    PIDX_debug_output(file, i, lvi, (lvi + lvc), file->idx->io_type);
  }
  PIDX_delete_timming_buffers1(time, vgc, file->idx->variable_count);


  ret = PIDX_io_finalize(file->io);
  if (ret != PIDX_success)
    return PIDX_err_io;


  int j = 0, p = 0;
  for (i = file->local_group_index; i < file->local_group_index + file->local_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    for (j = var_grp->local_variable_index; j < var_grp->local_variable_index + var_grp->local_variable_count; j++)
    {
      for(p = 0; p < var_grp->variable[j]->sim_patch_count; p++)
      {
        free(var_grp->variable[j]->sim_patch[p]);
        var_grp->variable[j]->sim_patch[p] = 0;
      }
    }
    var_grp->local_variable_index = var_grp->variable_index_tracker;
    var_grp->local_variable_count = 0;
  }

  file->local_group_index = file->idx->group_index_tracker;
  file->local_group_count = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_close(PIDX_file file)
{
  int ret;
  int i = 0;
  int j = 0;
  int k = 0;
  file->write_on_close = 1;

  ret = PIDX_flush(file);
  if (ret != PIDX_success)
    return PIDX_err_close;

  PIDX_time time = file->idx_d->time;
  time->sim_end = PIDX_get_time();

  for (i = 0; i < file->idx->variable_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    for (j = 0; j < var_grp->variable_count; j++)
    {
      for (k = 0; k < var_grp->variable[j]->sim_patch_count; k++)
      {
        free(var_grp->variable[j]->sim_patch[k]);
        var_grp->variable[j]->sim_patch[k] = 0;
      }
      free(var_grp->variable[j]);
      var_grp->variable[j] = 0;
    }
  }

  file->idx->variable_count = 0;

  for (i = 0; i < 16; i++)
  {
    free(file->idx->variable_grp[i]);
    file->idx->variable_grp[i] = 0;
  }

  PIDX_dump_state_finalize(file);

  free(file->idx);
  free(file->idx_d->restructured_grid);
  free(file->idx_d->time);
  free(file->idx_d);
  free(file->idx_dbg);
  free(file->idx_c);
  free(file->idx_cache);

  free(file);

  return PIDX_success;
}


static int approx_maxh(PIDX_file file)
{
  int maxh = log2(getPowerOf2(file->idx->bounds[0])) + log2(getPowerOf2(file->idx->bounds[1])) + log2(getPowerOf2(file->idx->bounds[2])) + 1;

  //fprintf(stderr, "[%d %d %d] mh - bpb + bpf %d - %d + %d\n", (int)log2(getPowerOf2(file->idx->bounds[0])), file->idx->bounds[1], file->idx->bounds[2], maxh, file->idx->bits_per_block, + (int)log2(file->idx->blocks_per_file));
  int lc = maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file));
  if (lc < 1)
    return 1;
  else
    return maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file));
}


static void PIDX_debug_output(PIDX_file file, int gi, int svi, int evi, int io_type)
{
  int i = 0;
  pidx_global_variable++;

#if DETAIL_OUTPUT
  if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
  {
#if 0
    if (file->idx->io_type == PIDX_io_type::PIDX_RAW_IO)
      fprintf(stderr, "PIDX_RAW_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_IDX_IO)
      fprintf(stderr, "PIDX_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_LOCAL_PARTITION_IDX_IO)
      fprintf(stderr, "PIDX_LOCAL_PARTITION_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_io_type::PIDX_GLOBAL_PARTITION_IDX_IO)
      fprintf(stderr, "PIDX_GLOBAL_PARTITION_IDX_IO %s\n", file->idx->filename);
#endif

    if (file->idx->io_type != PIDX_RAW_IO)
    {
      fprintf(stderr, "[%d : %d %d] [%d %d %d : %d]\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, (int) file->idx->bounds[0], (int) file->idx->bounds[1], (int) file->idx->bounds[2], file->idx->variable_count);
      fprintf(stderr, "Box set by user (PIDX_USER_RST_BOX) %d %d %d\n", (int)file->idx_d->restructured_grid->patch_size[0], (int)file->idx_d->restructured_grid->patch_size[1], (int)file->idx_d->restructured_grid->patch_size[2]);

      fprintf(stderr, "Compression Bit rate set to %f\n", file->idx->compression_bit_rate);

      if (file->idx->endian == PIDX_LITTLE_ENDIAN)
        fprintf(stderr, "Little Endian | ");
      else if (file->idx->endian == PIDX_BIG_ENDIAN)
        fprintf(stderr, "Big Endian | ");

      if (file->idx->flip_endian == 1)
        fprintf(stderr, "Endian Flipping Done\n");
      if (file->idx->flip_endian == 0)
        fprintf(stderr, "Endian Flipping Not Done\n");

      fprintf(stderr, "Partition count %d = %d x %d x %d Partitio size = %d x %d x %d\n", file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2], file->idx_d->partition_count[0], file->idx_d->partition_count[1], file->idx_d->partition_count[2], file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2]);
      fprintf(stderr, "Comp = %d\n", file->idx->compression_type);
      fprintf(stderr, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);
      fprintf(stderr, "Partition level : maxh = %d : %d\n", file->idx_d->total_partiton_level, file->idx_d->maxh);
    }
  }
#endif

  PIDX_time time = file->idx_d->time;
  time->sim_end = MPI_Wtime();

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->idx_c->global_comm);

  if (io_type != PIDX_RAW_IO)
  {
    if (max_time == total_time)
    {
      int si = 0;
      double group_init = time->init_end - time->init_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;

      double pre_group_total = group_init + group_reg_box;
      double post_group_total = group_bitstring + group_block_layout + header_io;

#if DETAIL_OUTPUT
      fprintf(stderr, "\n[T %d R %d N %d G %d]\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, pidx_global_variable);
      fprintf(stderr, "PRE INIT          :[%.4f + %.4f] = %.4f\n\n", group_init, group_reg_box, pre_group_total);
      fprintf(stderr, "POST INIT         :[%.4f + %.4f + %.4f] = %.4f\n\n", group_bitstring, group_block_layout, header_io, post_group_total);
#endif

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0;
      double hz_init = 0, hz_meta_data = 0, hz_buffer = 0, hz = 0, hz_compress = 0, hz_buffer_free = 0, hz_cleanup = 0, hz_total = 0, hz_all = 0;
      double chunk_init = 0, chunk_meta = 0, chunk_buffer = 0, chunk = 0, chunk_buffer_free = 0, chunk_cleanup = 0, chunk_total = 0, chunk_all = 0;
      double compression_init = 0, compression = 0, compression_total = 0, compression_all = 0;
      double io = 0, io_all = 0;

      for (si = svi; si < evi; si++)
      {
        rst_init = time->rst_init_end[gi][si] - time->rst_init_start[gi][si];
        rst_meta_data_create = time->rst_meta_data_create_end[gi][si] - time->rst_meta_data_create_start[gi][si];
        rst_meta_data_io = time->rst_meta_data_io_end[gi][si] - time->rst_meta_data_io_start[gi][si];
        rst_buffer = time->rst_buffer_end[gi][si] - time->rst_buffer_start[gi][si];
        rst_write_read = time->rst_write_read_end[gi][si] - time->rst_write_read_start[gi][si];
        rst_buff_agg = time->rst_buff_agg_end[gi][si] - time->rst_buff_agg_start[gi][si];
        rst_buff_agg_free = time->rst_buff_agg_free_end[gi][si] - time->rst_buff_agg_free_start[gi][si];
        rst_buff_agg_io = time->rst_buff_agg_io_end[gi][si] - time->rst_buff_agg_io_start[gi][si];
        rst_cleanup = time->rst_cleanup_end[gi][si] - time->rst_cleanup_start[gi][si];
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
      double w_all = 0;

      for (si = svi; si < evi; si++)
      {
        if (file->idx->variable_grp[gi]->variable_tracker[si] == 1)
        {
#if 0
          rst_init = time->rst_init_end[gi][si] - time->rst_init_start[gi][si];
          rst_meta_data_create = time->rst_meta_data_create_end[gi][si] - time->rst_meta_data_create_start[gi][si];
          rst_meta_data_io = time->rst_meta_data_io_end[gi][si] - time->rst_meta_data_io_start[gi][si];
          rst_buffer = time->rst_buffer_end[gi][si] - time->rst_buffer_start[gi][si];
          rst_write_read = time->rst_write_read_end[gi][si] - time->rst_write_read_start[gi][si];
          rst_buff_agg = time->rst_buff_agg_end[gi][si] - time->rst_buff_agg_start[gi][si];
          rst_buff_agg_free = time->rst_buff_agg_free_end[gi][si] - time->rst_buff_agg_free_start[gi][si];
          rst_buff_agg_io = time->rst_buff_agg_io_end[gi][si] - time->rst_buff_agg_io_start[gi][si];
          rst_cleanup = time->rst_cleanup_end[gi][si] - time->rst_cleanup_start[gi][si];
          rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
          rst_all = rst_all + rst_total;
#endif


          hz_init = time->hz_init_end[gi][si] - time->hz_init_start[gi][si];
          hz_meta_data = time->hz_meta_end[gi][si] - time->hz_meta_start[gi][si];
          hz_buffer = time->hz_buffer_end[gi][si] - time->hz_buffer_start[gi][si];
          hz = time->hz_end[gi][si] - time->hz_start[gi][si];
          hz_compress = time->hz_compress_end[gi][si] - time->hz_compress_start[gi][si];
          hz_buffer_free = time->hz_buffer_free_end[gi][si] - time->hz_buffer_free_start[gi][si];
          hz_cleanup = time->hz_cleanup_end[gi][si] - time->hz_cleanup_start[gi][si];
          hz_total = hz_init + hz_meta_data + hz_buffer + hz + hz_compress + hz_buffer_free + hz_cleanup;
          hz_all = hz_all + hz_total;


          chunk_init = time->chunk_init_end[gi][si] - time->chunk_init_start[gi][si];
          chunk_meta = time->chunk_meta_end[gi][si] - time->chunk_meta_start[gi][si];
          chunk_buffer = time->chunk_buffer_end[gi][si] - time->chunk_buffer_start[gi][si];
          chunk = time->chunk_end[gi][si] - time->chunk_start[gi][si];
          chunk_buffer_free = time->chunk_buffer_free_end[gi][si] - time->chunk_buffer_free_start[gi][si];
          chunk_cleanup = time->chunk_cleanup_end[gi][si] - time->chunk_cleanup_start[gi][si];
          chunk_total = chunk_init + chunk_meta + chunk_buffer + chunk + chunk_buffer_free + chunk_cleanup;
          chunk_all = chunk_all + chunk_total;


          compression_init = time->compression_init_end[gi][si] - time->compression_init_start[gi][si];
          compression = time->compression_end[gi][si] - time->compression_start[gi][si];
          compression_total = compression_init + compression;
          compression_all = compression_all + compression_total;

#if DETAIL_OUTPUT
          //fprintf(stderr, "RST           :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
          fprintf(stderr, "HZ            :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, hz_init, hz_meta_data, hz_buffer, hz_compress, hz, hz_buffer_free, hz_cleanup, hz_total);
          fprintf(stderr, "CHUNK         :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, chunk_init, chunk_meta, chunk_buffer, chunk, chunk_buffer_free, chunk_cleanup, chunk_total);
          fprintf(stderr, "CMP           :[%d] [%.4f + %.4f] = %.4f\n", si, compression_init, compression, compression_total);
#endif

          double hz_io = 0;
          for (i = file->idx->variable_grp[gi]->agg_level; i < file->idx->variable_grp[gi]->shared_layout_count + file->idx->variable_grp[gi]->nshared_layout_count ; i++)
          {
            hz_io = time->hz_io_end[gi][si][i] - time->hz_io_start[gi][si][i];
            hz_io_all = hz_io_all + hz_io;

#if DETAIL_OUTPUT
            fprintf(stderr, "[HZ I/O S %d %d]  :[%d] [%d] %.4f [%.4f] \n", file->idx->variable_grp[gi]->agg_level, file->idx->variable_grp[gi]->shared_end_layout_index, si, i, hz_io, hz_io_all);
#endif
          }
        }

        double w_comm_s_x_odd = 0, w_comp_s_x_odd = 0, w_comm_s_x_even = 0, w_comp_s_x_even;
        double w_comm_s_y_odd = 0, w_comp_s_y_odd = 0, w_comm_s_y_even = 0, w_comp_s_y_even;
        double w_comm_s_z_odd = 0, w_comp_s_z_odd = 0, w_comm_s_z_even = 0, w_comp_s_z_even;
        double w_comm_total = 0, w_comp_total = 0;
        double w_comp_s_x = 0, w_comp_s_y = 0, w_comp_s_z = 0;
        for (i = 0; i < file->idx_d->wavelet_levels ; i++)
        {
          if (file->idx_d->wavelet_imeplementation_type == WAVELET_STENCIL)
          {
            w_comm_s_x_odd = time->w_stencil_comm_x_odd_end[gi][si][i] - time->w_stencil_comm_x_odd_start[gi][si][i];
            w_comp_s_x_odd = time->w_stencil_comp_x_odd_end[gi][si][i] - time->w_stencil_comp_x_odd_start[gi][si][i];
            w_comm_s_x_even = time->w_stencil_comm_x_even_end[gi][si][i] - time->w_stencil_comm_x_even_start[gi][si][i];
            w_comp_s_x_even = time->w_stencil_comp_x_even_end[gi][si][i] - time->w_stencil_comp_x_even_start[gi][si][i];

            w_comm_s_y_odd = time->w_stencil_comm_y_odd_end[gi][si][i] - time->w_stencil_comm_y_odd_start[gi][si][i];
            w_comp_s_y_odd = time->w_stencil_comp_y_odd_end[gi][si][i] - time->w_stencil_comp_y_odd_start[gi][si][i];
            w_comm_s_y_even = time->w_stencil_comm_y_even_end[gi][si][i] - time->w_stencil_comm_y_even_start[gi][si][i];
            w_comp_s_y_even = time->w_stencil_comp_y_even_end[gi][si][i] - time->w_stencil_comp_y_even_start[gi][si][i];

            w_comm_s_z_odd = time->w_stencil_comm_z_odd_end[gi][si][i] - time->w_stencil_comm_z_odd_start[gi][si][i];
            w_comp_s_z_odd = time->w_stencil_comp_z_odd_end[gi][si][i] - time->w_stencil_comp_z_odd_start[gi][si][i];
            w_comm_s_z_even = time->w_stencil_comm_z_even_end[gi][si][i] - time->w_stencil_comm_z_even_start[gi][si][i];
            w_comp_s_z_even = time->w_stencil_comp_z_even_end[gi][si][i] - time->w_stencil_comp_z_even_start[gi][si][i];

            w_comm_total = w_comm_s_x_odd + w_comm_s_x_even + w_comm_s_y_odd + w_comm_s_y_even + w_comm_s_z_odd + w_comm_s_z_even;
            w_comp_total = w_comp_s_x_odd + w_comp_s_x_even + w_comp_s_y_odd + w_comp_s_y_even + w_comp_s_z_odd + w_comp_s_z_even;
            w_all = w_all + w_comm_total + w_comp_total;

#if DETAIL_OUTPUT
            fprintf(stderr, "[WAVELET STENCIL COMM][%d]   :[%d] [%f + %f] + [%f + %f] + [%f + %f] = %f [%f]\n", i, si, w_comm_s_x_odd, w_comm_s_x_even, w_comm_s_y_odd, w_comm_s_y_even, w_comm_s_z_odd, w_comm_s_z_even, w_comm_total, w_all);
            fprintf(stderr, "[WAVELET STENCIL COMP][%d]   :[%d] [%f + %f] + [%f + %f] + [%f + %f] = %f [%f]\n", i, si, w_comp_s_x_odd, w_comp_s_x_even, w_comp_s_y_odd, w_comp_s_y_even, w_comp_s_z_odd, w_comp_s_z_even, w_comp_total, w_all);
#endif
          }
          else
          {
            w_comp_s_x = time->w_rst_comp_x_end[gi][si][i] - time->w_rst_comp_x_start[gi][si][i];
            w_comp_s_y = time->w_rst_comp_y_end[gi][si][i] - time->w_rst_comp_y_start[gi][si][i];
            w_comp_s_z = time->w_rst_comp_z_end[gi][si][i] - time->w_rst_comp_z_start[gi][si][i];

            w_comp_total = w_comp_s_x + w_comp_s_y + w_comp_s_z;
            w_all = w_all + w_comm_total + w_comp_total;

#if DETAIL_OUTPUT
            fprintf(stderr, "[WAVELET RST COMP][%d]       :[%d] [%f + %f + %f] = %f [%f]\n", i, si, w_comp_s_x, w_comp_s_y, w_comp_s_z, w_comp_total, w_all);
#endif

          }
        }

        double agg_init = 0, agg_meta = 0, agg_buf = 0, agg = 0, agg_meta_cleanup = 0, agg_total = 0, agg_cmp = 0;
        for (i = file->idx->variable_grp[gi]->shared_start_layout_index; i < file->idx->variable_grp[gi]->agg_level ; i++)
        {
          agg_init = time->agg_init_end[gi][si][i] - time->agg_init_start[gi][si][i];
          agg_meta = time->agg_meta_end[gi][si][i] - time->agg_meta_start[gi][si][i];
          agg_buf = time->agg_buf_end[gi][si][i] - time->agg_buf_start[gi][si][i];
          agg = time->agg_end[gi][si][i] - time->agg_start[gi][si][i];
          agg_cmp = time->agg_compress_end[gi][si][i] - time->agg_compress_start[gi][si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[gi][si][i] - time->agg_meta_cleanup_start[gi][si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG S %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx->variable_grp[gi]->shared_start_layout_index, file->idx->variable_grp[gi]->shared_end_layout_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }

#if 0
        for (i = file->idx->variable_grp[gi]->shared_start_layout_index; i < file->idx->variable_grp[gi]->agg_l_shared ; i++)
        {
          agg_init = time->agg_init_end[gi][si][i] - time->agg_init_start[gi][si][i];
          agg_meta = time->agg_meta_end[gi][si][i] - time->agg_meta_start[gi][si][i];
          agg_buf = time->agg_buf_end[gi][si][i] - time->agg_buf_start[gi][si][i];
          agg = time->agg_end[gi][si][i] - time->agg_start[gi][si][i];
          agg_cmp = time->agg_compress_end[gi][si][i] - time->agg_compress_start[gi][si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[gi][si][i] - time->agg_meta_cleanup_start[gi][si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG S %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx->variable_grp[gi]->shared_start_layout_index, file->idx->variable_grp[gi]->shared_end_layout_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }

        for (i = file->idx->variable_grp[gi]->nshared_start_layout_index; i < file->idx->variable_grp[gi]->agg_l_nshared ; i++)
        {
          agg_init = time->agg_init_end[gi][si][i] - time->agg_init_start[gi][si][i];
          agg_meta = time->agg_meta_end[gi][si][i] - time->agg_meta_start[gi][si][i];
          agg_buf = time->agg_buf_end[gi][si][i] - time->agg_buf_start[gi][si][i];
          agg = time->agg_end[gi][si][i] - time->agg_start[gi][si][i];
          agg_cmp = time->agg_compress_end[gi][si][i] - time->agg_compress_start[gi][si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[gi][si][i] - time->agg_meta_cleanup_start[gi][si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_cmp + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          //fprintf(stderr, "[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);
          //fprintf(stderr, "[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta + agg_buf, agg,  agg_meta_cleanup, agg_total);

#if DETAIL_OUTPUT
          fprintf(stderr, "[AGG N %d %d]   :[%d] [%d] %f + %f + %f + %f + %f + %f = %f [%f]\n", file->idx->variable_grp[gi]->nshared_start_layout_index, file->idx->variable_grp[gi]->nshared_end_layout_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_cmp, agg_meta_cleanup, agg_total, agg_all);
#endif
        }
#endif

        io = time->io_end[gi][si] - time->io_start[gi][si];
        io_all = io_all + io;
#if DETAIL_OUTPUT
        fprintf(stderr, "IO [%d]        :[%d] %.4f\n", file->idx->variable_grp[gi]->variable_count, si, io);
#endif
      }

      grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all + w_all + hz_all + hz_io_all + chunk_all + compression_all + agg_all + io_all;

#if DETAIL_OUTPUT
      fprintf(stderr, "XIRPIWCCHHAI      :[%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f = %.4f] + %.4f [%.4f %.4f]\n", pre_group_total, rst_all, partition_time, post_group_total, w_all, chunk_all, compression_all, hz_all, hz_io_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);
#else
      fprintf(stderr, "[%s %d %d (%d %d %d)] IRPIWCCHHAI      :[%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f = %.4f] + %.4f [%.4f %.4f]\n", file->idx->filename, file->idx->current_time_step, (evi - svi), (int)file->idx->bounds[0], (int)file->idx->bounds[1], (int)file->idx->bounds[2], pre_group_total, rst_all, partition_time, post_group_total, w_all, chunk_all, compression_all, hz_all, hz_io_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);
#endif
    }
  }
  else
  {
    if (max_time == total_time)
    {
      int si = 0;
      double group_init = time->init_end - time->init_start;
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;
      double group_total = group_init + group_bitstring + group_reg_box + group_block_layout + header_io;

      //fprintf(stderr, "[T %d R %d N %d G %d] INIT          :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, pidx_global_variable, group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      //fprintf(stderr, "INIT      :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0, r1 = 0, r2 = 0, r3 = 0;
      //double grp_rst_hz_chunk_agg_io = group_total;

      for (si = svi; si < evi; si++)
      {
        if (file->idx->variable_grp[gi]->variable_tracker[si] == 1)
        {
          rst_init = time->rst_init_end[gi][si] - time->rst_init_start[gi][si];
          rst_meta_data_create = time->rst_meta_data_create_end[gi][si] - time->rst_meta_data_create_start[gi][si];
          rst_meta_data_io = time->rst_meta_data_io_end[gi][si] - time->rst_meta_data_io_start[gi][si];
          rst_buffer = time->rst_buffer_end[gi][si] - time->rst_buffer_start[gi][si];
          r1 = r1 + rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer;

          rst_write_read = time->rst_write_read_end[gi][si] - time->rst_write_read_start[gi][si];
          rst_buff_agg = time->rst_buff_agg_end[gi][si] - time->rst_buff_agg_start[gi][si];
          rst_buff_agg_free = time->rst_buff_agg_free_end[gi][si] - time->rst_buff_agg_free_start[gi][si];
          r2 = r2 + rst_write_read + rst_buff_agg + rst_buff_agg_free;

          rst_buff_agg_io = time->rst_buff_agg_io_end[gi][si] - time->rst_buff_agg_io_start[gi][si];
          rst_cleanup = time->rst_cleanup_end[gi][si] - time->rst_cleanup_start[gi][si];
          r3 = r3 + rst_buff_agg_io + rst_cleanup;

          rst_total = rst_init + rst_meta_data_create + rst_meta_data_io + rst_buffer + rst_write_read + rst_buff_agg + rst_buff_agg_free + rst_buff_agg_io + rst_cleanup;
          rst_all = rst_all + rst_total;

          //fprintf(stderr, "RST       :[%d] [%.4f + [MD] %.4f + %.4f + [B] %.4f + [R] %.4f + [MC] %.4f + %.4f + [IO] %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
        }
      }

      //grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all;
      fprintf(stderr, "[RAW] [%s] [%d %d %d : %d %d %d] [%d] [T %d R %d N %d V %d] : [%.4f + %.4f (%.4f = %.4f + %.4f + %.4f) = %.4f] + %.4f [%.4f %.4f]\n", file->idx->filename, (int)file->idx->bounds[0], (int)file->idx->bounds[1], (int)file->idx->bounds[2], (int)file->idx_d->restructured_grid->patch_size[0], (int)file->idx_d->restructured_grid->patch_size[1], (int)file->idx_d->restructured_grid->patch_size[2],  pidx_global_variable, file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, (evi - svi),
             group_total,
             rst_all,
             r1 + r2 + r3,
             r1, r2, r3,
             (group_total + rst_all),
             (time->SX - time->sim_start),
             (group_total + rst_all) + (time->SX - time->sim_start), max_time);
    }
  }
}


static PIDX_return_code PIDX_dump_state_finalize (PIDX_file file)
{
  if (file->idx_dbg->state_dump == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
  {
    fclose(file->idx_dbg->mpi_dump_fp);
    fclose(file->idx_dbg->local_dump_fp);
  }

  return PIDX_success;
}
