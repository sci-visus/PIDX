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

  file->io = PIDX_io_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg);
  if (file->io == NULL)
    return PIDX_err_flush;

  int block_layout_count = approx_maxh(file);

  PIDX_init_timming_buffers1(time, vgc, file->idx->variable_count, block_layout_count);
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
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
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


  for (i = file->local_group_index; i < file->local_group_index + file->local_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    /*
    for (j = var_grp->local_variable_index; j < var_grp->local_variable_index + var_grp->local_variable_count; j++)
    {
      for(p = 0; p < var_grp->variable[j]->sim_patch_count; p++)
      {
        free(var_grp->variable[j]->sim_patch[p]);
        var_grp->variable[j]->sim_patch[p] = 0;
      }
    }
    */
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

  free(file->idx);                  file->idx = 0;
  free(file->idx_d->time);          file->idx_d->time = 0;
  free(file->idx_d);                file->idx_d = 0;
  free(file->idx_dbg);              file->idx_dbg = 0;
  free(file->idx_c);                file->idx_c = 0;

  free(file);

  return PIDX_success;
}


static int approx_maxh(PIDX_file file)
{
  int maxh = log2(getPowerOf2(file->idx->bounds[0])) * log2(getPowerOf2(file->idx->bounds[1])) * log2(getPowerOf2(file->idx->bounds[2])) + 1;

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
  if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
  {
    if (file->idx->io_type == PIDX_RAW_IO)
      fprintf(stdout, "PIDX_RAW_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_IDX_IO)
      fprintf(stdout, "PIDX_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_LOCAL_PARTITION_IDX_IO)
      fprintf(stdout, "PIDX_LOCAL_PARTITION_IDX_IO %s\n", file->idx->filename);
    else if (file->idx->io_type == PIDX_GLOBAL_PARTITION_IDX_IO)
      fprintf(stdout, "PIDX_GLOBAL_PARTITION_IDX_IO %s\n", file->idx->filename);

    fprintf(stdout, "[%d : %d %d] [%d %d %d : %d]\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, (int) file->idx->bounds[0], (int) file->idx->bounds[1], (int) file->idx->bounds[2], file->idx->variable_count);

    if (file->idx->reg_box_set == PIDX_CLOSEST_POWER_TWO)
      fprintf(stdout, "Box set by user (PIDX_CLOSEST_POWER_TWO) %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);
    else if (file->idx->reg_box_set == PIDX_USER_RST_BOX)
      fprintf(stdout, "Box set by user (PIDX_USER_RST_BOX) %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);
    else if (file->idx->reg_box_set == PIDX_BOX_PER_PROCESS)
      fprintf(stdout, "Box set automatic for box per process case (PIDX_BOX_PER_PROCESS) %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);
    else if (file->idx->reg_box_set == PIDX_BOX_FROM_BITSTRING)
      fprintf(stdout, "Box set by bitstring (PIDX_BOX_FROM_BITSTRING) %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);

    if (file->idx->endian == 1)
      fprintf(stdout, "Little Endian | ");
    else if (file->idx->endian == 0)
      fprintf(stdout, "Big Endian | ");

    if (file->idx->flip_endian == 1)
      fprintf(stdout, "Endian Flipping Done\n");
    if (file->idx->flip_endian == 0)
      fprintf(stdout, "Endian Flipping Not Done\n");

    if (file->idx->io_type != PIDX_RAW_IO)
    {
      fprintf(stdout, "Partition count %d = %d x %d x %d Partitio size = %d x %d x %d\n", file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2], file->idx_d->partition_count[0], file->idx_d->partition_count[1], file->idx_d->partition_count[2], file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2]);
      fprintf(stdout, "Rst = %d Comp = %d\n", file->idx->enable_rst, file->idx->compression_type);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);
      fprintf(stdout, "Shared Block level : Partition level : maxh = %d : %d : %d\n", file->idx_d->shared_block_level, file->idx_d->total_partiton_level, file->idx_d->maxh);
    }
  }

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
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;
      double group_total = group_init + group_bitstring + group_reg_box + group_block_layout + header_io;
      printf("\n[T %d R %d N %d G %d] INIT          :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, pidx_global_variable, group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0;
      double hz_init = 0, hz_meta_data = 0, hz_buffer = 0, hz = 0, hz_buffer_free = 0, hz_cleanup = 0, hz_total = 0, hz_all = 0;
      double chunk_init = 0, chunk_meta = 0, chunk_buffer = 0, chunk = 0, chunk_buffer_free = 0, chunk_cleanup = 0, chunk_total = 0, chunk_all = 0;
      double compression_init = 0, compression = 0, compression_total = 0, compression_all = 0;
      double io = 0, io_all = 0;

      double grp_rst_hz_chunk_agg_io = group_total;
      double agg_all = 0;
      double hz_io_all = 0;
      for (si = svi; si < evi; si++)
      {
        if (file->idx->variable_grp[gi]->variable_tracker[si] == 1)
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


          hz_init = time->hz_init_end[gi][si] - time->hz_init_start[gi][si];
          hz_meta_data = time->hz_meta_end[gi][si] - time->hz_meta_start[gi][si];
          hz_buffer = time->hz_buffer_end[gi][si] - time->hz_buffer_start[gi][si];
          hz = time->hz_end[gi][si] - time->hz_start[gi][si];
          hz_buffer_free = time->hz_buffer_free_end[gi][si] - time->hz_buffer_free_start[gi][si];
          hz_cleanup = time->hz_cleanup_end[gi][si] - time->hz_cleanup_start[gi][si];
          hz_total = hz_init + hz_meta_data + hz_buffer + hz + hz_buffer_free + hz_cleanup;
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
          printf("RST           :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);
          printf("HZ            :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, hz_init, hz_meta_data, hz_buffer, hz, hz_buffer_free, hz_cleanup, hz_total);
          printf("CHUNK         :[%d] [%.4f + %.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", si, chunk_init, chunk_meta, chunk_buffer, chunk, chunk_buffer_free, chunk_cleanup, chunk_total);
          printf("CMP           :[%d] [%.4f + %.4f] = %.4f\n", si, compression_init, compression, compression_total);

          double hz_io = 0;
          for (i = file->idx->variable_grp[gi]->agg_l_shared; i < file->idx->variable_grp[gi]->shared_end_layout_index ; i++)
          {
            hz_io = time->hz_io_end[gi][si][i] - time->hz_io_start[gi][si][i];
            hz_io_all = hz_io_all + hz_io;

            printf("[HZ I/O S %d %d]  :[%d] [%d] %.4f [%.4f] \n", file->idx->variable_grp[gi]->agg_l_shared, file->idx->variable_grp[gi]->shared_end_layout_index, si, i, hz_io, hz_io_all);
          }

          for (i = file->idx->variable_grp[gi]->agg_l_nshared; i < file->idx->variable_grp[gi]->nshared_end_layout_index ; i++)
          {
            hz_io = time->hz_io_end[gi][si][i] - time->hz_io_start[gi][si][i];
            hz_io_all = hz_io_all + hz_io;

            printf("[HZ I/O N %d %d]  :[%d] [%d] %.4f [%.4f]\n", file->idx->variable_grp[gi]->agg_l_nshared, file->idx->variable_grp[gi]->nshared_end_layout_index, si, i, hz_io, hz_io_all);
          }

        }

        double agg_init = 0, agg_meta = 0, agg_buf = 0, agg = 0, agg_meta_cleanup = 0, agg_total = 0;
        for (i = file->idx->variable_grp[gi]->shared_start_layout_index; i < file->idx->variable_grp[gi]->agg_l_shared ; i++)
        {
          agg_init = time->agg_init_end[gi][si][i] - time->agg_init_start[gi][si][i];
          agg_meta = time->agg_meta_end[gi][si][i] - time->agg_meta_start[gi][si][i];
          agg_buf = time->agg_buf_end[gi][si][i] - time->agg_buf_start[gi][si][i];
          agg = time->agg_end[gi][si][i] - time->agg_start[gi][si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[gi][si][i] - time->agg_meta_cleanup_start[gi][si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          printf("[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);

          //printf("[AGG S %d %d]   :[%d] [%d] %f + %f + %f + %f + %f = %f [%f]\n", file->idx->variable_grp[gi]->shared_start_layout_index, file->idx->variable_grp[gi]->shared_end_layout_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total, agg_all);
        }

        for (i = file->idx->variable_grp[gi]->nshared_start_layout_index; i < file->idx->variable_grp[gi]->agg_l_nshared ; i++)
        {
          agg_init = time->agg_init_end[gi][si][i] - time->agg_init_start[gi][si][i];
          agg_meta = time->agg_meta_end[gi][si][i] - time->agg_meta_start[gi][si][i];
          agg_buf = time->agg_buf_end[gi][si][i] - time->agg_buf_start[gi][si][i];
          agg = time->agg_end[gi][si][i] - time->agg_start[gi][si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[gi][si][i] - time->agg_meta_cleanup_start[gi][si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          printf("[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total);
          //printf("[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + 0 = %f\n", si, i, agg_init + agg_meta + agg_buf, agg,  agg_meta_cleanup, agg_total);

          //printf("[AGG N %d %d]   :[%d] [%d] %.4f + %.4f + %.4f + %.4f + %.4f = %.4f [%.4f]\n", file->idx->variable_grp[gi]->nshared_start_layout_index, file->idx->variable_grp[gi]->nshared_end_layout_index, si, i, agg_init, agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total, agg_all);
        }

        io = time->io_end[gi][si] - time->io_start[gi][si];
        io_all = io_all + io;
        printf("IO [%d]        :[%d] %.4f\n", file->idx->variable_grp[gi]->variable_count, si, io);
      }

      grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all + hz_all + hz_io_all + chunk_all + compression_all + agg_all + io_all;

      printf("IRCCHHAI    :[%.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f + %.4f = %.4f] + %.4f [%.4f %.4f] \n", group_total, rst_all, chunk_all, compression_all, hz_all, hz_io_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);

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

      printf("[T %d R %d N %d G %d] INIT          :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", file->idx->current_time_step, file->idx_c->grank, file->idx_c->gnprocs, pidx_global_variable, group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      //printf("INIT      :[%.4f + %.4f + %.4f + %.4f + %.4f] = %.4f\n", group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0;
      double grp_rst_hz_chunk_agg_io = group_total;

      for (si = svi; si < evi; si++)
      {
        if (file->idx->variable_grp[gi]->variable_tracker[si] == 1)
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

          printf("RST       :[%d] [%.4f + [MD] %.4f + %.4f + [B] %.4f + [R] %.4f + [MC] %.4f + %.4f + [IO] %.4f + %.4f] = %.4f\n", si, rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);

        }
      }

      grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all;
      printf("TOT       :[%.4f + %.4f = %.4f] + %.4f [%.4f %.4f] \n\n", group_total, rst_all, grp_rst_hz_chunk_agg_io, (time->SX - time->sim_start), grp_rst_hz_chunk_agg_io + (time->SX - time->sim_start), max_time);
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
