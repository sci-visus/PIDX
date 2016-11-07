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

PIDX_return_code PIDX_flush(PIDX_file file)
{
  int i, j, p;
  int ret;
  if (file->idx->variable_count <= 0)
    return PIDX_err_variable;

  file->io = PIDX_io_init(file->idx, file->idx_d, file->idx_dbg);
  if (file->io == NULL)
    return PIDX_err_flush;

  ret = PIDX_io_set_communicator(file->io, file->comm);
  if (ret != PIDX_success)
    return PIDX_err_io;

  for (i = file->local_group_index; i < file->local_group_index + file->local_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    if (file->flags == MPI_MODE_CREATE)
      ret = PIDX_write(file->io, i, var_grp->local_variable_index, (var_grp->local_variable_index + var_grp->local_variable_count), file->idx->io_type);

    else if (file->flags == PIDX_MODE_RDONLY)
      ret = PIDX_read(file->io, i, var_grp->local_variable_index, (var_grp->local_variable_index + var_grp->local_variable_count), file->idx->io_type);

    if (ret != PIDX_success)
      return PIDX_err_io;
  }

  ret = PIDX_io_finalize(file->io);
  if (ret != PIDX_success)
    return PIDX_err_io;


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
  file->write_on_close = 1;

  ret = PIDX_flush(file);
  if (ret != PIDX_success)
    return PIDX_err_close;

  PIDX_time time = file->idx_d->time;
  time->sim_end = PIDX_get_time();

  int rank;
  int nprocs;
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
  if (rank == 0)
  {
    fprintf(stdout, "Time step %d File name %s\n", file->idx->current_time_step, file->idx->filename);
    fprintf(stdout, "Bitstring %s\n", file->idx->bitSequence);
    fprintf(stdout, "Global Data %lld %lld %lld Variables %d\n", (long long) file->idx->bounds[0], (long long) file->idx->bounds[1], (long long) file->idx->bounds[2], file->idx->variable_count);

    if (file->idx->reg_box_set == PIDX_USER_RST_BOX)
      fprintf(stdout, "Box set by user (PIDX_USER_RST_BOX)\n");
    else if (file->idx->reg_box_set == PIDX_BOX_PER_PROCESS)
      fprintf(stdout, "Box set automatic for box per process case (PIDX_BOX_PER_PROCESS)\n");
    else if (file->idx->reg_box_set == PIDX_BOX_FROM_BITSTRING)
      fprintf(stdout, "Box set by bitstring (PIDX_BOX_FROM_BITSTRING)\n");

    fprintf(stdout, "Restructuring Box Size %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);
    fprintf(stdout, "Partition count %d = %d x %d x %d\n", file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2], file->idx_d->partition_count[0], file->idx_d->partition_count[1], file->idx_d->partition_count[2]);
    fprintf(stdout, "Rst = %d Comp = %d\n", file->idx->enable_rst, file->idx->compression_type);
    fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);
    fprintf(stdout, "Shared Block level : Partition level : maxh = %d : %d : %d\n", file->idx_d->shared_block_level, file->idx_d->total_partiton_level, file->idx_d->maxh);
  }

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);

  if (max_time == total_time)
  {
    int si = 0, gi = 0;
    for (gi = 0; gi < file->idx->variable_group_count; gi++)
    {
      double group_init = time->init_end - time->init_start;
      double group_bitstring = time->bit_string_end - time->bit_string_start;
      double group_reg_box = time->set_reg_box_end - time->set_reg_box_start;
      double group_block_layout = time->layout_end - time->layout_start;
      double header_io = time->header_io_end - time->header_io_start;
      double group_total = group_init + group_bitstring + group_reg_box + group_block_layout + header_io;
      printf("GROUP INIT: [%f + %f + %f + %f + %f] = %f\n", group_init, group_bitstring, group_reg_box, group_block_layout, header_io, group_total);

      double rst_init = 0, rst_meta_data_create = 0, rst_meta_data_io = 0, rst_buffer = 0, rst_write_read = 0, rst_buff_agg = 0, rst_buff_agg_free = 0, rst_buff_agg_io = 0, rst_cleanup = 0, rst_total = 0, rst_all = 0;
      double hz_init = 0, hz_meta_data = 0, hz_buffer = 0, hz = 0, hz_buffer_free = 0, hz_cleanup = 0, hz_total = 0, hz_all = 0;
      double chunk_init = 0, chunk_meta = 0, chunk_buffer = 0, chunk = 0, chunk_buffer_free = 0, chunk_cleanup = 0, chunk_total = 0, chunk_all = 0;
      double compression_init = 0, compression = 0, compression_total = 0, compression_all = 0;
      double io = 0, io_all = 0;

      double grp_rst_hz_chunk_agg_io = group_total;
      double agg_all = 0;
      for (si = 0; si < file->idx->variable_grp[gi]->variable_count; si++)
      {
        if (file->idx->variable_grp[gi]->variable_tracker[si] == 1)
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
          printf("RST: [%f + %f + %f + %f + %f + %f + %f + %f + %f] = %f\n", rst_init, rst_meta_data_create, rst_meta_data_io, rst_buffer, rst_write_read, rst_buff_agg, rst_buff_agg_free, rst_buff_agg_io, rst_cleanup, rst_total);

          hz_init = time->hz_init_end[si] - time->hz_init_start[si];
          hz_meta_data = time->hz_meta_end[si] - time->hz_meta_start[si];
          hz_buffer = time->hz_buffer_end[si] - time->hz_buffer_start[si];
          hz = time->hz_end[si] - time->hz_start[si];
          hz_buffer_free = time->hz_buffer_free_end[si] - time->hz_buffer_free_start[si];
          hz_cleanup = time->hz_cleanup_end[si] - time->hz_cleanup_start[si];
          hz_total = hz_init + hz_meta_data + hz_buffer + hz + hz_buffer_free + hz_cleanup;
          hz_all = hz_all + hz_total;
          printf("HZ: [%f + %f + %f + %f + %f + %f] = %f\n", hz_init, hz_meta_data, hz_buffer, hz, hz_buffer_free, hz_cleanup, hz_total);

          chunk_init = time->chunk_init_end[si] - time->chunk_init_start[si];
          chunk_meta = time->chunk_meta_end[si] - time->chunk_meta_start[si];
          chunk_buffer = time->chunk_buffer_end[si] - time->chunk_buffer_start[si];
          chunk = time->chunk_end[si] - time->chunk_start[si];
          chunk_buffer_free = time->chunk_buffer_free_end[si] - time->chunk_buffer_free_start[si];
          chunk_cleanup = time->chunk_cleanup_end[si] - time->chunk_cleanup_start[si];
          chunk_total = chunk_init + chunk_meta + chunk_buffer + chunk + chunk_buffer_free + chunk_cleanup;
          chunk_all = chunk_all + chunk_total;
          printf("chunk: [%f + %f + %f + %f + %f + %f] = %f\n", chunk_init, chunk_meta, chunk_buffer, chunk, chunk_buffer_free, chunk_cleanup, chunk_total);

          compression_init = time->compression_init_end[si] - time->compression_init_start[si];
          compression = time->compression_end[si] - time->compression_start[si];
          compression_total = compression_init + compression;
          compression_all = compression_all + compression_total;
          printf("compression: [%f + %f] = %f\n", compression_init, compression, compression_total);
        }

        double agg_init = 0, agg_meta = 0, agg_buf = 0, agg = 0, agg_meta_cleanup = 0, agg_total;

        printf("SHARED %d %d\n", file->idx->variable_grp[gi]->shared_start_layout_index, file->idx->variable_grp[gi]->shared_end_layout_index);
        for (i = file->idx->variable_grp[gi]->shared_start_layout_index; i < file->idx->variable_grp[gi]->shared_end_layout_index ; i++)
        {
          agg_init = time->agg_init_end[si][i] - time->agg_init_start[si][i];
          agg_meta = time->agg_meta_end[si][i] - time->agg_meta_start[si][i];
          agg_buf = time->agg_buf_end[si][i] - time->agg_buf_start[si][i];
          agg = time->agg_end[si][i] - time->agg_start[si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[si][i] - time->agg_meta_cleanup_start[si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          printf("[S] [%d] [%d] %f + %f + %f + %f + %f = %f [%f]\n", si, i, agg_init, agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total, agg_all);
        }

        printf("NON SHARED %d %d\n", file->idx->variable_grp[gi]->nshared_start_layout_index, file->idx->variable_grp[gi]->nshared_end_layout_index);
        for (i = file->idx->variable_grp[gi]->nshared_start_layout_index; i < file->idx->variable_grp[gi]->nshared_end_layout_index ; i++)
        {
          agg_init = time->agg_init_end[si][i] - time->agg_init_start[si][i];
          agg_meta = time->agg_meta_end[si][i] - time->agg_meta_start[si][i];
          agg_buf = time->agg_buf_end[si][i] - time->agg_buf_start[si][i];
          agg = time->agg_end[si][i] - time->agg_start[si][i];
          agg_meta_cleanup = time->agg_meta_cleanup_end[si][i] - time->agg_meta_cleanup_start[si][i];
          agg_total = agg_init + agg_meta + agg_buf + agg + agg_meta_cleanup;
          agg_all = agg_all + agg_total;

          printf("[N] [%d] [%d] %f + %f + %f + %f + %f = %f [%f]\n", si, i, agg_init, agg_meta, agg_buf, agg, agg_meta_cleanup, agg_total, agg_all);
        }

        io = time->io_end[si] - time->io_start[si];
        io_all = io_all + io;
        printf("IO [%d]: %f\n", si, io);
      }

      grp_rst_hz_chunk_agg_io = grp_rst_hz_chunk_agg_io + rst_all + hz_all + chunk_all + compression_all + agg_all + io_all;

      printf("GROUP + RST + CHUNK + CMP + HZ + AGG + IO: [%f + %f + %f + %f + %f + %f + %f] = [%f %f] \n", group_total, rst_all, chunk_all, compression_all, hz_all, agg_all, io_all, grp_rst_hz_chunk_agg_io, max_time);
    }
  }

  PIDX_delete_timming_buffers2(file->idx_d->time, file->idx->variable_count);


  int j = 0;
  for (i = 0; i < file->idx->variable_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    for (j = 0; j < var_grp->variable_count; j++)
    {
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

  free(file->idx);                  file->idx = 0;
  free(file->idx_d->time);          file->idx_d->time = 0;
  free(file->idx_d);                file->idx_d = 0;
  free(file->idx_dbg);              file->idx_dbg = 0;

  free(file);

  return PIDX_success;
}
