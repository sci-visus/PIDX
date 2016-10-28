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

#if 0
  int rank;
  int nprocs;
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
  if (rank == 0)
  {
    fprintf(stdout, "Time step %d File name %s\n", file->idx->current_time_step, file->idx->filename);
    fprintf(stdout, "Bitstring %s\n", file->idx->bitSequence);
    fprintf(stdout, "Global Data %lld %lld %lld Variables %d\n", (long long) file->idx->bounds[0], (long long) file->idx->bounds[1], (long long) file->idx->bounds[2], file->idx->variable_count);

    fprintf(stdout, "Partition count %d = %d x %d x %d\n", file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2], file->idx_d->partition_count[0], file->idx_d->partition_count[1], file->idx_d->partition_count[2]);
    fprintf(stdout, "Rst = %d Comp = %d\n", file->idx->enable_rst, file->idx->compression_type);
    fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);
    fprintf(stdout, "Restructuring Box Size %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2]);
    fprintf(stdout, "Shared Block level : Partition level : maxh = %d : %d : %d\n", file->idx_d->shared_block_level, file->idx_d->total_partiton_level, file->idx_d->maxh);
  }
#endif

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  int rank = 0, nprocs = 1;

  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);

  if (max_time == total_time)
  {
    if (file->idx->io_type == PIDX_IDX_IO)
    {

    }
    else if (file->idx->io_type == PIDX_RAW_IO)
    {
      if (file->flags == PIDX_MODE_CREATE)
      {
        float component_time = (time->SX - time->sim_start) + (time->raw_write_rst_init_end  - time->raw_write_rst_init_start) + (time->raw_write_rst_end - time->raw_write_rst_start) + (time->raw_write_header_end - time->raw_write_header_start) + (time->raw_write_io_end - time->raw_write_io_start) + (time->raw_write_buffer_cleanup_end - time->raw_write_buffer_cleanup_start);

        printf("[%d (%d %d %d) %d] [%d %d] [%f %f] [B %f] + [M %f + R %f + F %f + I %f + C %f]\n", file->idx->current_time_step, (int) file->idx->bounds[0], (int) file->idx->bounds[1], (int) file->idx->bounds[2], file->idx->variable_count, rank, nprocs, max_time, component_time, (time->SX - time->sim_start), (time->raw_write_rst_init_end  - time->raw_write_rst_init_start), (time->raw_write_rst_end - time->raw_write_rst_start), (time->raw_write_header_end - time->raw_write_header_start), (time->raw_write_io_end - time->raw_write_io_start), (time->raw_write_buffer_cleanup_end - time->raw_write_buffer_cleanup_start));
      }
      else if (file->flags == PIDX_MODE_RDONLY)
      {

      }
    }
  }

#if 0
  //if (file->idx->io_type == PIDX_IDX_IO)
  {
    double total_time = time->sim_end - time->sim_start;
    double max_time = total_time;
    int rank = 0, nprocs = 1;

    MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm, &nprocs);

    if (max_time == total_time)
    {
      double total_computed_time = (time->file_create_time - time->sim_start) + (time->idx_init_end - time->idx_init_start) + (time->idx_rst_end - time->idx_rst_start) + (time->idx_partition_end - time->idx_partition_start) + (time->idx_bit_string_end - time->idx_bit_string_start) + (time->idx_hz_end - time->idx_hz_start) + (time->idx_comm_create_end - time->idx_comm_create_start) + (time->idx_layout_end - time->idx_layout_start) + (time->header_write_end - time->header_write_start) + (time->agg_buffer_end - time->agg_buffer_start) + (time->idx_agg_end - time->idx_agg_start) + (time->idx_io_end - time->idx_io_start) + (time->buffer_cleanup_end - time->buffer_cleanup_start);


      printf("[%d %d] Total Time %f = %f\n", rank, nprocs, max_time, total_computed_time);
      printf("[%f + %f + %f + %f + %f + %f + %f + %f + %f + %f + %f + %f + %f]\n", (time->file_create_time - time->sim_start), (time->idx_init_end - time->idx_init_start), (time->idx_rst_end - time->idx_rst_start), (time->idx_partition_end - time->idx_partition_start), (time->idx_bit_string_end - time->idx_bit_string_start), (time->idx_hz_end - time->idx_hz_start), (time->idx_comm_create_end - time->idx_comm_create_start), (time->idx_layout_end - time->idx_layout_start), (time->header_write_end - time->header_write_start), (time->agg_buffer_end - time->agg_buffer_start), (time->idx_agg_end - time->idx_agg_start), (time->idx_io_end - time->idx_io_start) ,(time->buffer_cleanup_end - time->buffer_cleanup_start) );

      /*
      //fprintf(stdout, "[P %d %d] Time Taken: %f Seconds [%f]\n", rank, nprocs, max_time, time->EX - time->SX);
      fprintf(stdout, "--------------------------------------------------------------------------------------------------------------------------\n");
      fprintf(stdout, "Init Time: %f Seconds\n", (time->file_create_time - time->sim_start));
      fprintf(stdout, "Partition time %f\n", time->partition_end_time - time->partition_start_time);
      fprintf(stdout, "[F Zero] Block layout creation time %f\n", time->populate_idx_end_time_f0 - time->populate_idx_start_time_f0);
      fprintf(stdout, "[Shared] Block layout creation time %f\n", time->populate_idx_end_time_s - time->populate_idx_start_time_s);
      fprintf(stdout, "[Non-Shared] Block layout creation time %f\n", time->populate_idx_end_time_ns - time->populate_idx_start_time_ns);

      double header_io_time = 0;
      for (var = 0; var < time->header_counter; var++)
      {
        header_io_time = header_io_time + (time->write_init_end[var] - time->write_init_start[var]);
        fprintf(stdout, "File Create time (+ header IO) %f\n", (time->write_init_end[var] - time->write_init_start[var]));
      }
      double stotal_time_ai = 0, stotal_time_bc = 0, stotal_time_a = 0, stotal_time_i = 0, stotal_time_pi = 0, stotal_time_m = 0;
      double ntotal_time_ai = 0, ntotal_time_bc = 0, ntotal_time_a = 0, ntotal_time_i = 0, ntotal_time_pi = 0, ntotal_time_m = 0;

      int p = 0;
      int g = 0;
      for (g = 0; g < file->idx->variable_group_count; g++)
      {
        PIDX_variable_group var_grp = file->idx->variable_grp[g];
        for (var = 0; var < file->idx->variable_count; var++)
        {
          PIDX_variable var1 = var_grp->variable[var];
          for (p = var_grp->shared_start_layout_index; p < var_grp->shared_end_layout_index; p++)
          {
            fprintf(stdout, "[S] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + %f = %f\n", var, p,
                   (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]),
                   (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]),
                   (time->agg_end[var][p] - time->agg_start[var][p]),
                   (time->io_end[var][p] - time->io_start[var][p]),
                   (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]),
                   (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]) + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]) + (time->agg_end[var][p] - time->agg_start[var][p]) + (time->io_end[var][p] - time->io_start[var][p]) + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]));

            stotal_time_bc = stotal_time_bc + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]);
            stotal_time_m = stotal_time_m + (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]);
            stotal_time_a = stotal_time_a + (time->agg_end[var][p] - time->agg_start[var][p]);
            stotal_time_i = stotal_time_i + (time->io_end[var][p] - time->io_start[var][p]);
            stotal_time_pi = stotal_time_pi + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]);
          }

          for (p = var_grp->nshared_start_layout_index; p < var_grp->nshared_end_layout_index; p++)
          {
            fprintf(stdout, "[N] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + %f = %f\n", var, p,
                   (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]),
                   (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]),
                   (time->agg_end[var][p] - time->agg_start[var][p]),
                   (time->io_end[var][p] - time->io_start[var][p]),
                   (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]),
                   (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]) + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]) + (time->agg_end[var][p] - time->agg_start[var][p]) + (time->io_end[var][p] - time->io_start[var][p]) + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]));

            ntotal_time_bc = ntotal_time_bc + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]);
            ntotal_time_m = ntotal_time_m + (time->agg_meta_end[var][p] - time->agg_meta_start[var][p]);
            ntotal_time_a = ntotal_time_a + (time->agg_end[var][p] - time->agg_start[var][p]);
            ntotal_time_i = ntotal_time_i + (time->io_end[var][p] - time->io_start[var][p]);
            ntotal_time_pi = ntotal_time_pi + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]);
          }
        }

        ntotal_time_ai = ntotal_time_m + ntotal_time_bc + ntotal_time_a + ntotal_time_i + ntotal_time_pi;
        stotal_time_ai = stotal_time_m + stotal_time_bc + stotal_time_a + stotal_time_i + stotal_time_pi;

        fprintf(stdout, "[ST] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + %f = %f\n", file->idx->variable_count, (var_grp->shared_end_layout_index - var_grp->shared_start_layout_index), stotal_time_m, stotal_time_bc, stotal_time_a, stotal_time_i, stotal_time_pi, stotal_time_ai);

        fprintf(stdout, "[NT] [%d %d] Agg meta + Agg Buf + Agg + AGG I/O + Per-Process I/O = %f + %f + %f + %f + %f = %f\n", file->idx->variable_count, (var_grp->nshared_end_layout_index - var_grp->nshared_start_layout_index), ntotal_time_m, ntotal_time_bc, ntotal_time_a, ntotal_time_i, ntotal_time_pi, ntotal_time_ai);

        fprintf(stdout, "HZ Time = %f\n", (time->hz_e_time - time->hz_s_time));
        fprintf(stdout, "Cleanup Time = %f\n", (time->buffer_cleanup_end - time->buffer_cleanup_start));

        fprintf(stdout, "PIDX Total Time = %f [%f + %f + (%f + %f + %f) + %f + %f + %f + %f + %f] [%f]\n", (time->file_create_time - time->sim_start) + (time->partition_end_time - time->partition_start_time) + (time->populate_idx_end_time_f0 - time->populate_idx_start_time_f0) + (time->populate_idx_end_time_s - time->populate_idx_start_time_s) + (time->populate_idx_end_time_ns - time->populate_idx_start_time_ns) + (time->hz_e_time - time->hz_s_time) + header_io_time + stotal_time_ai + ntotal_time_ai + (time->buffer_cleanup_end - time->buffer_cleanup_start),
                (time->file_create_time - time->sim_start),
                (time->partition_end_time - time->partition_start_time),
                (time->populate_idx_end_time_f0 - time->populate_idx_start_time_f0),
                (time->populate_idx_end_time_s - time->populate_idx_start_time_s),
                (time->populate_idx_end_time_ns - time->populate_idx_start_time_ns),
                (time->hz_e_time - time->hz_s_time),
                header_io_time,
                stotal_time_ai,
                ntotal_time_ai,
                (time->buffer_cleanup_end - time->buffer_cleanup_start),
                max_time);
        fprintf(stdout, "==========================================================================================================================\n");

      }
      */
    }
  }
#endif


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
