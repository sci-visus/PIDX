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
static char mpi_state_dump_dir_name[512];
static char local_state_dump_dir_name[512];
static PIDX_return_code dump_debug_data_init(PIDX_file file);

PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_chunk = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_compression(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_compress = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_hz(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_hz = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_agg = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_io(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_io = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_rst = debug_rst;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_hz = debug_hz;

  return PIDX_success;
}


PIDX_return_code PIDX_disable_rst(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_agg = 0;

  return PIDX_success;
}




PIDX_return_code PIDX_dump_state(PIDX_file file, int process_state)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->state_dump = process_state;

  char filename_skeleton[512];
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';

  sprintf(mpi_state_dump_dir_name, "%s_%d_global_state_dump", filename_skeleton, file->idx->current_time_step);
  sprintf(local_state_dump_dir_name, "%s_%d_local_state_dump", filename_skeleton, file->idx->current_time_step);

  if (dump_debug_data_init(file) != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d] PIDX debug failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  return PIDX_success;
}



static PIDX_return_code dump_debug_data_init(PIDX_file file)
{
  if (file->idx_dbg->state_dump == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
  {
    char mpi_file_name[1024];
    char local_file_name[1024];

    if (file->idx_c->grank == 0)
    {
#if 1
      char mkdir_line1[1024];
      char mkdir_line2[1024];
      sprintf(mkdir_line1, "mkdir -p %s", mpi_state_dump_dir_name);
      sprintf(mkdir_line2, "mkdir -p %s", local_state_dump_dir_name);
      system(mkdir_line1);
      system(mkdir_line2);
#else
      int ret;
      ret = mkdir(mpi_state_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
      if (ret != 0 && errno != EEXIST)
      {
        perror("mkdir");
        fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, mpi_state_dump_dir_name);
        return PIDX_err_io;
      }

      ret = mkdir(local_state_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
      if (ret != 0 && errno != EEXIST)
      {
        perror("mkdir");
        fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, local_state_dump_dir_name);
        return PIDX_err_io;
      }
#endif
    }

    MPI_Barrier(file->idx_c->global_comm);

    sprintf(mpi_file_name, "%s/rank_%d", mpi_state_dump_dir_name, file->idx_c->grank);
    file->idx_dbg->mpi_dump_fp = fopen(mpi_file_name, "a+");
    if (!file->idx_dbg->mpi_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, mpi_file_name);
      return PIDX_err_io;
    }

    sprintf(local_file_name, "%s/rank_%d", local_state_dump_dir_name, file->idx_c->grank);
    file->idx_dbg->local_dump_fp = fopen(local_file_name, "a+");
    if (!file->idx_dbg->local_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, local_file_name);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}
