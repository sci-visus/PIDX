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


PIDX_return_code PIDX_dump_rst_info(PIDX_file file, int dump_rst_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_d->dump_rst_info = dump_rst_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_d->rst_dump_dir_name, "%s_%d_rst_dump", filename_skeleton, file->idx->current_time_step);

  return PIDX_success;
}


PIDX_return_code PIDX_dump_agg_info(PIDX_file file, int dump_agg_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_d->dump_agg_info = dump_agg_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_d->agg_dump_dir_name, "%s_agg_dump", filename_skeleton);

  return PIDX_success;
}



PIDX_return_code PIDX_dump_io_info(PIDX_file file, int dump_io_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_d->dump_io_info = dump_io_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_d->io_dump_dir_name, "%s_io_dump", filename_skeleton);

  return PIDX_success;
}
