/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
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
static char mpi_debug_file_output_state_dir_name[512];
static char local_debug_file_output_state_dir_name[512];
static PIDX_return_code dump_debug_data_init(PIDX_file file);

PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_chunk = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_compression(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_compress = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_hz(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_hz = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_agg(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_agg = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_io(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_io = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_rst = debug_rst;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_hz = debug_hz;

  return PIDX_success;
}




PIDX_return_code PIDX_disable_agg(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->enable_agg = 0;

  return PIDX_success;
}




PIDX_return_code PIDX_dump_state(PIDX_file file, int process_state)
{
  if (!file)
    return PIDX_err_file;

  file->idx_dbg->debug_file_output_state = process_state;

  char filename_skeleton[512];
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';

  sprintf(mpi_debug_file_output_state_dir_name, "%s_%d_global_debug_file_output_state", filename_skeleton, file->idx->current_time_step);
  sprintf(local_debug_file_output_state_dir_name, "%s_%d_local_debug_file_output_state", filename_skeleton, file->idx->current_time_step);

  if (dump_debug_data_init(file) != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d] PIDX debug failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  return PIDX_success;
}



static PIDX_return_code dump_debug_data_init(PIDX_file file)
{
  if (file->idx_dbg->debug_file_output_state == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->debug_file_output_state == PIDX_NO_IO_AND_META_DATA_DUMP)
  {
    char local_file_name[1024];

    if (file->idx_c->simulation_rank == 0)
    {
      char mkdir_line2[1024];
      sprintf(mkdir_line2, "mkdir -p %s", local_debug_file_output_state_dir_name);
      system(mkdir_line2);
    }

    MPI_Barrier(file->idx_c->simulation_comm);

    sprintf(local_file_name, "%s/rank_%d", local_debug_file_output_state_dir_name, file->idx_c->simulation_rank);
    file->idx_dbg->debug_file_output_fp = fopen(local_file_name, "a+");
    if (!file->idx_dbg->debug_file_output_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, local_file_name);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}
