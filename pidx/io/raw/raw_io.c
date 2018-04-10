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
#include "../../PIDX_inc.h"

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code dump_process_extent(PIDX_io file);

/// Raw Write Steps
/********************************************************
*  Step 0: Setup Group related meta-data                *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform data Restructuring                   *
*  Step 3: Perform actual file IO                       *
*  Step 4: cleanup for Steps 1                          *
*********************************************************/

PIDX_return_code PIDX_raw_write(PIDX_io file, int gi, int svi, int evi)
{  
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 0
  time->set_reg_box_start = MPI_Wtime();
  if (set_rst_box_size_for_raw_write(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->set_reg_box_end = MPI_Wtime();

  ret = group_meta_data_init(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = dump_process_extent(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx_d->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
  {
    ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    ret = raw_restructure_setup(file, gi, si, ei, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 2: Perform data restructuring
    ret = raw_restructure(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Write out restructured data
    ret = raw_restructure_io(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4: Cleanup all buffers and ids
    ret = raw_restructure_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  return PIDX_success;
}



/// Raw Read Steps
/********************************************************
*  Step 0: Setup Group related meta-data                *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform actual file IO                       *
*  Step 3: Perform data Restructuring                   *
*  Step 4: cleanup for Steps 1                          *
*********************************************************/

PIDX_return_code PIDX_raw_read(PIDX_io file, int gi, int svi, int evi)
{
  int si = 0, ei = 0;
  PIDX_return_code ret;

  file->idx_d->variable_pipe_length = file->idx->variable_count;

  for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
  {
    ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    ret = raw_restructure_forced_read(file, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  return PIDX_success;
}


static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code dump_process_extent(PIDX_io file)
{
  int i, j, k;
  for (i = 0; i < file->idx->variable_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    for (j = 0; j < var_grp->variable_count; j++)
    {
      for (k = 0; k < var_grp->variable[j]->sim_patch_count; k++)
      {
        if (file->idx_dbg->state_dump == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
        {
          fprintf(file->idx_dbg->local_dump_fp, "[%d] [%d] %d %d %d %d %d %d\n", j, k, (int)var_grp->variable[j]->sim_patch[k]->offset[0], (int)var_grp->variable[j]->sim_patch[k]->offset[1], (int)var_grp->variable[j]->sim_patch[k]->offset[2], (int)var_grp->variable[j]->sim_patch[k]->size[0], (int)var_grp->variable[j]->sim_patch[k]->size[1], (int)var_grp->variable[j]->sim_patch[k]->size[2]);
          fflush(file->idx_dbg->local_dump_fp);
        }
      }
      if (file->idx_dbg->state_dump == PIDX_META_DATA_DUMP_ONLY || file->idx_dbg->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        fprintf(file->idx_dbg->local_dump_fp, "\n");
        fflush(file->idx_dbg->local_dump_fp);
      }
    }
  }

  return PIDX_success;
}
