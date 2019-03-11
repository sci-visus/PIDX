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

static int cvi = 0;

// Initialiazation and creation of buffers for restructuring phase
PIDX_return_code raw_restructure_setup(PIDX_io file, int svi, int evi, int mode)
{
  PIDX_time time = file->time;
  cvi = svi;

  // Initialize the restructuring phase
  time->rst_init_start[cvi] = PIDX_get_time();
  file->raw_rst_id = PIDX_raw_rst_init(file->idx, file->idx_c, file->idx_dbg, file->restructured_grid, svi, evi);
  time->rst_init_end[cvi] = PIDX_get_time();


  // Populates the relevant meta-data
  time->rst_meta_data_create_start[cvi] = PIDX_get_time();
  if (PIDX_raw_rst_meta_data_create(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->rst_meta_data_create_end[cvi] = PIDX_get_time();


  // Saving the metadata info needed for reading back the data.
  // Especially when number of cores is different from number of cores
  // used to create the dataset
  time->rst_meta_data_io_start[cvi] = PIDX_get_time();
  //if (file->idx->cached_ts == file->idx->current_time_step)
  //if (file->idx->current_time_step == 0)
  //{
    if (mode == PIDX_WRITE)
    {
      if (PIDX_raw_rst_meta_data_write(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
    }
  //}
  time->rst_meta_data_io_end[cvi] = PIDX_get_time();


  // Creating the buffers required for restructurig
  time->rst_buffer_start[cvi] = PIDX_get_time();
  if (PIDX_raw_rst_buf_create(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }


  // Aggregating the aligned small buffers after restructuring into one single buffer
  if (PIDX_raw_rst_aggregate_buf_create(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->rst_buffer_end[cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code raw_restructure(PIDX_io file, int mode)
{
  PIDX_time time = file->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Perform data restructuring
      time->rst_write_read_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_staged_write(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_write_read_end[cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        if (HELPER_raw_rst(file->raw_rst_id) != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }

      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_aggregate(file->raw_rst_id, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_end[cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_destroy(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
  }

  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_aggregate(file->raw_rst_id, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_end[cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        if (HELPER_raw_rst(file->raw_rst_id) != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }

      // Perform data restructuring
      time->rst_write_read_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_read(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_write_read_end[cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_destroy(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code raw_restructure_io(PIDX_io file, int mode)
{
  PIDX_time time = file->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Write out restructured data
      time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_aggregated_write(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
    }
  }
  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Read restructured data
      time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
      if (PIDX_raw_rst_buf_aggregated_read(file->raw_rst_id) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code raw_restructure_cleanup(PIDX_io file)
{
  PIDX_time time = file->time;

  // Destroy buffers allocated during restructuring phase
  time->rst_cleanup_start[cvi] = PIDX_get_time();
  if (PIDX_raw_rst_aggregate_buf_destroy(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  if (PIDX_raw_rst_meta_data_destroy(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  // Deleting the restructuring ID
  PIDX_raw_rst_finalize(file->raw_rst_id);
  time->rst_cleanup_end[cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code raw_restructure_forced_read(PIDX_io file, int svi, int evi)
{
  file->raw_rst_id = PIDX_raw_rst_init(file->idx, file->idx_c, file->idx_dbg, file->restructured_grid, svi, evi);

  if (PIDX_raw_rst_forced_raw_read(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  if (PIDX_raw_rst_finalize(file->raw_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}




PIDX_return_code set_rst_box_size_for_raw_write(PIDX_io file, int svi)
{
  PIDX_time time = file->time;
  time->set_reg_box_start = PIDX_get_time();

  if ((file->restructured_grid->patch_size[0] == -1 || file->restructured_grid->patch_size[1] == -1 || file->restructured_grid->patch_size[2] == -1)
        || (file->restructured_grid->patch_size[0] == 0 &&  file->restructured_grid->patch_size[1] == 0 && file->restructured_grid->patch_size[2] == 0))
  {
    fprintf(stderr,"Warning: restructuring box is not set, using default 32x32x32 size. [File %s Line %d]\n", __FILE__, __LINE__);
    file->restructured_grid->patch_size[0]=32;file->restructured_grid->patch_size[1]=32;file->restructured_grid->patch_size[2]=32;
  }

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}
