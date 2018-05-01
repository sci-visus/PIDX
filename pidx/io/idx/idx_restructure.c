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

#define ACTUAL_IO 1
#include "../../PIDX_inc.h"

static int cvi = 0;
static PIDX_return_code free_idx_rst_box(PIDX_io file);

// Initialiazation of metadata and creation of buffers for restructuring phase
PIDX_return_code idx_restructure_setup(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;
  cvi = svi;

  // Initialize the restructuring phase
  time->rst_init_start[cvi] = PIDX_get_time();
  file->idx_rst_id = PIDX_idx_rst_init(file->idx, file->idx_c, file->idx_dbg, file->restructured_grid, svi, evi);
  time->rst_init_end[cvi] = PIDX_get_time();


  // Populates the relevant meta-data
  time->rst_meta_data_create_start[cvi] = PIDX_get_time();
  if (PIDX_idx_rst_meta_data_create(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->rst_meta_data_create_end[cvi] = PIDX_get_time();


  // Creating the buffers required for restructurig
  time->rst_buffer_start[cvi] = PIDX_get_time();
  if (PIDX_idx_rst_buf_create(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }


  // Aggregating the aligned small buffers after restructuring into one single buffer
  if (PIDX_idx_rst_aggregate_buf_create(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->rst_buffer_end[cvi] = PIDX_get_time();


  return PIDX_success;
}



PIDX_return_code idx_restructure(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Perform data restructuring
      time->rst_write_read_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_staged_write(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        ret = HELPER_idx_rst(file->idx_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }

      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregate(file->idx_rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_destroy(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
  }

  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregate(file->idx_rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        ret = HELPER_idx_rst(file->idx_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }

      // Perform data restructuring
      time->rst_write_read_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_read(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_destroy(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code idx_restructure_io(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Write out restructured data
      time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregated_write(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
    }
  }
  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Read restructured data
      time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregated_read(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code idx_restructure_cleanup(PIDX_io file)
{
  PIDX_time time = file->time;

  // Destroy buffers allocated during restructuring phase
  time->rst_cleanup_start[cvi] = PIDX_get_time();
  if (PIDX_idx_rst_aggregate_buf_destroy(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  if (PIDX_idx_rst_meta_data_destroy(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  // Deleting the restructuring ID
  PIDX_idx_rst_finalize(file->idx_rst_id);
  free_idx_rst_box(file);


  time->rst_cleanup_end[cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code idx_restructure_forced_read(PIDX_io file, int svi, int evi)
{
  file->idx_rst_id = PIDX_idx_rst_init(file->idx, file->idx_c, file->idx_dbg, file->restructured_grid, svi, evi);

  if (PIDX_idx_rst_forced_raw_read(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  if (PIDX_idx_rst_finalize(file->idx_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}



PIDX_return_code idx_restructure_rst_comm_create(PIDX_io file, int svi)
{
  PIDX_variable var0 = file->idx->variable[svi];

  MPI_Comm_rank(file->idx_c->simulation_comm, &(file->idx_c->simulation_rank));

  // the processes that are holding the super patch are grouped into the communicator rst_comm
  MPI_Comm_split(file->idx_c->simulation_comm, var0->restructured_super_patch_count, file->idx_c->simulation_rank, &(file->idx_c->rst_comm));
  MPI_Comm_rank(file->idx_c->rst_comm, &(file->idx_c->rrank));
  MPI_Comm_size(file->idx_c->rst_comm, &(file->idx_c->rnprocs));

  // copying the restructuring comm to partition comm
  // this is just inititalizing the partition communicator, later on depending if partitioning is done or not
  // the restructuring comm can be split to partition comms as well
  file->idx_c->partition_comm = file->idx_c->rst_comm;
  MPI_Comm_rank(file->idx_c->partition_comm, &(file->idx_c->partition_rank));
  MPI_Comm_size(file->idx_c->partition_comm, &(file->idx_c->partition_nprocs));

  return PIDX_success;
}



PIDX_return_code free_restructured_communicators(PIDX_io file)
{
  if(MPI_Comm_free(&(file->idx_c->rst_comm)) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}


static PIDX_return_code free_idx_rst_box(PIDX_io file)
{
  uint64_t *rgp = file->restructured_grid->total_patch_count;
  uint64_t total_patch_count = rgp[0] * rgp[1] * rgp[2];

  for (uint64_t i = 0; i < total_patch_count; i++)
    free(file->restructured_grid->patch[i]);

  free(file->restructured_grid->patch);

  return PIDX_success;
}
