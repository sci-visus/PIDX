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


#include "../PIDX_inc.h"

PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_metadata_cache meta_data_cache, idx_blocks idx_b, PIDX_restructured_grid restructured_grid, PIDX_time time, int fs_block_size, int variable_index_tracker)
{
  //Creating the restructuring ID
  PIDX_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_dbg = idx_dbg;
  idx_io_id->idx_c = idx_c;
  idx_io_id->idx_b = idx_b;
  idx_io_id->meta_data_cache = meta_data_cache;
  idx_io_id->fs_block_size = fs_block_size;
  idx_io_id->variable_index_tracker = variable_index_tracker;
  idx_io_id->restructured_grid = restructured_grid;
  idx_io_id->time = time;

  return (idx_io_id);
}



PIDX_return_code PIDX_write(PIDX_io file, int svi, int evi, int MODE)
{
  // Four IO modes are supported
  // 1. idx io (no partitioning): does not create partitioning communicators
  // 2. local partitioned IDX io
  // 3. raw io: using restructuring phase
  // 4. particle io

  // current scheme of deciding between particle io and the others is not good
  // we need to come up with a better way of sperating the two

  PIDX_return_code ret = 0;
  file->time->SX = PIDX_get_time();

#if 0
  if (file->idx_c->simulation_nprocs == 1)
  {
    if (MODE == PIDX_IDX_IO || MODE == PIDX_LOCAL_PARTITION_IDX_IO)
      ret = PIDX_serial_idx_write(file, gi, svi, evi);
    else if (MODE == PIDX_RAW_IO)
      ret = PIDX_raw_write(file, gi, svi, evi);
  }
  else
#endif

  if (MODE == PIDX_IDX_IO)
    ret = PIDX_idx_write(file, svi, evi);

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
    ret = PIDX_local_partition_idx_write(file, svi, evi);

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_write(file, svi, evi);

  else if (MODE == PIDX_PARTICLE_IO)
    ret = PIDX_particle_file_per_process_write(file, svi, evi);

  else if (MODE == PIDX_RST_PARTICLE_IO)
    ret = PIDX_particle_rst_write(file, svi, evi);


  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->time->EX = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code PIDX_read(PIDX_io file, int svi, int evi, int MODE)
{
  // Four IO modes are supported
  // 1. idx io (no partitioning): does not create partitioning communicators
  // 2. local partitioned IDX io
  // 3. raw io: using restructuring phase
  // 4. particle io

  // current scheme of deciding between particle io and the others is not good
  // we need to come up with a better way of sperating the two

  PIDX_return_code ret = 0;
  file->time->SX = PIDX_get_time();


  if (MODE == PIDX_IDX_IO)
    ret = PIDX_idx_read(file, svi, evi);

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
  {
    ret = PIDX_local_partition_idx_generic_read(file, svi, evi);

    // Switch to this when you are sure that you are reading with the same number of
    // processes you used to write the data and also the per-process configuration is same
    //ret = PIDX_local_partition_idx_read(file, svi, evi);
  }

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_read(file, svi, evi);

  else if (MODE == PIDX_PARTICLE_IO)
    ret = PIDX_particle_vis_read(file, svi, evi);

  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->time->EX = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code PIDX_io_finalize(PIDX_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}
