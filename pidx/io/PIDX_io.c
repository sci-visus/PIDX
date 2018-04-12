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
/*
 * 1. One global dataset.
 * 2. Partitions in local index space.
 * 3. Partitions in global index space.
 * 4. Partition divided into non-shared, shared, file 0 (global index)
 * 5. Partition divided into non-shared, shared, file 0 (local index)
 */


#include "../PIDX_inc.h"



PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, idx_metadata_cache idx_cache)
{
  //Creating the restructuring ID
  PIDX_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_d = idx_derived_ptr;
  idx_io_id->idx_dbg = idx_dbg;
  idx_io_id->idx_c = idx_c;
  idx_io_id->idx_cache = idx_cache;

  return (idx_io_id);
}



PIDX_return_code PIDX_write(PIDX_io file, int gi, int svi, int evi, int MODE)
{
  PIDX_return_code ret = 0;

  file->idx_d->time->SX = PIDX_get_time();

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  if (var0->is_particle == 1)
    ret = PIDX_particle_write(file, gi, svi, evi);
  else
  {

#if 0
    if (file->idx_c->gnprocs == 1)
    {
      if (MODE == PIDX_IDX_IO || MODE == PIDX_LOCAL_PARTITION_IDX_IO || MODE == PIDX_GLOBAL_PARTITION_IDX_IO)
        ret = PIDX_serial_idx_write(file, gi, svi, evi);
      else if (MODE == PIDX_RAW_IO)
        ret = PIDX_raw_write(file, gi, svi, evi);
    }
    else
#endif
    {
      if (MODE == PIDX_IDX_IO)
        ret = PIDX_idx_write(file, gi, svi, evi);

      else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
        ret = PIDX_local_partition_idx_write(file, gi, svi, evi);

      else if (MODE == PIDX_RAW_IO)
        ret = PIDX_raw_write(file, gi, svi, evi);
    }
  }

  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  file->idx_d->time->EX = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code PIDX_read(PIDX_io file, int gi, int svi, int evi, int MODE)
{
  PIDX_return_code ret = 0;

  file->idx_d->time->SX = PIDX_get_time();

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  if (var0->is_particle == 1)
    ret = PIDX_particle_restart_read(file, gi, svi, evi);
  else
  {
    if (MODE == PIDX_IDX_IO)
      ret = PIDX_idx_read(file, gi, svi, evi);

    else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
      ret = PIDX_parallel_local_partition_idx_read(file, gi, svi, evi);

    else if (MODE == PIDX_RAW_IO)
      ret = PIDX_raw_read(file, gi, svi, evi);
  }


  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  file->idx_d->time->EX = PIDX_get_time();

  return PIDX_success;
}




PIDX_return_code PIDX_io_finalize(PIDX_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}

