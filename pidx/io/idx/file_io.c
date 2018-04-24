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



PIDX_return_code file_io(PIDX_io file, int svi, int mode)
{
  int ret = 0;
  assert(file->idx_b->file0_agg_group_from_index == 0);
  PIDX_time time = file->time;

  time->io_start[svi] = PIDX_get_time();
  for (uint32_t j = file->idx_b->file0_agg_group_from_index; j < file->idx_b->agg_level; j++)
  {
    Agg_buffer temp_agg = file->idx->agg_buffer[svi][j];
    PIDX_block_layout temp_layout = file->idx_b->block_layout_by_agg_group[j];

    file->io_id[svi][j] = PIDX_file_io_init(file->idx, file->idx_c, file->fs_block_size, svi, svi);

    if (file->idx_dbg->debug_do_io == 1)
    {
      if (mode == PIDX_WRITE)
        ret = PIDX_file_io_blocking_write(file->io_id[svi][j], temp_agg, temp_layout, file->idx->filename_template_partition);
      else
        ret = PIDX_file_io_blocking_read(file->io_id[svi][j], temp_agg, temp_layout, file->idx->filename_template_partition);

      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
    PIDX_file_io_finalize(file->io_id[svi][j]);
  }
  time->io_end[svi] = PIDX_get_time();

  return PIDX_success;
}
