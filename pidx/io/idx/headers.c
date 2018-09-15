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
#include "timming.h"


static PIDX_return_code write_idx_headers_layout(PIDX_io file, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl);

PIDX_return_code write_global_idx(PIDX_io file, int start_var_index, int end_var_index, int mode)
{
  // populate the filename template
  generate_file_name_template(file->idx->maxh, file->idx->bits_per_block, file->idx->filename, file->idx->current_time_step, file->idx->filename_template);

  if (mode != PIDX_WRITE)
    return PIDX_success;

  // if file io is turned on (debugging mode)
  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_c, file->idx_b, file->restructured_grid, file->fs_block_size, start_var_index, end_var_index);

    // write the .idx file
    if (PIDX_header_io_global_idx_write (file->header_io_id, file->idx->filename) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    PIDX_header_io_finalize(file->header_io_id);
  }

  return PIDX_success;
}


PIDX_return_code write_headers(PIDX_io file, int start_var_index, int end_var_index, int mode)
{
  // When using non-partitioned idx IO there is no need to create different IDX file per partitions
  // so we use one single .idx file.
  // Note: This strcpy is a little hacky, but works.
  // The .idx file will be written exactly as a single partition file
  if (file->idx->io_type == PIDX_IDX_IO)
    strcpy(file->idx->filename_partition,file->idx->filename);
  
  generate_file_name_template(file->idx->maxh, file->idx->bits_per_block, file->idx->filename_partition, file->idx->current_time_step, file->idx->filename_template_partition);

  //fprintf(stderr, "(maxh %d) FN %s FNT %s\n",file->idx->maxh, file->idx->filename_partition, file->idx->filename_template_partition);

  if (mode == PIDX_READ)
    return PIDX_success;

  if (write_idx_headers_layout(file, start_var_index, end_var_index, file->idx->filename_partition, file->idx->filename_template_partition, file->idx_b->block_layout) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



static PIDX_return_code write_idx_headers_layout(PIDX_io file, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl)
{
  if (file->idx_dbg->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_c, file->idx_b, file->restructured_grid, file->fs_block_size, start_var_index, end_var_index);

    if (PIDX_header_io_idx_file_create(file->header_io_id, bl, filename_template) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    /* STEP 2 */
    if (file->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      if (PIDX_header_io_idx_file_write(file->header_io_id, bl, filename_template,  0) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_header;
      }
    }

    if (file->variable_index_tracker == file->idx->variable_count)
    {
      if (PIDX_header_io_idx_file_write(file->header_io_id, bl, filename_template, 1) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_header;
      }
    }

    /* STEP 3 */
    if (PIDX_header_io_partition_idx_write(file->header_io_id, filename) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    if (PIDX_header_io_finalize(file->header_io_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }
  }

  return PIDX_success;
}



//
PIDX_return_code raw_headers_create_folder_structure(PIDX_io file, int start_var_index, int end_var_index, char* filename)
{
  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_c, file->idx_b, file->restructured_grid, file->fs_block_size, start_var_index, end_var_index);

    if (PIDX_header_io_raw_dir_create(file->header_io_id, filename) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    if (PIDX_header_io_finalize(file->header_io_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }
  }

  return PIDX_success;
}



PIDX_return_code raw_headers_create_idx_file(PIDX_io file, int start_var_index, int end_var_index, char* filename)
{
  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_c, file->idx_b, file->restructured_grid, file->fs_block_size, start_var_index, end_var_index);

    if (PIDX_header_io_raw_idx_write (file->header_io_id, filename) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    if (PIDX_header_io_finalize(file->header_io_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }
  }

  return PIDX_success;
}
