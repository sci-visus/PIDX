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
static int levi = 0;

static PIDX_return_code hz_init(PIDX_io file, int svi, int evi);
static PIDX_return_code chunk_init(PIDX_io file, int svi, int evi);
static PIDX_return_code compression_init(PIDX_io file, int svi, int evi);
static PIDX_return_code meta_data_create(PIDX_io file);
static PIDX_return_code buffer_create(PIDX_io file);
static PIDX_return_code compress_and_encode(PIDX_io file);
static PIDX_return_code encode_and_uncompress(PIDX_io file);
static PIDX_return_code hz_cleanup(PIDX_io file);
static PIDX_return_code chunk_cleanup(PIDX_io file);



PIDX_return_code hz_encode_setup(PIDX_io file, int svi, int evi)
{
  cvi = svi;
  levi = evi;

  // Init
  if ( hz_init(file, svi, evi) || chunk_init(file, svi, evi) || compression_init(file, svi, evi) != PIDX_success )
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  // Meta data
  if (meta_data_create(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  // Buffer create
  if (buffer_create(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}



PIDX_return_code hz_encode(PIDX_io file, int mode)
{
  if (mode == PIDX_READ)
  {
    // encode, decompress, unchunk
    if (encode_and_uncompress(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
  }
  if (mode == PIDX_WRITE)
  {
    // chunk, compress and encode
    if (compress_and_encode(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
  }

  return PIDX_success;
}



PIDX_return_code hz_io(PIDX_io file, int mode)
{
  if (file->idx_dbg->debug_do_io == 1)
  {
    for (int j = file->idx_b->agg_level; j < file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count; j++)
    {
      file->time->hz_io_start[cvi][j] = MPI_Wtime();
      if (PIDX_file_io_per_process(file->hz_id, file->idx_b->block_layout_by_agg_group[j], mode) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      file->time->hz_io_end[cvi][j] = MPI_Wtime();
    }
  }

  return PIDX_success;
}



PIDX_return_code hz_encode_cleanup(PIDX_io file)
{
  // meta data and buffer cleanup
  if (hz_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  // meta data and buffer cleanup
  if (chunk_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}


static PIDX_return_code hz_init(PIDX_io file, int svi, int evi)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->hz_init_start[cvi] = PIDX_get_time();
  // Create the HZ encoding ID
  file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_c, file->idx_dbg, file->meta_data_cache, file->fs_block_size, svi, evi);

  // resolution for HZ encoding
  ret = PIDX_hz_encode_set_resolution(file->hz_id, file->idx_b->reduced_resolution_factor);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->hz_init_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code chunk_init(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;

  time->chunk_init_start[cvi] = PIDX_get_time();
  // Create the chunking ID
  file->chunk_id = PIDX_chunk_init(file->idx, file->idx_c, svi, evi);
  time->chunk_init_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code compression_init(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;

  time->compression_init_start[cvi] = PIDX_get_time();
  // Create the compression ID
  file->comp_id = PIDX_compression_init(file->idx, file->idx_c, svi, evi);
  time->compression_init_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code meta_data_create(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->chunk_meta_start[cvi] = PIDX_get_time();
  // metadata for chunking phase
  ret = PIDX_chunk_meta_data_create(file->chunk_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->chunk_meta_end[cvi] = PIDX_get_time();


  time->hz_meta_start[cvi] = PIDX_get_time();
  // metadata for hz encoding phase
  ret = PIDX_hz_encode_meta_data_create(file->hz_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  time->hz_meta_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code buffer_create(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->chunk_buffer_start[cvi] = PIDX_get_time();
  // Creating the buffers required for chunking
  ret = PIDX_chunk_buf_create(file->chunk_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_chunk;
  }
  time->chunk_buffer_end[cvi] = PIDX_get_time();

  time->hz_buffer_start[cvi] = PIDX_get_time();
  // Creating the buffers required for HZ encoding
  ret = PIDX_hz_encode_buf_create(file->hz_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  time->hz_buffer_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code compress_and_encode(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->chunk_start[cvi] = PIDX_get_time();
  // Perform Chunking
  if (file->idx_dbg->debug_do_chunk == 1)
  {
    ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
  }
  time->chunk_end[cvi] = PIDX_get_time();


  time->compression_start[cvi] = PIDX_get_time();
  // Perform Compression
  if (file->idx_dbg->debug_do_compress == 1)
  {
    ret = PIDX_compression(file->comp_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_compress;
    }
  }
  time->compression_end[cvi] = PIDX_get_time();


  // Perform HZ encoding
  if (file->idx_dbg->debug_do_hz == 1)
  {
    time->hz_start[cvi] = PIDX_get_time();
    ret = PIDX_hz_encode_write(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    // Verify the HZ encoding
    if (file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        //return PIDX_err_hz;
      }
    }
    time->hz_end[cvi] = PIDX_get_time();
  }

//

  time->chunk_buffer_free_start[cvi] = PIDX_get_time();
  // Destroy buffers allocated during chunking phase
  if (PIDX_chunk_buf_destroy(file->chunk_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_chunk;
  }
  time->chunk_buffer_free_end[cvi] = PIDX_get_time();

  return PIDX_success;
}

static PIDX_return_code encode_and_uncompress(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->hz_start[cvi] = PIDX_get_time();
  // Verify the HZ encoding
  if (file->idx_dbg->debug_hz == 1)
  {
    ret = HELPER_Hz_encode(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      //return PIDX_err_hz;
    }
  }

  // Perform HZ encoding
  if (file->idx_dbg->debug_do_hz == 1)
  {
    ret = PIDX_hz_encode_read(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }
  }
  time->hz_end[cvi] = PIDX_get_time();


  time->compression_start[cvi] = PIDX_get_time();
  // Perform Compression
  if (file->idx_dbg->debug_do_compress == 1)
  {
    ret = PIDX_decompression(file->comp_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_compress;
    }
  }
  time->compression_end[cvi] = PIDX_get_time();


  time->chunk_start[cvi] = PIDX_get_time();
  // Perform Chunking
  if (file->idx_dbg->debug_do_chunk == 1)
  {
    ret = PIDX_chunk(file->chunk_id, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
  }
  time->chunk_end[cvi] = PIDX_get_time();


  time->chunk_buffer_free_start[cvi] = PIDX_get_time();
  // Destroy buffers allocated during chunking phase
  ret = PIDX_chunk_buf_destroy(file->chunk_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_chunk;
  }
  time->chunk_buffer_free_end[cvi] = PIDX_get_time();

  return PIDX_success;
}

static PIDX_return_code hz_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->hz_buffer_free_start[cvi] = PIDX_get_time();
  // Destroy buffers allocated during HZ encoding phase
  ret = PIDX_hz_encode_buf_destroy(file->hz_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  time->hz_buffer_free_end[cvi] = PIDX_get_time();


  time->hz_cleanup_start[cvi] = PIDX_get_time();
  ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  PIDX_hz_encode_finalize(file->hz_id);
  time->hz_cleanup_end[cvi] = PIDX_get_time();

  return PIDX_success;
}


static PIDX_return_code chunk_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  time->chunk_cleanup_start[cvi] = PIDX_get_time();
  ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  PIDX_compression_finalize(file->comp_id);

  PIDX_chunk_finalize(file->chunk_id);
  time->chunk_cleanup_end[cvi] = PIDX_get_time();

  return PIDX_success;
}
