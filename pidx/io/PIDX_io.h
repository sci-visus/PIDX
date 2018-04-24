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
#ifndef __PIDX_IO_H
#define __PIDX_IO_H

///
/// \brief The PIDX_io_descriptor struct
/// Contains descriptor of all io phases
///
struct PIDX_io_descriptor
{
  // FS info
  int fs_block_size;                                ///< file system block size which is queryed once at the beginning

  int variable_index_tracker;                       ///< tracking upto which variable io has been done (used for flushing)

  // Different IO phases
  PIDX_header_io_id header_io_id;                   ///< Creates the file hierarchy and populates the raw header
  // only one of the three is activated at a time
  PIDX_particles_rst_id particles_rst_id;           ///< Multi-patch restructuring phase id for particles
  PIDX_raw_rst_id raw_rst_id;                       ///< Multi-patch restructuring phase id for raw io
  PIDX_idx_rst_id idx_rst_id;                       ///< Multi-patch restructuring phase id for idx io
  PIDX_chunk_id chunk_id;                           ///< Block restructuring id (prepration for compression)
  PIDX_comp_id comp_id;                             ///< Compression (lossy and lossless) id
  PIDX_hz_encode_id hz_id;                          ///< HZ encoding phase id
  PIDX_agg_id** agg_id;                             ///< Aggregation phase
  PIDX_file_io_id** io_id;                          ///< File io

  // IDX related
  idx_dataset idx;                                  ///< Contains all IDX related info
  idx_blocks idx_b;                                 ///< idx block related
  idx_comm idx_c;                                   ///< MPI related
  idx_debug idx_dbg;                                ///< Flags for debugging

  // for caching HZ indices
  PIDX_metadata_cache meta_data_cache;              ///< enables caching across time steps

  // for restructuring and partitioning
  PIDX_restructured_grid restructured_grid;         ///< contains information of the restructured grid

  // Timming
  PIDX_time time;                                   ///< For detailed time profiling of all phases
};
typedef struct PIDX_io_descriptor* PIDX_io;


///
/// \brief PIDX_io_init
/// \param idx_meta_data
/// \param idx_c
/// \param idx_dbg
/// \param meta_data_cache
/// \param idx_b
/// \param restructured_grid
/// \param time
/// \param fs_block_size
/// \param variable_index_tracker
/// \return
///
PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_metadata_cache meta_data_cache, idx_blocks idx_b, PIDX_restructured_grid restructured_grid, PIDX_time time, int fs_block_size, int variable_index_tracker);


///
/// \brief PIDX_write
/// \param file
/// \param start_var_index
/// \param end_var_index
/// \param MODE
/// \return
///
PIDX_return_code PIDX_write(PIDX_io file, int start_var_index, int end_var_index, int MODE);


///
/// \brief PIDX_read
/// \param file
/// \param svi
/// \param evi
/// \param MODE
/// \return
///
PIDX_return_code PIDX_read(PIDX_io file, int svi, int evi, int MODE);


///
/// \brief PIDX_io_finalize
/// \param file
/// \return
///
PIDX_return_code PIDX_io_finalize(PIDX_io file);

#endif
