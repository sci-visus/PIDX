/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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


#include "PIDX.h"

///
/// \brief The PIDX_file_descriptor struct is the PIDX File descriptor
/// (equivalent to the descriptor returned by) POSIX or any other IO framework
///
struct PIDX_file_descriptor
{
  int flags;                                    ///< idx file open and create mode

  // file system info
  int fs_block_size;                            ///< file system block size which is queryed once at the beginning

  // flush related
  int variable_index_tracker;                   ///< tracking upto which variable io has been done (used for flushing)
  int local_variable_index;                     ///< starting index of variable that needs to be written out before a flush
  int local_variable_count;                     ///< total number of variables that is written out in a flush

  // IDX related
  idx_dataset idx;                              ///< Contains all IDX related info
  idx_blocks idx_b;                             ///< idx block related
  idx_comm idx_c;                               ///< MPI related
  idx_debug idx_dbg;                            ///< Flags for debugging

  // IO phases
  PIDX_io io;                                   ///< this descriptor contains pointers to descriptors to all other sub-phases like restructuring, HZ encoding and aggregation

  // Timming
  PIDX_time time;                               ///< For detailed time profiling of all phases

  // for caching HZ indices
  PIDX_metadata_cache meta_data_cache;          ///< enables caching across time steps

  // for restructuring and partitioning
  PIDX_restructured_grid restructured_grid;     ///< contains information of the restructured grid
};
