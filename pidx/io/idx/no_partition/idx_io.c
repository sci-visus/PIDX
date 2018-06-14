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

#include "../../../PIDX_inc.h"


// IDX Write Steps
PIDX_return_code PIDX_idx_write(PIDX_io file, int svi, int evi)
{
  // Steps 1-4 are global in nature that they apply to the overall idx file
  // Steps 5-14 apply to variables of the idx file
  // Steps 15-16 are cleanup steps applying on the entire idx file

  // Step 1: Compute the restructuring box (grid) size
  // We need to identify the restructuring box size (super patch first because we directly use that to populate the
  // bitstring, which is first written out to the .idx file and used in almost every phase of IO
  if (set_rst_box_size_for_write(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2: Setting the stage for restructuring (meta data)
  if (idx_restructure_setup(file, svi, evi - 1) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 3: Perform data restructuring
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 4: Group the processes holding the super patch into new communicator rst_comm
  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // proceed only if a process holds a superpatch, others just wait at completion of io
  PIDX_variable var0 = file->idx->variable[svi];
  if (var0->restructured_super_patch_count == 1)
  {
    // Step 5: populate idx blocks and agg related buffers
    if (populate_block_layout_and_buffers(file, svi, evi, PIDX_WRITE, PIDX_IDX_IO) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // variable_pipe_length is computed based on the configuration of the run. If there are enough number of processes
    // then all the variables are aggregated at once and then variable_pipe_length is (variable_count - 1), otherwise
    // the variables are worked in smaller packs.
    for (uint32_t si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      uint32_t ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 6: Setup HZ buffers
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 7: Perform HZ encoding
      if (hz_encode(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 8: This is not performed by default, it only happens when aggregation is
      // turned off or when there are a limited number of aggregators
      if (hz_io(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 9: Setup aggregation data buffers
      if (aggregation_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 10: Performs data aggregation
      if (aggregation(file, si, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 11: Performs actual file io
      if (file_io(file, si, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 12: free aggregation buffers
      if (aggregation_cleanup(file, si) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 13: Cleanup hz buffers and ids
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 14: free block layout and agg related buffer
    if (destroy_block_layout_and_buffers(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 15: free the restructured communicator
  if (free_restructured_communicators(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 16: cleanup restructuring buffers
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



// IDX Read Steps
PIDX_return_code PIDX_idx_read(PIDX_io file, int svi, int evi)
{
  // Step 1: Setup restructuring buffers
  if (set_rst_box_size_for_read(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2: setup the restructuring buffers
  if (idx_restructure_setup(file, svi, evi - 1) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 3: Create restructuring communicator, discarding the processes that are not going to hold any super patch
  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  // proceed only if a process holds a superpatch, others just wait at completion of io
  PIDX_variable var0 = file->idx->variable[svi];
  if (var0->restructured_super_patch_count == 1)
  {
    // Step 4: populate idx block layout based on the super patches
    if (populate_block_layout_and_buffers(file, svi, evi, PIDX_READ, PIDX_IDX_IO) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (uint32_t si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      uint32_t ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 5: Setup HZ buffers
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 6: this again is an optional phase, which gets activated when aggregation is turned off or when we do not have enough aggregators
      if (hz_io(file, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 7: Setting for file io phase by creating aggregation buffers
      if (aggregation_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 8: Performs actual file io
      if (file_io(file, si, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 9: Scatter data from aggregators to all processes (reverse aggregation)
      if (aggregation(file, si, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 10: free aggregation buffers
      if (aggregation_cleanup(file, si) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 11: Perform reverse HZ encoding, converting data from idx layout to row/column major application layout
      if (hz_encode(file, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 12: Cleanup hz buffers and ids
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 13: free block layout structure
    if (destroy_block_layout_and_buffers(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 14: free the restructured communicator
  if (free_restructured_communicators(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 15: Perform reverse of data restructuring, moving from super patches to simulation patches
  if (idx_restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 16: Cleanup of restructuring (super patch) buffers
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}
