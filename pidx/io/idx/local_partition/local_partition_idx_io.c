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


/// local Partitioned IDX Write Steps
PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int svi, int evi)
{

#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 0] : var range %d %d\n", svi, evi);
#endif


  // Step 1: Compute the restructuring box (grid) size
  // We need to identify the restructuring box size (super patch first because we directly use that to populate the
  // bitstring, which is first written out to the .idx file and used in almost every phase of IO
  if (set_rst_box_size_for_write(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 1] : rst box size %lld %lld %lld\n", (unsigned long long)file->restructured_grid->patch_size[0], (unsigned long long)file->restructured_grid->patch_size[1], (unsigned long long)file->restructured_grid->patch_size[2]);
#endif


  // Step 2: Setting the stage for restructuring (meta data)
  if (idx_restructure_setup(file, svi, evi - 1) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 2] : Setting up the restructuring phase\n");
#endif


  // Step 3: Perform data restructuring
  // restructuring is done keeping partitioning in mind
  // restructuring phase ensures that a process is not holding more that one super patch
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 3] : Restructuring phase\n");
#endif


  // Step 4: Group the processes holding the super patch into new communicator rst_comm
  // This step is intended to reduce MPI overhead by having a comm with fewer processes
  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 4] : Creating restructuring communicator\n");
#endif


  // proceed only if a process holds a superpatch, others just wait at completion of io
  // essential for partitioning
  PIDX_variable var0 = file->idx->variable[svi];
  if (var0->restructured_super_patch_count == 1)
  {
    // Step 5:  Perform partitioning
    // Create new communicator
    if (partition(file, svi) != PIDX_success)
    {
      fprintf(stderr,"Error [File %s Line %d]: partition\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 5] : Creating partitions (count %d %d %d) (size %d %d %d) (offset %d %d %d)\n", (int)file->idx->partition_count[0], (int)file->idx->partition_count[1], (int)file->idx->partition_count[2], (int)file->idx->partition_size[0], (int)file->idx->partition_size[1], (int)file->idx->partition_size[2], (int)file->idx->partition_offset[0], (int)file->idx->partition_offset[1], (int)file->idx->partition_offset[2]);
#endif
    

    // Step 6: Populate the bit string of the global dataset
    if (populate_bit_string(file, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 6] : bounds %lld %lld %lld bitstring %s\n", (unsigned long long)file->idx->bounds[0], (unsigned long long)file->idx->bounds[1], (unsigned long long)file->idx->bounds[2], file->idx->bitSequence);
#endif


    // Step 7: write the global .idx file
    if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 7] : Writing global idx file\n");
#endif


    // Step 8: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 8] : Adjusting the offset of boxes\n");
#endif


    // Step 9:  Populate block layout for the partitions (this is a local step, from now onwards the processes work within its partition only)
    if (populate_block_layout_and_buffers(file, svi, evi, PIDX_WRITE, PIDX_LOCAL_PARTITION_IDX_IO) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 9] : Populating blocks maxh %d max file count %d\n", file->idx->maxh, file->idx->max_file_count);
#endif


    // Iterate through the variables, variable_pipe_length at a point
    // variable_pipe_length is computed based on the configuration of the run. If there are enough number of processes
    // then all the variables are aggregated at once and then variable_pipe_length is (variable_count - 1), otherwise
    // the variables are worked in smaller packs.
    for (uint32_t si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      uint32_t ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 10:  Setup HZ encoding Phase
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 10]: Setting up hz phase %d %d\n", si, ei);
#endif

      // Step 11: Perform HZ encoding
      if (hz_encode(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 11]: Performing hz phase\n");
#endif


      // Step 12: This is not performed by default, it only happens when aggregation is
      // turned off or when there are a limited number of aggregators
      if (hz_io(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 12]: Performing hz io\n");
#endif


      // Setup 13: Setup aggregation buffers
      if (aggregation_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 13]: Setting up aggregation phase\n");
#endif


      // Setup 14: Performs data aggregation
      if (aggregation(file, si, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 14]: Performing aggregation\n");
#endif


      // Setup 15: Performs actual file io
      if (file_io(file, si, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 15]: Performing file io\n");
#endif


      // Step 16: free aggregation buffers
      if (aggregation_cleanup(file, si) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 16]: Cleaning up aggregation phase\n");
#endif


      // Step 17: Cleanup HZ buffers
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
#if DEBUG_OUTPUT
      if (file->idx_c->simulation_rank == 0)
        fprintf(stderr, "[Local partition idx io 17]: Cleaning up hz phase\n");
#endif
    }

    // Step 18: free block layout and agg related buffer
    if (destroy_block_layout_and_buffers(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 18]: Cleaning up block layout\n");
#endif

    // Step 19: free partitioning buffer
    if (partition_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if DEBUG_OUTPUT
    if (file->idx_c->simulation_rank == 0)
      fprintf(stderr, "[Local partition idx io 18]: Cleaning up partition\n");
#endif
  }

  // Step 20: free the restructured communicator
  if (free_restructured_communicators(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 18]: Freeing up restructuring communicators\n");
#endif

  // Step 21: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if DEBUG_OUTPUT
  if (file->idx_c->simulation_rank == 0)
    fprintf(stderr, "[Local partition idx io 18]: Cleaning restructuring phase\n");
#endif


  return PIDX_success;
}




/// local Partitioned IDX Read Steps
PIDX_return_code PIDX_local_partition_idx_read(PIDX_io file, int svi, int evi)
{
  // This function is called for restart styke reads only when the number of processes
  // needed to read the data is same as the number of processes used to create the dataset


  // Step 1: Compute the restructuring box (grid) size
  // We need to identify the restructuring box size (super patch first because we directly use that to populate the
  // bitstring, which is first written out to the .idx file and used in almost every phase of IO
  // This is the same as the bix used to write the data as number of processes is same
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

  // Step 4: Group the processes holding the super patch into new communicator rst_comm
  // This step is intended to reduce MPI overhead by having a comm with fewer processes
  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // proceed only if a process holds a superpatch, others just wait at completion of io
  // essential for partitioning
  PIDX_variable var0 = file->idx->variable[svi];
  if (var0->restructured_super_patch_count == 1)
  {
    // Step 5:  Perform partitioning
    // Create new communicator
    if (partition(file, svi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 6: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 7:  Populate block layout for the partitions (this is a local step, from now onwards the processes work within its partition only)
    if (populate_block_layout_and_buffers(file, svi, evi, PIDX_READ, PIDX_LOCAL_PARTITION_IDX_IO) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (uint32_t si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      uint32_t ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 8:  Setup HZ encoding Phase
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 9: This is not performed by default, it only happens when aggregation is
      // turned off or when there are a limited number of aggregators
      if (hz_io(file, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 10: Setup aggregation buffers
      if (aggregation_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 11: Performs actual file io
      if (file_io(file, si, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 12: Performs data aggregation
      if (aggregation(file, si, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 13: free aggregation buffers
      if (aggregation_cleanup(file, si) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 14: Perform HZ encoding
      if (hz_encode(file, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 15: Cleanup hz buffers
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 16: Cleanup the group and IDX related meta-data
    if (destroy_block_layout_and_buffers(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 17: Partition cleanup
    if (partition_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 18: re-adjust the box size
    if (re_adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 19:  Restrucure
  if (idx_restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 20: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}




