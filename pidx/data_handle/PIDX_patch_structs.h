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



#ifndef __PIDX_PATCH_STRUCTS_H
#define __PIDX_PATCH_STRUCTS_H


/// Struct to store the row/column major chunk of data given by application
struct PIDX_patch_struct
{
  uint64_t offset[PIDX_MAX_DIMENSIONS];               ///< logical offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension) in the 3D global space
  uint64_t size[PIDX_MAX_DIMENSIONS];                 ///< logical size (extents) in each of the dimensions for the data chunk

  double physical_offset[PIDX_MAX_DIMENSIONS];        ///< physical offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension) in the 3D global space
  double physical_size[PIDX_MAX_DIMENSIONS];          ///< physical size (extents) in each of the dimensions for the data chunk

  // TODO WILL: The particles needing to modify inputs and buffer sizes
  // which were previously thought to be fixed makes this struct really now
  // a dual-mode pain to deal with. Do something better.
  uint64_t particle_count;
  uint64_t *read_particle_count;

  // TODO WILL: Do we want some other way of differentiating the particle/grid patches?
  union {
    /// The data buffer for grid variables (particle_count = 0)
    /// (or a redundant case with variable->is_particle == 1)
    unsigned char* buffer;
    unsigned char** tbuffer;
    /// The data buffer for particle variables (particle_count != 0). This
    /// buffer is allocated by PIDX. TODO: Should the user free it? Or
    /// should we have a PIDX_free_variable? Does such a function already exist?
    unsigned char** read_particle_buffer;
  };
  uint64_t read_particle_buffer_capacity;
};
typedef struct PIDX_patch_struct* PIDX_patch;


#endif
