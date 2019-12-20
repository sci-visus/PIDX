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



#ifndef __PIDX_RESTRUCTURED_GRID_STRUCTS_H
#define __PIDX_RESTRUCTURED_GRID_STRUCTS_H


/// Struct to store the patch info of the restructured grid
struct PIDX_Ndim_empty_patch_struct
{
  int rank;                                             ///< rank associated with the patch
  int is_boundary_patch;                                ///< 1 if the patch is at the boundary (non-pwer two dataset) 0 otherwise

  uint64_t offset[PIDX_MAX_DIMENSIONS];                 ///< offset of the data chunk (of PIDX_MAX_DIMENSIONS dimension)
  uint64_t size[PIDX_MAX_DIMENSIONS];                   ///< size (extents) in each of the dimensions for the data chunk

  double physical_offset[PIDX_MAX_DIMENSIONS];          ///< physical extents (for particles)
  double physical_size[PIDX_MAX_DIMENSIONS];            ///< physical extents (for particles)
};
typedef struct PIDX_Ndim_empty_patch_struct* Ndim_empty_patch;



/// Struct to store the restructured grid
struct PIDX_grid_struct
{
  double physical_patch_size[PIDX_MAX_DIMENSIONS];      ///< restructured grid size (physical), used for particle simulations
  uint64_t patch_size[PIDX_MAX_DIMENSIONS];             ///< restructured grid size (logical)

  uint64_t total_patch_count[PIDX_MAX_DIMENSIONS];      ///< total number of patches forming the restructured super patch

  Ndim_empty_patch* patch;                              ///< patch contained in the super patch
};
typedef struct PIDX_grid_struct* PIDX_restructured_grid;

#endif
