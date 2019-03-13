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

/**
 * \file PIDX_memory_layout_data_structs.h
 *
 * \author Sidharth Kumar
 *
 * Generic data structs-
 * -- N dimensional patch
 * -- Set of N dimension patches
 * -- HZ encoded buffer
 * -- Aggregation buffer
 */

#ifndef __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H
#define __PIDX_MEMORY_LAYOUT_DATA_STRUCTS_H

// data struct to store data for an application (simulation) N-dimensinal (1/2/3) patch
#include "PIDX_patch_structs.h"


// data struct to store a super patch which is also in application layout (n-dimension)
// a super patch is rectilinear block of data composed of alligned (non-overlapping) patchess
#include "PIDX_super_patch_structs.h"


// data struct to store the restructured grid (details like extents of patches and ranks associated with the patches)
#include "PIDX_restructured_grid_structs.h"


// data struct to store encoded hz buffer (hz buffer is 2D buffer (#of levels x #of samples in the level)
#include "PIDX_hz_buffer_structs.h"


// data struct to store aggregated buffer
#include "PIDX_aggregation_buffer_structs.h"


#endif
