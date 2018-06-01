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

#ifndef __PIDX_VARIABLE_STRUCTS_H
#define __PIDX_VARIABLE_STRUCTS_H


struct PIDX_variable_struct
{
  // General Info
  char var_name[1024];                                       ///< Variable name
  int vps;                                                   ///< values per sample, Vector(3), scalar(1), or n
  int bpv;                                                   ///< Number of bits each need
  PIDX_data_type type_name;                                  ///< Name of the type uint8, bob
  PIDX_data_layout data_layout;                              ///< Row major or column major


  // buffer (before, after HZ encoding phase)
  int sim_patch_count;                                       ///< Number of patches/blocks that the simulation feeds to PIDX (application layout)
  PIDX_patch sim_patch[1024];                                ///< Pointer to the patches


  // buffer after restructuring
  int32_t restructured_super_patch_count;                    ///< Number of super patches after restructuring, the way restructuring is setup now, this can only be 1 (a process has a super patch) or a 0 (a process does not have a super patch)
  PIDX_super_patch restructured_super_patch;                 ///< Pointer to the super patch


  // buffer for chunked data (only used with zfp compression)
  PIDX_super_patch chunked_super_patch;                      ///< Pointer to the super patch formed after block restructuring


  // buffer to hold the HZ encoded data
  HZ_buffer hz_buffer;                                       ///< HZ encoded buffer of the super patche


  // this is used only in raw io mode. With raw io, it is possible for a process to hold more than one super patch.
  int raw_io_restructured_super_patch_count;                ///< number of super patch after restructuring, can be greater than equal to 0
  PIDX_super_patch* raw_io_restructured_super_patch;        ///< pointer to the restructured super patches
};
typedef struct PIDX_variable_struct* PIDX_variable;

#endif
