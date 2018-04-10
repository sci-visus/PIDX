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
#ifndef __PIDX_INC_H
#define __PIDX_INC_H

#include "PIDX_define.h"

#include <time.h>
//#include <byteswap.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <fcntl.h>
#include <stdint.h>

#if defined _MSC_VER
  #include "PIDX_windows_define.h"
#else
  #include <unistd.h>
  #include <arpa/inet.h>
#endif

#include <mpi.h>

#if PIDX_HAVE_ZFP
  #include <zfp.h>
#endif

#if PIDX_HAVE_PMT
  #include <pidx_insitu.h>
#endif

#if PIDX_HAVE_VTK
  #include <renderer.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define PIDX_MIN(a,b) (((a)<(b))?(a):(b))


#include "./utils/PIDX_error_codes.h"
#include "./utils/PIDX_point.h"
#include "./utils/PIDX_utils.h"
#include "./utils/PIDX_file_name.h"
#include "./utils/PIDX_file_access_modes.h"
#include "./utils/PIDX_buffer.h"

#include "./comm/PIDX_comm.h"

#include "./metadata/PIDX_metadata_cache.h"

#include "./data_handle/PIDX_data_layout.h"
#include "./data_handle/PIDX_blocks.h"
#include "./data_handle/PIDX_idx_data_structs.h"

#include "./core/PIDX_header/PIDX_header_io.h"
#include "./core/PIDX_raw_rst/PIDX_raw_rst.h"
#include "./core/PIDX_idx_rst/PIDX_idx_rst.h"
#include "./core/PIDX_particles_rst/PIDX_particles_rst.h"
#include "./core/PIDX_hz/PIDX_hz_encode.h"
#include "./core/PIDX_block_rst/PIDX_block_restructure.h"
#include "./core/PIDX_cmp/PIDX_compression.h"
#include "./core/PIDX_agg/PIDX_agg.h"
#include "./core/PIDX_file_io/PIDX_file_io.h"
#include "./core/PIDX_wavelet/PIDX_wavelet.h"


#include "./io/PIDX_io.h"
#include "./io/idx/particle/particle_io.h"
#include "./io/idx/serial/serial_idx_io.h"
#include "./io/idx/no_partition/idx_io.h"
#include "./io/idx/local_partition/local_partition_idx_io.h"
#include "./io/raw/raw_io.h"


#include "./io/idx/io_setup.h"
#include "./io/idx/bit_string.h"
#include "./io/idx/idx_restructure.h"
#include "./io/idx/particles_restructure.h"
#include "./io/raw/raw_restructure.h"
#include "./io/idx/partition.h"
#include "./io/idx/local_buffer.h"
#include "./io/idx/headers.h"
#include "./io/idx/rst_blocks.h"
#include "./io/idx/sim_blocks.h"
#include "./io/idx/hz_buffer.h"
#include "./io/idx/agg_io.h"
#include "./io/idx/timming.h"

#ifdef __cplusplus
}
#endif

#endif
