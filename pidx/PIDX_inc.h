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
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <fcntl.h>
#include <unistd.h>

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

#include "./comm/PIDX_comm.h"

#include "./meta_data_cache/PIDX_meta_data_cache.h"

#include "./data_handle/PIDX_data_layout.h"
#include "./data_handle/PIDX_blocks.h"
#include "./data_handle/PIDX_idx_data_structs.h"

#include "./core/PIDX_header/PIDX_header_io.h"
#include "./core/PIDX_raw_rst/PIDX_raw_rst.h"
#include "./core/PIDX_idx_rst/PIDX_idx_rst.h"
#include "./core/PIDX_hz/PIDX_hz_encode.h"
#include "./core/PIDX_block_rst/PIDX_block_restructure.h"
#include "./core/PIDX_cmp/PIDX_compression.h"
#include "./core/PIDX_agg/PIDX_agg.h"
#include "./core/PIDX_file_io/PIDX_file_io.h"
#include "./core/PIDX_wavelet/PIDX_wavelet.h"


#include "./io/PIDX_io.h"
#include "./io/idx/serial/serial_idx_io.h"
#include "./io/idx/no_partition/idx_io.h"
#include "./io/idx/local_partition/local_partition_idx_io.h"
#include "./io/idx/global_partition/global_partition_idx_io.h"
#include "./io/raw/raw_io.h"


#include "./io/idx/io_setup.h"
#include "./io/idx/bit_string.h"
#include "./io/idx/idx_restructure.h"
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
