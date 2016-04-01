#ifndef __PIDX_INC_H
#define __PIDX_INC_H

#define PIDX_MAX_DIMENSIONS 5

#define PIDX_HAVE_MPI 1
#define PIDX_HAVE_ZFP 0
#define PIDX_HAVE_PNETCDF 0
#define PIDX_HAVE_NETCDF 0
#define PIDX_HAVE_HDF5 0
#define PIDX_HAVE_NVISUSIO 0

//#include "PIDX_config.h"

#define SIMULATE_IO 0
#define PIDX_MAX_TEMPLATE_DEPTH 6

#ifndef __cplusplus
#  define _XOPEN_SOURCE 600
#endif

#ifdef BGQ
  #define _XOPEN_SOURCE 600
#ifndef _GNU_SOURCE
    #define _GNU_SOURCE
#endif
#endif

//#define pmin(x, y) ((x) < (y) ? (x) : (y))
//#ifndef __cplusplus
//  #define max(x, y) ((x) > (y) ? (x) : (y))
//#endif

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
//#include <sys/types.h>
//#include <time.h>
//#include <stdint.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#else
  #include <sys/time.h>
#endif

#if PIDX_HAVE_ZFP
  #include <zfp.h>  
#endif


#if defined(BGL) || defined(BGP) || defined(BGQ)

#include <mpi.h>
#include <math.h>

#ifdef BGL 

#include <bglpersonality.h>
#include <rts.h>

#define   get_personality                rts_get_personality
#define   get_processor_id               rts_get_processor_id
#define   Personality                    BGLPersonality
#define   Personality_getLocationString  BGLPersonality_getLocationString
#define   Personality_numIONodes         BGLPersonality_numIONodes
#define   Personality_numPsets           BGLPersonality_numPsets
#define   Personality_numNodesInPset     BGLPersonality_numNodesInPset
#define   Personality_rankInPset         BGLPersonality_rankInPset
#define   Personality_psetNum            BGLPersonality_psetNum

#endif
#ifdef BGP

#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>

#define   get_personality                Kernel_GetPersonality
#define   get_processor_id               Kernel_PhysicalProcessorID
#define   Personality                    _BGP_Personality_t
#define   Personality_getLocationString  BGP_Personality_getLocationString
#define   Personality_numIONodes         BGP_Personality_numIONodes
#define   Personality_numNodesInPset     BGP_Personality_psetSize
#define   Personality_rankInPset         BGP_Personality_rankInPset
#define   Personality_psetNum            BGP_Personality_psetNum

#endif

#ifdef BGQ

#include <kernel/process.h>
#include <kernel/location.h>
#include <firmware/include/personality.h>
#include <mpix.h>

#define   get_personality                Kernel_GetPersonality
#define   get_processor_id               Kernel_PhysicalProcessorID
#define   Personality                    Personality_t

#endif
#endif

enum IO_MODE {PIDX_READ, PIDX_WRITE};

#ifdef __cplusplus
extern "C" {
#endif

#define PIDX_NO_COMPRESSION 0
#define PIDX_CHUNKING_ONLY 1
#define PIDX_CHUNKING_ZFP 2

#define PIDX_row_major                           0
#define PIDX_column_major                        1


#include "./utils/PIDX_error_codes.h"
#include "./utils/PIDX_point.h"
#include "./utils/PIDX_utils.h"
#include "./utils/PIDX_file_name.h"
#include "./utils/PIDX_file_access_modes.h"

#include "./comm/PIDX_comm.h"


#include "./data_handle/PIDX_data_layout.h"
#include "./data_handle/PIDX_data_types.h"
#include "./data_handle/PIDX_blocks.h"
#include "./data_handle/PIDX_idx_data_structs.h"


#include "./topology/PIDX_topology.h"

#include "./core/PIDX_header/PIDX_header_io.h"
#include "./core/PIDX_rst/PIDX_rst.h"
#include "./core/PIDX_multi_patch_rst/PIDX_multi_patch_rst.h"
#include "./core/PIDX_hz/PIDX_hz_encode.h"
#include "./core/PIDX_block_rst/PIDX_block_restructure.h"
#include "./core/PIDX_cmp/PIDX_compression.h"
#include "./core/PIDX_agg/PIDX_agg.h"
#include "./core/PIDX_file_io/PIDX_file_io.h"



#ifdef __cplusplus
}
#endif

#endif
