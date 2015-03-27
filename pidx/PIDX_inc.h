#ifndef __PIDX_INC_H
#define __PIDX_INC_H

#include "PIDX_config.h"

#ifdef __bg__
  #define _XOPEN_SOURCE 600
  #define _GNU_SOURCE
#endif

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#else
  #include <sys/time.h>
#endif

#ifdef PIDX_HAVE_LOSSY_ZFP
  #include "zfp.h"
  #include "fpzip.h"
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


#ifdef __cplusplus
extern "C" {
#endif

#include "PIDX_data_layout.h"
#include "PIDX_data_types.h"
#include "PIDX_error_codes.h"
#include "PIDX_file_access_modes.h"

#include "PIDX_blocks.h"
#include "PIDX_idx_data_structs.h"
#include "PIDX_comm.h"
#include "PIDX_utils.h"
#include "PIDX_point.h"
#include "PIDX_file_name.h"
#include "PIDX_topology.h"

#include "PIDX_header_io.h"
#include "PIDX_rst.h"
#include "PIDX_hz_encode.h"
#include "PIDX_block_restructure.h"
#include "PIDX_compression.h"
#include "PIDX_agg.h"
#include "PIDX_io.h"

#ifdef __cplusplus
}
#endif

#endif
