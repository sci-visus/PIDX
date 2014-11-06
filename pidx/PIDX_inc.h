#ifndef __PIDX_INC_H
#define __PIDX_INC_H

#include "PIDX_config.h"

#ifdef BLUEGENEP
  #define _XOPEN_SOURCE 600
#endif

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

#include "PIDX_error_codes.h"
#include "PIDX_file_access_modes.h"
#include "PIDX_data_types.h"
#include "PIDX_data_layout.h"

#include "Generic_data_structs.h"
#include "PIDX_comm.h"
#include "PIDX_blocks.h"
#include "PIDX_data_structs.h"
#include "PIDX_utils.h"
#include "PIDX_point.h"

#include "PIDX_header_io.h"
#include "PIDX_rst.h"
#include "PIDX_hz_encode.h"
#include "PIDX_agg.h"
#include "PIDX_io.h"

#endif