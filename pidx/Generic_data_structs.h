/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/
 
 
#define _XOPEN_SOURCE 600
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

#include "PIDX_comm.h"

#define PIDX_MAX_DIMENSIONS 5

#ifndef __GENERIC_DATA_STRUCTS_H
#define __GENERIC_DATA_STRUCTS_H

struct PIDX_Ndim_buffer_struct
{
  int offset[PIDX_MAX_DIMENSIONS];
  int count[PIDX_MAX_DIMENSIONS];
  unsigned char* buffer;
};
typedef struct PIDX_Ndim_buffer_struct* Ndim_buffer;

struct PIDX_Ndim_buffer_group_struct
{
  int type;
  int count;
  Ndim_buffer *block;
  int power_two_offset[PIDX_MAX_DIMENSIONS];
  int power_two_count[PIDX_MAX_DIMENSIONS];
};
typedef struct PIDX_Ndim_buffer_group_struct* Ndim_buffer_group; 

struct PIDX_HZ_buffer_struct
{
  int HZ_level_from;
  int HZ_level_to;
  unsigned char** buffer;
  int* samples_per_level;
  long long *allign_start_hz;
  long long *allign_end_hz;
  int **allign_offset;
  int **allign_count;
  long long* buffer_index;
};
typedef struct PIDX_HZ_buffer_struct* HZ_buffer;

struct PIDX_HZ_Agg_buffer_struct
{
  unsigned char* buffer;
  int buffer_size;
  
  int file_number;
  int var_number;
  int sample_number;
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

#endif
