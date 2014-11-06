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
 
#ifndef __GENERIC_DATA_STRUCTS_H
#define __GENERIC_DATA_STRUCTS_H

struct PIDX_Ndim_buffer_struct
{
  long long offset[PIDX_MAX_DIMENSIONS];
  long long count[PIDX_MAX_DIMENSIONS];
  unsigned char* buffer;
};
typedef struct PIDX_Ndim_buffer_struct* Ndim_buffer;

struct PIDX_Ndim_buffer_group_struct
{
  int type;
  int count;
  Ndim_buffer *block;
  long power_two_offset[PIDX_MAX_DIMENSIONS];
  long power_two_count[PIDX_MAX_DIMENSIONS];
};
typedef struct PIDX_Ndim_buffer_group_struct* Ndim_buffer_group; 

struct PIDX_HZ_buffer_struct
{
  /// starting HZ level
  int HZ_level_from;
  
  /// ending HZ level
  int HZ_level_to;
  
  /// number of samples in the hz levels (#level = HZ_level_from - HZ_level_to + 1)
  long long *samples_per_level;
  
  /// Starting HZ index at of the data at all the HZ levels
  long long *start_hz_index;
  
  /// Ending HZ index at of the data at all the HZ levels
  long long *end_hz_index;
  
  /// HZ indices of the data (used only when no restructuring phsae is used)
  long long *buffer_index;
  
  // data buffer at all the HZ levels
  unsigned char** buffer;
};
typedef struct PIDX_HZ_buffer_struct* HZ_buffer;

struct PIDX_HZ_Agg_buffer_struct
{
  int file_number;
  int var_number;
  int sample_number;
  
  unsigned char* buffer;
  int buffer_size;
};
typedef struct PIDX_HZ_Agg_buffer_struct* Agg_buffer;

#endif