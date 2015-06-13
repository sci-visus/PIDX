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

#include "PIDX_inc.h"

PIDX_data_type INT8          = "1*int8";
PIDX_data_type INT8_GA       = "2*int8";
PIDX_data_type INT8_RGB      = "3*int8";
PIDX_data_type INT8_RGBA     = "4*int8";

PIDX_data_type UINT8         = "1*uint8";
PIDX_data_type UINT8_GA      = "2*uint8";
PIDX_data_type UINT8_RGB     = "3*uint8";
PIDX_data_type UINT8_RGBA    = "4*uint8";

PIDX_data_type INT16         = "1*int16";
PIDX_data_type INT16_GA      = "2*int16";
PIDX_data_type INT16_RGB     = "3*int16";
PIDX_data_type INT16_RGBA    = "4*int16";

PIDX_data_type UINT16        = "1*uint16";
PIDX_data_type UINT16_GA     = "2*uint16";
PIDX_data_type UINT16_RGB    = "3*uint16";
PIDX_data_type UINT16_RGBA   = "4*uint16";

PIDX_data_type INT32         = "1*int32";
PIDX_data_type INT32_GA      = "2*int32";
PIDX_data_type INT32_RGB     = "3*int32";
PIDX_data_type INT32_RGBA    = "4*int32";

PIDX_data_type UINT32        = "1*uint32";
PIDX_data_type UINT32_GA     = "2*uint32";
PIDX_data_type UINT32_RGB    = "3*uint32";
PIDX_data_type UINT32_RGBA   = "4*uint32";

PIDX_data_type INT64         = "1*int64";
PIDX_data_type INT64_GA      = "2*int64";
PIDX_data_type INT64_RGB     = "3*int64";
PIDX_data_type INT64_RGBA    = "4*int64";

PIDX_data_type UINT64        = "1*uint64";
PIDX_data_type UINT64_GA     = "2*uint64";
PIDX_data_type UINT64_RGB    = "3*uint64";
PIDX_data_type UINT64_RGBA   = "4*uint64";

PIDX_data_type FLOAT32       = "1*float32";
PIDX_data_type FLOAT32_GA    = "2*float32";
PIDX_data_type FLOAT32_RGB   = "3*float32";
PIDX_data_type FLOAT32_RGBA  = "4*float32";

PIDX_data_type FLOAT64       = "1*float64";
PIDX_data_type FLOAT64_GA    = "2*float64";
PIDX_data_type FLOAT64_RGB   = "3*float64";
PIDX_data_type FLOAT64_RGBA  = "4*float64";


PIDX_return_code PIDX_default_bits_per_datatype(PIDX_data_type type, int* bits)
{
  if (strcmp(type, INT8) == 0)
    *bits = 8;
  else if (strcmp(type, INT8_GA) == 0)
    *bits = 16;
  else if (strcmp(type, INT8_RGB) == 0)
    *bits = 24;
  else if (strcmp(type, INT8_RGBA) == 0)
    *bits = 32;

  if (strcmp(type, UINT8) == 0)
    *bits = 8;
  else if (strcmp(type, UINT8_GA) == 0)
    *bits = 16;
  else if (strcmp(type, UINT8_RGB) == 0)
    *bits = 24;
  else if (strcmp(type, UINT8_RGBA) == 0)
    *bits = 32;

  if (strcmp(type, INT16) == 0)
    *bits = 16;
  else if (strcmp(type, INT16_GA) == 0)
    *bits = 32;
  else if (strcmp(type, INT16_RGB) == 0)
    *bits = 48;
  else if (strcmp(type, INT16_RGBA) == 0)
    *bits = 64;

  if (strcmp(type, UINT16) == 0)
    *bits = 16;
  else if (strcmp(type, UINT16_GA) == 0)
    *bits = 32;
  else if (strcmp(type, UINT16_RGB) == 0)
    *bits = 48;
  else if (strcmp(type, UINT16_RGBA) == 0)
    *bits = 64;

  if (strcmp(type, INT32) == 0)
    *bits = 32;
  else if (strcmp(type, INT32_GA) == 0)
    *bits = 64;
  else if (strcmp(type, INT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(type, INT32_RGBA) == 0)
    *bits = 128;

  if (strcmp(type, UINT32) == 0)
    *bits = 32;
  else if (strcmp(type, UINT32_GA) == 0)
    *bits = 64;
  else if (strcmp(type, UINT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(type, UINT32_RGBA) == 0)
    *bits = 128;

  if (strcmp(type, INT64) == 0)
    *bits = 64;
  else if (strcmp(type, INT64_GA) == 0)
    *bits = 128;
  else if (strcmp(type, INT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(type, INT64_RGBA) == 0)
    *bits = 256;

  if (strcmp(type, UINT64) == 0)
    *bits = 64;
  else if (strcmp(type, UINT64_GA) == 0)
    *bits = 128;
  else if (strcmp(type, UINT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(type, UINT64_RGBA) == 0)
    *bits = 256;

  if (strcmp(type, FLOAT32) == 0)
    *bits = 32;
  else if (strcmp(type, FLOAT32_GA) == 0)
    *bits = 64;
  else if (strcmp(type, FLOAT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(type, FLOAT32_RGBA) == 0)
    *bits = 128;

  if (strcmp(type, FLOAT64) == 0)
    *bits = 64;
  else if (strcmp(type, FLOAT64_GA) == 0)
    *bits = 128;
  else if (strcmp(type, FLOAT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(type, FLOAT64_RGBA) == 0)
    *bits = 256;
  else
    *bits = 0;

  return PIDX_success;
}
