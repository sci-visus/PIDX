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

PIDX_type INT8          = "int8";
PIDX_type INT8_GA       = "int8[2]";
PIDX_type INT8_RGB      = "int8[3]";
PIDX_type INT8_RGBA     = "int8[4]";

PIDX_type UINT8         = "uint8";
PIDX_type UINT8_GA      = "uint8[2]";
PIDX_type UINT8_RGB     = "uint8[3]";
PIDX_type UINT8_RGBA    = "uint8[4]";

PIDX_type INT16         = "int16";
PIDX_type INT16_GA      = "int16[2]";
PIDX_type INT16_RGB     = "int16[3]";
PIDX_type INT16_RGBA    = "int16[4]";

PIDX_type UINT16        = "uint16";
PIDX_type UINT16_GA     = "uint16[2]";
PIDX_type UINT16_RGB    = "uint16[3]";
PIDX_type UINT16_RGBA   = "uint16[4]";

PIDX_type INT32         = "int32";
PIDX_type INT32_GA      = "int32[2]";
PIDX_type INT32_RGB     = "int32[3]";
PIDX_type INT32_RGBA    = "int32[4]";

PIDX_type UINT32        = "uint32";
PIDX_type UINT32_GA     = "uint32[2]";
PIDX_type UINT32_RGB    = "uint32[3]";
PIDX_type UINT32_RGBA   = "uint32[4]";

PIDX_type INT64         = "int64";
PIDX_type INT64_GA      = "int64[2]";
PIDX_type INT64_RGB     = "int64[3]";
PIDX_type INT64_RGBA    = "int64[4]";

PIDX_type UINT64        = "uint64";
PIDX_type UINT64_GA     = "uint64[2]";
PIDX_type UINT64_RGB    = "uint64[3]";
PIDX_type UINT64_RGBA   = "uint64[4]";

PIDX_type FLOAT32       = "float32";
PIDX_type FLOAT32_GA    = "float32[2]";
PIDX_type FLOAT32_RGB   = "float32[3]";
PIDX_type FLOAT32_RGBA  = "float32[4]";

PIDX_type FLOAT64       = "float64";
PIDX_type FLOAT64_GA    = "float64[2]";
PIDX_type FLOAT64_RGB   = "float64[3]";
PIDX_type FLOAT64_RGBA  = "float64[4]";

PIDX_return_code PIDX_default_bytes_per_datatype(PIDX_type type, int* bytes)
{
  if (strcmp(type, INT8) == 0)
    *bytes = 8;
  else if (strcmp(type, INT8_GA) == 0)
    *bytes = 16;
  else if (strcmp(type, INT8_RGB) == 0)
    *bytes = 24;
  else if (strcmp(type, INT8_RGBA) == 0)
    *bytes = 32;

  if (strcmp(type, UINT8) == 0)
    *bytes = 8;
  else if (strcmp(type, UINT8_GA) == 0)
    *bytes = 16;
  else if (strcmp(type, UINT8_RGB) == 0)
    *bytes = 24;
  else if (strcmp(type, UINT8_RGBA) == 0)
    *bytes = 32;


  if (strcmp(type, INT16) == 0)
    *bytes = 16;
  else if (strcmp(type, INT16_GA) == 0)
    *bytes = 32;
  else if (strcmp(type, INT16_RGB) == 0)
    *bytes = 48;
  else if (strcmp(type, INT16_RGBA) == 0)
    *bytes = 64;

  if (strcmp(type, UINT16) == 0)
    *bytes = 16;
  else if (strcmp(type, UINT16_GA) == 0)
    *bytes = 32;
  else if (strcmp(type, UINT16_RGB) == 0)
    *bytes = 48;
  else if (strcmp(type, UINT16_RGBA) == 0)
    *bytes = 64;


  if (strcmp(type, INT32) == 0)
    *bytes = 32;
  else if (strcmp(type, INT32_GA) == 0)
    *bytes = 64;
  else if (strcmp(type, INT32_RGB) == 0)
    *bytes = 96;
  else if (strcmp(type, INT32_RGBA) == 0)
    *bytes = 128;

  if (strcmp(type, UINT32) == 0)
    *bytes = 32;
  else if (strcmp(type, UINT32_GA) == 0)
    *bytes = 64;
  else if (strcmp(type, UINT32_RGB) == 0)
    *bytes = 96;
  else if (strcmp(type, UINT32_RGBA) == 0)
    *bytes = 128;

  if (strcmp(type, INT64) == 0)
    *bytes = 64;
  else if (strcmp(type, INT64_GA) == 0)
    *bytes = 128;
  else if (strcmp(type, INT64_RGB) == 0)
    *bytes = 192;
  else if (strcmp(type, INT64_RGBA) == 0)
    *bytes = 256;

  if (strcmp(type, UINT64) == 0)
    *bytes = 64;
  else if (strcmp(type, UINT64_GA) == 0)
    *bytes = 128;
  else if (strcmp(type, UINT64_RGB) == 0)
    *bytes = 192;
  else if (strcmp(type, UINT64_RGBA) == 0)
    *bytes = 256;


  if (strcmp(type, FLOAT32) == 0)
    *bytes = 32;
  else if (strcmp(type, FLOAT32_GA) == 0)
    *bytes = 64;
  else if (strcmp(type, FLOAT32_RGB) == 0)
    *bytes = 96;
  else if (strcmp(type, FLOAT32_RGBA) == 0)
    *bytes = 128;

  if (strcmp(type, FLOAT64) == 0)
    *bytes = 64;
  else if (strcmp(type, FLOAT64_GA) == 0)
    *bytes = 128;
  else if (strcmp(type, FLOAT64_RGB) == 0)
    *bytes = 192;
  else if (strcmp(type, FLOAT64_RGBA) == 0)
    *bytes = 256;
  else
    *bytes = 0;

  return PIDX_success;
}
