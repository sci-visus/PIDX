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

#include <ctype.h>

#include "PIDX_file_handler.h"

struct pidx_dtype PIDX_DType = {
        .INT8       = "1*int8",
        .INT8_GA    = "2*int8",
        .INT8_RGB   = "3*int8",
        .INT8_RGBA  = "4*int8",

        .UINT8       = "1*uint8",
        .UINT8_GA    = "2*uint8",
        .UINT8_RGB   = "3*uint8",
        .UINT8_RGBA  = "4*uint8",

        .INT16       = "1*int16",
        .INT16_GA    = "2*int16",
        .INT16_RGB   = "3*int16",
        .INT16_RGBA  = "4*int16",

        .UINT16       = "1*uint16",
        .UINT16_GA    = "2*uint16",
        .UINT16_RGB   = "3*uint16",
        .UINT16_RGBA  = "4*uint16",

        .INT32       = "1*int32",
        .INT32_GA    = "2*int32",
        .INT32_RGB   = "3*int32",
        .INT32_RGBA  = "4*int32",

        .UINT32       = "1*uint32",
        .UINT32_GA    = "2*uint32",
        .UINT32_RGB   = "3*uint32",
        .UINT32_RGBA  = "4*uint32",

        .INT64       = "1*int64",
        .INT64_GA    = "2*int64",
        .INT64_RGB   = "3*int64",
        .INT64_RGBA  = "4*int64",

        .UINT64       = "1*uint64",
        .UINT64_GA    = "2*uint64",
        .UINT64_RGB   = "3*uint64",
        .UINT64_RGBA  = "4*uint64",

        .FLOAT32       = "1*float32",
        .FLOAT32_GA    = "2*float32",
        .FLOAT32_RGB   = "3*float32",
        .FLOAT32_RGBA  = "4*float32",
        .FLOAT32_7STENCIL  = "7*float32",

        .FLOAT64       = "1*float64",
        .FLOAT64_GA    = "2*float64",
        .FLOAT64_RGB   = "3*float64",
        .FLOAT64_RGBA  = "4*float64",
        .FLOAT64_7STENCIL  = "7*float64",

        .FLOAT64_9TENSOR  = "9*float64",
        .FLOAT32_9TENSOR  = "9*float32",
        .INT64_9TENSOR  = "9*int64",
        .INT32_9TENSOR  = "9*int32"
};

PIDX_return_code PIDX_default_bits_per_datatype(PIDX_data_type type, int* bits)
{
  char full_typename[64];
  strncpy(full_typename, type, 64);
  // if typename does not start with the number of components change to "1*dataype"
  // TODO handle also case with datatype written as datatype[n_comp]
  if(!isdigit(full_typename[0])){
    sprintf(full_typename, "1*%s", type);
  }

  if (strcmp(full_typename, PIDX_DType.INT8) == 0)
    *bits = 8;
  else if (strcmp(full_typename, PIDX_DType.INT8_GA) == 0)
    *bits = 16;
  else if (strcmp(full_typename, PIDX_DType.INT8_RGB) == 0)
    *bits = 24;
  else if (strcmp(full_typename, PIDX_DType.INT8_RGBA) == 0)
    *bits = 32;

  else if (strcmp(full_typename, PIDX_DType.UINT8) == 0)
    *bits = 8;
  else if (strcmp(full_typename, PIDX_DType.UINT8_GA) == 0)
    *bits = 16;
  else if (strcmp(full_typename, PIDX_DType.UINT8_RGB) == 0)
    *bits = 24;
  else if (strcmp(full_typename, PIDX_DType.UINT8_RGBA) == 0)
    *bits = 32;

  else if (strcmp(full_typename, PIDX_DType.INT16) == 0)
    *bits = 16;
  else if (strcmp(full_typename, PIDX_DType.INT16_GA) == 0)
    *bits = 32;
  else if (strcmp(full_typename, PIDX_DType.INT16_RGB) == 0)
    *bits = 48;
  else if (strcmp(full_typename, PIDX_DType.INT16_RGBA) == 0)
    *bits = 64;

  else if (strcmp(full_typename, PIDX_DType.UINT16) == 0)
    *bits = 16;
  else if (strcmp(full_typename, PIDX_DType.UINT16_GA) == 0)
    *bits = 32;
  else if (strcmp(full_typename, PIDX_DType.UINT16_RGB) == 0)
    *bits = 48;
  else if (strcmp(full_typename, PIDX_DType.UINT16_RGBA) == 0)
    *bits = 64;

  else if (strcmp(full_typename, PIDX_DType.INT32) == 0)
    *bits = 32;
  else if (strcmp(full_typename, PIDX_DType.INT32_GA) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.INT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(full_typename, PIDX_DType.INT32_RGBA) == 0)
    *bits = 128;

  else if (strcmp(full_typename, PIDX_DType.UINT32) == 0)
    *bits = 32;
  else if (strcmp(full_typename, PIDX_DType.UINT32_GA) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.UINT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(full_typename, PIDX_DType.UINT32_RGBA) == 0)
    *bits = 128;

  else if (strcmp(full_typename, PIDX_DType.INT64) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.INT64_GA) == 0)
    *bits = 128;
  else if (strcmp(full_typename, PIDX_DType.INT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(full_typename, PIDX_DType.INT64_RGBA) == 0)
    *bits = 256;

  else if (strcmp(full_typename, PIDX_DType.UINT64) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.UINT64_GA) == 0)
    *bits = 128;
  else if (strcmp(full_typename, PIDX_DType.UINT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(full_typename, PIDX_DType.UINT64_RGBA) == 0)
    *bits = 256;

  else if (strcmp(full_typename, PIDX_DType.FLOAT32) == 0)
    *bits = 32;
  else if (strcmp(full_typename, PIDX_DType.FLOAT32_GA) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.FLOAT32_RGB) == 0)
    *bits = 96;
  else if (strcmp(full_typename, PIDX_DType.FLOAT32_RGBA) == 0)
    *bits = 128;
  else if (strcmp(full_typename, PIDX_DType.FLOAT32_7STENCIL) == 0)
    *bits = 224;
  else if (strcmp(full_typename, PIDX_DType.FLOAT32_9TENSOR) == 0)
    *bits = 288;

  else if (strcmp(full_typename, PIDX_DType.FLOAT64) == 0)
    *bits = 64;
  else if (strcmp(full_typename, PIDX_DType.FLOAT64_GA) == 0)
    *bits = 128;
  else if (strcmp(full_typename, PIDX_DType.FLOAT64_RGB) == 0)
    *bits = 192;
  else if (strcmp(full_typename, PIDX_DType.FLOAT64_RGBA) == 0)
    *bits = 256;
  else if (strcmp(full_typename, PIDX_DType.FLOAT64_7STENCIL) == 0)
    *bits = 448;
  else if (strcmp(full_typename, PIDX_DType.FLOAT64_9TENSOR) == 0)
    *bits = 576;
  else
    *bits = 0;

  return PIDX_success;
}

PIDX_return_code PIDX_values_per_datatype(PIDX_data_type type, int* values, int* bits)
{
    char full_typename[64];
    strncpy(full_typename, type, 64);
    if(!isdigit(full_typename[0])){
      sprintf(full_typename, "1*%s", type);
    }

    if (strcmp(full_typename, PIDX_DType.INT8) == 0)
    {
      *bits = 8;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.INT8_GA) == 0)
    {
      *bits = 8;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.INT8_RGB) == 0)
    {
      *bits = 8;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.INT8_RGBA) == 0)
    {
      *bits = 8;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.UINT8) == 0)
    {
      *bits = 8;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT8_GA) == 0)
    {
      *bits = 8;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT8_RGB) == 0)
    {
      *bits = 8;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT8_RGBA) == 0)
    {
      *bits = 8;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.INT16) == 0)
    {
      *bits = 16;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.INT16_GA) == 0)
    {
      *bits = 16;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.INT16_RGB) == 0)
    {
      *bits = 16;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.INT16_RGBA) == 0)
    {
      *bits = 16;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.UINT16) == 0)
    {
      *bits = 16;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT16_GA) == 0)
    {
      *bits = 16;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT16_RGB) == 0)
    {
      *bits = 16;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT16_RGBA) == 0)
    {
      *bits = 16;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.INT32) == 0)
    {
      *bits = 32;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.INT32_GA) == 0)
    {
      *bits = 32;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.INT32_RGB) == 0)
    {
      *bits = 32;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.INT32_RGBA) == 0)
    {
      *bits = 32;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.UINT32) == 0)
    {
      *bits = 32;
      *values = 1;
    }
    else if (strcmp(type, PIDX_DType.UINT32_GA) == 0)
    {
      *bits = 32;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT32_RGB) == 0)
    {
      *bits = 32;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT32_RGBA) == 0)
    {
      *bits = 32;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.INT64) == 0)
    {
      *bits = 64;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.INT64_GA) == 0)
    {
      *bits = 64;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.INT64_RGB) == 0)
    {
      *bits = 64;
      *values = 3;
    }
    else if (strcmp(full_typename, PIDX_DType.INT64_RGBA) == 0)
    {
      *bits = 64;
      *values = 4;
    }

    else if (strcmp(full_typename, PIDX_DType.UINT64) == 0)
    {
      *bits = 64;
      *values = 1;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT64_GA) == 0)
    {
      *bits = 64;
      *values = 2;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT64_RGB) == 0)
    {
      *values = 3;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.UINT64_RGBA) == 0)
    {
      *values = 4;
      *bits = 64;
    }

    else if (strcmp(full_typename, PIDX_DType.FLOAT32) == 0)
    {
      *values = 1;
      *bits = 32;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT32_GA) == 0)
    {
      *values = 2;
      *bits = 32;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT32_RGB) == 0)
    {
      *values = 3;
      *bits = 32;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT32_RGBA) == 0)
    {
      *values = 4;
      *bits = 32;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT32_7STENCIL) == 0)
    {
      *values = 7;
      *bits = 32;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT32_9TENSOR) == 0)
    {
      *values = 9;
      *bits = 32;
    }

    else if (strcmp(full_typename, PIDX_DType.FLOAT64) == 0)
    {
      *values = 1;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT64_GA) == 0)
    {
      *values = 2;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT64_RGB) == 0)
    {
      *values = 3;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT64_RGBA) == 0)
    {
      *values = 4;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT64_7STENCIL) == 0)
    {
      *values = 7;
      *bits = 64;
    }
    else if (strcmp(full_typename, PIDX_DType.FLOAT64_9TENSOR) == 0)
    {
      *values = 9;
      *bits = 64;
    }
    else
      *values = 0;

    return PIDX_success;
}
