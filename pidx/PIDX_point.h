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
#ifndef __PIDX_POINT_H
#define __PIDX_POINT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "PIDX_error_codes.h"
#include "Generic_data_structs.h"

extern const int PIDX_default_bits_per_block;
extern const int PIDX_default_blocks_per_file;

typedef int PIDX_point[5];

/// Utility functions to set or get the dimensions of an offset or a box (defined as points)

PIDX_return_code PIDX_create_point(PIDX_point* point);
PIDX_return_code PIDX_delete_point(PIDX_point* point);
PIDX_return_code PIDX_set_point_1D(int  x, PIDX_point point);
PIDX_return_code PIDX_get_point_1D(int* x, PIDX_point point);
PIDX_return_code PIDX_set_point_2D(int  x,int  y, PIDX_point point);
PIDX_return_code PIDX_get_point_2D(int* x,int* y, PIDX_point point);
PIDX_return_code PIDX_set_point_3D(int  x,int  y,int  z, PIDX_point point);
PIDX_return_code PIDX_get_point_3D(int* x,int* y,int* z, PIDX_point point);
PIDX_return_code PIDX_set_point_4D(int  x,int  y,int  z,int  u, PIDX_point point);
PIDX_return_code PIDX_get_point_4D(int* x,int* y,int* z,int* u, PIDX_point point);
PIDX_return_code PIDX_set_point_5D(int  x,int  y,int  z,int  u,int  v, PIDX_point point);
PIDX_return_code PIDX_get_point_5D(int* x,int* y,int* z,int* u,int* v, PIDX_point point);

#ifdef __cplusplus
}
#endif

#endif
