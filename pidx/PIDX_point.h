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

/**
 * \file PIDX_point.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Functions and datatype for point datatype
 *
 */

#ifndef __PIDX_POINT_H
#define __PIDX_POINT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "PIDX_error_codes.h"
#include "Generic_data_structs.h"

//TODO: this should be longs
typedef long long PIDX_point[5];

/// Utility functions to set or get the dimensions of an offset or a box (defined as points)
PIDX_return_code PIDX_set_point_1D(long long  x, PIDX_point point);
PIDX_return_code PIDX_get_point_1D(long long* x, PIDX_point point);
PIDX_return_code PIDX_set_point_2D(long long  x,long long  y, PIDX_point point);
PIDX_return_code PIDX_get_point_2D(long long* x,long long* y, PIDX_point point);
PIDX_return_code PIDX_set_point_3D(long long  x,long long  y,long long  z, PIDX_point point);
PIDX_return_code PIDX_get_point_3D(long long* x,long long* y,long long* z, PIDX_point point);
PIDX_return_code PIDX_set_point_4D(long long  x,long long  y,long long  z,long long  u, PIDX_point point);
PIDX_return_code PIDX_get_point_4D(long long* x,long long* y,long long* z,long long* u, PIDX_point point);
PIDX_return_code PIDX_set_point_5D(long long  x,long long  y,long long  z,long long  u,long long  v, PIDX_point point);
PIDX_return_code PIDX_get_point_5D(long long* x,long long* y,long long* z,long long* u,long long* v, PIDX_point point);

PIDX_return_code PIDX_inner_product(PIDX_point point, long long *inner_product);

#ifdef __cplusplus
}
#endif

#endif