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

typedef unsigned long long PIDX_point[PIDX_MAX_DIMENSIONS];

/// Utility functions to set or get the dimensions of an offset or a box (defined as points)
PIDX_return_code PIDX_set_point_5D(PIDX_point point, unsigned long long  x,unsigned long long  y,unsigned long long  z,unsigned long long  u,unsigned long long  v);
PIDX_return_code PIDX_get_point_5D(unsigned long long* x,unsigned long long* y,unsigned long long* z,unsigned long long* u,unsigned long long* v, PIDX_point point);

/// PIDX_inner_product
PIDX_return_code PIDX_inner_product(unsigned long long *inner_product, PIDX_point point);


#endif
