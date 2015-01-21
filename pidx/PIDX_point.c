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

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_point_5D(PIDX_point point, int64_t  x, int64_t  y, int64_t  z, int64_t  u, int64_t  v)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = z;
  point[3] = u;
  point[4] = u;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_point_5D(int64_t* x, int64_t* y, int64_t* z, int64_t* u, int64_t* v, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  *z = point[2];
  *u = point[3];
  *v = point[4];
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_inner_product(int64_t *inner_product, PIDX_point point)
{
  *inner_product = point[0] * point[1] * point[2] * point[3] * point[4];
  //safe_add, result=a+b overflow happens when (a+b)>MAX ---> b>MAX-a
  //if ((b>0?+b:-b)>(NumericLimits<T>::highest()-(a>0?+a:-a))) return false;
  //TODO: ensure there was no overflow
  return PIDX_success;
}
