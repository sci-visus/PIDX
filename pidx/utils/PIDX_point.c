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

#include "../PIDX_inc.h"

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_point(PIDX_point point, unsigned long long  x, unsigned long long  y, unsigned long long  z)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = z;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_point(unsigned long long* x, unsigned long long* y, unsigned long long* z, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  *z = point[2];
  
  return PIDX_success;
}


PIDX_return_code PIDX_set_physical_point(PIDX_physical_point point, double  x, double  y, double  z)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = z;

  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_physical_point(double* x, double* y, double* z, PIDX_physical_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  *z = point[2];

  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_inner_product(unsigned long long *inner_product, PIDX_point point)
{
  *inner_product = point[0] * point[1] * point[2];
  //safe_add, result=a+b overflow happens when (a+b)>MAX ---> b>MAX-a
  //if ((b>0?+b:-b)>(NumericLimits<T>::highest()-(a>0?+a:-a))) return false;
  //TODO: ensure there was no overflow
  return PIDX_success;
}
