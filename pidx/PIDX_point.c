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

#include "PIDX_point.h"

PIDX_return_code PIDX_create_point(PIDX_point* point)
{
  //*point = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
  //memset(*point, 0, sizeof(int) * PIDX_MAX_DIMENSIONS);
  
  return PIDX_success;
}

PIDX_return_code PIDX_delete_point(PIDX_point* point)
{
  //if(*point == NULL)
  //  return PIDX_err_point;
  
  //free(*point);
  //*point = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_point_1D(int x, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = 0;
  point[2] = 0;
  point[3] = 0;
  point[4] = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_point_1D(int* x, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_point_2D(int  x, int  y, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = 0;
  point[3] = 0;
  point[4] = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_point_2D(int* x, int* y, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_point_3D(int  x, int  y, int  z, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = z;
  point[3] = 0;
  point[4] = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_point_3D(int* x, int* y, int* z, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  *z = point[2];
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_point_4D(int  x, int  y, int  z, int u, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  point[0] = x;
  point[1] = y;
  point[2] = z;
  point[3] = u;
  point[4] = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_point_4D(int* x, int* y, int* z, int* u, PIDX_point point)
{
  if(point == NULL)
    return PIDX_err_point;

  *x = point[0];
  *y = point[1];
  *z = point[2];
  *u = point[3];
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_point_5D(int  x,int  y,int  z,int  u,int  v, PIDX_point point)
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

PIDX_return_code PIDX_get_point_5D(int* x, int* y, int* z, int* u, int* v, PIDX_point point)
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