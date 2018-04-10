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
