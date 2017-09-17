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


PIDX_return_code PIDX_create_meta_data_cache(PIDX_meta_data_cache* cache)
{

  *cache = malloc(sizeof (*(*cache)));
  memset(*cache, 0, sizeof (*(*cache)));

  return PIDX_success;
}


PIDX_return_code PIDX_free_meta_data_cache(PIDX_meta_data_cache cache)
{
  free(cache);
  return PIDX_success;
}
