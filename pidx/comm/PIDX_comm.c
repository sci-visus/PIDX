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

PIDX_return_code PIDX_create_access(PIDX_access* access)
{
  *access = (PIDX_access)malloc(sizeof (*(*access)));
  memset(*access, 0, sizeof (*(*access)));
  
  (*access)->comm = MPI_COMM_NULL;
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_mpi_access(PIDX_access access, MPI_Comm comm)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->comm = comm;

  return PIDX_success;
}

PIDX_return_code PIDX_close_access(PIDX_access access)
{
  if(access == NULL)
    return PIDX_err_access;
  
  free(access);
  
  return PIDX_success;
}
