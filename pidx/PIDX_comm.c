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

#include "PIDX_comm.h"

PIDX_return_code PIDX_create_access(PIDX_access* access)
{
  *access = malloc(sizeof (*(*access)));
  memset(*access, 0, sizeof (*(*access)));
  
  (*access)->parallel = 0;
  
#if PIDX_HAVE_MPI
  (*access)->comm = NULL;
#endif
  
  return PIDX_success;
}

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_set_mpi_access(PIDX_access access, MPI_Comm comm)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->parallel = 1;
  
  MPI_Comm_dup(comm, &(access->comm));
  
  return PIDX_success;
}
#endif

PIDX_return_code PIDX_close_access(PIDX_access* access)
{
  if(*access == NULL)
    return PIDX_err_access;
  
#if PIDX_HAVE_MPI
  MPI_Comm_free(&((*access)->comm));
#endif  
  
  free(*access);
  *access = 0;
  
  return PIDX_success;
}
