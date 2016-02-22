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

PIDX_return_code PIDX_create_access(PIDX_access* access)
{
  *access = (PIDX_access)malloc(sizeof (*(*access)));
  memset(*access, 0, sizeof (*(*access)));
  
  (*access)->parallel = 0;

  (*access)->idx_count[0] = 1;
  (*access)->idx_count[1] = 1;
  (*access)->idx_count[2] = 1;

#if PIDX_HAVE_MPI
  (*access)->comm = MPI_COMM_NULL;
#endif
  
  return PIDX_success;
}

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_set_mpi_access(PIDX_access access, MPI_Comm comm)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->parallel = 1;

  int ret;
  ret = MPI_Comm_dup(comm, &(access->comm));
  if (ret != MPI_SUCCESS)
    return PIDX_err_access;
  
  return PIDX_success;
}


PIDX_return_code PIDX_set_idx_count(PIDX_access access, int idx_count_x, int idx_count_y, int idx_count_z)
{
  if(access == NULL)
    return PIDX_err_access;

  access->idx_count[0] = idx_count_x;
  access->idx_count[1] = idx_count_y;
  access->idx_count[2] = idx_count_z;

  return PIDX_success;
}


PIDX_return_code PIDX_set_process_extent(PIDX_access access, int sub_div_x, int sub_div_y, int sub_div_z)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->sub_div[0] = sub_div_x;
  access->sub_div[1] = sub_div_y;
  access->sub_div[2] = sub_div_z;
  
  return PIDX_success;
}

/*
PIDX_return_code PIDX_enable_topology_aware_io(PIDX_access access, int topology_io)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->topology_aware_io = topology_io;
  
  return PIDX_success;
}
*/

PIDX_return_code PIDX_set_process_rank_decomposition(PIDX_access access, int rank_x, int rank_y, int rank_z)
{
  if(access == NULL)
    return PIDX_err_access;
  
  access->rank_component[0] = rank_x;
  access->rank_component[1] = rank_y;
  access->rank_component[2] = rank_z;
  
  return PIDX_success;
}
#endif

PIDX_return_code PIDX_close_access(PIDX_access access)
{
  if(access == NULL)
    return PIDX_err_access;
  
#if PIDX_HAVE_MPI
  if (access->comm != MPI_COMM_NULL)
  {
    if (access->parallel == 1)
      MPI_Comm_free(&(access->comm));
    else
      fprintf(stderr, "PIDX ERROR: Trying to free a NULL communicator. Application must specify MPI communicator using PIDX_set_mpi_access.\n");
  }

#endif  
  
  free(access);
  
  return PIDX_success;
}
