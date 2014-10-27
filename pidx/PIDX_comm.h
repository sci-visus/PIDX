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
 * \file PIDX_comm.h
 *
 * \mainpage
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 * \author Giorgio Scorzelli
 * \author Valerio Pascucci
 * \date   10/09/14
 *
 * PIDX is an I/O library that enables HPC applications to write distributed 
 * multi-dimensional data directly into a hierarchical multi-resolution 
 * data format (IDX) with minimal overhead. 
 *
 */
 
#ifndef __PIDX_COMM_H
#define __PIDX_COMM_H 

#include "PIDX_error_codes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if PIDX_HAVE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct PIDX_access_struct
{
  int parallel;
  
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
};
typedef struct PIDX_access_struct* PIDX_access;

PIDX_return_code PIDX_create_access(PIDX_access* access);
PIDX_return_code PIDX_close_access(PIDX_access* access);

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_set_mpi_access(PIDX_access access, MPI_Comm comm);
#endif

#ifdef __cplusplus
} //extern C
#endif

#endif
