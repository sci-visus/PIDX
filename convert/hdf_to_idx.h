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

#ifndef _HDF_TO_IDX_H
#define _HDF_TO_IDX_H

#include "PIDX.h"

#if PIDX_OPTION_HDF5

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "hdf5.h"

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

#define H5FILE_NAME     "/home/sid/data/Post2D_h2air_flame_4.0000E-03/Solution.h5"

int main(int argc, char **argv);

#endif
#endif