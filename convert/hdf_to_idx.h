#ifndef _HDF_TO_IDX_H
#define _HDF_TO_IDX_H

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

#define H5FILE_NAME     "/home/sid/research/KAUST_data/Post2D_h2air_flame_0.0000E+00/Solution.h5"

int main(int argc, char **argv);

#endif