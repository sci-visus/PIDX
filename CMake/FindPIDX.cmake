##
## BSD 3-Clause License
## 
## Copyright (c) 2010-2018 ViSUS L.L.C., 
## Scientific Computing and Imaging Institute of the University of Utah
## 
## ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
## University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
##  
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## * Redistributions of source code must retain the above copyright notice, this
## list of conditions and the following disclaimer.
## 
## * Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
## 
## * Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## For additional information about this project contact: pascucci@acm.org
## For support: support@visus.net
## 
##


# This module finds if PIDX is installed, and sets the following variables
# indicating where it is:
#
# PIDX_FOUND               - system has PIDX
# PIDX_INCLUDE_DIR         - path to PIDX include directory
# PIDX_LIBRARIES           - all PIDX libraries
#
# Execute cmake with "-DPIDX_DIR=/path/to/pidx" to help find the library.
#

FIND_PATH(PIDX_INCLUDE_DIR NAMES PIDX.h PATHS ${PIDX_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(HDF5_INCLUDE_DIR NAMES PIDX.h)

IF (PIDX_INCLUDE_DIR)

  # Read pidx_config.h to see if we need to link MPI
  FILE(STRINGS "${PIDX_INCLUDE_DIR}/PIDX_config.h" config)
  LIST(FIND config "#define PIDX_HAVE_MPI 1" have_mpi)
  IF (have_mpi GREATER 0)
    FIND_PACKAGE(MPI REQUIRED)
    IF (MPI_C_FOUND)
      SET(PIDX_INCLUDE_DIR ${PIDX_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
      SET(PIDX_MPI_LIBS ${MPI_C_LIBRARIES})
    ELSE()
      MESSAGE("PIDX was configured with MPI support, but we could not locate MPI on your system!")
    ENDIF()
  ENDIF()
    
  SET(PIDX_INCLUDE_DIRS "${PIDX_INCLUDE_DIR}")

  FIND_LIBRARY(PIDX_LIBRARY     pidx    PATHS ${PIDX_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(HDF5_LIBRARY     pidx    )

  SET(PIDX_LIBRARIES 
    ${PIDX_LIBRARY}
    ${PIDX_MPI_LIBS}
    )

  SET(PIDX_FOUND "YES" CACHE BOOL "PIDX library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using PIDX_INCLUDE_DIR  = ") 
    FOREACH(inc ${PIDX_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found PIDX_LIBRARIES    = ")
    FOREACH(lib ${PIDX_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (PIDX_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "PIDX library not found. Try setting PIDX_DIR")
  ELSE()
    MESSAGE("WARNING: PIDX library not found. Try setting PIDX_DIR")
  ENDIF()
ENDIF ()
                         
