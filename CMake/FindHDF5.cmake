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


# This module finds if HDF5 is installed, and sets the following variables
# indicating where it is:
#
# HDF5_FOUND               - system has HDF5
# HDF5_INCLUDE_DIR         - path to HDF5 include directory
# HDF5_LIBRARIES           - all HDF5 libraries
#
# Execute cmake with "-DHDF5_DIR=/path/to/hdf5" to help find the library.
#

FIND_PATH(HDF5_INCLUDE_DIR NAMES hdf5.h PATHS ${HDF5_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(HDF5_INCLUDE_DIR NAMES hdf5.h)

IF (HDF5_INCLUDE_DIR)
  FIND_PACKAGE(MPI REQUIRED)
  IF (MPI_C_FOUND)
    SET(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
    SET(HDF5_MPI_LIBS ${MPI_C_LIBRARIES})
  ELSE()
    MESSAGE("HDF5 requires MPI, but we could not locate MPI on your system!")
  ENDIF()
    
  SET(HDF5_INCLUDE_DIRS "${HDF5_INCLUDE_DIR}")

  FIND_LIBRARY(HDF5_LIBRARY     hdf5 hdf5-shared    PATHS ${HDF5_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(HDF5_LIBRARY     hdf5    )

  SET(HDF5_LIBRARIES 
    ${HDF5_LIBRARY}
    ${HDF5_MPI_LIBS}
    )

  SET(HDF5_FOUND "YES" CACHE BOOL "HDF5 library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using HDF5_INCLUDE_DIR  = ") 
    FOREACH(inc ${HDF5_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found HDF5_LIBRARIES    = ")
    FOREACH(lib ${HDF5_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (HDF5_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "HDF5 library not found. Try setting HDF5_DIR")
  ELSE()
    MESSAGE("WARNING: HDF5 library not found. Try setting HDF5_DIR")
  ENDIF()
ENDIF ()

