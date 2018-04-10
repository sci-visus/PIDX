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


# This module finds if PNETCDF is installed, and sets the following variables
# indicating where it is:
#
# PNETCDF_FOUND               - system has PNETCDF
# PNETCDF_INCLUDE_DIR         - path to PNETCDF include directory
# PNETCDF_LIBRARIES           - all PNETCDF libraries
#
# Execute cmake with "-DPNETCDF_DIR=/path/to/PnetCDF" to help find the library.
#

FIND_PATH(PNETCDF_INCLUDE_DIR NAMES pnetcdf.h PATHS ${PNETCDF_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(PNETCDF_INCLUDE_DIR NAMES pnetcdf.h)

IF (PNETCDF_INCLUDE_DIR)
  FIND_PACKAGE(MPI REQUIRED)
  IF (MPI_C_FOUND)
    SET(PNETCDF_INCLUDE_DIR ${PNETCDF_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
    SET(PNETCDF_MPI_LIBS ${MPI_C_LIBRARIES})
  ELSE()
    MESSAGE("PNETCDF requires MPI, but we could not locate MPI on your system!")
  ENDIF()
    
  SET(PNETCDF_INCLUDE_DIRS "${PNETCDF_INCLUDE_DIR}"             CACHE INTERNAL "" FORCE)

  FIND_LIBRARY(PNETCDF_LIBRARY     pnetcdf    PATHS ${PNETCDF_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(PNETCDF_LIBRARY     pnetcdf    )

  SET(PNETCDF_LIBRARIES 
    ${PNETCDF_LIBRARY}
    ${PNETCDF_MPI_LIBS}
    )

  SET(PNETCDF_FOUND "YES" CACHE BOOL "PNETCDF library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using PNETCDF_INCLUDE_DIR  = ") 
    FOREACH(inc ${PNETCDF_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found PNETCDF_LIBRARIES    = ")
    FOREACH(lib ${PNETCDF_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (PNETCDF_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "PNETCDF library not found. Try setting PNETCDF_DIR")
  ELSE()
    MESSAGE("WARNING: PNETCDF library not found. Try setting PNETCDF_DIR")
  ENDIF()
ENDIF ()

