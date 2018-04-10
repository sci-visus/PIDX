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


# This module finds if NETCDF is installed, and sets the following variables
# indicating where it is:
#
# NETCDF_FOUND               - system has NETCDF
# NETCDF_INCLUDE_DIR         - path to NETCDF include directory
# NETCDF_LIBRARIES           - all NETCDF libraries
#
# Execute cmake with "-DNETCDF_DIR=/path/to/NetCDF" to help find the library.
#

FIND_PATH(NETCDF_INCLUDE_DIR NAMES netcdf_par.h PATHS ${NETCDF_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(NETCDF_INCLUDE_DIR NAMES netcdf_par.h)

IF (NETCDF_INCLUDE_DIR)
  FIND_PACKAGE(MPI REQUIRED)
  IF (MPI_C_FOUND)
    SET(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
    SET(NETCDF_MPI_LIBS ${MPI_C_LIBRARIES})
  ELSE()
    MESSAGE("NETCDF requires MPI, but we could not locate MPI on your system!")
  ENDIF()
    
  FIND_PACKAGE(HDF5 REQUIRED)
  IF (HDF5_FOUND)
    SET(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR} ${HDF5_INCLUDE_PATH})
    SET(NETCDF_HDF5_LIBS ${HDF5_LIBRARIES})
  ELSE()
    MESSAGE("NETCDF requires HDF5, but we could not locate MPI on your system!")
  ENDIF()
    
  SET(NETCDF_INCLUDE_DIRS "${NETCDF_INCLUDE_DIR}"             CACHE INTERNAL "" FORCE)

  FIND_LIBRARY(NETCDF_LIBRARY     netcdf    PATHS ${NETCDF_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(NETCDF_LIBRARY     netcdf    )

  SET(NETCDF_LIBRARIES 
    ${NETCDF_LIBRARY}
    ${NETCDF_MPI_LIBS}
    ${NETCDF_HDF5_LIBS}
    )

  SET(NETCDF_FOUND "YES" CACHE BOOL "NETCDF library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using NETCDF_INCLUDE_DIR  = ") 
    FOREACH(inc ${NETCDF_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found NETCDF_LIBRARIES    = ")
    FOREACH(lib ${NETCDF_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (NETCDF_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "NETCDF library not found. Try setting NETCDF_DIR and ensure NetCDF is compiled with parallel support")
  ELSE()
    MESSAGE("WARNING: NETCDF library not found. Try setting NETCDF_DIR and ensure NetCDF is compiled with parallel support")
  ENDIF()
ENDIF ()

