#/*****************************************************
# **  PIDX Parallel I/O Library                      **
# **  Copyright (c) 2010-2014 University of Utah     **
# **  Scientific Computing and Imaging Institute     **
# **  72 S Central Campus Drive, Room 3750           **
# **  Salt Lake City, UT 84112                       **
# **                                                 **
# **  PIDX is licensed under the Creative Commons    **
# **  Attribution-NonCommercial-NoDerivatives 4.0    **
# **  International License. See LICENSE.md.         **
# **                                                 **
# **  For information about this project see:        **
# **  http://www.cedmav.com/pidx                     **
# **  or contact: pascucci@sci.utah.edu              **
# **  For support: PIDX-support@visus.net            **
# **                                                 **
# *****************************************************/


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

