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

