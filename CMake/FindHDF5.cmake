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

