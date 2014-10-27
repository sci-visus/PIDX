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

IF (PIDX_INCLUDE_DIR)

  # Read pidx_config.h to see if we need to link MPI
  FILE(STRINGS "${PIDX_INCLUDE_DIR}/pidx_config.h" config)
  LIST(FIND config "PIDX_HAVE_MPI 1" idx)
  IF (IDX GREATER 0)
    FIND_PACKAGE(MPI REQUIRED)
    IF (MPI_C_FOUND)
      SET(PIDX_INCLUDE_DIR ${PIDX_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
      SET(PIDX_MPI_LIBS ${MPI_C_LIBRARIES})
    ELSE()
      MESSAGE("PIDX was configured with MPI support, but we could not locate MPI on your system!")
    ENDIF()
  ENDIF()
    
  FIND_LIBRARY(PIDX_LIB     pidx    ${PIDX_DIR}/lib)

  SET(PIDX_LIBRARIES 
    ${PIDX_LIB}
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
  MESSAGE("ERROR PIDX library not found on the system")
ENDIF ()
                         
