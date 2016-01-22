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


# This module finds if ZFP is installed, and sets the following variables
# indicating where it is:
#
# ZFP_FOUND               - system has ZFP
# ZFP_INCLUDE_DIR         - path to ZFP include directory
# ZFP_LIBRARIES           - all ZFP libraries
#
# Execute cmake with "-DZFP_DIR=/path/to/zfp" to help find the library.
#

FIND_PATH(ZFP_INCLUDE_DIR NAMES zfp.h PATHS ${ZFP_DIR}/inc ${ZFP_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(ZFP_INCLUDE_DIR NAMES zfp.h)

IF (ZFP_INCLUDE_DIR)

  SET(ZFP_INCLUDE_DIRS "${ZFP_INCLUDE_DIR}")

  FIND_LIBRARY(ZFP_LIBRARY     zfp    PATHS ${ZFP_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(ZFP_LIBRARY     zfp    )

  SET(ZFP_LIBRARIES 
    ${ZFP_LIBRARY}
    stdc++
    )

  SET(ZFP_FOUND "YES" CACHE BOOL "ZFP library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using ZFP_INCLUDE_DIR  = ") 
    FOREACH(inc ${ZFP_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found ZFP_LIBRARIES    = ")
    FOREACH(lib ${ZFP_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (ZFP_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ZFP library not found. Try setting ZFP_DIR")
  ELSE()
    MESSAGE("WARNING: ZFP library not found. Try setting ZFP_DIR")
  ENDIF()
ENDIF ()

