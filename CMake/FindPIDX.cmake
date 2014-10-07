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

   FIND_PACKAGE(ZLIB REQUIRED)
   FIND_PACKAGE(MPI REQUIRED)

   FIND_LIBRARY(PIDX_LIB     pidx    ${PIDX_DIR}/lib)

   SET(PIDX_LIBRARIES 
       ${PIDX_LIB}
       ${ZLIB_LIBRARY}
       ${MPI_CXX_LIBRARIES}
   )

   SET(PIDX_FOUND "YES" CACHE BOOL "PIDX library found" FORCE)

   IF (CMAKE_VERBOSE_MAKEFILE)
      MESSAGE("Using PIDX_INCLUDE_DIR  = " ${PIDX_INCLUDE_DIR}) 
      MESSAGE("Found PIDX_LIBRARIES    = " ${PIDX_LIBRARIES}) 
      FOREACH(lib ${PIDX_LIBRARIES})
         MESSAGE("is: " ${lib})
      ENDFOREACH()
   ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
   MESSAGE("ERROR PIDX library not found on the system")
ENDIF ()
                         
