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

IF (WIN32)
    MESSAGE(ERROR, "Unsupported platform. PIDX is designed for Unix.")
ELSE()
  SITE_NAME(HOSTNAME)
	MESSAGE("Configuring PIDX for ${CMAKE_SYSTEM_NAME} (${CMAKE_SYSTEM})")

  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
  endif()

  IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    ADD_DEFINITIONS(-D_DEBUG=1)
  ENDIF()

  IF (MPI_CXX_FOUND OR MPI_C_FOUND)
     MESSAGE("MPI_ENABLED ${ENABLE_MPI}")
     ADD_DEFINITIONS(-DPIDX_HAVE_MPI)
  ENDIF()

  # ///////////////////////////////////////////////
  # HOST-SPECIFIC CONFIGURATIONS
  # ///////////////////////////////////////////////

  IF (0) # NOTE: you can use regexp here, e.g. IF (${HOSTNAME} MATCHES "hopper.*")
    MESSAGE(ERROR "UNKNOWN SYSTEM")  
  ELSE()
    ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE)
    SET(CMAKE_C_FLAGS    ${CMAKE_C_FLAGS}   "-fPIC -Wall")
    SET(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-fPIC -Wall")
  ENDIF()

ENDIF()
