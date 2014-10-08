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

IF (PIDX_BUILD_TEST)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  FILE(GLOB PIDXTEST_SOURCES *.h *.c)


  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})


  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  PIDX_ADD_EXECUTABLE(pidxtest "${PIDXTEST_SOURCES}")
  ADD_DEPENDENCIES(pidxtest pidx)
  TARGET_LINK_LIBRARIES(pidxtest pidx ${ZLIB_LIBRARY} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})


ENDIF()
