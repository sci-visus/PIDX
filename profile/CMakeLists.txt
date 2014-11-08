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

IF (PIDX_BUILD_PROFILE)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  SET(IDXCOMPARE_SOURCES idx-compare.c)
  SET(DUMPHEADER_SOURCES dump-header.c)


  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  
  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

  PIDX_ADD_EXECUTABLE(idxcompare "${IDXCOMPARE_SOURCES}")
  ADD_DEPENDENCIES(idxcompare pidx)
  TARGET_LINK_LIBRARIES(idxcompare pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(dumpheader "${DUMPHEADER_SOURCES}")
  ADD_DEPENDENCIES(dumpheader pidx)
  TARGET_LINK_LIBRARIES(dumpheader pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

ENDIF()
