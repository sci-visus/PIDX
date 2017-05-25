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

IF (PIDX_BUILD_TOOLS)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  SET(IDXCOMPARE_SOURCES idx-compare.c)
  SET(DUMPHEADER_SOURCES idx-dump-header.c)
  SET(IDXVERIFY_SOURCES idx-verify.c)
  SET(IDXMINMAX_SOURCES idx-minmax.c)

  SET(TOOLS_LINK_LIBS pidx ${PIDX_LINK_LIBS})
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring utility with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring profiling with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  
  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(idxcompare "${IDXCOMPARE_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxdumpheader "${DUMPHEADER_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxverify "${IDXVERIFY_SOURCES}")

  PIDX_ADD_CEXECUTABLE(minmax "${IDXMINMAX_SOURCES}")
  
  TARGET_LINK_LIBRARIES(idxverify m ${TOOLS_LINK_LIBS})
  TARGET_LINK_LIBRARIES(minmax ${TOOLS_LINK_LIBS})

ENDIF()
