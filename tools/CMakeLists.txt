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
  # SET(HDF5VERIFY_SOURCES hdf5-verify.c)

  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  #INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  #INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  
  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(idxcompare "${IDXCOMPARE_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxdumpheader "${DUMPHEADER_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxverify "${IDXVERIFY_SOURCES}")
  
  IF (ZFP_FOUND)
    MESSAGE("Configuring profiling with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(PROFILE_LINK_LIBS ${PROFILE_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  TARGET_LINK_LIBRARIES(idxverify m ${PROFILE_LINK_LIBS})

ENDIF()
