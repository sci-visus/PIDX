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


IF (PIDX_BUILD_UTILITY)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(UTILITY_LINK_LIBS pidx)
  IF (MPI_C_FOUND)
    MESSAGE("Configuring utility with MPI support")
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    SET(UTILITY_LINK_LIBS ${UTILITY_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring utility with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(UTILITY_LINK_LIBS ${UTILITY_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////


  PIDX_ADD_CEXECUTABLE(min_max "min_max.c")
  TARGET_LINK_LIBRARIES(min_max ${UTILITY_LINK_LIBS})

ENDIF()
