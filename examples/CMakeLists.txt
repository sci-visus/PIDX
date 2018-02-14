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

IF (PIDX_BUILD_EXAMPLES)
  
  # ////////////////////////////////////////
  # includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(EXAMPLES_LINK_LIBS pidx ${PIDX_LINK_LIBS}) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring examples with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring examples with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(idx_write "basic-io/idx_write.c")
  TARGET_LINK_LIBRARIES(idx_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_read "basic-io/idx_read.c")
  TARGET_LINK_LIBRARIES(idx_read ${EXAMPLES_LINK_LIBS})


ENDIF()

IF (PIDX_BUILD_CONVERT)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(EXAMPLES_LINK_LIBS pidx ${PIDX_LINK_LIBS}) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring examples with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring examples with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(raw_to_idx "convert/raw-to-idx.c")
  TARGET_LINK_LIBRARIES(raw_to_idx ${EXAMPLES_LINK_LIBS})

ENDIF()

