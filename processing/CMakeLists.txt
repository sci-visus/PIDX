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

IF (PIDX_BUILD_PROCESSING)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(PROCESSING_LINK_LIBS pidx) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring PROCESSINGs with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(PROCESSING_LINK_LIBS ${PROCESSING_LINK_LIBS} ${MPI_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring PROCESSINGs with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(PROCESSING_LINK_LIBS ${PROCESSING_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  IF (PMT_FOUND)
    MESSAGE("Configuring pidx with PMT support")
    INCLUDE_DIRECTORIES(${PMT_INSITU_INCLUDE_DIR})
    SET(PROCESSING_LINK_LIBS ${PROCESSING_LINK_LIBS} ${PMT_INSITU_LINK_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////


  SET(MERGE_SOURCES merge.c)
  PIDX_ADD_CEXECUTABLE(merge "${MERGE_SOURCES}")
  TARGET_LINK_LIBRARIES(merge ${PROCESSING_LINK_LIBS})

  SET(GLOBAL_IDX_MERGE_SOURCES global_idx_merge.c)
  PIDX_ADD_CEXECUTABLE(global_idx_merge "${GLOBAL_IDX_MERGE_SOURCES}")
  TARGET_LINK_LIBRARIES(global_idx_merge ${PROCESSING_LINK_LIBS})

ENDIF()
