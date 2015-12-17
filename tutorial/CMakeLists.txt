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

IF (PIDX_BUILD_TUTORIAL)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  #FILE(GLOB EXAMPLES_SOURCES *.h *.c *.cpp)
  
  # ////////////////////////////////////////
  # includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(TUTORIAL_LINK_LIBS pidx) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring tutorials with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${MPI_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring tutorials with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  SET(CHECKPOINT_SOURCES 1_PIDX_Basics_Dumps/checkpoint.c)
  PIDX_ADD_EXECUTABLE(checkpoint "${CHECKPOINT_SOURCES}")
  TARGET_LINK_LIBRARIES(checkpoint ${TUTORIAL_LINK_LIBS})

  SET(RESTART_SOURCES 1_PIDX_Basics_Dumps/restart.c)
  PIDX_ADD_EXECUTABLE(restart "${RESTART_SOURCES}")
  TARGET_LINK_LIBRARIES(restart ${TUTORIAL_LINK_LIBS})

  SET(CHECKPOINT_MULTI_BUFFER_SOURCES 1_PIDX_Basics_Dumps/checkpoint_multi_buffer.c)
  PIDX_ADD_EXECUTABLE(checkpoint_multi_buffer "${CHECKPOINT_MULTI_BUFFER_SOURCES}")
  TARGET_LINK_LIBRARIES(checkpoint_multi_buffer ${TUTORIAL_LINK_LIBS})

  SET(ROI_SOURCES 2_PIDX_Advanced_ROI/roi_writes.c)
  PIDX_ADD_EXECUTABLE(roi_writes "${ROI_SOURCES}")
  TARGET_LINK_LIBRARIES(roi_writes ${TUTORIAL_LINK_LIBS})

  SET(REDUCED_RESOLUTION_WRITE_SOURCES 2_PIDX_Advanced_ROI/reduced_resolution_writes.c)
  PIDX_ADD_EXECUTABLE(reduced_resolution_writes "${REDUCED_RESOLUTION_WRITE_SOURCES}")
  TARGET_LINK_LIBRARIES(reduced_resolution_writes ${TUTORIAL_LINK_LIBS})

  SET(ADAPTIVE_ROI_SOURCES 2_PIDX_Advanced_ROI/adaptive_roi_writes.c)
  PIDX_ADD_EXECUTABLE(adaptive_roi_writes "${ADAPTIVE_ROI_SOURCES}")
  TARGET_LINK_LIBRARIES(adaptive_roi_writes ${TUTORIAL_LINK_LIBS})

  SET(RAW_TO_IDX_SOURCES 3_PIDX_Utility_Converters/raw_to_idx.c)
  PIDX_ADD_EXECUTABLE(raw_to_idx "${RAW_TO_IDX_SOURCES}")
  TARGET_LINK_LIBRARIES(raw_to_idx ${TUTORIAL_LINK_LIBS})
  
  SET(FILE_PER_PROCESS_SOURCES 3_PIDX_Utility_Converters/file_per_process_to_idx.c)
  PIDX_ADD_EXECUTABLE(file_per_process_to_idx "${FILE_PER_PROCESS_SOURCES}")
  TARGET_LINK_LIBRARIES(file_per_process_to_idx ${TUTORIAL_LINK_LIBS})
  
  SET(TEMPLATE_CONVETER_SOURCES 3_PIDX_Utility_Converters/template_converter.c)
  PIDX_ADD_EXECUTABLE(template_converter "${TEMPLATE_CONVETER_SOURCES}")
  TARGET_LINK_LIBRARIES(template_converter ${TUTORIAL_LINK_LIBS})

ENDIF()
