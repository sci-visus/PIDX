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

IF (PIDX_OPTION_HDF5)

  SET(HDF_LINK_LIBS ${HDF5_LIBRARIES} ${MPI_C_LIBRARIES})
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_PATH})
  ADD_EXECUTABLE(hdf5_write "HDF5/test-HDF5-writer.c")
  TARGET_LINK_LIBRARIES(hdf5_write ${HDF_LINK_LIBS})

ENDIF()

IF (PIDX_BUILD_TUTORIAL)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(TUTORIAL_LINK_LIBS pidx)
  IF (MPI_C_FOUND)
    MESSAGE("Configuring tutorials with MPI support")
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring tutorials with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////


  PIDX_ADD_CEXECUTABLE(multi_buffer_single_variable_read "basic_io/read/multi_buffer_single_variable_read.c")
  TARGET_LINK_LIBRARIES(multi_buffer_single_variable_read ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_single_variable_read "basic_io/read/single_buffer_single_variable_read.c")
  TARGET_LINK_LIBRARIES(single_buffer_single_variable_read ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(multi_buffer_all_variable_read "basic_io/read/multi_buffer_all_variable_read.c")
  TARGET_LINK_LIBRARIES(multi_buffer_all_variable_read ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_all_variable_read "basic_io/read/single_buffer_all_variable_read.c")
  TARGET_LINK_LIBRARIES(single_buffer_all_variable_read ${TUTORIAL_LINK_LIBS})



  PIDX_ADD_CEXECUTABLE(multi_buffer_idx_write "basic_io/write/idx/multi_buffer_idx_write.c")
  TARGET_LINK_LIBRARIES(multi_buffer_idx_write ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_idx_write "basic_io/write/idx/single_buffer_idx_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_idx_write ${TUTORIAL_LINK_LIBS})




  PIDX_ADD_CEXECUTABLE(multi_buffer_raw_write "basic_io/write/raw/multi_buffer_raw_write.c")
  TARGET_LINK_LIBRARIES(multi_buffer_raw_write ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_raw_write "basic_io/write/raw/single_buffer_raw_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_raw_write ${TUTORIAL_LINK_LIBS})





  PIDX_ADD_CEXECUTABLE(multi_buffer_global_partitioned_idx_write "basic_io/write/global_partitioned_idx/multi_buffer_global_partitioned_idx_write.c")
  TARGET_LINK_LIBRARIES(multi_buffer_global_partitioned_idx_write ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_global_partitioned_idx_write "basic_io/write/global_partitioned_idx/single_buffer_global_partitioned_idx_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_global_partitioned_idx_write ${TUTORIAL_LINK_LIBS})



  PIDX_ADD_CEXECUTABLE(multi_buffer_local_partitioned_idx_write "basic_io/write/local_partitioned_idx/multi_buffer_local_partitioned_idx_write.c")
  TARGET_LINK_LIBRARIES(multi_buffer_local_partitioned_idx_write ${TUTORIAL_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_local_partitioned_idx_write "basic_io/write/local_partitioned_idx/single_buffer_local_partitioned_idx_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_local_partitioned_idx_write ${TUTORIAL_LINK_LIBS})


ENDIF()
