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
  # includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(PIDXTEST_LINK_LIBS pidx) 
  IF (MPI_C_FOUND)
    MESSAGE("Configuring pidxtest with MPI support")
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    SET(PIDXTEST_LINK_LIBS ${PIDXTEST_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()


  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  PIDX_ADD_EXECUTABLE(pidxtest "${PIDXTEST_SOURCES}")
  ADD_DEPENDENCIES(pidxtest pidx)
  TARGET_LINK_LIBRARIES(pidxtest ${PIDXTEST_LINK_LIBS})
  if (PIDX_OPTION_HDF5)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    ADD_DEPENDENCIES(pidxtest ${HDF5_LIBRARIES})
    TARGET_LINK_LIBRARIES(pidxtest ${HDF5_LIBRARIES})
  ENDIF()




ENDIF()
