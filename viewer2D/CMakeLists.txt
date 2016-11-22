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
# ////////////////////////////////////////
# sources
# ////////////////////////////////////////

IF (PIDX_BUILD_VIEWER)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(TUTORIAL_LINK_LIBS pidx) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring tutorials with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${MPI_LIBRARIES})
  ENDIF()

  IF (OPENGL_FOUND)
    SET(PIDX_VIEWER2D_SOURCES PIDX_slice_viewer.c)
    PIDX_ADD_CEXECUTABLE(slice_viewer "${PIDX_VIEWER2D_SOURCES}")
    TARGET_LINK_LIBRARIES(slice_viewer ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${TUTORIAL_LINK_LIBS})
  ENDIF()

ENDIF()
