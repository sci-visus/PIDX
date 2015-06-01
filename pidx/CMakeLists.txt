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

FILE(GLOB PIDX_SOURCES *.h *.c)


# ////////////////////////////////////////
# includes
# ////////////////////////////////////////

IF (MPI_C_FOUND)
  MESSAGE("Configuring pidx with MPI support")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${MPI_C_LIBRARIES})
ENDIF()

IF (PIDX_HAVE_ZFP)
  MESSAGE("Configuring pidx with ZFP support")
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${ZFP_DIR}/lib/libzfp.a)
ENDIF()

# ////////////////////////////////////////
# library
# ////////////////////////////////////////

PIDX_ADD_LIBRARY(pidx "${PIDX_SOURCES}")
IF (PIDX_HAVE_ZFP)
  TARGET_LINK_LIBRARIES(pidx m stdc++ ${PIDX_LINK_LIBS})
ELSE()
  TARGET_LINK_LIBRARIES(pidx m ${PIDX_LINK_LIBS})
ENDIF()
SET_TARGET_PROPERTIES(pidx PROPERTIES LINKER_LANGUAGE CXX)
