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

ADD_SUBDIRECTORY(PIDX_generic_rst)
#ADD_SUBDIRECTORY(PIDX_wavelet_rst)
#ADD_SUBDIRECTORY(PIDX_rst)
ADD_SUBDIRECTORY(PIDX_multi_patch_rst)
ADD_SUBDIRECTORY(PIDX_block_rst)
ADD_SUBDIRECTORY(PIDX_cmp)
ADD_SUBDIRECTORY(PIDX_hz)
ADD_SUBDIRECTORY(PIDX_in_situ_interface)
#ADD_SUBDIRECTORY(PIDX_in_transit_interface)
ADD_SUBDIRECTORY(PIDX_agg)
#ADD_SUBDIRECTORY(PIDX_shared_block_agg)
ADD_SUBDIRECTORY(PIDX_file_io)
ADD_SUBDIRECTORY(PIDX_header)
ADD_SUBDIRECTORY(PIDX_wavelet)


FILE(GLOB PIDX_CORE_SOURCES *.h *.c)


# ////////////////////////////////////////
# includes
# ////////////////////////////////////////

IF (MPI_C_FOUND)
  #MESSAGE("Configuring pidx with MPI support")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(PIDX_CORE_LINK_LIBS ${PIDX_CORE_LINK_LIBS} ${MPI_C_LIBRARIES})
ENDIF()

IF (ZFP_FOUND)
  MESSAGE("Configuring pidx with ZFP support")
  INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
  SET(PIDX_CORE_LINK_LIBS ${PIDX_CORE_LINK_LIBS} ${ZFP_LIBRARIES})
ENDIF()

IF (PMT_FOUND)
  MESSAGE("Configuring pidx with PMT support")
  INCLUDE_DIRECTORIES(${PMT_INSITU_INCLUDE_DIR})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${PMT_INSITU_LINK_LIBRARIES})
ENDIF()

IF (VTK_RENDERER_FOUND)
  INCLUDE_DIRECTORIES(${VTK_RENDERER_INCLUDE_DIR})
ENDIF()

# ////////////////////////////////////////
# library
# ////////////////////////////////////////

PIDX_ADD_LIBRARY(pidx_core "${PIDX_CORE_SOURCES}")
TARGET_LINK_LIBRARIES(pidx_core m pidx_in_situ_interface pidx_utils pidx_data_handle pidx_header pidx_multi_patch_rst pidx_generic_rst pidx_in_situ_interface pidx_wavelet pidx_cmp pidx_hz pidx_agg pidx_file_io ${PIDX_CORE_LINK_LIBS})
SET_TARGET_PROPERTIES(pidx_core PROPERTIES LINKER_LANGUAGE CXX)
