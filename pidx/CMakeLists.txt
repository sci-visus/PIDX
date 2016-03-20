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


ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(comm)
ADD_SUBDIRECTORY(data_handle)
ADD_SUBDIRECTORY(topology)
ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(io)

FILE(GLOB PIDX_SOURCES *.h *.c ./utils/*.h ./utils/*.c ./comm/*.h ./comm/*.c ./data_handle/*.h ./data_handle/*.c ./core/PIDX_agg/*.h ./core/PIDX_agg/*.c ./core/PIDX_block_rst/*.h ./core/PIDX_block_rst/*.c ./core/PIDX_cmp/*.h ./core/PIDX_cmp/*.c ./core/PIDX_file_io/*.h ./core/PIDX_file_io/*.c ./core/PIDX_header/*.h ./core/PIDX_header/*.c ./core/PIDX_hz/*.h ./core/PIDX_hz/*.c ./core/PIDX_rst/*.h ./core/PIDX_rst/*.c ./io/*.h ./io/*.c ./io/PIDX_idx_io/*.h ./io/PIDX_idx_io/*.c ./io/PIDX_partitioned_idx_io/*.h ./io/PIDX_partitioned_idx_io/*.c ./io/PIDX_partition_merge_idx_io/*.h ./io/PIDX_partition_merge_idx_io/*.c ./io/PIDX_raw_io/*.h ./io/PIDX_raw_io/*.c
)


# ////////////////////////////////////////
# includes
# ////////////////////////////////////////

IF (MPI_C_FOUND)
  MESSAGE("Configuring pidx with MPI support")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${MPI_C_LIBRARIES})
ENDIF()

IF (ZFP_FOUND)
  MESSAGE("Configuring pidx with ZFP support")
  INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${ZFP_LIBRARIES})
ENDIF()

# ////////////////////////////////////////
# library
# ////////////////////////////////////////
MESSAGE("XXXX")
MESSAGE("${PIDX_SOURCES}")
MESSAGE("XXXX")
MESSAGE("${PIDX_UTILS_SOURCES}")
PIDX_ADD_LIBRARY(pidx "${PIDX_SOURCES}")
TARGET_LINK_LIBRARIES(pidx m pidx_comm pidx_data_handle pidx_topology pidx_utils pidx_core pidx_io ${PIDX_LINK_LIBS})


SET_TARGET_PROPERTIES(pidx PROPERTIES LINKER_LANGUAGE CXX)
