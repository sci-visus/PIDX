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

INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})


# ////////////////////////////////////////
# library
# ////////////////////////////////////////

PIDX_ADD_LIBRARY(pidx "${PIDX_SOURCES}")
ADD_DEPENDENCIES(pidx ${ZLIB_LIBRARY})
TARGET_LINK_LIBRARIES(pidx m ${MPI_C_LIBRARIES})

