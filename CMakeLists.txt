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


cmake_minimum_required(VERSION 2.6)

PROJECT(PIDX)


# ///////////////////////////////////////////////
# Cmake configurations
# ///////////////////////////////////////////////

SET(CMAKE_NUM_PROCS 8)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Use a verbose makefile")
OPTION(BUILD_SHARED_LIBS "Build shared libraries." FALSE)
OPTION(ENABLE_MPI "Enable MPI." TRUE)


# ///////////////////////////////////////////////
# external packages
# ///////////////////////////////////////////////

FIND_PACKAGE(ZLIB REQUIRED)

IF (ENABLE_MPI)
   FIND_PACKAGE(MPI REQUIRED)
ENDIF()


# ///////////////////////////////////////////////
# platform configuration
# ///////////////////////////////////////////////

INCLUDE(PIDXPlatform)
INCLUDE(PIDXMacros)


# ///////////////////////////////////////////////
# PIDX_GIT_REVISION
# ///////////////////////////////////////////////

FIND_PROGRAM(GIT_CMD git)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GIT REQUIRED_VARS GIT_CMD)
IF (GIT_CMD)
	EXECUTE_PROCESS(COMMAND ${GIT_CMD} rev-list -n 1 HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE PIDX_GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(${result} TEMP_VAR)
ELSE()
  SET(PIDX_GIT_REVISION "")
ENDIF()
ADD_DEFINITIONS(-DPIDX_GIT_REVISION=\"${PIDX_GIT_REVISION}\")
MESSAGE("PIDX_GIT_REVISION ${PIDX_GIT_REVISION}")


# ///////////////////////////////////////////////
# PIDX_BUILD_XXXX
# ///////////////////////////////////////////////

OPTION(PIDX_BUILD_TEST "Build pidxtest" TRUE)
MESSAGE("PIDX_BUILD_TEST ${PIDX_BUILD_TEST}")

OPTION(PIDX_BUILD_PROFILE "Build pidxprofile" FALSE)
MESSAGE("PIDX_BUILD_PROFILE ${PIDX_BUILD_PROFILE}")

OPTION(PIDX_BUILD_TOOLS "Build pidx tools" TRUE)
MESSAGE("PIDX_BUILD_TOOLS ${PIDX_BUILD_TOOLS}")


# ///////////////////////////////////////////////
# PIDX_OPTION_XXXX
# ///////////////////////////////////////////////

OPTION(PIDX_OPTION_PNETCDF "Test Pnetcdf" FALSE)
IF (PIDX_OPTION_PNETCDF)
  ADD_DEFINITIONS(-DPIDX_OPTION_PNETCDF=1)
  FIND_PACKAGE(Pnetcdf REQUIRED)
ENDIF()
MESSAGE("PIDX_OPTION_PNETCDF ${PIDX_OPTION_PNETCDF}")


# ///////////////////////////////////////////////
# subdirectories
# ///////////////////////////////////////////////

ADD_SUBDIRECTORY(pidx)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(profile)
ADD_SUBDIRECTORY(tools)


# ///////////////////////////////////////////////
# installation
# ///////////////////////////////////////////////

MESSAGE("CMAKE_INSTALL_PREFIX is " ${CMAKE_INSTALL_PREFIX})
INSTALL(DIRECTORY pidx/ DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(FILES     "${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindPIDX.cmake" DESTINATION .)

