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


cmake_minimum_required(VERSION 2.8)

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

IF (ENABLE_MPI)
   FIND_PACKAGE(MPI REQUIRED)
   IF (MPI_C_FOUND)
     SET(PIDX_HAVE_MPI 1)
   ENDIF()
ENDIF()


# ///////////////////////////////////////////////
# platform configuration
# ///////////////////////////////////////////////

INCLUDE(PIDXPlatform)
INCLUDE(PIDXMacros)
PIDX_SET_COMPILER_OPTIONS()
PIDX_SET_MACHINE_SPECIFIC_OPTIONS()


# ///////////////////////////////////////////////
# PIDX_GIT_REVISION
# ///////////////////////////////////////////////

FIND_PROGRAM(GIT_CMD git)
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GIT REQUIRED_VARS GIT_CMD)
IF (GIT_CMD)
	EXECUTE_PROCESS(COMMAND ${GIT_CMD} rev-list -n 1 HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE PIDX_GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(${result} TEMP_VAR)
ELSE()
  SET(PIDX_GIT_REVISION "")
ENDIF()
MESSAGE("PIDX_GIT_REVISION ${PIDX_GIT_REVISION}")


# ///////////////////////////////////////////////
# PIDX_BUILD_XXXX
# ///////////////////////////////////////////////

OPTION(PIDX_BUILD_TUTORIAL "Build tutorial" FALSE)
MESSAGE("PIDX_BUILD_TUTORIAL ${PIDX_BUILD_TUTORIAL}")

OPTION(PIDX_BUILD_UTILITY "Build utility" FALSE)
MESSAGE("PIDX_BUILD_TUTORIAL ${PIDX_BUILD_UTILITY}")

OPTION(PIDX_BUILD_EXAMPLE "Build example" FALSE)
MESSAGE("PIDX_BUILD_EXAMPLE ${PIDX_BUILD_EXAMPLE}")

OPTION(PIDX_BUILD_CONVERT "Build pidx convert tools" FALSE)
MESSAGE("PIDX_BUILD_CONVERT ${PIDX_BUILD_CONVERT}")

OPTION(PIDX_BUILD_PROFILE "pidx profile tools" FALSE)
MESSAGE("PIDX_BUILD_PROFILE ${PIDX_BUILD_PROFILE}")

OPTION(PIDX_BUILD_PROCESSING "Build processing" TRUE)
MESSAGE("PIDX_BUILD_PROCESSING ${PIDX_BUILD_PROCESSING}")

OPTION(PIDX_BUILD_IN_SITU_MERGE "Build in situ merge" FALSE)
MESSAGE("PIDX_BUILD_IN_SITU_MERGE ${PIDX_BUILD_IN_SITU_MERGE}")
IF (PIDX_BUILD_IN_SITU_MERGE)
  set(DATAFLOW_DIR "${CMAKE_PREFIX_PATH}../../PIDX_external/dataflow/build" CACHE STRING "DataFlow installation directory")
  find_library(DATAFLOW_LIB NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
  include_directories(${DATAFLOW_DIR}/include)

  set(PMT_SRC_DIR "${CMAKE_PREFIX_PATH}../../PIDX_external/parallelmergetree/src" CACHE STRING "PMT src directory")

  file(GLOB PMT_SRC ${PMT_SRC_DIR}/*.cpp)
  include_directories(${PMT_SRC_DIR}/)

  add_executable(your_app your_app.cpp ${PMT_SRC})
  target_link_libraries(your_app ${DATAFLOW_LIB})
ENDIF()

# ///////////////////////////////////////////////
# PIDX_OPTION_XXXX
# ///////////////////////////////////////////////

#SET(PIDX_OPTION_ZFP 0)
OPTION(PIDX_OPTION_ZFP "Enable ZFP." FALSE)
MESSAGE("PIDX_OPTION_ZFP ${PIDX_OPTION_ZFP}")
IF (PIDX_OPTION_ZFP)
  SET(PIDX_HAVE_ZFP 1)
ENDIF()


MESSAGE("PIDX_BUILD_VIEWER ${PIDX_BUILD_VIEWER}")
OPTION(PIDX_BUILD_VIEWER "Enable OpenGL" FALSE)
IF (PIDX_BUILD_VIEWER)
   find_package(OpenGL REQUIRED)
   find_package(GLUT REQUIRED)
   IF (OPENGL_FOUND)
     SET(PIDX_HAVE_OPENGL 1)
   ENDIF()
ENDIF()



OPTION(PIDX_OPTION_NVISUSIO "Test nvisusio library (hint: set CMAKE_MODULE_PATH to path of ViSUS installation)." FALSE)
MESSAGE("PIDX_OPTION_NVISUSIO ${PIDX_OPTION_NVISUSIO}")
IF (PIDX_OPTION_NVISUSIO)
  FIND_PACKAGE(NVISUSIO REQUIRED)
  IF (NVISUSIO_FOUND)
    SET(PIDX_HAVE_NVISUSIO 1)
  ENDIF()
ENDIF()

MESSAGE("PIDX_OPTION_PNETCDF ${PIDX_OPTION_PNETCDF}")
OPTION(PIDX_OPTION_PNETCDF "Test Pnetcdf (netcdf-1)" FALSE)
IF (PIDX_OPTION_PNETCDF)
  FIND_PACKAGE(PNETCDF REQUIRED)
  IF (PNETCDF_FOUND)
    SET(PIDX_HAVE_PNETCDF 1)
  ENDIF()
ENDIF()

MESSAGE("PIDX_OPTION_NETCDF ${PIDX_OPTION_NETCDF}")
OPTION(PIDX_OPTION_NETCDF "Test parallel netcdf-4" FALSE)
IF (PIDX_OPTION_NETCDF)
  FIND_PACKAGE(NETCDF REQUIRED)
  IF (NETCDF_FOUND)
    SET(PIDX_HAVE_NETCDF 1)
  ENDIF()
ENDIF()

MESSAGE("PIDX_OPTION_HDF5 ${PIDX_OPTION_HDF5}")
OPTION(PIDX_OPTION_HDF5 "Test HDF5" FALSE)
IF (PIDX_OPTION_HDF5)
  FIND_PACKAGE(HDF5 REQUIRED)
  IF (HDF5_FOUND)
    SET(PIDX_HAVE_HDF5 1)
  ENDIF()
ENDIF()

# ///////////////////////////////////////////////
# configuration
# ///////////////////////////////////////////////

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/PIDX_config.h.in" "${PROJECT_BINARY_DIR}/PIDX_config.h" )
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

# ///////////////////////////////////////////////
# subdirectories
# ///////////////////////////////////////////////

ADD_SUBDIRECTORY(external)
ADD_SUBDIRECTORY(pidx)
ADD_SUBDIRECTORY(tutorial)
ADD_SUBDIRECTORY(utility)
ADD_SUBDIRECTORY(viewer2D)
ADD_SUBDIRECTORY(profile)
ADD_SUBDIRECTORY(processing)
# ADD_SUBDIRECTORY(examples)

# ///////////////////////////////////////////////
# installation
# ///////////////////////////////////////////////

MESSAGE("CMAKE_INSTALL_PREFIX is " ${CMAKE_INSTALL_PREFIX})
INSTALL(DIRECTORY pidx/ DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(FILES     "${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindPIDX.cmake" DESTINATION .)
INSTALL(FILES     "${PROJECT_BINARY_DIR}/PIDX_config.h" DESTINATION include)
