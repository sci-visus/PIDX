##
## BSD 3-Clause License
## 
## Copyright (c) 2010-2018 ViSUS L.L.C., 
## Scientific Computing and Imaging Institute of the University of Utah
## 
## ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
## University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
##  
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## * Redistributions of source code must retain the above copyright notice, this
## list of conditions and the following disclaimer.
## 
## * Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
## 
## * Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## For additional information about this project contact: pascucci@acm.org
## For support: support@visus.net
## 
##


cmake_minimum_required(VERSION 2.8.4)

PROJECT(PIDX)


# ///////////////////////////////////////////////
# Cmake configurations
# ///////////////////////////////////////////////

SET(CMAKE_NUM_PROCS 8)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Use a verbose makefile")
OPTION(BUILD_SHARED_LIBS "Build shared libraries." FALSE)
OPTION(ENABLE_MPI "Enable MPI." TRUE)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "" FORCE)
ENDIF()

SET(ZFP_FOUND 1)

# ///////////////////////////////////////////////
# external packages
# ///////////////////////////////////////////////

IF (ENABLE_MPI)
   FIND_PACKAGE(MPI REQUIRED)
   IF (MPI_C_FOUND)
     SET(PIDX_HAVE_MPI 1)
   ENDIF ()
ENDIF ()


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
ENDIF ()
MESSAGE("PIDX_GIT_REVISION ${PIDX_GIT_REVISION}")


# ///////////////////////////////////////////////
# PIDX_BUILD_XXXX
# ///////////////////////////////////////////////

OPTION(PIDX_BUILD_EXAMPLES "Build examples" FALSE)
MESSAGE("PIDX_BUILD_EXAMPLES ${PIDX_BUILD_EXAMPLES}")

OPTION(PIDX_BUILD_TOOLS "Build pidx profile and testing tools" FALSE)
MESSAGE("PIDX_BUILD_TOOLS ${PIDX_BUILD_TOOLS}")

#OPTION(PIDX_BUILD_CONVERT "Build pidx convert tools" FALSE)
#MESSAGE("PIDX_BUILD_CONVERT ${PIDX_BUILD_CONVERT}")

#OPTION(PIDX_BUILD_PROCESSING "Build processing" FALSE)
#MESSAGE("PIDX_BUILD_PROCESSING ${PIDX_BUILD_PROCESSING}")

#OPTION(PIDX_BUILD_IN_SITU_MERGE "Build in situ merge" FALSE)
#MESSAGE("PIDX_BUILD_IN_SITU_MERGE ${PIDX_BUILD_IN_SITU_MERGE}")

#OPTION(PIDX_VTK_RENDERER "Build in vtk renderer" FALSE)
#MESSAGE("PIDX_VTK_RENDERER ${PIDX_VTK_RENDERER}")


#IF (PIDX_BUILD_IN_SITU_MERGE)
#  set(DATAFLOW_DIR "${CMAKE_PREFIX_PATH}../../PIDX_external/dataflow/build" CACHE STRING "DataFlow installation directory")
#  find_library(DATAFLOW_LIB NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
#  include_directories(${DATAFLOW_DIR}/include)

#  set(PMT_SRC_DIR "${CMAKE_PREFIX_PATH}../../PIDX_external/parallelmergetree/src" CACHE STRING "PMT src directory")

#  file(GLOB PMT_SRC ${PMT_SRC_DIR}/*.cpp)
#  include_directories(${PMT_SRC_DIR}/)

#  add_executable(your_app your_app.cpp ${PMT_SRC})
#  target_link_libraries(your_app ${DATAFLOW_LIB})
#ENDIF ()

# ///////////////////////////////////////////////
# PIDX_OPTION_XXXX
# ///////////////////////////////////////////////

SET(PIDX_HAVE_ZFP 1)  #ZFP is embedded and installed with PIDX


#MESSAGE("PIDX_BUILD_VIEWER ${PIDX_BUILD_VIEWER}")
#OPTION(PIDX_BUILD_VIEWER "Enable OpenGL" FALSE)
#IF (PIDX_BUILD_VIEWER)
#   find_package(OpenGL REQUIRED)
#   find_package(GLUT REQUIRED)
#   IF (OPENGL_FOUND)
#     SET(PIDX_HAVE_OPENGL 1)
#   ENDIF ()
#ENDIF ()

#OPTION(PIDX_OPTION_NVISUSIO "Test nvisusio library (hint: set CMAKE_MODULE_PATH to path of ViSUS installation)." FALSE)
#MESSAGE("PIDX_OPTION_NVISUSIO ${PIDX_OPTION_NVISUSIO}")
#IF (PIDX_OPTION_NVISUSIO)
#  FIND_PACKAGE(NVISUSIO REQUIRED)
#  IF (NVISUSIO_FOUND)
#    SET(PIDX_HAVE_NVISUSIO 1)
#  ENDIF ()
#ENDIF ()

#MESSAGE("PIDX_OPTION_PNETCDF ${PIDX_OPTION_PNETCDF}")
#OPTION(PIDX_OPTION_PNETCDF "Test Pnetcdf (netcdf-1)" FALSE)
#IF (PIDX_OPTION_PNETCDF)
#  FIND_PACKAGE(PNETCDF REQUIRED)
#  IF (PNETCDF_FOUND)
#    SET(PIDX_HAVE_PNETCDF 1)
#  ENDIF ()
#ENDIF ()

#MESSAGE("PIDX_OPTION_NETCDF ${PIDX_OPTION_NETCDF}")
#OPTION(PIDX_OPTION_NETCDF "Test parallel netcdf-4" FALSE)
#IF (PIDX_OPTION_NETCDF)
#  FIND_PACKAGE(NETCDF REQUIRED)
#  IF (NETCDF_FOUND)
#    SET(PIDX_HAVE_NETCDF 1)
#  ENDIF ()
#ENDIF ()

#MESSAGE("PIDX_OPTION_HDF5 ${PIDX_OPTION_HDF5}")
#OPTION(PIDX_OPTION_HDF5 "Test HDF5" FALSE)
#IF (PIDX_OPTION_HDF5)
#  FIND_PACKAGE(HDF5 REQUIRED)
#  IF (HDF5_FOUND)
#    SET(PIDX_HAVE_HDF5 1)
#  ENDIF ()
#ENDIF ()


#MESSAGE("PIDX_OPTION_VTK ${PIDX_OPTION_VTK}")
#OPTION(PIDX_OPTION_VTK "Build with vtk" FALSE)
#IF (PIDX_OPTION_VTK)
#  set(VTK_DIR "/home/sidharth/research/PIDX/external/VTK-7.1.0/build/")
#  find_package(VTK REQUIRED)
#  include(${VTK_USE_FILE})
#  #add_executable(Test Test.cxx)
#  #target_link_libraries(Test ${VTK_LIBRARIES})
#ENDIF ()

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
ADD_SUBDIRECTORY(viewer2D)
ADD_SUBDIRECTORY(tools)
ADD_SUBDIRECTORY(processing)
ADD_SUBDIRECTORY(examples)

# ///////////////////////////////////////////////
# installation
# ///////////////////////////////////////////////

MESSAGE("CMAKE_INSTALL_PREFIX is " ${CMAKE_INSTALL_PREFIX})
INSTALL(DIRECTORY pidx/ DESTINATION include FILES_MATCHING PATTERN "*.h")
#INSTALL(DIRECTORY external/ DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(FILES     "${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindPIDX.cmake" DESTINATION .)
INSTALL(FILES     "${PROJECT_BINARY_DIR}/PIDX_config.h" DESTINATION include)
