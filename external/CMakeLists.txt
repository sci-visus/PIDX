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

MESSAGE("Preparing external ZFP library.")
#IF (NOT ZFP_DIR)
#  SET(ZFP_DIR ${CMAKE_SOURCE_DIR}/external/zfp)
#ENDIF ()

include(ExternalProject)

ExternalProject_Add(zfp_llnl
        URL https://github.com/LLNL/zfp/releases/download/0.5.3/zfp-0.5.3.tar.gz
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/zfp_llnl
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}"
        )

ExternalProject_Get_Property(zfp_llnl install_dir)

SET(ZFP_INCLUDE_DIR  "${install_dir}/include" PARENT_SCOPE)
SET(ZFP_LIB_DIR "${install_dir}/lib" PARENT_SCOPE)

INSTALL(DIRECTORY "${install_dir}/include/" DESTINATION include/)
INSTALL(DIRECTORY "${install_dir}/lib/" DESTINATION lib/)

IF (PIDX_BUILD_IN_SITU_MERGE)
  MESSAGE("Compiling included mergetree library.")
  IF (NOT MERGE_TREE_DIR)
    SET(MERGE_TREE_DIR ${CMAKE_SOURCE_DIR}/external/parallelmergetree)
  ENDIF ()

#  set(DATAFLOW_DIR ${CMAKE_SOURCE_DIR}/external/dataflow/build)
#  find_library(DATAFLOW_LIBRARIES NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
# find_library(DATAFLOW_LIB NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
#  include_directories(${DATAFLOW_DIR}/include)

#  set(PMT_SRC_DIR ${CMAKE_SOURCE_DIR}/external/parallelmergetree/src)

#  file(GLOB PMT_SRC ${PMT_SRC_DIR}/*.cpp)
#  include_directories(${PMT_SRC_DIR}/)

  ADD_SUBDIRECTORY(dataflow)
  ADD_SUBDIRECTORY(parallelmergetree)
ENDIF ()

IF (PIDX_VTK_RENDERER)
  MESSAGE("Compiling vtk renderer")
  ADD_SUBDIRECTORY(vtkRendered)
ENDIF ()
