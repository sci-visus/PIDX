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

IF (PIDX_BUILD_TOOLS)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  SET(IDXCOMPARE_SOURCES idx-compare.c)
  SET(DUMPHEADER_SOURCES idx-dump-header.c)
  SET(IDXVERIFY_SOURCES idx-verify.c)
  SET(IDXMINMAX_SOURCES idx-minmax.c)

  SET(TOOLS_LINK_LIBS pidx ${PIDX_LINK_LIBS})
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring utility with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring profiling with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  
  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(idxcompare "${IDXCOMPARE_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxdumpheader "${DUMPHEADER_SOURCES}")
  PIDX_ADD_CEXECUTABLE(idxverify "${IDXVERIFY_SOURCES}")

  PIDX_ADD_CEXECUTABLE(minmax "${IDXMINMAX_SOURCES}")
  
  TARGET_LINK_LIBRARIES(idxverify m ${TOOLS_LINK_LIBS})
  TARGET_LINK_LIBRARIES(minmax ${TOOLS_LINK_LIBS})

ENDIF()
