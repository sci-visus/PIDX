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


# This module finds if ZFP is installed, and sets the following variables
# indicating where it is:
#
# ZFP_FOUND               - system has ZFP
# ZFP_INCLUDE_DIR         - path to ZFP include directory
# ZFP_LIBRARIES           - all ZFP libraries
#
# Execute cmake with "-DZFP_DIR=/path/to/zfp" to help find the library.
#

FIND_PATH(ZFP_INCLUDE_DIR NAMES zfp.h PATHS ${ZFP_DIR}/inc ${ZFP_DIR}/include NO_DEFAULT_PATH)
FIND_PATH(ZFP_INCLUDE_DIR NAMES zfp.h)

IF (ZFP_INCLUDE_DIR)

  SET(ZFP_INCLUDE_DIRS "${ZFP_INCLUDE_DIR}")

  FIND_LIBRARY(ZFP_LIBRARY     zfp    PATHS ${ZFP_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(ZFP_LIBRARY     zfp    )

  SET(ZFP_LIBRARIES 
    ${ZFP_LIBRARY}
    stdc++
    )

  SET(ZFP_FOUND "YES" CACHE BOOL "ZFP library found" FORCE)

  IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE("Using ZFP_INCLUDE_DIR  = ") 
    FOREACH(inc ${ZFP_INCLUDE_DIR})
      MESSAGE("  " ${inc})
    ENDFOREACH()
    MESSAGE("Found ZFP_LIBRARIES    = ")
    FOREACH(lib ${ZFP_LIBRARIES})
      MESSAGE("  " ${lib})
    ENDFOREACH()
  ENDIF (CMAKE_VERBOSE_MAKEFILE)

ELSE ()
  IF (ZFP_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ZFP library not found. Try setting ZFP_DIR")
  ELSE()
    MESSAGE("WARNING: ZFP library not found. Try setting ZFP_DIR")
  ENDIF()
ENDIF ()

