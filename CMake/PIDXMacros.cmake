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

# See https://cmake.org/Wiki/CMake_Useful_Variables

###########################################
# PIDX_SET_COMPILER_OPTIONS
###########################################

MACRO(PIDX_SET_COMPILER_OPTIONS)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
  IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    message("Using clang compiler")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=implicit-function-declaration")
  ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    message("Using gcc compiler")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration")
  ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    message("Using intel compiler")
  ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    message("Using msvc compiler")
    message(WARNING "PIDX is not yet supported on Windows")
  ENDIF ()
ENDMACRO()

###########################################
# PIDX_SET_MACHINE_SPECIFIC_OPTIONS
###########################################

MACRO(PIDX_SET_MACHINE_SPECIFIC_OPTIONS)
  IF ($ENV{HOSTNAME} MATCHES "mira")
    #mira-specIFic options
  ELSEIF ($ENV{HOSTNAME} MATCHES "vulcan")
    #vulcan-specIFic options
  ENDIF ()
ENDMACRO()

###########################################
# PIDX_ADD_EXECUTABLE
###########################################

MACRO(PIDX_ADD_CEXECUTABLE targetname files)

  ADD_EXECUTABLE(${targetname} "MACOSX_BUNDLE" ${files})
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE C)
  INSTALL(TARGETS ${targetname} 
                  RUNTIME DESTINATION bin
                  BUNDLE DESTINATION  bin)
ENDMACRO()

MACRO(PIDX_ADD_CXXEXECUTABLE targetname files)

  ADD_EXECUTABLE(${targetname} "MACOSX_BUNDLE" ${files})
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE CXX)
  INSTALL(TARGETS ${targetname} 
                  RUNTIME DESTINATION bin
                  BUNDLE DESTINATION  bin)
ENDMACRO()

###########################################
# PIDX_ADD_LIBRARY
###########################################

MACRO(PIDX_ADD_LIBRARY targetname files)
  #MESSAGE("Adding library " ${targetname} " from sources: " ${files})                        
  IF (${BUILD_SHARED_LIBS})
    ADD_LIBRARY(${targetname} SHARED ${files})
  ELSE()
    ADD_LIBRARY(${targetname} STATIC ${files})
  ENDIF ()
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE C)
  INSTALL(TARGETS ${targetname} 
                  ARCHIVE DESTINATION lib
                  LIBRARY DESTINATION lib)
ENDMACRO()
