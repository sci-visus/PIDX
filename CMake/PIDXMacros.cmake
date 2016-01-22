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

# See https://cmake.org/Wiki/CMake_Useful_Variables

###########################################
# PIDX_SET_COMPILER_OPTIONS
###########################################

MACRO(PIDX_SET_COMPILER_OPTIONS)
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
  ENDIF()
ENDMACRO()

###########################################
# PIDX_SET_MACHINE_SPECIFIC_OPTIONS
###########################################

MACRO(PIDX_SET_MACHINE_SPECIFIC_OPTIONS)
  IF ($ENV{HOSTNAME} MATCHES "mira")
    #mira-specIFic options
  ELSEIF($ENV{HOSTNAME} MATCHES "vulcan")
    #vulcan-specIFic options
  ENDIF()
ENDMACRO()

###########################################
# PIDX_ADD_EXECUTABLE
###########################################

MACRO(PIDX_ADD_EXECUTABLE targetname files)
  #MESSAGE("Adding executable " ${targetname} " from sources: " ${files})                        
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
  ENDIF()
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE C)
  INSTALL(TARGETS ${targetname} 
                  ARCHIVE DESTINATION lib
                  LIBRARY DESTINATION lib)
ENDMACRO()
