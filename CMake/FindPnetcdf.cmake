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
# From https://github.com/CESM-Development/CMake_Fortran_utils/blob/master/FindPnetcdf.cmake
# TODO: not yet working for PIDX

include(FindPackageHandleStandardArgs)

FIND_PATH(PNETCDF_INCLUDE_DIR 
          pnetcdf.h
          HINTS ${PNETCDF_DIR}/include)


IF (${PREFER_SHARED})
  FIND_LIBRARY(PNETCDF_LIBRARY 
               NAMES pnetcdf
               HINTS ${PNETCDF_DIR}/lib)

ELSE ()
  FIND_LIBRARY(PNETCDF_LIBRARY 
               NAMES libpnetcdf.a pnetcdf
               HINTS ${PNETCDF_DIR}/lib)
ENDIF ()

find_file( PNETCDFTEST NAMES TryPnetcdf_mod.f90  PATHS ${CMAKE_MODULE_PATH} NO_DEFAULT_PATH)
get_filename_component( CMAKE_TEST_PATH ${PNETCDFTEST} PATH)

TRY_COMPILE(PNETCDF_MOD  ${CMAKE_CURRENT_BINARY_DIR}/tryPnetcdf_mod
                          ${CMAKE_TEST_PATH}/TryPnetcdf_mod.f90
			  COMPILE_DEFINITIONS -I${PNETCDF_INCLUDE_DIR}
			  CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${PNETCDF_LIBRARIES}"
                           OUTPUT_VARIABLE Pnet_OUT)

if(NOT PNETCDF_MOD)
  TRY_COMPILE(PNETCDF_INC  ${CMAKE_CURRENT_BINARY_DIR}/tryPnetcdf_inc
                          ${CMAKE_TEST_PATH}/TryPnetcdf_inc.f90
			  COMPILE_DEFINITIONS -I${PNETCDF_INCLUDE_DIR}
			  CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${PNETCDF_LIBRARIES}"
                           OUTPUT_VARIABLE Pnet_OUT)
endif()



SET(PNETCDF_LIBRARIES ${PNETCDF_LIBRARY} )
SET(PNETCDF_INCLUDE_DIRS ${PNETCDF_INCLUDE_DIR} )

# Handle QUIETLY and REQUIRED.
find_package_handle_standard_args(pnetcdf DEFAULT_MSG
  PNETCDF_LIBRARY PNETCDF_INCLUDE_DIR  )

mark_as_advanced(PNETCDF_INCLUDE_DIR PNETCDF_LIBRARY PNETCDF_INC PNETCDF_MOD)