#####################################################
#   PIDX Parallel I/O Library                      ##
#   Copyright (c) 2010-2014 University of Utah     ##
#   Scientific Computing and Imaging Institute     ##
#   72 S Central Campus Drive, Room 3750           ##
#   Salt Lake City, UT 84112                       ##
#                                                  ##
#   PIDX is licensed under the Creative Commons    ##
#   Attribution-NonCommercial-NoDerivatives 4.0    ##
#   International License. See LICENSE.md.         ##
#                                                  ##
#   For information about this project see:        ##
#   http://www.cedmav.com/pidx                     ##
#   or contact: pascucci@sci.utah.edu              ##
#   For support: PIDX-support@visus.net            ##
#                                                  ##
#####################################################

IF (PIDX_BUILD_CONVERT)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  
  SET(HDF_TO_IDX_SOURCES hdf_to_idx.c)


  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  
  if (PIDX_OPTION_HDF5)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

    PIDX_ADD_EXECUTABLE(hdf5_to_idx "${HDF_TO_IDX_SOURCES}")
    ADD_DEPENDENCIES(hdf5_to_idx pidx)
    TARGET_LINK_LIBRARIES(hdf5_to_idx pidx)

    IF (MPI_C_FOUND)
      TARGET_LINK_LIBRARIES(hdf5_to_idx ${MPI_C_LIBRARIES})
    ENDIF()

    ADD_DEPENDENCIES(hdf5_to_idx ${HDF5_LIBRARIES})
    TARGET_LINK_LIBRARIES(hdf5_to_idx ${HDF5_LIBRARIES})
  ENDIF()

ENDIF()

