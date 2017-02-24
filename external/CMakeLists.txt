

IF (PIDX_OPTION_ZFP)
  MESSAGE("Compiling included zfp library.")
  IF (NOT ZFP_DIR)
    SET(ZFP_DIR ${CMAKE_SOURCE_DIR}/external/zfp)
  ENDIF()
  ADD_SUBDIRECTORY(zfp)
ENDIF()

IF (PIDX_BUILD_IN_SITU_MERGE)
  MESSAGE("Compiling included mergetree library.")
  IF (NOT MERGE_TREE_DIR)
    SET(MERGE_TREE_DIR ${CMAKE_SOURCE_DIR}/external/parallelmergetree)
  ENDIF()

#  set(DATAFLOW_DIR ${CMAKE_SOURCE_DIR}/external/dataflow/build)
#  find_library(DATAFLOW_LIBRARIES NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
# find_library(DATAFLOW_LIB NAMES dataflow HINTS ${DATAFLOW_DIR}/lib)
#  include_directories(${DATAFLOW_DIR}/include)

#  set(PMT_SRC_DIR ${CMAKE_SOURCE_DIR}/external/parallelmergetree/src)

#  file(GLOB PMT_SRC ${PMT_SRC_DIR}/*.cpp)
#  include_directories(${PMT_SRC_DIR}/)

  ADD_SUBDIRECTORY(dataflow)
  ADD_SUBDIRECTORY(parallelmergetree)
ENDIF()



IF (PIDX_VTK_RENDERER)
  MESSAGE("Compiling vtk renderer")
  ADD_SUBDIRECTORY(vtkRendered)
ENDIF()
