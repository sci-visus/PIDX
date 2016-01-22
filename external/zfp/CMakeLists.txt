INCLUDE_DIRECTORIES(./inc)
FILE(GLOB ZFP_SRC ./inc/*.h ./src/*.cpp)
IF (${BUILD_SHARED_LIBS})
  add_definitions(-DZFP_DLL)
ENDIF()
PIDX_ADD_LIBRARY(zfp "${ZFP_SRC}")

#set vars as in FindZFP.cmake
SET(ZFP_FOUND 1                                        CACHE BOOL "" FORCE)
SET(ZFP_INCLUDE_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/inc" CACHE STRING "" FORCE)
SET(ZFP_INCLUDE_DIRS "${ZFP_INCLUDE_DIR}"              CACHE STRING "" FORCE)
SET(ZFP_LIBRARY      zfp                               CACHE STRING "" FORCE)
SET(ZFP_LIBRARIES    stdc++ ${ZFP_LIBRARY}             CACHE STRING "" FORCE)
