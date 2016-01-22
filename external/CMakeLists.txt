

IF (PIDX_OPTION_ZFP)
  MESSAGE("Compiling included zfp library.")
  IF (NOT ZFP_DIR)
    SET(ZFP_DIR ${CMAKE_SOURCE_DIR}/external/zfp)
  ENDIF()
  ADD_SUBDIRECTORY(zfp)
ENDIF()
