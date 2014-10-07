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

IF (PIDX_BUILD_PROFILE)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  SET(FORTRANREADER_SOURCES fortran-reader.c)
  SET(FORTRANRST_SOURCES fortran-rst.c)
  SET(FORTRANTOIDX_SOURCES fortran-to-idx.c)
  SET(AGGREGATIONBENCHMARK_SOURCES aggregation-benchmark.c)
  SET(AGGREGATIONVERIFYBGQ_SOURCES aggregation-verify-bgq.c)
  SET(AGGREGATIONVERIFY_SOURCES aggregation-verify.c)
  SET(WRITEPNETCDF_SOURCES pnetcdf-write-standard.c)
  SET(READPNETCDF_SOURCES pnetcdf-read-standard.c)
  SET(IDXCOMPARECLUBED_SOURCES idx-compare-clubed-variables.c)
  SET(IDXCOMPARE_SOURCES idx-compare.c)
  SET(IDXVERIFY_SOURCES idx-verify.c)
  SET(KBENCHMARK1_SOURCES k-benchmark-1.c)
  SET(KBENCHMARK2_SOURCES k-benchmark-2.c)
  SET(PARALLETESTREAD_SOURCES parallel-test-read.c)
  SET(RESTRUCTUREDRANKORDER_SOURCES restructured_rank_order.c)
  SET(RSTPHASE_SOURCES rst-phase.c)
  SET(S3DVELOCITY2IDX_SOURCES s3d-velocity-2idx.c)
  SET(SHAREDMEMBENCHMARK1_SOURCES shared-memory-benchmark1.c)
  SET(SPARSEDATAWRITER_SOURCES sparse-data-writer.c)


  # ////////////////////////////////////////
  # includes
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

  if (PIDX_OPTION_PNETCDF)
    INCLUDE_DIRECTORIES(${PNETCDF_INCLUDE_DIR})
  ENDIF()


  # ////////////////////////////////////////
  # executables
  # ////////////////////////////////////////

  PIDX_ADD_EXECUTABLE(fortranreader "${FORTRANREADER_SOURCES}")
  ADD_DEPENDENCIES(fortranreader pidx)
  TARGET_LINK_LIBRARIES(fortranreader pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(fortranrst "${FORTRANRST_SOURCES}")
  ADD_DEPENDENCIES(fortranrst pidx)
  TARGET_LINK_LIBRARIES(fortranrst pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(fortrantoidx "${FORTRANTOIDX_SOURCES}")
  ADD_DEPENDENCIES(fortrantoidx pidx)
  TARGET_LINK_LIBRARIES(fortrantoidx pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(idxcompareclubed "${IDXCOMPARECLUBED_SOURCES}")
  ADD_DEPENDENCIES(idxcompareclubed pidx)
  TARGET_LINK_LIBRARIES(idxcompareclubed pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(idxcompare "${IDXCOMPARE_SOURCES}")
  ADD_DEPENDENCIES(idxcompare pidx)
  TARGET_LINK_LIBRARIES(idxcompare pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(idxverify "${IDXVERIFY_SOURCES}")
  ADD_DEPENDENCIES(idxverify pidx)
  TARGET_LINK_LIBRARIES(idxverify pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(kbenchmark1 "${KBENCHMARK1_SOURCES}")
  ADD_DEPENDENCIES(kbenchmark1 pidx)
  TARGET_LINK_LIBRARIES(kbenchmark1 pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(kbenchmark2 "${KBENCHMARK2_SOURCES}")
  ADD_DEPENDENCIES(kbenchmark2 pidx)
  TARGET_LINK_LIBRARIES(kbenchmark2 pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(paralletestread "${PARALLETESTREAD_SOURCES}")
  ADD_DEPENDENCIES(paralletestread pidx)
  TARGET_LINK_LIBRARIES(paralletestread pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(restructuredrankorder "${RESTRUCTUREDRANKORDER_SOURCES}")
  ADD_DEPENDENCIES(restructuredrankorder pidx)
  TARGET_LINK_LIBRARIES(restructuredrankorder pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(rstphase "${RSTPHASE_SOURCES}")
  ADD_DEPENDENCIES(rstphase pidx)
  TARGET_LINK_LIBRARIES(rstphase pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(s3dvelocity2idx "${S3DVELOCITY2IDX_SOURCES}")
  ADD_DEPENDENCIES(s3dvelocity2idx pidx)
  TARGET_LINK_LIBRARIES(s3dvelocity2idx pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(sharedmembenchmark1 "${SHAREDMEMBENCHMARK1_SOURCES}")
  ADD_DEPENDENCIES(sharedmembenchmark1 pidx)
  TARGET_LINK_LIBRARIES(sharedmembenchmark1 pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(sparsedatawriter "${SPARSEDATAWRITER_SOURCES}")
  ADD_DEPENDENCIES(sparsedatawriter pidx)
  TARGET_LINK_LIBRARIES(sparsedatawriter pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(aggregation-verify "${AGGREGATIONVERIFY_SOURCES}")
  ADD_DEPENDENCIES(aggregation-verify pidx)
  TARGET_LINK_LIBRARIES(aggregation-verify pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(aggregation-verify-bgq "${AGGREGATIONVERIFYBGQ_SOURCES}")
  ADD_DEPENDENCIES(aggregation-verify-bgq pidx)
  TARGET_LINK_LIBRARIES(aggregation-verify-bgq pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  PIDX_ADD_EXECUTABLE(aggregation-benchmark "${AGGREGATIONBENCHMARK_SOURCES}")
  ADD_DEPENDENCIES(aggregation-benchmark pidx)
  TARGET_LINK_LIBRARIES(aggregation-benchmark pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})

  if (PIDX_OPTION_PNETCDF)
    PIDX_ADD_EXECUTABLE(writepnetcdf "${WRITEPNETCDF_SOURCES}")
    ADD_DEPENDENCIES(writepnetcdf pidx)
    TARGET_LINK_LIBRARIES(writepnetcdf pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})
    ADD_DEPENDENCIES(writepnetcdf ${PNETCDF_LIBRARY})
    TARGET_LINK_LIBRARIES(writepnetcdf ${PNETCDF_LIBRARY})

    PIDX_ADD_EXECUTABLE(readpnetcdf "${READPNETCDF_SOURCES}")
    ADD_DEPENDENCIES(readpnetcdf pidx)
    TARGET_LINK_LIBRARIES(readpnetcdf pidx ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})
    ADD_DEPENDENCIES(readpnetcdf ${PNETCDF_LIBRARY})
    TARGET_LINK_LIBRARIES(readpnetcdf ${PNETCDF_LIBRARY})
  ENDIF()

ENDIF()

