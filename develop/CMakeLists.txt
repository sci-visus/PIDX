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

IF (PIDX_BUILD_DEVELOP)

  # ////////////////////////////////////////
  # sources
  # ////////////////////////////////////////

  #FILE(GLOB EXAMPLES_SOURCES *.h *.c *.cpp)
  
  # ////////////////////////////////////////
  # includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(EXAMPLES_LINK_LIBS pidx) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring examples with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${MPI_CXX_LIBRARIES})
  ENDIF ()

  IF (ZFP_FOUND)
    MESSAGE("Configuring examples with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF ()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

#  SET(S3D_CHECKPOINT_SOURCES checkpoint-restart/PIDX-s3d-checkpoint.c)
#  PIDX_ADD_CEXECUTABLE(pidx-s3d-checkpoint "${S3D_CHECKPOINT_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-s3d-checkpoint ${EXAMPLES_LINK_LIBS})

#  SET(S3D_RESTART_SOURCES checkpoint-restart/PIDX-s3d-restart.c)
#  PIDX_ADD_CEXECUTABLE(pidx-s3d-restart "${S3D_RESTART_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-s3d-restart ${EXAMPLES_LINK_LIBS})

#  SET(UINTAH_CHECKPOINT_SOURCES checkpoint-restart/PIDX-uintah-checkpoint.c)
#  PIDX_ADD_CEXECUTABLE(pidx-uintah-checkpoint "${UINTAH_CHECKPOINT_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-uintah-checkpoint ${EXAMPLES_LINK_LIBS})

#  SET(UINTAH_RESTART_SOURCES checkpoint-restart/PIDX-uintah-restart.c)
#  PIDX_ADD_CEXECUTABLE(pidx-uintah-restart "${UINTAH_RESTART_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-uintah-restart ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_ROI_SOURCES ROI/PIDX-roi.c)
#  PIDX_ADD_CEXECUTABLE(pidx-roi "${PIDX_ROI_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-roi ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_REDUCED_RESOLUTION_WRITE_SOURCES ROI/PIDX-reduced-resolution-write.c)
#  PIDX_ADD_CEXECUTABLE(pidx-reduced-resolution-write "${PIDX_REDUCED_RESOLUTION_WRITE_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-reduced-resolution-write ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_REDUCED_RESOLUTION_READ_SOURCES ROI/PIDX-reduced-resolution-read.c)
#  PIDX_ADD_CEXECUTABLE(pidx-reduced-resolution-read "${PIDX_REDUCED_RESOLUTION_READ_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-reduced-resolution-read ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_ADAPTIVE_ROI_SOURCES ROI/PIDX-adaptive-roi.c)
#  PIDX_ADD_CEXECUTABLE(pidx-adaptive-roi "${PIDX_ADAPTIVE_ROI_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-adaptive-roi ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_WRITE_COMPRESSED_SOURCES compression/PIDX-write-compressed.c)
#  PIDX_ADD_CEXECUTABLE(pidx-write-compressed "${PIDX_WRITE_COMPRESSED_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-write-compressed ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_READ_COMPRESSED_SOURCES compression/PIDX-read-compressed.c)
#  PIDX_ADD_CEXECUTABLE(pidx-read-compressed "${PIDX_READ_COMPRESSED_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-read-compressed ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_CHANGE_FILE_COUNT_SOURCES file-variations/PIDX-change-file-count.c)
#  PIDX_ADD_CEXECUTABLE(pidx-change-file-count "${PIDX_CHANGE_FILE_COUNT_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-change-file-count ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_WRITE_BY_FIELD_SOURCES field-data-access/PIDX-write-by-field.c)
#  PIDX_ADD_CEXECUTABLE(pidx-write-by-field "${PIDX_WRITE_BY_FIELD_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-write-by-field ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_READ_BY_FIELD_SOURCES field-data-access/PIDX-read-by-field.c)
#  PIDX_ADD_CEXECUTABLE(pidx-read-by-field "${PIDX_READ_BY_FIELD_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-read-by-field ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_PARAMETER_STUDY_SOURCES parameter-study/PIDX-parameter-study.c)
#  PIDX_ADD_CEXECUTABLE(pidx-parameter-study "${PIDX_PARAMETER_STUDY_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-parameter-study ${EXAMPLES_LINK_LIBS})

#  SET(PIDX_WRITE_SIMULATE_SOURCES parameter-study/PIDX-write-simulate.c)
#  PIDX_ADD_CEXECUTABLE(pidx-write-simulate "${PIDX_WRITE_SIMULATE_SOURCES}")
#  TARGET_LINK_LIBRARIES(pidx-write-simulate ${EXAMPLES_LINK_LIBS})

#  IF (PIDX_BUILD_CONVERT)

#    SET(PNETCDF_TO_IDX_SOURCES convert/pnetcdf-to-idx.c)  
#    SET(HDF5_TO_IDX_SOURCES convert/hdf5-to-idx.c)
#    SET(FILE_PER_PROCESS_SOURCES convert/file-per-process-to-idx.c)
#    SET(RAW_TO_IDX_SOURCES convert/raw-to-idx.c)

#    PIDX_ADD_EXECUTABLE(file-per-process-to-idx "${FILE_PER_PROCESS_SOURCES}")
#    TARGET_LINK_LIBRARIES(file-per-process-to-idx ${EXAMPLES_LINK_LIBS})

#    PIDX_ADD_EXECUTABLE(raw-to-idx "${RAW_TO_IDX_SOURCES}")
#    TARGET_LINK_LIBRARIES(raw-to-idx ${EXAMPLES_LINK_LIBS})

#    if (PIDX_OPTION_HDF5)
#      if (HDF5_FOUND)
#        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
#        PIDX_ADD_EXECUTABLE(hdf5_to_idx "${HDF5_TO_IDX_SOURCES}")
#        TARGET_LINK_LIBRARIES(hdf5_to_idx pidx ${HDF5_LIBRARIES})
#      else()
#         message(ERROR ": HDF5 library required to compile ${HDF5_TO_IDX_SOURCES}")
#      endif ()
          
#    ENDIF ()

#    if (PIDX_OPTION_PNETCDF)
#      if (PNETCDF_FOUND)
#        INCLUDE_DIRECTORIES(${PNETCDF_INCLUDE_DIR})

#        PIDX_ADD_EXECUTABLE(pnetcdf_to_idx "${PNETCDF_TO_IDX_SOURCES}")
#        TARGET_LINK_LIBRARIES(pnetcdf_to_idx pidx)  

#        TARGET_LINK_LIBRARIES(pnetcdf_to_idx ${PNETCDF_LIBRARIES})
#      else()
#         message(ERROR ": PnetCDF library required to compile ${PNETCDF_TO_IDX_SOURCES}")
#      endif ()
#    ENDIF ()  

#  ENDIF ()


  PIDX_ADD_CEXECUTABLE(single_buffer_single_variable_read "basic-io/read/single_buffer_single_variable_read.c")
  TARGET_LINK_LIBRARIES(single_buffer_single_variable_read ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_raw_write "basic-io/write/raw/single_buffer_raw_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_raw_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_idx_write "basic-io/write/idx/single_buffer_idx_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_idx_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(multi_buffer_idx_write "basic-io/write/idx/multi_buffer_idx_write.c")
  TARGET_LINK_LIBRARIES(multi_buffer_idx_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(single_buffer_local_partitioned_idx_write "basic-io/write/local_partitioned_idx/single_buffer_local_partitioned_idx_write.c")
  TARGET_LINK_LIBRARIES(single_buffer_local_partitioned_idx_write ${EXAMPLES_LINK_LIBS})




  #PIDX_ADD_CEXECUTABLE(multi_buffer_local_partitioned_idx_write "basic-io/write/local_partitioned_idx/multi_buffer_local_partitioned_idx_write.c")
  #TARGET_LINK_LIBRARIES(multi_buffer_local_partitioned_idx_write ${EXAMPLES_LINK_LIBS})

  #PIDX_ADD_CEXECUTABLE(single_buffer_global_partitioned_idx_write "basic-io/write/global_partitioned_idx/single_buffer_global_partitioned_idx_write.c")
  #TARGET_LINK_LIBRARIES(single_buffer_global_partitioned_idx_write ${EXAMPLES_LINK_LIBS})

  #PIDX_ADD_CEXECUTABLE(multi_buffer_global_partitioned_idx_write "basic-io/write/global_partitioned_idx/multi_buffer_global_partitioned_idx_write.c")
  #TARGET_LINK_LIBRARIES(multi_buffer_global_partitioned_idx_write ${EXAMPLES_LINK_LIBS})


ENDIF ()
