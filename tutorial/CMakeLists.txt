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

IF (PIDX_OPTION_HDF5)

#  SET(HDF_LINK_LIBS ${HDF5_LIBRARIES} ${MPI_C_LIBRARIES})
#  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
#  ADD_EXECUTABLE(hdf5_write "HDF5/test-HDF5-writer.c")
#  TARGET_LINK_LIBRARIES(hdf5_write ${HDF5_C_LIBRARIES})

ENDIF()

IF (PIDX_BUILD_TUTORIAL)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(TUTORIAL_LINK_LIBS pidx)
  IF (MPI_C_FOUND)
    MESSAGE("Configuring tutorials with MPI support")
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF()

  IF (ZFP_FOUND)
    MESSAGE("Configuring tutorials with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF()

  IF (PMT_FOUND)
    MESSAGE("Configuring pidx with PMT XX support")
    include_directories(${DATAFLOW_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${PMT_INSITU_INCLUDE_DIR})
#    SET(TUTORIAL_LINK_LIBS ${TUTORIAL_LINK_LIBS} ${DATAFLOW_LIBRARIES} ${PMT_INSITU_LINK_LIBRARIES})
  ENDIF()

  IF (VTK_RENDERER_FOUND)
    INCLUDE_DIRECTORIES(${VTK_RENDERER_INCLUDE_DIR})
  ENDIF()

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////


#  PIDX_ADD_CEXECUTABLE(multi_buffer_single_variable_read "basic_io/read/multi_buffer_single_variable_read.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_single_variable_read ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_single_variable_read "basic_io/read/single_buffer_single_variable_read.c")
#  TARGET_LINK_LIBRARIES(single_buffer_single_variable_read ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(multi_buffer_all_variable_read "basic_io/read/multi_buffer_all_variable_read.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_all_variable_read ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_all_variable_read "basic_io/read/single_buffer_all_variable_read.c")
#  TARGET_LINK_LIBRARIES(single_buffer_all_variable_read ${TUTORIAL_LINK_LIBS})



#  PIDX_ADD_CEXECUTABLE(multi_buffer_idx_write "basic_io/write/idx/multi_buffer_idx_write.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_idx_write ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_idx_write "basic_io/write/idx/single_buffer_idx_write.c")
#  TARGET_LINK_LIBRARIES(single_buffer_idx_write ${TUTORIAL_LINK_LIBS})




#  PIDX_ADD_CEXECUTABLE(multi_buffer_raw_write "basic_io/write/raw/multi_buffer_raw_write.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_raw_write ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_raw_write "basic_io/write/raw/single_buffer_raw_write.c")
#  TARGET_LINK_LIBRARIES(single_buffer_raw_write ${TUTORIAL_LINK_LIBS})





#  PIDX_ADD_CEXECUTABLE(multi_buffer_global_partitioned_idx_write "basic_io/write/global_partitioned_idx/multi_buffer_global_partitioned_idx_write.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_global_partitioned_idx_write ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_global_partitioned_idx_write "basic_io/write/global_partitioned_idx/single_buffer_global_partitioned_idx_write.c")
#  TARGET_LINK_LIBRARIES(single_buffer_global_partitioned_idx_write ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(two_var "basic_io/write/global_partitioned_idx/two_var.c")
#  TARGET_LINK_LIBRARIES(two_var ${TUTORIAL_LINK_LIBS})



#  PIDX_ADD_CEXECUTABLE(multi_buffer_local_partitioned_idx_write "basic_io/write/local_partitioned_idx/multi_buffer_local_partitioned_idx_write.c")
#  TARGET_LINK_LIBRARIES(multi_buffer_local_partitioned_idx_write ${TUTORIAL_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(single_buffer_local_partitioned_idx_write "basic_io/write/local_partitioned_idx/single_buffer_local_partitioned_idx_write.c")
#  TARGET_LINK_LIBRARIES(single_buffer_local_partitioned_idx_write ${TUTORIAL_LINK_LIBS})


ENDIF()
