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


ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(comm)
ADD_SUBDIRECTORY(data_handle)
ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(io)
ADD_SUBDIRECTORY(metadata)

FILE(GLOB PIDX_SOURCES *.h *.c ./utils/*.h ./utils/*.c ./comm/*.h ./comm/*.c ./meta_data_access/*.h ./meta_data_access/*.c ./data_handle/*.h ./data_handle/*.c ./core/PIDX_in_transit_interface/*.h ./core/PIDX_in_transit_interface/*.c ./core/PIDX_agg/*.h ./core/PIDX_agg/*.c ./core/PIDX_block_rst/*.h ./core/PIDX_block_rst/*.c ./core/PIDX_cmp/*.h ./core/PIDX_cmp/*.c ./core/PIDX_file_io/*.h ./core/PIDX_file_io/*.c ./core/PIDX_header/*.h ./core/PIDX_header/*.c ./core/PIDX_hz/*.h ./core/PIDX_hz/*.c ./core/PIDX_idx_rst/*.h ./core/PIDX_idx_rst/*.c ./core/PIDX_raw_rst/*.c ./core/PIDX_raw_rst/*.h ./core/PIDX_wavelet/*.h ./core/PIDX_wavelet/*.c ./io/*.h ./io/*.c ./io/raw/*.c ./io/raw/*.h ./io/idx/*.c ./io/idx/*.h ./io/idx/local_partition/*.c ./io/idx/no_partition/*.c ./io/idx/serial/*.c ./io/idx/particle/*.c metadata/*.h metadata/*.c)

#FILE(GLOB PIDX_SOURCES *.h *.c ./utils/*.h ./utils/*.c ./comm/*.h ./comm/*.c ./data_handle/*.h ./data_handle/*.c ./core/PIDX_agg/PIDX_global_agg/*.h ./core/PIDX_agg/PIDX_global_agg/*.c ./core/PIDX_agg/PIDX_local_agg/*.h ./core/PIDX_agg/PIDX_local_agg/*.c ./core/PIDX_agg/*.h ./core/PIDX_agg/*.c ./core/PIDX_block_rst/*.h ./core/PIDX_block_rst/*.c ./core/PIDX_cmp/*.h ./core/PIDX_cmp/*.c ./core/PIDX_file_io/*.h ./core/PIDX_file_io/*.c ./core/PIDX_header/*.h ./core/PIDX_header/*.c ./core/PIDX_hz/*.h ./core/PIDX_hz/*.c ./core/PIDX_rst/*.h ./core/PIDX_rst/*.c ./io/*.h ./io/*.c ./io/PIDX_global_idx_io/*.h ./io/PIDX_global_idx_io/*.c ./io/PIDX_io/*.h ./io/PIDX_io/*.c ./io/PIDX_idx_io/*.h ./io/PIDX_idx_io/*.c ./io/PIDX_multi_patch_idx_io/*.h ./io/PIDX_multi_patch_idx_io/*.c ./io/PIDX_partitioned_idx_io/*.h ./io/PIDX_partitioned_idx_io/*.c ./io/PIDX_partition_merge_idx_io/*.h ./io/PIDX_partition_merge_idx_io/*.c ./io/PIDX_raw_io/*.h ./io/PIDX_raw_io/*.c)

#FILE(GLOB PIDX_SOURCES *.h *.c)

MESSAGE("using OS_SPECIFIC_LIBS: ${OS_SPECIFIC_LIBS}")

# ////////////////////////////////////////
# includes
# ////////////////////////////////////////

IF (MPI_C_FOUND)
  MESSAGE("Configuring pidx with MPI support")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${OS_SPECIFIC_LIBS} ${MPI_C_LIBRARIES})
ENDIF()

IF (ZFP_FOUND)
  MESSAGE("Configuring pidx with ZFP support")
  INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${ZFP_LIBRARIES})
ENDIF()

IF (PMT_FOUND)
  MESSAGE("Configuring pidx with PMT support YY" ${PMT_INSITU_LINK_LIBRARIES} "ZZ")
  INCLUDE_DIRECTORIES(${PMT_INSITU_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${DATAFLOW_INCLUDE_DIR})
  SET(PIDX_LINK_LIBS ${PIDX_LINK_LIBS} ${PMT_INSITU_LINK_LIBRARIES} ${DATAFLOW_LIBRARY})
ENDIF()

IF (VTK_RENDERER_FOUND)
  INCLUDE_DIRECTORIES(${VTK_RENDERER_INCLUDE_DIR})
ENDIF()

# ////////////////////////////////////////
# library
# ////////////////////////////////////////
PIDX_ADD_LIBRARY(pidx "${PIDX_SOURCES}")
TARGET_LINK_LIBRARIES(pidx ${OS_SPECIFIC_LIBS} pidx_comm pidx_metadata pidx_data_handle pidx_utils pidx_core pidx_io zfp ${PIDX_LINK_LIBS})


SET_TARGET_PROPERTIES(pidx PROPERTIES LINKER_LANGUAGE CXX)
