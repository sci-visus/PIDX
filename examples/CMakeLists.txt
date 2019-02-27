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

IF (PIDX_BUILD_EXAMPLES)
  
  # ////////////////////////////////////////
  # includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pidx)
  SET(EXAMPLES_LINK_LIBS pidx ${PIDX_LINK_LIBS}) 
  IF (MPI_CXX_FOUND)
    MESSAGE("Configuring examples with MPI support")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${MPI_C_LIBRARIES})
  ENDIF ()

  IF (ZFP_FOUND)
    #MESSAGE("Configuring examples with ZFP support")
    INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
    link_directories(${ZFP_LIB_DIR})
    SET(EXAMPLES_LINK_LIBS ${EXAMPLES_LINK_LIBS} ${ZFP_LIBRARIES})
  ENDIF ()

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/examples)
  SET(EXAMPLES_UTILS pidx_examples_utils.h)

  # ////////////////////////////////////////
  # executable
  # ////////////////////////////////////////

  PIDX_ADD_CEXECUTABLE(idx_write "grids/idx_write.c" ${EXAMPLES_UTILS})
  TARGET_LINK_LIBRARIES(idx_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(raw_write "grids/raw_write.c" ${EXAMPLES_UTILS})
  TARGET_LINK_LIBRARIES(raw_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_read "grids/idx_read.c" ${EXAMPLES_UTILS})
  TARGET_LINK_LIBRARIES(idx_read ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(particle_write "particles/particle_write.c")
  TARGET_LINK_LIBRARIES(particle_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(particle_non_uniform_write "particles/particle_non_uniform_write.c")
  TARGET_LINK_LIBRARIES(particle_non_uniform_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(particle_gradient_write "particles/particle_gradient_write.c")
  TARGET_LINK_LIBRARIES(particle_gradient_write ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(particle_read "particles/particle_read.c")
  TARGET_LINK_LIBRARIES(particle_read ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_write_partitioned "grids/idx_write_partitioned.c" ${EXAMPLES_UTILS})
  TARGET_LINK_LIBRARIES(idx_write_partitioned ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_write_compressed "grids/idx_write_compressed.c")
  TARGET_LINK_LIBRARIES(idx_write_compressed ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_write_multibuffer "grids/idx_write_multibuffer.c")
  TARGET_LINK_LIBRARIES(idx_write_multibuffer ${EXAMPLES_LINK_LIBS})

  PIDX_ADD_CEXECUTABLE(idx_checkpoint_restart "grids/idx_checkpoint_restart.cpp")
  TARGET_LINK_LIBRARIES(idx_checkpoint_restart ${EXAMPLES_LINK_LIBS})

# Experimental reduced resolution
#  PIDX_ADD_CEXECUTABLE(idx_write_reduce_resolution "experimental/reduced_resolution/idx_write_reduced_resolution.c")
#  TARGET_LINK_LIBRARIES(idx_write_reduce_resolution ${EXAMPLES_LINK_LIBS})

#  PIDX_ADD_CEXECUTABLE(idx_read_reduce_resolution "experimental/reduced_resolution/idx_read_reduced_resolution.c")
#  TARGET_LINK_LIBRARIES(idx_read_reduce_resolution ${EXAMPLES_LINK_LIBS})

ENDIF ()

