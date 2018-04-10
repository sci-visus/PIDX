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
add_library(zfp
            bitstream.c decode1f.c decode1d.c encode1f.c encode1d.c decode2f.c
            decode2d.c encode2f.c encode2d.c decode3f.c decode3d.c encode3f.c
            encode3d.c zfp.c)

set(ZFP_HEADERS
    ${CMAKE_SOURCE_DIR}/inc/bitstream.h
    ${CMAKE_SOURCE_DIR}/inc/macros.h
    ${CMAKE_SOURCE_DIR}/inc/system.h
    ${CMAKE_SOURCE_DIR}/inc/types.h
    ${CMAKE_SOURCE_DIR}/inc/zfp.h)

set_target_properties(zfp PROPERTIES
    PUBLIC_HEADER "${ZFP_HEADERS}"
    POSITION_INDEPENDENT_CODE ON)
target_include_directories(zfp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include/zfp>)

if(MSVC)
    target_compile_definitions(zfp PUBLIC _CRT_SECURE_NO_WARNINGS)
    target_compile_options(zfp PUBLIC /W4)
elseif(UNIX OR APPLE OR MINGW)
    target_compile_options(zfp PUBLIC -Wall -Wno-unused-fuction)
endif()

option(BIT_STREAM_STRIDED "enable strided access for progressive zfp streams" FALSE)
option(TESTZFP_MEDIUM_ARRAYS "run regression test with medium sized arrays" FALSE)
option(TESTZFP_LARGE_ARRAYS "run regression test with large arrays" FALSE)
option(WITHOUT_COMPRESSION "run diffusion without compression" FALSE)
option(ALIGNED_ALLOC "use aligned memory allocation" FALSE)
option(CACHE_TWOWAY "use two-way skew-associative cache" FALSE)
option(CACHE_FAST_HASH "use faster but more collision prone hash function" FALSE)
option(CACHE_PROFILE "count cache misses" FALSE)
option(HAVE_C99_MATH "for recent versions of MSVC" FALSE)
set(BIT_STREAM_WORD_TYPE uint64 CACHE STRING "use smaller bit stream word type for finer rate granularity")
set_property(CACHE BIT_STREAM_WORD_TYPE PROPERTY STRINGS uint8 uint16 uint32 uint64)

if(BIT_STREAM_STRIDED)
    target_compile_definitions(zfp PUBLIC BIT_STREAM_STRIDED)
endif()
if(TESTZFP_MEDIUM_ARRAYS)
    target_compile_definitions(zfp PUBLIC TESTZFP_MEDIUM_ARRAYS)
endif()
if(TESTZFP_LARGE_ARRAYS)
    target_compile_definitions(zfp PUBLIC TESTZFP_LARGE_ARRAYS)
endif()
if(WITHOUT_COMPRESSION)
    target_compile_definitions(zfp PUBLIC WITHOUT_COMPRESSION)
endif()
if(ALIGNED_ALLOC)
    target_compile_definitions(zfp PUBLIC ALIGNED_ALLOC)
endif()
if(CACHE_TWOWAY)
    target_compile_definitions(zfp PUBLIC CACHE_TWOWAY)
endif()
if(CACHE_FAST_HASH)
    target_compile_definitions(zfp PUBLIC CACHE_FAST_HASH)
endif()
if(CACHE_PROFILE)
    target_compile_definitions(zfp PUBLIC CACHE_PROFILE)
endif()
if(HAVE_C99_MATH)
    target_compile_definitions(zfp PUBLIC HAVE_C99_MATH)
endif()
target_compile_definitions(zfp PUBLIC BIT_STREAM_WORD_TYPE=${BIT_STREAM_WORD_TYPE})

install(TARGETS zfp
    EXPORT zfp
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)
install(EXPORT zfp DESTINATION lib)

export(EXPORT zfp FILE ${CMAKE_BINARY_DIR}/zfp.cmake)
