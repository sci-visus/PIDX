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
