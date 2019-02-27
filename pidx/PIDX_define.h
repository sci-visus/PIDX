/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */
#ifndef __PIDX_DEFINE_H
#define __PIDX_DEFINE_H

#ifdef __cplusplus
extern "C" {
#endif


#define PIDX_CURR_METADATA_VERSION "6.1"
  
#define PIDX_MAX_DIMENSIONS 3
#define MULTI_BOX 0
#define PIDX_STRING_SIZE 512

#define PIDX_HAVE_ZFP 1
#define PIDX_HAVE_PMT 0
#define PIDX_HAVE_VTK 0
#define PIDX_HAVE_PNETCDF 0
#define PIDX_HAVE_NETCDF 0
#define PIDX_HAVE_HDF5 0
#define PIDX_HAVE_NVISUSIO 0

#define PIDX_MAX_TEMPLATE_DEPTH 6

#define DETAIL_OUTPUT 0

#define DEBUG_OUTPUT 0

#ifndef __cplusplus
#  define _XOPEN_SOURCE 600
#endif

#ifdef BGQ
  #define _XOPEN_SOURCE 600
#ifndef _GNU_SOURCE
    #define _GNU_SOURCE
#endif
#endif

// PIDX write and read mode
// PIDX_READ - Read only mode
// PIDX_WRITE - Write only mode
enum IO_READ_WRITE {PIDX_READ, PIDX_WRITE};

// If a block (application order or HZ order is on the boundary or not, a block is on the boundary
// if it does not align with a power in to dimension block.
enum boundary_type {power_two_block = 1,
                   non_power_two_block = 2};

// No process dumps any meta data info
#define PIDX_NO_META_DATA_DUMP             0

// Every process writes the MPI related meta data into a seperate file while continuing actual MPI and IO calls
#define PIDX_META_DATA_DUMP_ONLY             1

// Every process writes the MPI related meta data into a seperate file while preventing any actual MPI and IO call
#define PIDX_NO_IO_AND_META_DATA_DUMP        2


#define PIDX_NO_COMPRESSION 0
#define PIDX_CHUNKING_ONLY 1
#define PIDX_CHUNKING_ZFP 2

// Data in buffer is in row order
#define PIDX_row_major                           0

// Data in buffer is in column order
#define PIDX_column_major                        1


enum PIDX_io_type {
  PIDX_IDX_IO=0,                    /// Writes data in IDX format
  PIDX_LOCAL_PARTITION_IDX_IO=1,    /// Writes data in partitioned space with local indexing
  PIDX_RAW_IO=2,                    /// Writes data in raw format

  PIDX_PARTICLE_IO=3,              /// Writes particle data using file per process io
  PIDX_RST_PARTICLE_IO=4           /// Writes particles data after restructuring
};

enum PIDX_endian_type{
  PIDX_BIG_ENDIAN=0,                     /// Use big endianess
  PIDX_LITTLE_ENDIAN=1                   /// Use little endianess
};
  
// Calls merge tree analysis code (in-situ mode)
#define PIDX_MERGE_TREE_ANALYSIS                      5


#define PIDX_default_bits_per_block              15
#define PIDX_default_blocks_per_file             256

#define PIDX_FILE_PATH_LENGTH                    1024

/// Create the file if it does not exist.
#define PIDX_MODE_CREATE              1

/// Error creating a file that already exists.
#define PIDX_MODE_EXCL               64

#define PIDX_MODE_RDONLY              2  /* ADIO_RDONLY */
#define PIDX_MODE_WRONLY              4  /* ADIO_WRONLY  */
#define PIDX_MODE_RDWR                8  /* ADIO_RDWR  */
#define PIDX_MODE_DELETE_ON_CLOSE    16  /* ADIO_DELETE_ON_CLOSE */
#define PIDX_MODE_UNIQUE_OPEN        32  /* ADIO_UNIQUE_OPEN */

#define PIDX_MODE_APPEND            128  /* ADIO_APPEND */
#define PIDX_MODE_SEQUENTIAL        256  /* ADIO_SEQUENTIAL */


/// IDX specifies generic types using simple strings consisting of an unambiguous data type and
/// C array syntax, e.g. "float32[3]".  In the PIDX library, we declare types using strings so
/// users can easily override the provided defaults.

typedef unsigned int PIDX_data_layout;
typedef char PIDX_data_type[512];

// PLEASE NOTE: these are example types, not a complete list of possible IDX types

struct pidx_dtype {
    PIDX_data_type INT8;
    PIDX_data_type INT8_GA;
    PIDX_data_type INT8_RGB;
    PIDX_data_type INT8_RGBA;

    PIDX_data_type UINT8;
    PIDX_data_type UINT8_GA;
    PIDX_data_type UINT8_RGB;
    PIDX_data_type UINT8_RGBA;

    PIDX_data_type INT16;
    PIDX_data_type INT16_GA;
    PIDX_data_type INT16_RGB;
    PIDX_data_type INT16_RGBA;

    PIDX_data_type UINT16;
    PIDX_data_type UINT16_GA;
    PIDX_data_type UINT16_RGB;
    PIDX_data_type UINT16_RGBA;

    PIDX_data_type INT32;
    PIDX_data_type INT32_GA;
    PIDX_data_type INT32_RGB;
    PIDX_data_type INT32_RGBA;

    PIDX_data_type UINT32;
    PIDX_data_type UINT32_GA;
    PIDX_data_type UINT32_RGB;
    PIDX_data_type UINT32_RGBA;

    PIDX_data_type INT64;
    PIDX_data_type INT64_GA;
    PIDX_data_type INT64_RGB;
    PIDX_data_type INT64_RGBA;

    PIDX_data_type UINT64;
    PIDX_data_type UINT64_GA;
    PIDX_data_type UINT64_RGB;
    PIDX_data_type UINT64_RGBA;

    PIDX_data_type FLOAT32;
    PIDX_data_type FLOAT32_GA;
    PIDX_data_type FLOAT32_RGB;
    PIDX_data_type FLOAT32_RGBA;
    PIDX_data_type FLOAT32_7STENCIL;
    PIDX_data_type FLOAT32_9TENSOR;

    PIDX_data_type FLOAT64;
    PIDX_data_type FLOAT64_GA;
    PIDX_data_type FLOAT64_RGB;
    PIDX_data_type FLOAT64_RGBA;
    PIDX_data_type FLOAT64_7STENCIL;
    PIDX_data_type FLOAT64_9TENSOR;

    PIDX_data_type INT64_9TENSOR;
    PIDX_data_type INT32_9TENSOR;
};

extern struct pidx_dtype PIDX_DType;

#ifdef __cplusplus
}
#endif

#endif
