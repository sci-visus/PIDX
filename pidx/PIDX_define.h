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





// Aggregation modes
// AGG_SETUP - Only sets the aggregation buffers
// AGG_PERFORM - Performs aggregation phase
// AGG_SETUP_AND_PERFORM - Sets aggregation buffers and also performs aggregation
enum AGG_MODES {AGG_SETUP, AGG_PERFORM, AGG_SETUP_AND_PERFORM};


enum WAVELET_MODES {WAVELET_STENCIL, WAVELET_RST};

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
  PIDX_GLOBAL_PARTITION_IDX_IO=1,   /// Writes data in partitioned space with global indexing
  PIDX_LOCAL_PARTITION_IDX_IO=2,    /// Writes data in partitioned space with local indexing
  PIDX_RAW_IO=3,                    /// Writes data in raw format
  PIDX_WAVELET_IO=4                 /// Calls wavelet computation code
};

enum PIDX_endian_type{
  PIDX_BIG_ENDIAN=0,                     /// Use big endianess
  PIDX_LITTLE_ENDIAN=1                   /// Use little endianess
};
  
// Calls merge tree analysis code (in-situ mode)
#define PIDX_MERGE_TREE_ANALYSIS                      5

#define PIDX_default_bits_per_block              15
#define PIDX_default_blocks_per_file             256


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

//double pa1, pa2, pa3, pa4, pa5, pa6;

/// IDX specifies generic types using simple strings consisting of an unambiguous data type and
/// C array syntax, e.g. "float32[3]".  In the PIDX library, we declare types using strings so
/// users can easily override the provided defaults.

typedef char PIDX_data_type[512];

// PLEASE NOTE: these are example types, not a complete list of possible IDX types

extern PIDX_data_type INT8;
extern PIDX_data_type INT8_GA;
extern PIDX_data_type INT8_RGB;
extern PIDX_data_type INT8_RGBA;

extern PIDX_data_type UINT8;
extern PIDX_data_type UINT8_GA;
extern PIDX_data_type UINT8_RGB;
extern PIDX_data_type UINT8_RGBA;

extern PIDX_data_type INT16;
extern PIDX_data_type INT16_GA;
extern PIDX_data_type INT16_RGB;
extern PIDX_data_type INT16_RGBA;

extern PIDX_data_type UINT16;
extern PIDX_data_type UINT16_GA;
extern PIDX_data_type UINT16_RGB;
extern PIDX_data_type UINT16_RGBA;

extern PIDX_data_type INT32;
extern PIDX_data_type INT32_GA;
extern PIDX_data_type INT32_RGB;
extern PIDX_data_type INT32_RGBA;

extern PIDX_data_type UINT32;
extern PIDX_data_type UINT32_GA;
extern PIDX_data_type UINT32_RGB;
extern PIDX_data_type UINT32_RGBA;

extern PIDX_data_type INT64;
extern PIDX_data_type INT64_GA;
extern PIDX_data_type INT64_RGB;
extern PIDX_data_type INT64_RGBA;

extern PIDX_data_type UINT64;
extern PIDX_data_type UINT64_GA;
extern PIDX_data_type UINT64_RGB;
extern PIDX_data_type UINT64_RGBA;

extern PIDX_data_type FLOAT32;
extern PIDX_data_type FLOAT32_GA;
extern PIDX_data_type FLOAT32_RGB;
extern PIDX_data_type FLOAT32_RGBA;
extern PIDX_data_type FLOAT32_9TENSOR;

extern PIDX_data_type FLOAT64;
extern PIDX_data_type FLOAT64_GA;
extern PIDX_data_type FLOAT64_RGB;
extern PIDX_data_type FLOAT64_RGBA;
extern PIDX_data_type FLOAT64_7STENCIL;
extern PIDX_data_type FLOAT64_9TENSOR;

#ifdef __cplusplus
}
#endif

#endif
