/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

/**
 * \file PIDX_typedefs.h
 *
 * PIDX common definitions and constants.
 *
 */

#ifndef __PIDX_TYPEDEFS_H
#define __PIDX_TYPEDEFS_H



/////////////////////////////////////////////////
// IDX DATA TYPES
/////////////////////////////////////////////////

/// IDX specifies generic types using simple strings consisting of an unambiguous data type and
/// C array syntax, e.g. "float32[3]".  In the PIDX library, we declare types using strings so
/// users can easily override the provided defaults.

typedef const char* PIDX_type;

// PLEASE NOTE: these are example types, not a complete list of possible IDX types

#define PIDX_FLOAT32 "float32[1]"
#define PIDX_FLOAT64 "float64[1]"
#define PIDX_UINT8   "uint8[1]"
#define PIDX_UINT16  "uint16[1]"
#define PIDX_UINT32  "uint32[1]"
#define PIDX_UINT64  "uint64[1]"
#define PIDX_INT8    "int8[1]"
#define PIDX_INT16   "int16[1]"
#define PIDX_INT32   "int32[1]"
#define PIDX_INT64   "int64[1]"
#define PIDX_POINT3  "float32[3]"



/////////////////////////////////////////////////
// ERROR CODES
/////////////////////////////////////////////////

typedef int PIDX_return_code;

extern PIDX_return_code PIDX_success;
extern PIDX_return_code PIDX_err_unsupported_flags;
extern PIDX_return_code PIDX_err_file_exists;
extern PIDX_return_code PIDX_err_name;
extern PIDX_return_code PIDX_err_box;
extern PIDX_return_code PIDX_err_file;
extern PIDX_return_code PIDX_err_time;
extern PIDX_return_code PIDX_err_block;
extern PIDX_return_code PIDX_err_comm;
extern PIDX_return_code PIDX_err_count;
extern PIDX_return_code PIDX_err_size;
extern PIDX_return_code PIDX_err_offset;
extern PIDX_return_code PIDX_err_type;
extern PIDX_return_code PIDX_err_variable;
extern PIDX_return_code PIDX_err_not_implemented;
extern PIDX_return_code PIDX_err_point;
extern PIDX_return_code PIDX_err_access;



/////////////////////////////////////////////////
// FILE ACCESS MODES
/////////////////////////////////////////////////

typedef const unsigned PIDX_flags;

extern PIDX_flags PIDX_file_excl;       // Error creating a file that already exists.
extern PIDX_flags PIDX_file_trunc;      // Create the file if it does not exist.
extern PIDX_flags PIDX_file_rdwr;       // Read only.
extern PIDX_flags PIDX_file_rdonly;     // Reading and writing.
extern PIDX_flags PIDX_file_wronly;     // Write only. 



/////////////////////////////////////////////////
// DATA LAYOUT
/////////////////////////////////////////////////

typedef unsigned int PIDX_data_layout;

/// Data layout options 
/// \param PIDX_row_major row major layout
/// \param PIDX_column_major column major layout
///
extern PIDX_data_layout PIDX_row_major;
extern PIDX_data_layout PIDX_column_major;

#endif

