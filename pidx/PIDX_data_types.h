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
 * \file PIDX_data_types.h
 *
 * PIDX common definitions and constants.
 *
 */

#ifndef __PIDX_DATA_TYPES_H
#define __PIDX_DATA_TYPES_H



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


#endif
