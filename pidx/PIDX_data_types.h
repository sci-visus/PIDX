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

extern PIDX_data_type FLOAT64;
extern PIDX_data_type FLOAT64_GA;
extern PIDX_data_type FLOAT64_RGB;
extern PIDX_data_type FLOAT64_RGBA;

PIDX_return_code PIDX_default_bits_per_datatype(PIDX_data_type type, int* bits);

#endif
