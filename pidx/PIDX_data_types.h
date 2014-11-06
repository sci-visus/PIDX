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


#ifndef __PIDX_DATA_TYPES_H
#define __PIDX_DATA_TYPES_H


/// String describing a type. PIDX can support types like Float32 (e.g. for pressure) , 3*float64 (velocity vector field)
/// and more general structures. We use a string here so that a user who knows what he/she is doing can put an arbitrary
/// string, say "BOB". Obviously there will be need for a tool that understands the type "BOB".
typedef const char* PIDX_type;

#define PIDX_FLOAT    "1*float32"
#define PIDX_DOUBLE   "1*float64"
#define PIDX_INT      "1*int32"
#define PIDX_LONG     "1*int64"
#define PIDX_2DFLOAT  "2*float32"

#endif
