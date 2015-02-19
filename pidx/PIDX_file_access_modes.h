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
 * \file PIDX_file_access_modes.h
 *
 * PIDX common definitions and constants.
 *
 */

#ifndef __PIDX_FILE_ACCESS_MODES_H
#define __PIDX_FILE_ACCESS_MODES_H


/////////////////////////////////////////////////
// FILE ACCESS MODES
/////////////////////////////////////////////////

typedef const unsigned PIDX_flags;

extern PIDX_flags PIDX_file_excl;       // Error creating a file that already exists.
extern PIDX_flags PIDX_file_trunc;      // Create the file if it does not exist.
extern PIDX_flags PIDX_file_rdwr;       // Read only.
extern PIDX_flags PIDX_file_rdonly;     // Reading and writing.
extern PIDX_flags PIDX_file_wronly;     // Write only. 

#endif
