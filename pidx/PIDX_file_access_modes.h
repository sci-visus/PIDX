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
 * \author Sidharth Kumar
 * PIDX file access modes
 *
 */

#ifndef __PIDX_FILE_ACCESS_MODES_H
#define __PIDX_FILE_ACCESS_MODES_H


typedef const uint8_t PIDX_flags;

/// Error creating a file that already exists.
extern PIDX_flags PIDX_file_excl;

/// Create the file if it does not exist.
extern PIDX_flags PIDX_file_trunc;

/// Read only.
extern PIDX_flags PIDX_file_rdwr;

/// Reading and writing.
extern PIDX_flags PIDX_file_rdonly;

/// Write only.
extern PIDX_flags PIDX_file_wronly;

#endif
