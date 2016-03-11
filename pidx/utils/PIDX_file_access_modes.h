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


typedef const unsigned int PIDX_flags;
/*
/// Error creating a file that already exists.
extern PIDX_flags PIDX_file_excl;

/// Create the file if it does not exist.
extern PIDX_flags PIDX_MODE_CREATE;

*/
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


#endif
