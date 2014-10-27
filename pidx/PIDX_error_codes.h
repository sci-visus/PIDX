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
#ifndef __PIDX_ERROR_CODES_H
#define __PIDX_ERROR_CODES_H

typedef unsigned int PIDX_return_code;

extern PIDX_return_code PIDX_success;
extern PIDX_return_code PIDX_err_unsopperted_flags;
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

#endif
