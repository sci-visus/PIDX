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
#include "PIDX_error_codes.h"


/////////////////////////////////////////////////
// ERROR CODES
/////////////////////////////////////////////////

PIDX_return_code PIDX_success                           = 0;
PIDX_return_code PIDX_err_unsupported_flags             = 1;
PIDX_return_code PIDX_err_file_exists                   = 2;
PIDX_return_code PIDX_err_name                          = 3;
PIDX_return_code PIDX_err_box                           = 4;
PIDX_return_code PIDX_err_file                          = 5;
PIDX_return_code PIDX_err_time                          = 6;
PIDX_return_code PIDX_err_block                         = 7;
PIDX_return_code PIDX_err_comm                          = 8;
PIDX_return_code PIDX_err_count                         = 9;
PIDX_return_code PIDX_err_size                          = 10;
PIDX_return_code PIDX_err_offset                        = 11;
PIDX_return_code PIDX_err_type                          = 12;
PIDX_return_code PIDX_err_variable                      = 13;
PIDX_return_code PIDX_err_not_implemented               = 14;
PIDX_return_code PIDX_err_point                         = 15;
PIDX_return_code PIDX_err_access                        = 16;
PIDX_return_code PIDX_err_id                            = 17;
PIDX_return_code PIDX_err_mpi                           = 18;
PIDX_return_code PIDX_err_unsupported_compression_type  = 19;
PIDX_return_code PIDX_err_rst                           = 20;
PIDX_return_code PIDX_err_compress                      = 21;
PIDX_return_code PIDX_err_hz                            = 22;
PIDX_return_code PIDX_err_agg                           = 23;
PIDX_return_code PIDX_err_io                            = 24;
PIDX_return_code PIDX_err_chunk                         = 25;
PIDX_return_code PIDX_err_close                         = 26;
PIDX_return_code PIDX_err_flush                         = 27;
