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
PIDX_return_code PIDX_err_name                          = 4;
PIDX_return_code PIDX_err_box                           = 8;
PIDX_return_code PIDX_err_file                          = 16;
PIDX_return_code PIDX_err_time                          = 32;
PIDX_return_code PIDX_err_block                         = 64;
PIDX_return_code PIDX_err_comm                          = 128;
PIDX_return_code PIDX_err_count                         = 256;
PIDX_return_code PIDX_err_size                          = 512;
PIDX_return_code PIDX_err_offset                        = 1024;
PIDX_return_code PIDX_err_type                          = 2048;
PIDX_return_code PIDX_err_variable                      = 4096;
PIDX_return_code PIDX_err_not_implemented               = 8192;
PIDX_return_code PIDX_err_point                         = 16384;
PIDX_return_code PIDX_err_access                        = 32768;
PIDX_return_code PIDX_err_id                            = 65536;
PIDX_return_code PIDX_err_mpi                           = 131072;
PIDX_return_code PIDX_err_unsupported_compression_type  = 262144;
PIDX_return_code PIDX_err_rst                           = 262144*2;
PIDX_return_code PIDX_err_compress= 262144*4;
PIDX_return_code PIDX_err_hz= 262144*8;
PIDX_return_code PIDX_err_agg= 262144*16;
PIDX_return_code PIDX_err_io= 262144*32;
PIDX_return_code PIDX_err_chunk= 262144*64;
