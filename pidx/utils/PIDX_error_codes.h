/*
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */

/**
 * \file PIDX_error_codes.h
 *
 * PIDX common definitions and constants.
 *
 */

#ifndef __PIDX_ERROR_CODES_H
#define __PIDX_ERROR_CODES_H


/////////////////////////////////////////////////
// ERROR CODES
/////////////////////////////////////////////////

typedef int PIDX_return_code;

extern PIDX_return_code PIDX_success;
extern PIDX_return_code PIDX_err_id;
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
extern PIDX_return_code PIDX_err_mpi;
extern PIDX_return_code PIDX_err_rst;
extern PIDX_return_code PIDX_err_chunk;
extern PIDX_return_code PIDX_err_compress;
extern PIDX_return_code PIDX_err_hz;
extern PIDX_return_code PIDX_err_agg;
extern PIDX_return_code PIDX_err_io;
extern PIDX_return_code PIDX_err_unsupported_compression_type;
extern PIDX_return_code PIDX_err_close;
extern PIDX_return_code PIDX_err_flush;
extern PIDX_return_code PIDX_err_header;
extern PIDX_return_code PIDX_err_wavelet;
extern PIDX_return_code PIDX_err_metadata;

#endif
