/*
 * BSD 3-Clause License
 * 
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
#include "../PIDX_inc.h"


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
PIDX_return_code PIDX_err_header                        = 28;
PIDX_return_code PIDX_err_wavelet                       = 29;
PIDX_return_code PIDX_err_metadata                      = 30;
