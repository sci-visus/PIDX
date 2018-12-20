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
#ifndef __PIDX_WINDOWS_DEFINE_H
#define __PIDX_WINDOWS_DEFINE_H

/** Macros to swap bytes in a multi-byte value, to convert from big-endian data
to little-endian data and vice versa. These are taken from the Boost library.*/
#ifndef __has_builtin
#define __has_builtin(x) 0 // Compatibility with non-clang compilers
#endif
#define htonl(x) _byteswap_ulong(x)
#define ntohl(x) htonl(x)

#define inline __inline
#define snprintf _snprintf

#include <io.h>
#include <stdio.h>
#include <string.h>
#include <basetsd.h>

typedef long int __uint64_t;
typedef long int __off64_t;
//typedef __uint64_t uint64_t;

//#define uint64_t uint64_t
#define PATH_MAX 256
#define mkdir(x,y) mkdir(x) // This returns error if directory exists already
#define PIDX_HAVE_MPI 1  // TODO: remove this later cause it is going to be deprecated

inline int pwrite(int fd, const void *buf, uint64_t nbytes, uint64_t offset)
{
  long ret = _lseek(fd, offset, SEEK_SET);

  if (ret == -1) {
    return(-1);
  }
  return(_write(fd, buf, nbytes));
}

inline int pread(int fd, void *buf, uint64_t nbytes, uint64_t offset)
{
  if (_lseek(fd, offset, SEEK_SET) != offset) {
    return -1;
  }
  return read(fd, buf, nbytes);
}

#endif
