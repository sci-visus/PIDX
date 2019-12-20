/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

PIDX_buffer PIDX_buffer_create_empty()
{
  PIDX_buffer b = { NULL, 0, 0 };
  return b;
}
PIDX_buffer PIDX_buffer_create_with_capacity(uint64_t capacity)
{
  PIDX_buffer b = { malloc(capacity), 0, capacity };
  return b;
}
void PIDX_buffer_free(PIDX_buffer *b)
{
  free(b->buffer);
  b->size = 0;
  b->capacity = 0;
}
void PIDX_buffer_append(PIDX_buffer *b, const unsigned char *data, const uint64_t size)
{
  if (b->capacity - b->size < size) {
    b->capacity = Max2ab(b->capacity + size, b->capacity * 1.5);
    b->buffer = realloc(b->buffer, b->capacity);
  }
  memcpy(b->buffer + b->size, data, size);
  b->size += size;
}
void PIDX_buffer_resize(PIDX_buffer *b, const uint64_t size)
{
  if (b->capacity < size) {
    b->capacity = size;
    b->buffer = realloc(b->buffer, b->capacity);
  }
  b->size = size;
}

