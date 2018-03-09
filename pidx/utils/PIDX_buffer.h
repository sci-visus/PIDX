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
 
#ifndef __PIDX_BUFFER_H
#define __PIDX_BUFFER_H

// A resizeable buffer, similar to std::vector
typedef struct PIDX_buffer {
  char *buffer;
  size_t size;
  size_t capacity;
} PIDX_buffer;

PIDX_buffer PIDX_buffer_create_empty();
PIDX_buffer PIDX_buffer_create_with_capacity(size_t capacity);
void PIDX_buffer_free(PIDX_buffer *b);
void PIDX_buffer_append(PIDX_buffer *b, const char *data, const size_t size);
void PIDX_buffer_resize(PIDX_buffer *b, const size_t size);

#endif

