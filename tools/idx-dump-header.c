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

#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <inttypes.h>

/**
 * \file idx-dump-header.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Dumps header data
 *
 */

int main(int argc, char **argv) 
{
  int j, var;
  int fd1;
  int blocks_per_file = 0, variable_count = 0;
  uint64_t binheader_count;
  uint32_t *binheader1;
  uint64_t data_length1;
  uint64_t data_offset1;
  
  if ( argv[1] == NULL || argv[2] == NULL || argv[3] == NULL )
  {
    fprintf(stderr, "Missing arguments\n");
    return 0;
  }
  
  fd1 = open(argv[1], O_RDONLY); 
  if (fd1 < 0)
  {
    fprintf(stderr, "Invalid input binary file %s\n", argv[1]);
    return(-1);
  }
  
  blocks_per_file = atoi(argv[2]);
  variable_count = atoi(argv[3]);
  
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader1 = malloc(sizeof(*binheader1)*(binheader_count));
  
  uint64_t rc = 0;
  rc = read(fd1, binheader1, (sizeof(*binheader1) * binheader_count));
  if (rc != (sizeof(*binheader1) * binheader_count))
  {
    fprintf(stderr, "Error reading header 1\n");
    exit(0);
  }
  
  for (var = 0; var < variable_count; var++ )
  {
    for (j = 0 ; j < blocks_per_file; j++)
    {
      data_length1 = 0;
      data_offset1 = 0;
      data_length1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 14]);
      data_offset1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 12]);
      fprintf(stderr, "[%d] Block %d: Offset %lld Count %lld\n", var, j, (unsigned long long)data_offset1, (unsigned long long)data_length1);

      double *data_buffer = malloc(data_length1);
      rc = pread(fd1, data_buffer, data_length1, data_offset1);
      //if (rc != data_length1)
      //{
      //  fprintf(stderr, "Error reading header 2 %d %d\n", (int)rc, (int)data_length1);
      //  exit(0);
      //}
      int k = 0;
      for (k = 0; k < data_length1/sizeof(double); k++)
        fprintf(stderr, "Value at %d = %f\n", k, data_buffer[k]);
      
    }
  }
  
  free(binheader1);
  binheader1 = 0;
  close(fd1);
  
  return(0);
}
