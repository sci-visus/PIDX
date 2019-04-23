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
 * \file idx-compare.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Compares samples between two idx binary file, assumes that
 * the data is in double format
 *
 */

int main(int argc, char **argv) 
{
  int i, j, var;
  int fd1, fd2;
  int blocks_per_file = 0, variable_count = 0;
  uint64_t binheader_count;
  uint32_t *binheader1, *binheader2;
  uint64_t data_length1, data_length2;
  uint64_t data_offset1, data_offset2;
  double *data_buffer1, *data_buffer2;
  uint64_t *zero_count = 0, *non_zero_equal_count = 0, *non_zero_unequal_count = 0;

  if ( argv[1] == NULL || argv[2] == NULL || argv[3] == NULL || argv[4] == NULL )
  {
    fprintf(stderr, "Missing arguments\n");
    return 0;
  }
  
  fd1 = open(argv[1], O_RDONLY | O_BINARY);
  if (fd1 < 0)
  {
    fprintf(stderr, "Invalid input binary file %s\n", argv[1]);
    return(-1);
  }
  
  fd2 = open (argv[2], O_RDONLY | O_BINARY);
  if (fd2 < 0)
  {
    fprintf(stderr, "Invalid input binary file %s\n", argv[1]);
    return(-1);
  }
  
  blocks_per_file = atoi(argv[3]);
  variable_count = atoi(argv[4]);
  
  zero_count = malloc(sizeof(uint64_t) * variable_count);
  non_zero_equal_count = malloc(sizeof(uint64_t) * variable_count);
  non_zero_unequal_count = malloc(sizeof(uint64_t) * variable_count);
  memset(zero_count, 0, sizeof(uint64_t) * variable_count);
  memset(non_zero_equal_count, 0, sizeof(uint64_t) * variable_count);
  memset(non_zero_unequal_count, 0, sizeof(uint64_t) * variable_count);
  
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader1 = malloc(sizeof(*binheader1)*(binheader_count));
  binheader2 = malloc(sizeof(*binheader2)*(binheader_count));
  
  uint64_t rc = 0;
  rc = read(fd1, binheader1, (sizeof(*binheader1) * binheader_count));
  if (rc != (sizeof(*binheader1) * binheader_count))
  {
    //fprintf(stderr, "Error reading header %d %d\n", rc, (sizeof(*binheader1) * binheader_count));
    exit(0);
  }

  rc = read(fd2, binheader2, (sizeof(*binheader2) * binheader_count));
  if (rc != (sizeof(*binheader2) * binheader_count))
  {
    //fprintf(stderr, "Error reading header %d %d\n", rc, (sizeof(*binheader2) * binheader_count));
    exit(0);
  }
  
  for (var = 0; var < variable_count; var++ )
  {
    for (j = 0 ; j < blocks_per_file; j++)
    {
      data_length1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 14]);
      data_offset1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 12]);
      data_buffer1 = (double*)malloc(data_length1);
      memset(data_buffer1, 0, data_length1);
      assert(data_buffer1 != NULL);
      
      data_length2 = ntohl(binheader2[(j + blocks_per_file * var) * 10 + 14]);
      data_offset2 = ntohl(binheader2[(j + blocks_per_file * var) * 10 + 12]);
      data_buffer2 = (double*)malloc(data_length2);
      memset(data_buffer2, 0, data_length2);
      assert(data_buffer2 != NULL);

      //fprintf(stderr, "%ld %ld :: %ld %ld\n", data_length1, data_offset1, data_length2, data_offset2);
      
      int ret2;
      ret2 = pread(fd1, data_buffer1, data_length1, data_offset1);
      //fprintf(stderr, "I O %d %d\n", ret2, data_length);
      assert(ret2 == data_length1);

      ret2 = pread(fd2, data_buffer2, data_length2, data_offset2);
      //fprintf(stderr, "I O %d %d\n", ret2, data_length2);
      assert(ret2 == data_length2);
      
      if (data_length1 != ret2)
      {
        fprintf(stderr, "Different block length for block %d (%ld %ld)\n", j, data_length1, data_length2);
        return -1;
      }
      
      for (i = 0; i < data_length1 / sizeof(double); i++)
      {
        if (data_buffer1[i] == 0 && data_buffer2[i] == 0)
          zero_count[var]++;
        else
        {
          if (data_buffer1[i] == data_buffer2[i])
          {
            non_zero_equal_count[var]++;
            //fprintf(stderr, "[R] Values in Block %d at index %d = %f %f\n", j, i, data_buffer1[i], data_buffer2[i]);
          }
          else
          {
            non_zero_unequal_count[var]++;
            //fprintf(stderr, "[W] Values in Block %d at index %d = %f %f\n", j, i, data_buffer1[i], data_buffer2[i]);
          }
        }
      }
      free(data_buffer1);
      free(data_buffer2);
    }
  }
  
  free(binheader1);
  binheader1 = 0;
  free(binheader2);
  binheader2 = 0;
  close(fd1);
  close(fd2);

  for (var = 0; var < variable_count; var++ )
    fprintf(stderr, "[%d] Zero Equal %lld Non-Zero Equal %lld Non-Zero Non Equal %lld\n", var, (long long)zero_count[var], (long long)non_zero_equal_count[var], (long long)non_zero_unequal_count[var]);
  
  return(0);
}
