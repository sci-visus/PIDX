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
  size_t binheader_count;
  uint32_t *binheader1, *binheader2;
  size_t data_length1, data_length2;
  off_t data_offset1, data_offset2;
  double *data_buffer1, *data_buffer2;
  long long *zero_count = 0, *non_zero_equal_count = 0, *non_zero_unequal_count = 0;

  if ( argv[1] == NULL || argv[2] == NULL || argv[3] == NULL || argv[4] == NULL )
  {
    fprintf(stderr, "Missing arguments\n");
    return 0;
  }
  
  fd1 = open(argv[1], O_RDONLY); 
  if(fd1 < 0)
  {
    fprintf(stderr, "Invalid input binary file %s\n", argv[1]);
    return(-1);
  }
  
  fd2 = open (argv[2], O_RDONLY); 
  if (fd2 < 0)
  {
    fprintf(stderr, "Invalid input binary file %s\n", argv[1]);
    return(-1);
  }
  
  blocks_per_file = atoi(argv[3]);
  variable_count = atoi(argv[4]);
  
  zero_count = malloc(sizeof(long long) * variable_count);
  non_zero_equal_count = malloc(sizeof(long long) * variable_count);
  non_zero_unequal_count = malloc(sizeof(long long) * variable_count);
  memset(zero_count, 0, sizeof(long long) * variable_count);
  memset(non_zero_equal_count, 0, sizeof(long long) * variable_count);
  memset(non_zero_unequal_count, 0, sizeof(long long) * variable_count);
  
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader1 = malloc(sizeof(*binheader1)*(binheader_count));
  binheader2 = malloc(sizeof(*binheader2)*(binheader_count));
  
  read(fd1, binheader1, (sizeof(*binheader1) * binheader_count));
  read(fd2, binheader2, (sizeof(*binheader2) * binheader_count));
  
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

      pread(fd1, data_buffer1, data_length1, data_offset1);
      pread(fd2, data_buffer2, data_length2, data_offset2);
      
      if(data_length1 != data_length2)
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
	    non_zero_equal_count[var]++;
	  else
	  {
	    non_zero_unequal_count[var]++;
	    printf("Values in Block %d at index %d = %f %f\n", j, i, data_buffer1[i], data_buffer2[i]);
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
    printf("[%d] Zero Equal %lld Non-Zero Equal %lld Non-Zero Non Equal %lld\n", var, zero_count[var], non_zero_equal_count[var], non_zero_unequal_count[var]);
  
  return(0);
}