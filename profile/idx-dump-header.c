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
  int i, j, var;
  int fd1;
  int blocks_per_file = 0, variable_count = 0;
  size_t binheader_count;
  uint32_t *binheader1;
  size_t data_length1;
  off_t data_offset1;
  
  if ( argv[1] == NULL || argv[2] == NULL || argv[3] == NULL )
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
  
  blocks_per_file = atoi(argv[2]);
  variable_count = atoi(argv[3]);
  
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader1 = malloc(sizeof(*binheader1)*(binheader_count));
  
  read(fd1, binheader1, (sizeof(*binheader1) * binheader_count));
  
  for (var = 0; var < variable_count; var++ )
  {
    for (j = 0 ; j < blocks_per_file; j++)
    {
      data_length1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 14]);
      data_offset1 = ntohl(binheader1[(j + blocks_per_file * var) * 10 + 12]);
      printf("[%d] Block %d: Offset %ld Count %ld\n", var, j, data_offset1, data_length1);
      
    }
  }
  
  free(binheader1);
  binheader1 = 0;
  close(fd1);
  
  return(0);
}