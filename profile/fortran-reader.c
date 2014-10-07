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

#include <string.h>
#include <strings.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <arpa/inet.h>
 #include <stdio.h>
#include <math.h>

int main(int argc, char **argv) 
{
      int i, j, ret = 0, val = 0;
      FILE* idx_file;
      int fd;
      int data_size;
      int data_offset, counter = 0;
      double* data_buffer = NULL;
      
      fd = open("/home/sid/Research/FINAL-PIDX/Working-PIDX/parallel_idx/Working_copy/working_copy/source/PIDX-lib/d-1.0000E-09/field.00000", O_RDONLY); 
      if(fd < 0)
      {
	  fprintf(stderr, "[File : %s] [Line : %d] open\n",  __FILE__, __LINE__);
      } 
      
      data_offset = 0;
      data_size = 22*36*22 * 16 * sizeof(double);
      
      data_buffer = (double*)malloc(data_size);
      assert(data_buffer != NULL);
      memset(data_buffer, 0, data_size);
      
      ret = pread(fd, data_buffer, data_size, data_offset);
      printf("ReadLines = %d\n", ret);
      //assert(ret == data_size);
      
      for(val = 0; val < ret / sizeof(double) ; val++)
      {
	  //printf("Value at %d is %f\n", val, data_buffer[val]);
	  if(data_buffer[val] == (double)0)
	    counter++;
      }
      printf("Eleemnts with 0 value = %d\n", counter);
      free(data_buffer);
      data_buffer = 0;
      
}
