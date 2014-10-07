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

/* displays samples from a specific block in an idx binary file, assumes that
 * the data is in double format
 */
int main(int argc, char **argv) 
{
    char buffer[256];
    int blocksperfile = -1;
    int ret;
    int fd, fd2, fd3, fd4, fd5;
    uint32_t* binheader;
    uint32_t* binheader1;
    int binheader_count;
    int target_block = 0;
    double* data_buffer = NULL;
    size_t data_size, data_size1;
    off_t data_offset, data_offset1;
    int i, j, k, count = 0;
    int number_of_directory;
    size_t length, length1;
  
    
    fd = open(argv[1], O_RDONLY); 
    if(fd < 0)
    {
        perror("open");
        return(-1);
    }
    printf("File names : %s\n", argv[1]);
    blocksperfile = 8;
    binheader_count = 10 + 10*blocksperfile*1;
    
    
    binheader = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count));
    ret = read(fd, binheader, (sizeof(*binheader) * binheader_count));
    
    int cnt = 0, xer = 0;
    int countE = 0;
    int countB = 0;
    int countDD = 0;
    int no_elements1 = 0, no_elements1A = 0;
    int no_elements2 = 0, no_elements2A = 0;
    int no_elements3 = 0, no_elements3A = 0;
    int no_elements4 = 0, no_elements4A = 0;
    int no_elements_v1 = 0;
    int no_elements_v2 = 0;
    int no_elements_v3 = 0;
    int no_elements_v4 = 0;
    
    for(j = 0 ; j < blocksperfile; j++)
    {
	length = ntohl(binheader[j*10 + 14]);
	data_offset = ntohl(binheader[j*10 + 12]);
	
	//printf("%d %d\n", length, length1);
	if( length == 32768*8*1 )
	{
		data_buffer = (double*)malloc(length);
		assert(data_buffer != NULL);
		
		
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i< length/(sizeof(double)); i++)
		{
		    if((int)data_buffer[i] == (int)101325)
		      countE++;
		    
		    //printf("VAR : [%d] : %16.16f  \n", i, (double)data_buffer[i]);
 			  
		}
		free(data_buffer);
	}  
    }
    close(fd);
    
    
    //assert(no_elements1 == no_elements2);
    //assert(no_elements2 == no_elements3);
    //assert(no_elements3 == no_elements4);
    //printf("%d : Total Number of elements = [%d %d] [%d %d] [%d %d] [%d %d]\n", xer,  no_elements1, no_elements1A, no_elements2, no_elements2A, no_elements3, no_elements3A, no_elements4, no_elements4A);
    printf("Equal number of elements for variable 1 is %d\n", countE);
    //printf("Equal number of elements for variable 2 is %d\n", no_elements_v2);
    //printf("Equal number of elements for variable 3 is %d\n", no_elements_v3);
    //printf("Equal number of elements for variable 4 is %d\n", no_elements_v4);
    return(0);
} 
