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
    double* data_buffer1 = NULL;
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
    printf("File names : %s %s %s %s : %s\n", argv[2], argv[3], argv[4], argv[5], argv[1]);
    blocksperfile = 32;
    int blocksperfile1 = 32;
    int binheader_count1=0;
    binheader_count = 10 + 10*blocksperfile*1;
    binheader_count1 = 10 + 10*blocksperfile1*1;
    
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
    
    fd2 = open(argv[2], O_RDONLY); 
    binheader1 = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count1));
    ret = read(fd2, binheader1, (sizeof(*binheader1) * (binheader_count1)));
    for(j = 0 ; j < blocksperfile; j++)
    {
	length = ntohl(binheader[j*10 + 14]);
	data_offset = ntohl(binheader[j*10 + 12]);
	
	length1 = ntohl(binheader1[j*10 + 14]);
	data_offset1 = ntohl(binheader1[j*10 + 12]);
 	//printf("%d %d\n", length, length1);
	if( length == 32768*8*1 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (double*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd2, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i<(length/(sizeof(double)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements1++;
		    if(data_buffer1[1*i] !=0 )
			no_elements1A++;
		    
		    if( (data_buffer[1*i] == (double)data_buffer1[1*i]))
		    {
 		      if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
 		      {
 		      }
 		      else
		      {
			//printf("%f :: %f\n", data_buffer[0], data_buffer1[0]);
			no_elements_v1++;
		      }
		    }
		    else
		    {
 			  if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
 			  {
 			  }
 			  else
 			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				printf("VAR 1 : bad : %f %f \n", (double)data_buffer[1*i], data_buffer1[1*i]);
 			  }
		    }  
		    cnt++;
		}
		free(data_buffer);
		free(data_buffer1);
	}  
    }
    //free(binheader1);
    //binheader1 = 0;
    //close(fd2);
    
    /*
    //fd3 = open(argv[3], O_RDONLY); 
    //binheader1 = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count1));
    //ret = read(fd3, binheader1, (sizeof(*binheader1) * (binheader_count1)));
    for(j = 0 ; j < blocksperfile; j++)
    {
	length = ntohl(binheader[(j + blocksperfile)*10 + 14]);
	data_offset = ntohl(binheader[(j + blocksperfile)*10 + 12]);
	
	length1 = ntohl(binheader1[(j + blocksperfile)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j + blocksperfile)*10 + 12]);
	printf("%d %d\n", data_offset, data_offset1);
	
	if( length == 32768*8*1 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (double*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd2, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i<(length/(sizeof(double)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements2++;
		    if(data_buffer1[1*i] !=0 )
			no_elements2A++;
		    //printf("%f :: %f\n", data_buffer[0], data_buffer1[0]);
		    if( (data_buffer[1*i] == (double)data_buffer1[1*i]))
		    {
		      if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
 		      {
 		      }
 		      else
			no_elements_v2++;
		    }
		    else
		    {
 			  if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
 			  {
 			  }
 			  else
 			  {
				xer++;
				printf("VAR 2 : bad : %f %f \n", (double)data_buffer[1*i], data_buffer1[1*i]);
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 2 : bad : %d %d \n", (int)data_buffer[1*i], (int)data_buffer1[1*i]);
 			  }
		    }
		    cnt++;
		}
		
		free(data_buffer);
		free(data_buffer1);
	}
    }
    //free(binheader1);
    //binheader1 = 0;
    //close(fd3);
    
    //fd4 = open(argv[4], O_RDONLY); 
    //binheader1 = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count1));
    //ret = read(fd4, binheader1, (sizeof(*binheader1) * (binheader_count1)));
    for(j = 0 ; j < blocksperfile; j++)
    {
	length = ntohl(binheader[(j + blocksperfile*2)*10 + 14]);
	data_offset = ntohl(binheader[(j + blocksperfile*2)*10 + 12]);
	
	length1 = ntohl(binheader1[(j + blocksperfile*2)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j + blocksperfile*2)*10 + 12]);
	//printf("%d %d\n", length, length1);
	
	if( length == 32768*8 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (double*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);
  
		ret = pread(fd2, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		for(i=0; i<(length/(sizeof(double)*1)); i++)
		{
		    if(data_buffer[i] !=0 )
			no_elements3++;
		    if(data_buffer1[1*i] !=0 )
			no_elements3A++;
		    
		    if( (data_buffer[i] == data_buffer1[1*i]) )
		    {
			if(data_buffer[i] == 0 && data_buffer1[1*i] == 0)
			{
			}
			else
			  no_elements_v3++;
		    }
		    else
		    {
			  if(data_buffer[i] == 0 && data_buffer1[1*i] == 0)
		          {
				
			  }
			  else
			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 3 : bad %d : %f\n", cnt, data_buffer[i], data_buffer1[3*i]);
			  }
		    }
		    
		    cnt++;
		}
		free(data_buffer);
		free(data_buffer1);
	}  
    }
    //free(binheader1);
    //binheader1 = 0;
    //close(fd4);
    
    //fd5 = open(argv[5], O_RDONLY); 
    //binheader1 = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count1));
    //ret = read(fd5, binheader1, (sizeof(*binheader1) * (binheader_count1)));
    for(j = 0 ; j < blocksperfile; j++)
    {
	length = ntohl(binheader[(j + blocksperfile*3)*10 + 14]);
	data_offset = ntohl(binheader[(j + blocksperfile*3)*10 + 12]);
	
	length1 = ntohl(binheader1[(j + blocksperfile*3)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j + blocksperfile*3)*10 + 12]);
	//printf("%d %d\n", length, length1);
	
	if( length == 32768*8*1 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (double*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd2, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i<(length/(sizeof(double)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements4++;
		    if(data_buffer1[1*i] !=0 )
			no_elements4A++;
		    
		    if( (data_buffer[1*i] == (double)data_buffer1[1*i]))
		    {
			if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
			{
			}
			else
			  no_elements_v4++;
		    }
		    else
		    {
 			  if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
 			  {
 				
 			  }
 			  else
 			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 4 : bad : %d %d \n", data_buffer[1*i], (int)data_buffer1[11*i]/);
				//printf("VAR 4 : bad : %d %d \n", (int)data_buffer[1*i], (int)data_buffer1[11*i]);
 			  }
		    }
		    
		    cnt++;
		}
		
		free(data_buffer);
		free(data_buffer1);
	}
    }
    free(binheader1);
    binheader1 = 0;
    close(fd2);
    */
    
    //assert(no_elements1 == no_elements2);
    //assert(no_elements2 == no_elements3);
    //assert(no_elements3 == no_elements4);
    printf("%d : Total Number of elements = [%d %d] [%d %d] [%d %d] [%d %d]\n", xer,  no_elements1, no_elements1A, no_elements2, no_elements2A, no_elements3, no_elements3A, no_elements4, no_elements4A);
    printf("Equal number of elements for variable 1 is %d\n", no_elements_v1);
    printf("Equal number of elements for variable 2 is %d\n", no_elements_v2);
    printf("Equal number of elements for variable 3 is %d\n", no_elements_v3);
    printf("Equal number of elements for variable 4 is %d\n", no_elements_v4);
    return(0);
} 
