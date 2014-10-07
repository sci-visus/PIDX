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
    int fd;
    int fd1;
    uint32_t* binheader;
    uint32_t* binheader1;
    int binheader_count;
    int target_block = 0;
    int* data_buffer = NULL;
    double* data_buffer1 = NULL;
    size_t data_size, data_size1;
    off_t data_offset, data_offset1;
    int i, j, k, count = 0;
    int number_of_directory;
    size_t length, length1;
  
    fd1 = open(argv[2], O_RDONLY); 
    fd = open(argv[1], O_RDONLY); 
    if(fd < 0)
    {
        perror("open");
        return(-1);
    }
		
    blocksperfile = 256;
    int blocksperfile1 = 256*4;
    int binheader_count1=0;
    binheader_count = 10 + 10*blocksperfile;
    binheader_count1 = 10 + 10*blocksperfile1;
    
    binheader = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count));
    binheader1 = (uint32_t*)malloc(sizeof(*binheader)*(binheader_count1));
    
    ret = read(fd, binheader, (sizeof(*binheader) * binheader_count));
    ret = read(fd1, binheader1, (sizeof(*binheader1) * (binheader_count1)));
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
    
    
    for(j = 0 ; j < 256; j++)
    {
	length = ntohl(binheader[j*10 + 14]);
	data_offset = ntohl(binheader[j*10 + 12]);
	
	length1 = ntohl(binheader1[j*10 + 14]);
	data_offset1 = ntohl(binheader1[j*10 + 12]);
	//printf("%d %d\n", length, length1);
	if( length == 32768*4*1 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (int*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd1, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i<(length/(sizeof(int)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements1++;
		    if(data_buffer1[1*i] !=0 )
			no_elements1A++;
		    
		    if( (data_buffer[1*i] + 100 == (int)data_buffer1[1*i]) )
		    {
		      if(data_buffer[1*i] == 0 && data_buffer1[1*i] == 0)
		      {
		      }
		      else
			no_elements_v1++;
		    }
		    else
		    {
			  if(data_buffer[1*i] == 0 && (int)data_buffer1[1*i] == 0)
			  {
				
			  }
			  else
			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 1 : bad %d : %f %f \n", cnt, data_buffer[1*i], data_buffer1[1*i]);
			  }
		    }
		    
		    cnt++;
		}
		
		free(data_buffer);
		free(data_buffer1);
	}  
    }
    
    
    for(j = 0 ; j < 256; j++)
    {
	length = ntohl(binheader[(j)*10 + 14]);
	data_offset = ntohl(binheader[(j)*10 + 12]);
	
	length1 = ntohl(binheader1[(j+256)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j+256)*10 + 12]);
	//printf("%d %d\n", length, length1);
	
	if( length == 32768*4*1 && length1 == 32768*8*1 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (int*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd1, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		for(i=0; i<(length/(sizeof(int)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements2++;
		    if(data_buffer1[1*i] !=0 )
			no_elements2A++;
		    
		    if( (data_buffer[1*i] + 100 == (int)data_buffer1[1*i]) )
		    {
		      if(data_buffer[1*i] == 0 && (int)data_buffer1[1*i] == 0)
		      {
		      }
		      else
		      {
			no_elements_v2++;
			
		      }
		    }
		    else
		    {
			  if(data_buffer[1*i] == 0 && (int)data_buffer1[1*i] == 0)
			  {
				
			  }
			  else
			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 2 : bad %d : %f %f \n", cnt, data_buffer[1*i], data_buffer1[1*i]);
				
			  }
		    }
		    
		    cnt++;
		}
		
		free(data_buffer);
		free(data_buffer1);
	}  
    }
    
    
    for(j = 0 ; j < 256; j++)
    {
	length = ntohl(binheader[(j)*10 + 14]);
	data_offset = ntohl(binheader[(j)*10 + 12]);
	
	length1 = ntohl(binheader1[(j+512)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j+512)*10 + 12]);
	
	//printf("%d %d\n", length, length1);
	
	if( length == 32768*4 && length1 == 32768*8*3 )
	{
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (int*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);
  
		ret = pread(fd1, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		for(i=0; i<(length/(sizeof(int)*1)); i++)
		{
		    if(data_buffer[i] !=0 )
			no_elements3++;
		    if(data_buffer1[3*i] !=0 )
			no_elements3A++;
		    
		    if( (data_buffer[i] + 100 == data_buffer1[3*i]) && (data_buffer[i] + 101 == data_buffer1[3*i + 1]) && (data_buffer[i] + 102 == data_buffer1[3*i + 2]))
		    {
			if(data_buffer[i] == 0 && data_buffer1[3*i] == 0)
			{
			}
			else
			  no_elements_v3++;
		    }
		    else
		    {
			  if(data_buffer[i] == 0 && data_buffer1[3*i] == 0)
		          {
				
			  }
			  else
			  {
				xer++;
				//printf("VAR 3 : bad %d : %d\n", data_buffer[i], (int)data_buffer1[3*i]);
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
    
    
    for(j = 0 ; j < 256; j++)
    {
	length = ntohl(binheader[(j)*10 + 14]);
	data_offset = ntohl(binheader[(j)*10 + 12]);
	
	length1 = ntohl(binheader1[(j+512+256)*10 + 14]);
	data_offset1 = ntohl(binheader1[(j+512 + 256)*10 + 12]);
	//printf("%d %d\n", length, length1);
	
	if( length == 32768*4*1 && length1 == 32768*8*11 )
	{
		
		data_buffer1 = (double*)malloc(length1);
		data_buffer = (int*)malloc(length);
		assert(data_buffer != NULL);
		assert(data_buffer1 != NULL);

		ret = pread(fd1, data_buffer1, length1, data_offset1);
		ret = pread(fd, data_buffer, length, data_offset);
		
		
		for(i=0; i<(length/(sizeof(int)*1)); i++)
		{
		    if(data_buffer[1*i] !=0 )
			no_elements4++;
		    if(data_buffer1[11*i] !=0 )
			no_elements4A++;
		    
		    if( (data_buffer[1*i] + 100 == (int)data_buffer1[11*i])  && (data_buffer[1*i] + 101 == (int)data_buffer1[11*i + 1]) && (data_buffer[1*i] + 102 == (int)data_buffer1[11*i + 2]) && (data_buffer[1*i] + 103 == (int)data_buffer1[11*i + 3]) && (data_buffer[1*i] + 104 == (int)data_buffer1[11*i + 4]) && (data_buffer[1*i] + 105 == (int)data_buffer1[11*i + 5]) && (data_buffer[1*i] + 106 == (int)data_buffer1[11*i + 6]) && (data_buffer[1*i] + 107 == (int)data_buffer1[11*i + 7]) && (data_buffer[1*i] + 108 == (int)data_buffer1[11*i + 8]) && (data_buffer[1*i] + 109 == (int)data_buffer1[11*i + 9]) && (data_buffer[1*i] + 110 == (int)data_buffer1[11*i + 10]))
		    {
			if(data_buffer[1*i] == 0 && data_buffer1[11*i] == 0)
			{
			}
			else
			  no_elements_v4++;
		    }
		    else
		    {
 			  if(data_buffer[1*i] == 0 && data_buffer1[11*i] == 0)
 			  {
 				
 			  }
 			  else
 			  {
				xer++;
				//if(data_buffer[i] - data_buffer1[i] > 0.0001 || data_buffer[i] - data_buffer1[i] < 0.0001)
				//printf("VAR 4 : bad : %d %d \n", data_buffer[1*i], (int)data_buffer1[11*i]/);
 			  }
		    }
		    
		    cnt++;
		}
		
		free(data_buffer);
		free(data_buffer1);
	}
    }
    
    close(fd);
    close(fd1);
    //assert(no_elements1 == no_elements2);
    //assert(no_elements2 == no_elements3);
    //assert(no_elements3 == no_elements4);
    printf("%d : Total Number of elements = [%d %d] [%d %d] [%d %d] [%d %d]\n", xer,  no_elements1, no_elements1A, no_elements2, no_elements2A, no_elements3, no_elements3A, no_elements4, no_elements4A);
    printf("Equal number of elements for variable 1 is %d\n", no_elements_v1);
    printf("Equal number of elements for variable 2 is %d\n", no_elements_v2);
    printf("Equal number of elements for variable 3 is %d\n", no_elements_v3);
    printf("Equal number of elements for variable 4  is %d\n", no_elements_v4);
    return(0);
} 
