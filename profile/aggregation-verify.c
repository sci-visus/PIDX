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

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

static char dir_name[512] = {0};
static int variable_count = 0;
static int nfiles = 0;

static int parse_args(int argc, char **argv);
static void usage(void);

int main(int argc, char **argv)
{
      int i, j = 0;
      char file_name[512] = {0};
      FILE* fp;
      double *buffer;
      int counter = 0;
      size_t result;
      parse_args(argc, argv);
	      
      for(i = 0 ; i < nfiles ; i++)
      {
	  sprintf(file_name,"%s/file_%d.bin", dir_name, i);
	  printf("File Name %s\n", file_name);
	  fp = fopen(file_name, "rb");
	  if(!fp)
	  {
	      //
	  }
	  buffer = (double*)malloc(sizeof(double) * 512 * 32768 * variable_count);
	  assert(buffer);
	  memset(buffer, 0 , (sizeof(double) * 512 * 32768 * variable_count));
	  
	  result = fread(buffer, sizeof(double), 512 * 32768 * variable_count, fp);
	  assert(result == 512 * 32768 * variable_count);
	  
	  for( j = 0 ; j < 512 * 32768 * variable_count ; j++)
	  {
		if(buffer[j] != 17.85)
		{
		  counter++;
		  printf("j = %d : %f\n", j, buffer[j]);
		}
	  }
	  printf("Number of unequal elements in file %d are %d\n", i, counter);
	  fclose(fp);
      }
}

static int parse_args(int argc, char **argv)
{
    char flags[] = "n:v:f:";
    int one_opt = 0, i =0;
    
    while((one_opt = getopt(argc, argv, flags)) != EOF)
    {
        /* postpone error checking for after while loop */
        switch(one_opt)
        {
	    case('n'):
                sscanf(optarg, "%d", &nfiles);
                break;
	    case('v'):
                sscanf(optarg, "%d", &variable_count);
                break;
            case('f'): 
		sprintf(dir_name, "%s", optarg);
		break;
            case('?'):
                return(-1);
        }
    }
    /* need positive dimensions */
    
    return(0);
}

/* prints usage instructions */
static void usage(void)
{
    printf("Usage: idx-verify -f Filename.idx\n");
    printf("  -f: IDX Filename\n");
    printf("\n");
    return;
} 
