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

static int nx = 0, ny = 0, nz = 0;

static int parse_args(int argc, char **argv);
static void usage(void);

int main(int argc, char **argv)
{
      int i = 0, j = 0, k = 0;
      parse_args(argc, argv);      
      
      FILE* fp;
	
       fp = fopen("MPICH_RANK_ORDER", "w");
       if (!fp) {
 	  exit(0);	  
       }
      for(k = 0 ; k < nz ; k = k + 2)
      {
	  for(j = 0 ; j < ny ; j = j + 2)
	  {
	      for(i = 0 ; i < nx ; i = i + 2)
	      {
		  fprintf(fp, "%d,", (nx *ny * k) + (nx *  j) + i);
		  fprintf(fp, "%d,", (nx *ny * k) + (nx *  j) + i+1);
		  fprintf(fp, "%d,", (nx *ny * k) + (nx *  (j+1)) + i);
		  fprintf(fp, "%d,", (nx *ny * k) + (nx *  (j+1)) + i+1);
		  
		  fprintf(fp, "%d,", (nx *ny * (k+1)) + (nx *  j) + i);
		  fprintf(fp, "%d,", (nx *ny * (k+1)) + (nx *  j) + i+1);
		  fprintf(fp, "%d,", (nx *ny * (k+1)) + (nx *  (j+1)) + i);
		  fprintf(fp, "%d,", (nx *ny * (k+1)) + (nx *  (j+1)) + i+1);
	      }
	  }
      }
      fclose(fp);
}

static int parse_args(int argc, char **argv)
{
    char flags[] = "x:y:z:";
    int one_opt = 0, i =0;
    
    while((one_opt = getopt(argc, argv, flags)) != EOF)
    {
        /* postpone error checking for after while loop */
        switch(one_opt)
        {
	    case('x'):
                sscanf(optarg, "%d", &nx);
                break;
	    case('y'):
                sscanf(optarg, "%d", &ny);
                break;
            case('z'): 
		sscanf(optarg, "%d", &nz);
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
