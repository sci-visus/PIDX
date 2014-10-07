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

#define _XOPEN_SOURCE 600
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
 #include <mpi.h>
#include <stdio.h>
#include <math.h>


static int parse_args(int argc, char **argv);
static void usage(void);

static int* extents = 0;
static int* local_extents = 0;
static char** dir_path;
static char output_file_temp[512] = {0};
static int number_of_variables;
static int time_step;
static int MODE = 0;

int main(int argc, char **argv) 
{
  
    int ret = 0, j = 0, i = 0;
    int rank, nprocs;
    uint32_t* headers;
    MPI_Init(&argc, &argv);
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      
    /*Rank 0 parses the command Line Arguments*/
    extents = (int*)malloc(3 * sizeof(int));
    assert(extents);
    memset(extents, 0, 3 * sizeof(int));
    
    local_extents = (int*)malloc(3 * sizeof(int));
    assert(local_extents);
    memset(local_extents, 0, 3 * sizeof(int));
    
    
    if(rank == 0)
    {
	ret = parse_args(argc, argv);
	if(ret < 0)
	{
	    usage();
	    MPI_Abort(MPI_COMM_WORLD, -1);
	}
    }
    
    MPI_Bcast(&extents[0], 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&extents[1], 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&extents[2], 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(&local_extents[0], 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&local_extents[1], 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&local_extents[2], 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&number_of_variables, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&output_file_temp, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    struct stat stat_buf;
    ret = stat(output_file_temp, &stat_buf);
    if(ret != 0)
    {
	  fprintf(stderr, "Error A \n");
	  MPI_Abort(MPI_COMM_WORLD, -1);
    }
    uint32_t data_offset = 0;
    uint32_t header_offset = 0;
    int fs_block_size = stat_buf.st_blksize;
    int total_header_size = (10+(10 * 512 /*block_per_file*/))*sizeof(*headers) * number_of_variables;
    int start_fs_block = total_header_size / fs_block_size;
    header_offset += start_fs_block * fs_block_size;
    
    //printf("Data Offset : %jd [%d : %d : %d]\n", (intmax_t)data_offset, fs_block_size, total_header_size, start_fs_block);
    
    int spv[2] = {1,3};
    double mode0A[2] = {0,0};
    double mode0B[2] = {0,0};
    double mode0C[2] = {0,0};
    double mode0D[2] = {0,0};
    
    double mode1A = 0;
    double mode1B = 0;
    double mode1C = 0;
    double mode1D = 0;
    
    int global_offset[2] = {0, (extents[0] * extents[1] * extents[2] * sizeof(double))};
    printf("XXXX  :  [%d : %d]\n", global_offset[0], global_offset[1]);
    printf("YYYY  :  [%d : %d]\n", spv[0], spv[1]);
    MPI_Status status;
    MPI_File fh;

   
    
    if(MODE == 0) //No aggregation
    {
	double **local_buffer = (double**)malloc(sizeof(double*) * number_of_variables);
	for(i = 0 ; i < number_of_variables ; i++)
	{
	      local_buffer[i] = (double*)malloc(sizeof(double) * local_extents[0] * local_extents[1] * local_extents[2] * spv[i]);
	      assert(local_buffer[i]);
	      memset(local_buffer[i], 0, sizeof(double) * local_extents[0] * local_extents[1] * local_extents[2] * spv[i]);
	}
	
	for(i = 0 ; i < number_of_variables ; i ++)
	{
	      data_offset = header_offset + (rank * local_extents[0] * local_extents[1] * local_extents[2] * spv[i] * sizeof(double)) + global_offset[i];
	      if(rank == 0)
	      printf("A [%d] [%d] : OFFSET %jd COUNT %d\n", rank, i, (intmax_t)data_offset, (local_extents[0] * local_extents[1] * local_extents[2] * spv[i]));
	      
	      mode0A[i] = MPI_Wtime();
	      
	      ret = MPI_File_open(MPI_COMM_SELF, output_file_temp, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	      if(ret != MPI_SUCCESS)
	      {
		    fprintf(stderr, "Error 1 \n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
	      }
	      
	      mode0B[i] = MPI_Wtime();
	      
	      ret = MPI_File_read_at(fh, data_offset, local_buffer[i], (local_extents[0] * local_extents[1] * local_extents[2] * spv[i]), MPI_DOUBLE, &status);
	      if(ret != MPI_SUCCESS)
	      {
		    fprintf(stderr, "Error 2 \n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
	      }
	      if(rank == 0)
	      {
		  int t = 0;
		  for(t = 0 ; t < 20 ; t++)
		    printf("[%d] value at %d = %f\n", i, t, local_buffer[i][t]);
	      }
	      
	      mode0C[i] = MPI_Wtime();
	      
	      MPI_File_close(&fh);
	      
	      mode0D[i] = MPI_Wtime();
	}
	
		
	for(i = 0 ; i < number_of_variables ; i++)
	{
	    //printf("[%d] Time_0 : [%f] [%f %f %f]\n", rank, (mode0D[i] - mode0A[i]), (mode0B[i] - mode0A[i]), (mode0C[i] - mode0B[i]), (mode0D[i] - mode0C[i]));
	}
	
	int m=0, n=0, b=0, block_limit = 0, block_negative_offset = 0;
	long long initial_offset = 0;
	
	if(rank == 0)
	{
	    ret = MPI_File_open(MPI_COMM_SELF, "test.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	    if(ret != MPI_SUCCESS)
	    {
		
	    }
	    
	    
	    ret = MPI_File_close(&fh);
	    if(ret != MPI_SUCCESS)
	    {
		
	    }
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	
	ret = MPI_File_open(MPI_COMM_SELF, "test.bin", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if(ret != MPI_SUCCESS)
	{
	      
	}
		
	for( n = 0 ; n < number_of_variables ; n++ )
	{
	    data_offset = 32768 + (rank * local_extents[0] * local_extents[1] * local_extents[2] * spv[n] * sizeof(double)) + global_offset[n];
	    printf("B [%d] [%d] : OFFSET %jd COUNT %d\n", rank, n, (intmax_t)data_offset, (local_extents[0] * local_extents[1] * local_extents[2] * spv[n]));
	    
	    ret = MPI_File_write_at(fh, data_offset , local_buffer[n], local_extents[0] * local_extents[1] * local_extents[2] * spv[n], MPI_DOUBLE, &status);
	    if(ret != MPI_SUCCESS)
	    {
		return(-1);
	    }
	}
	MPI_File_close(&fh);
	 
	for(i = 0 ; i < number_of_variables ; i++)
	{
	    free(local_buffer[i]);
	    local_buffer[i] = 0;
	}
	free(local_buffer);
	local_buffer = 0;
	
    }
    else //With Aggregation
    {
	  if(rank == 0 || rank == nprocs / 4 || rank == (2 * nprocs) / 4 || rank == (3 * nprocs) / 4)
	  {
	    
		double *agg_buffer = (double*)malloc(sizeof(double) * 512 * 32768);
		
		agg_buffer = (double*)malloc(sizeof(double) * 512 * 32768);
		assert(agg_buffer);
		memset(agg_buffer, 0, sizeof(double) * 512 * 32768);
		if(rank == 0)
		    data_offset = data_offset + 0;
		if(rank ==  nprocs / 4)
		    data_offset = data_offset + 512 * 32768 * sizeof(double);
		if(rank ==  2 * nprocs / 4)
		    data_offset = data_offset + 512 * 32768 * 2 * sizeof(double);
		if(rank ==  3 * nprocs / 4)
		    data_offset = data_offset + 512 * 32768 * 3 * sizeof(double);
		
		printf("[%d] [%d] : OFFSET %jd \n", rank, i, (intmax_t)data_offset);
		
		
		mode1A = MPI_Wtime();
		ret = MPI_File_open(MPI_COMM_SELF, output_file_temp, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		if(ret != MPI_SUCCESS)
		{
		      fprintf(stderr, "Error 3 \n");
		      MPI_Abort(MPI_COMM_WORLD, -1);
		}
		mode1B = MPI_Wtime();
		ret = MPI_File_read_at(fh, data_offset, agg_buffer, (512 * 32768), MPI_DOUBLE, &status);
		if(ret != MPI_SUCCESS)
		{
		    fprintf(stderr, "Error 4 \n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
		}
		mode1C = MPI_Wtime();
		MPI_File_close(&fh);
		mode1D = MPI_Wtime();
		
		printf("[%d] Time_1 : [%f] [%f %f %f]\n", rank, (mode0D - mode0A), (mode0B - mode0A), (mode0C - mode0B), (mode0D - mode0C));
	  }
    }
    
    
    
    MPI_Finalize();
}

static int parse_args(int argc, char **argv)
{
    char flags[] = "g:l:f:v:";
    int one_opt = 0, i =0;
    
    while((one_opt = getopt(argc, argv, flags)) != EOF)
    {
        switch(one_opt)
        {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
                break;
	    case('l'):
                sscanf(optarg, "%dx%dx%d", &local_extents[0], &local_extents[1], &local_extents[2]);
                break;
	    case('f'): 
		sprintf(output_file_temp, "%s", optarg);
		break;
	    case('v'):
                sscanf(optarg, "%d", &number_of_variables);
		break;
		
	    /*
	    case('t'):
                sscanf(optarg, "%d", &time_step);
                break;
	    case('m'):
                sscanf(optarg, "%d", &MODE);
		break;
	    */


	    case('?'):
                return(-1);
        }
    }
    return(0);
}

static void usage(void)
{
    printf("Usage: test-PIDX -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("  -f: IDX Filename\n");
    printf("  -t: number of timesteps\n");
    printf("\n");
    return;
} 
