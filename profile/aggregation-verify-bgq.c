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
#include "mpi.h"

int main(int argc, char **argv)
{
    MPI_Win win;
    int i, counter = 0, ts;
    int nprocs, rank;
    int retval=0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
    if(nprocs % 16 != 0)
    {
      fprintf(stdout, "Program Works only when nprocs %% 16 = 0 && nprocs = 16 * 2^i (where i = 0,1,2,3 ...)\ni.e. nprocs = 16, 32, 64, 128, 256 and so on\n");
      MPI_Abort(MPI_COMM_WORLD, retval);
    }
    
    
    double* agg_buffer;
    double* local_buffer;
    
    int time_step_count = 10;	//Number of Time steps
    int local_resolution = 32*32*32;
    int variable_count = 16;
    int proc_buffer_size = local_resolution * variable_count;
    int agg_count = 16;
    int agg_buffer_size = (proc_buffer_size * nprocs)/agg_count;
    
    MPI_Info info;
        
    for(ts = 0 ; ts < time_step_count ; ts++) {
      	counter = 0;
	
	MPI_Info_create(&info);
	MPI_Info_set(info, "no_locks", "1");
	
	local_buffer = (double*) malloc(sizeof(double) * (proc_buffer_size));	//Per Process Buffer
	for(i = 0 ; i  < (proc_buffer_size) ; i++)
	    local_buffer[i] = 88;	//randomly choose 88 (to verify in the end)
	
	if(rank % (nprocs / agg_count) == 0) {		//Chosing 16 aggregators uniformly placed across all ranks. nproc has to be a multiple of 16
	    agg_buffer = (double*) malloc(sizeof(double) * agg_buffer_size);	//aggregator buffer
	    
	    retval = MPI_Win_create(agg_buffer, agg_buffer_size, sizeof(double), info, MPI_COMM_WORLD, &(win));	//if aggregator then create window of agg_buffer_size
	    if (MPI_SUCCESS != retval) {
		fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, retval);
	    }
	}
	else {
	    retval = MPI_Win_create(0, 0, 1, info, MPI_COMM_WORLD, &(win));	//if not aggregator then create window of size 0
	    if (MPI_SUCCESS != retval) {
		fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, retval);
	    }
	}
	    
	retval = MPI_Win_fence(0, win);		//First Fence
	if (MPI_SUCCESS != retval) {
	    fprintf(stderr, "Error Creating MPI_Win_fence\n");
	    MPI_Abort(MPI_COMM_WORLD, retval);
	}
	
	/*
	 * Each process sends data to all 16 aggregators.
	 * Since each process (32x32x32) x 16 data samples, each aggregator gets (32x32x32)
	 * Following PIDX mode of communication, the 32x32x32 block of data that every process sends to an aggregator is actually sent in smaller chunks.
	 * The size of chunk uncreases with powers of two, for instance 32x32x32 (32768), is actually sent as 16 data packets of sizes:
	 * 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 1 
	 */
	int send_count = 1;
	for(i = 0 ; i < agg_count ; i++) {
	    send_count = 1;
	    
	    if(rank == 0 && ts == 0)
	        printf("Rank: [%d] :: Target Rank %d Target Offset %d Count %d\n", rank, (i)*(nprocs/agg_count), (proc_buffer_size/agg_count)*(rank)+(send_count-1), (send_count));
		    
	    //First data packet of size 1
	    assert(((i*nprocs)/agg_count) % (nprocs / agg_count) == 0);
	    assert(((i*(proc_buffer_size/agg_count) + (send_count - 1)) + (send_count)) <= proc_buffer_size);
	    assert(((proc_buffer_size/agg_count)*(rank) + (send_count - 1)) <= agg_buffer_size);
	    
	    retval = MPI_Put(local_buffer+(i*(proc_buffer_size/agg_count) + (send_count - 1)), (send_count), MPI_DOUBLE, (i)*(nprocs/agg_count), (proc_buffer_size/agg_count)*(rank) + (send_count-1), (send_count), MPI_DOUBLE, win);
	    if (MPI_SUCCESS != retval) {
		fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] MPI_Put\n", rank, __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, -1);
	    }
	    
	    while(send_count != (proc_buffer_size/agg_count)) {
	      //remaining data packets of size 1 to 32x32x16
	      assert(((i*nprocs)/agg_count) % (nprocs / agg_count) == 0);
	      assert(((i*(proc_buffer_size/agg_count) + (send_count)) + send_count) <= proc_buffer_size);
	      assert(((proc_buffer_size/agg_count)*(rank) + (send_count) + send_count) <= agg_buffer_size);
	      
	      retval = MPI_Put(local_buffer+(i*(proc_buffer_size/agg_count) + (send_count)), (send_count), MPI_DOUBLE, (i)*(nprocs/agg_count), (proc_buffer_size/agg_count)*(rank) + (send_count), (send_count), MPI_DOUBLE, win);
	      if (MPI_SUCCESS != retval) {
		  fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] MPI_Put\n", rank, __FILE__, __LINE__);
		  MPI_Abort(MPI_COMM_WORLD, -1);
	      }
	      if(rank == 0 && ts == 0)
		printf("Rank: [%d] :: Target Rank %d Target Offset %d Count %d\n", rank, (i)*(nprocs/agg_count), (proc_buffer_size/agg_count)*(rank)+(send_count), (send_count));
	      send_count = send_count * 2;
	    }
	    
        }
	// printf("\n");  //if this line is uncommented then wierdly the code works fine for 256 processes    
	retval = MPI_Win_fence(0, win);	//second fence
	if (MPI_SUCCESS != retval) {
	    fprintf(stderr, "Error With MPI_Fence\n");
	    MPI_Abort(MPI_COMM_WORLD, retval);
	}

	retval = MPI_Win_free(&(win));	//freeing window
	if (MPI_SUCCESS != retval) {
	    fprintf(stderr, "Error with freeing window\n");
	    MPI_Abort(MPI_COMM_WORLD, retval);
	}

	//verification code
	int total_volume = 0;
	if(rank % (nprocs/agg_count) == 0) {	//if aggregator check count all elements that have 88
	    for(i = 0 ; i  < agg_buffer_size ; i++)
	      if(agg_buffer[i] == (double)88)
		counter++;
	    /*  
	    if(counter != agg_buffer_size)
		printf("TS %d : [%d] RMA Failure: %d (%f %f %f %f)\n", ts, rank, counter, agg_buffer[0], agg_buffer[1], agg_buffer[2], agg_buffer[3]);
	    else
		printf("TS %d : [%d] RMA Pass: %d\n", ts, rank, counter);
	    */
	    free(agg_buffer);
	    agg_buffer = 0;
	}
	MPI_Allreduce(&counter, &total_volume, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	//count samples with value 88 across all aggregators
	if(total_volume != (proc_buffer_size)*nprocs) {
	    fprintf(stdout, "Aggregation Failed (count recorded %d)\n", total_volume);
	    MPI_Abort(MPI_COMM_WORLD, -1);
	}
	else {
	    if(rank == 0)
		fprintf(stdout, "Time Step [%d] Aggregation Passed (recorded count %d) (actual count %d)\n", ts, total_volume, (32*32*32*16*nprocs));
	}
	   
	free(local_buffer);
	local_buffer = 0;
	MPI_Info_free(&info);
    }
        
    MPI_Finalize();
    return (0);
}
