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
#include "PIDX_rst.h"
#include "mpi.h"


#define PIDX_MAX_DIMENSIONS 5

static int parse_args(int argc, char **argv);
static void usage(void);
static void print_error(char *error_message, char* file, int line);

/* global dimensions of 3D volume */
static int extents[5] = {0, 0, 0, 0, 0};

/* per-process dimensions of each sub-block */
static int count_local[5] = {0, 0, 0, 0, 0};

/* Number of time-steps */
static int time_step = 0;

int main(int argc, char **argv) 
{
    int i = 0, j = 0, k = 0;
    int spv = 0; /*samples per variable*/
    int ts, vc; /*time step counter and variable counter*/
    int nprocs, rank; /*process count and rank*/
    int slice, ret;
    int sub_div[5], offset_local[5];
    
    const int *gextent; /*Global Extensions of the dataset (64 64 64 0 0)*/
    const int *count; /*Local extents of each process*/
    const int *offset; /*Local counts of each process*/

    int number_of_variables = 1;
    int* sample_per_variable_buffer;

    double** write_data;

    int dimension = 3;
    int send_c = 0, send_o = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start_time=0;
    double end_time=0;
    
    int local_rank_count[5] = {1,1,1,1,1};
    int local_rank_offset[5] = {0,0,0,0,0,};
    
    int* expanded_box_dimension;

    int *set_expanded_box_dimension;
    double* expanded_box;
    
    int holding_data = 0;
    PIDX_rst_id rst_id;
    PIDX_Ndim_buffer*** rst_output_buffer;
    PIDX_Ndim_buffer** in_buf;
    int buffer_size;
    int rst_output_buffer_count = 0, var = 0;
    
    FILE* fp;
    
    /*Rank 0 parses the command Line Arguments*/
    if (rank == 0) {
	ret = parse_args(argc, argv);
        if (ret < 0) {
            usage();
            print_error("ret error", __FILE__, __LINE__);
        }
        if (count_local[0] == 0 || count_local[1] == 0 || count_local[2] == 0) {
            usage();
            print_error("Local Dimension cannot be 0!!!!!!!!!\n", __FILE__, __LINE__);
        } else {
            if ((extents[0] / count_local[0]) * (extents[1] / count_local[1]) * (extents[2] / count_local[2]) != nprocs) {
                usage();
                print_error("Wrong Number of Processes\n", __FILE__, __LINE__);
            }
        }
    }
    extents[3] = 1;
    extents[4] = 1;

#ifdef MPI
    /*   The command line arguments are shared by all processes  */
    MPI_Bcast(extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    /* Calculating every process's offset and count */
    gextent = extents;
    sub_div[0] = (extents[0] / count_local[0]);
    sub_div[1] = (extents[1] / count_local[1]);
    sub_div[2] = (extents[2] / count_local[2]);
    offset_local[2] = (rank / (sub_div[0] * sub_div[1])) * count_local[2];
    slice = rank % (sub_div[0] * sub_div[1]);
    offset_local[1] = (slice / sub_div[0]) * count_local[1];
    offset_local[0] = (slice % sub_div[0]) * count_local[0];

    offset_local[3] = 0;
    offset_local[4] = 0;
    count_local[3] = 1;
    count_local[4] = 1;

    offset = offset_local;
    count = count_local;
    assert(offset[0] < extents[0]);
    assert(offset[1] < extents[1]);
    assert(offset[2] < extents[2]);
    assert(offset[3] < extents[3]);
    assert(offset[4] < extents[4]);

    sample_per_variable_buffer = (int*) malloc(sizeof (int) * number_of_variables);
    if (!sample_per_variable_buffer)
        print_error("Error Allocating {sample_per_variable_buffer} buffer", __FILE__, __LINE__);

    sample_per_variable_buffer[0] = 1;
    //sample_per_variable_buffer[1] = 4;
    //sample_per_variable_buffer[2] = 4;
    //sample_per_variable_buffer[3] = 4;
    double *rst_init_start, *rst_init_end, **rst_start, **rst_end, **internal_rst_start, **internal_rst_end, *rst_cleanup_start, *rst_cleanup_end;

    rst_init_start = (double*) malloc(sizeof (double) * time_step);
    rst_init_end = (double*) malloc(sizeof (double) * time_step);
    rst_cleanup_start = (double*) malloc(sizeof (double) * time_step);
    rst_cleanup_end = (double*) malloc(sizeof (double) * time_step);
    memset(rst_init_start, 0, (sizeof (double) * time_step));
    memset(rst_init_end, 0, (sizeof (double) * time_step));
    memset(rst_cleanup_start, 0, (sizeof (double) * time_step));
    memset(rst_cleanup_end, 0, (sizeof (double) * time_step));

    
    rst_start = (double**) malloc(sizeof (double*) * time_step);
    rst_end = (double**) malloc(sizeof (double*) * time_step);
    internal_rst_start = (double**) malloc(sizeof (double*) * time_step);
    internal_rst_end = (double**) malloc(sizeof (double*) * time_step);
    memset(rst_start, 0, (sizeof (double*) * time_step));
    memset(rst_end, 0, (sizeof (double*) * time_step));
    memset(internal_rst_start, 0, (sizeof (double*) * time_step));
    memset(internal_rst_end, 0, (sizeof (double*) * time_step));

    for (i = 0; i < time_step; i++) 
    {
        rst_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
        rst_end[i] = (double*) malloc(sizeof (double) * number_of_variables);
        internal_rst_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
        internal_rst_end[i] = (double*) malloc(sizeof (double) * number_of_variables);

        memset(rst_start[i], 0, (sizeof (double) * number_of_variables));
        memset(rst_end[i], 0, (sizeof (double) * number_of_variables));
        memset(internal_rst_start[i], 0, (sizeof (double) * number_of_variables));
        memset(internal_rst_end[i], 0, (sizeof (double) * number_of_variables));
    }
    
    expanded_box_dimension = (int*)malloc(sizeof(int) * 5);
    set_expanded_box_dimension = (int*)malloc(sizeof(int) * 5);
    
    if(rank == 0)
    {
	FILE* fp;
	fp = fopen("input", "r");
	fscanf(fp, "%d %d %d", set_expanded_box_dimension, set_expanded_box_dimension+1, set_expanded_box_dimension+2);
	fclose(fp);
	printf("READ:  %d %d %d\n", set_expanded_box_dimension[0], set_expanded_box_dimension[1], set_expanded_box_dimension[2]);
    }
    set_expanded_box_dimension[3] = 1;
    set_expanded_box_dimension[4] = 1;
    
    MPI_Bcast(set_expanded_box_dimension, 5, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank == 0)
    	fp = fopen("output", "w");
	
    for (ts = 0; ts < time_step; ts++) 
    {
	//WRITE BUFFER
	write_data = (double**) malloc(sizeof (double*) * number_of_variables);
	if (!write_data)
	    print_error("Error Allocating Buffer", __FILE__, __LINE__);
	memset(write_data, 0, sizeof (double*) * number_of_variables);

	for (vc = 0; vc < number_of_variables; vc++) 
	{
	    write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * count[3] * count[4] * sample_per_variable_buffer[vc]);
	    if (!write_data[vc])
		print_error("Error Allocating Buffer", __FILE__, __LINE__);

	    for (k = 0; k < count[2]; k++)
		for (j = 0; j < count[1]; j++)
		    for (i = 0; i < count[0]; i++) 
		    {
			long long index = (long long) (count[0] * count[1] * k) + (count[0] * j) + i;
			for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
			    write_data[vc][index * sample_per_variable_buffer[vc] + spv] = 100 + (ts) +spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i);
		    }
	}
	
	start_time = MPI_Wtime();
	rst_init_start[ts] = MPI_Wtime();
	in_buf = (PIDX_Ndim_buffer**) malloc(number_of_variables * sizeof (PIDX_Ndim_buffer*));
	for (i = 0; i < number_of_variables; i++)
	    in_buf[i] = (PIDX_Ndim_buffer*) malloc(sizeof (PIDX_Ndim_buffer));


	rst_output_buffer = (PIDX_Ndim_buffer***) malloc(sizeof (*rst_output_buffer) * number_of_variables);
	memset(rst_output_buffer, 0, sizeof (*rst_output_buffer) * number_of_variables);
	rst_output_buffer_count = 0;


	rst_id = PIDX_rst_init(MPI_COMM_WORLD, dimension, (int*) gextent, count, offset, 1, set_expanded_box_dimension, &rst_output_buffer_count);
	rst_init_end[ts] = MPI_Wtime();

	

	for (var = 0; var < number_of_variables; var++) 
	{
	    rst_start[ts][var] = MPI_Wtime();
	    PIDX_rst_buf_init(in_buf[var], 3, offset, count, write_data[var], sample_per_variable_buffer[var], MPI_DOUBLE, "var_name", var);
	    rst_output_buffer[var] = (PIDX_Ndim_buffer**) malloc((rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
	    memset(rst_output_buffer[var], 0, (rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
	    PIDX_rst_restructure(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count));
	    PIDX_rst_restructure_IO(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count), 1);
	    //HELPER_rst(rst_output_buffer[var], rst_id, rst_output_buffer_count, sample_per_variable_buffer[var]);
	    rst_end[ts][var] = MPI_Wtime();

	    internal_rst_start[ts][var] = MPI_Wtime();
	    if (rst_output_buffer_count != 0) 
	    {
		buffer_size = 0;
		expanded_box_dimension = PIDX_rst_get_box_dimension(rst_id);
		if (dimension == 1)
		    buffer_size = expanded_box_dimension[0];
		else if (dimension == 2)
		    buffer_size = expanded_box_dimension[0] * expanded_box_dimension[1];
		else if (dimension == 3)
		    buffer_size = expanded_box_dimension[0] * expanded_box_dimension[1] * expanded_box_dimension[2];

		expanded_box = (double*) malloc(sizeof (double) * buffer_size * sample_per_variable_buffer[var]);
		memset(expanded_box, 0, (sizeof (double) * buffer_size * sample_per_variable_buffer[var]));


		int k1, j1, i1, r, index = 0, recv_o = 0;
		for (r = 0; r < rst_output_buffer_count; r++) 
		{
		    for (k1 = rst_output_buffer[var][r]->lower_bounds[2]; k1 < rst_output_buffer[var][r]->upper_bounds[2]; k1++)
			for (j1 = rst_output_buffer[var][r]->lower_bounds[1]; j1 < rst_output_buffer[var][r]->upper_bounds[1]; j1++)
			    for (i1 = rst_output_buffer[var][r]->lower_bounds[0]; i1 < rst_output_buffer[var][r]->upper_bounds[0]; i1 = i1 + rst_output_buffer[var][r]->upper_bounds[0] - rst_output_buffer[var][r]->lower_bounds[0]) {

				index = ((rst_output_buffer[var][r]->upper_bounds[0] - rst_output_buffer[var][r]->lower_bounds[0])* (rst_output_buffer[var][r]->upper_bounds[1] - rst_output_buffer[var][r]->lower_bounds[1]) * (k1 - rst_output_buffer[var][r]->lower_bounds[2])) +
					((rst_output_buffer[var][r]->upper_bounds[0] - rst_output_buffer[var][r]->lower_bounds[0]) * (j1 - rst_output_buffer[var][r]->lower_bounds[1])) +
					(i1 - rst_output_buffer[var][r]->lower_bounds[0]);

				send_o = index * sample_per_variable_buffer[var];
				send_c = (rst_output_buffer[var][r]->upper_bounds[0] - rst_output_buffer[var][r]->lower_bounds[0]) * sample_per_variable_buffer[var];


				recv_o = ((rst_output_buffer[var][r]->regular_upper_bounds[0] - rst_output_buffer[var][r]->regular_lower_bounds[0]) * (rst_output_buffer[var][r]->regular_upper_bounds[1] - rst_output_buffer[var][r]->regular_lower_bounds[1]) * (k1 - rst_output_buffer[var][r]->regular_lower_bounds[2])) +
					((rst_output_buffer[var][r]->regular_upper_bounds[0] - rst_output_buffer[var][r]->regular_lower_bounds[0])* (j1 - rst_output_buffer[var][r]->regular_lower_bounds[1])) +
					(i1 - rst_output_buffer[var][r]->regular_lower_bounds[0]);

				memcpy(expanded_box + (recv_o * sample_per_variable_buffer[var]), rst_output_buffer[var][r]->buffer + send_o, send_c * sizeof (double));
			    }
		}
	    }
	    internal_rst_end[ts][var] = MPI_Wtime();
	    
	    
	    if (rst_output_buffer_count != 0) {
		holding_data = 1;
		local_rank_count[0] = rst_output_buffer[var][0]->regular_upper_bounds[0] - rst_output_buffer[var][0]->regular_lower_bounds[0];
		local_rank_count[1] = rst_output_buffer[var][0]->regular_upper_bounds[1] - rst_output_buffer[var][0]->regular_lower_bounds[1];
		local_rank_count[2] = rst_output_buffer[var][0]->regular_upper_bounds[2] - rst_output_buffer[var][0]->regular_lower_bounds[2];
		local_rank_count[3] = rst_output_buffer[var][0]->regular_upper_bounds[3] - rst_output_buffer[var][0]->regular_lower_bounds[3];
		local_rank_count[4] = rst_output_buffer[var][0]->regular_upper_bounds[4] - rst_output_buffer[var][0]->regular_lower_bounds[4];
		
		local_rank_offset[0] = rst_output_buffer[var][0]->regular_lower_bounds[0];
		local_rank_offset[1] = rst_output_buffer[var][0]->regular_lower_bounds[1];
		local_rank_offset[2] = rst_output_buffer[var][0]->regular_lower_bounds[2];
		local_rank_offset[3] = rst_output_buffer[var][0]->regular_lower_bounds[3];
		local_rank_offset[4] = rst_output_buffer[var][0]->regular_lower_bounds[4];
	    }
	
	    if (rst_output_buffer_count != 0) 
	    {
		free(expanded_box);
		expanded_box = 0;
	    }

	    rst_cleanup_start[ts] = MPI_Wtime();
	    ret = PIDX_rst_buf_destroy(in_buf[var]);
	    for (j = 0; j < rst_output_buffer_count; j++) {
		ret = PIDX_rst_buf_destroy(rst_output_buffer[var][j]);
		if (ret == -1) {
		    fprintf(stderr, "PIDX : [%d] Error in PIDX_rst_buf_destroy\n", rank);
		    MPI_Abort(MPI_COMM_WORLD, -1);
		}
	    }
	    if (rst_output_buffer_count != 0) {
		free(rst_output_buffer[var]);
		rst_output_buffer[var] = 0;
	    }
	    
	    int *rank_status;
	    int *rank_offset;
	    int *rank_count;
	    
	    rank_status = (int*) malloc(sizeof (int) * nprocs);
	    if (!rank_status) PIDX_rst_print_error("Memory : rank_status", __FILE__, __LINE__);
	    memset(rank_status, 0, (sizeof (int) * nprocs));
	    
	    rank_offset = (int*) malloc(sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS);
	    if (!rank_offset) PIDX_rst_print_error("Memory : rank_offset", __FILE__, __LINE__);
	    memset(rank_offset, 0, (sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS));
	    
	    rank_count = (int*) malloc(sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS);
	    if (!rank_count) PIDX_rst_print_error("Memory : rank_count", __FILE__, __LINE__);
	    memset(rank_count, 0, (sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS));
	    
	    ret = MPI_Gather( &holding_data, 1, MPI_INT, rank_status, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Allgather : rank_status", __FILE__, __LINE__);

 	    ret = MPI_Gather( local_rank_offset, PIDX_MAX_DIMENSIONS, MPI_INT, rank_offset, PIDX_MAX_DIMENSIONS, MPI_INT, 0, MPI_COMM_WORLD);
 	    if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Allgather : rank_offset", __FILE__, __LINE__);
	    
 	    ret = MPI_Gather( local_rank_count, PIDX_MAX_DIMENSIONS, MPI_INT, rank_count, PIDX_MAX_DIMENSIONS, MPI_INT, 0, MPI_COMM_WORLD);
 	    if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Allgather : rank_count", __FILE__, __LINE__);
	    
  	    int i1;
 	    if(rank == 0)
	    {
#if 0
 		for(i1 = 0 ; i1 < nprocs ; i1++)
 		    printf("Rank %d : %d :: %d %d %d :: %d %d %d\n", i1, rank_status[i1], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 0], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 1], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 2],
		      rank_count[PIDX_MAX_DIMENSIONS * i1 + 0], rank_count[PIDX_MAX_DIMENSIONS * i1 + 1], rank_count[PIDX_MAX_DIMENSIONS * i1 + 2]);
#endif		
		for(i1 = 0 ; i1 < nprocs ; i1++)
 		    fprintf(fp, "%d %d %d %d %d %d %d %d %d %d\n",ts, var, i1, rank_status[i1], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 0], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 1], rank_offset[PIDX_MAX_DIMENSIONS * i1 + 2],
		      rank_count[PIDX_MAX_DIMENSIONS * i1 + 0], rank_count[PIDX_MAX_DIMENSIONS * i1 + 1], rank_count[PIDX_MAX_DIMENSIONS * i1 + 2]);
		
		
	    }
	    
	    rst_cleanup_end[ts] = MPI_Wtime();
	}
	free(rst_output_buffer);
	rst_output_buffer = 0;

	end_time = MPI_Wtime();


	for (vc = 0; vc < number_of_variables; vc++) {
	    free(write_data[vc]);
	    write_data[vc] = 0;
	}
	free(write_data);
	write_data = 0;
	
	double max_time = 0, total_time = 0;
	total_time = end_time - start_time;
	
	long long total_data = 0;
	int sample_sum = 0;
	for (var = 0; var < number_of_variables; var++) {
	    sample_sum = sample_sum + sample_per_variable_buffer[var];
	}
	total_data = (gextent[0] * gextent[1] * gextent[2] * sample_sum * sizeof (double));
	MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if (total_time == max_time) {
	    printf("[Rank %d Nproc %d] [Time Step %d] [Total Time %f] [Total Volume %d x %d x %d x %d] [Throughput %f MiB/sec]\n", rank, nprocs, ts, max_time, gextent[0], gextent[1], gextent[2], sample_sum, (float) total_data / (1024 * 1024 * max_time));
	    printf("Restructuring Init Time: %f\n", (rst_init_end[ts] - rst_init_start[ts]));
	    for (var = 0; var < number_of_variables; var++) 
	    {
		  printf("[%d] Retructuring Time :: Internal Restructuring Time [%f %f]\n", var, (rst_end[ts][var] - rst_start[ts][var]), (internal_rst_end[ts][var] - internal_rst_end[ts][var]));
		
	    }
	    printf("Cleanup Time: %f\n", (rst_cleanup_end[ts] - rst_cleanup_start[ts]));
	    
	    printf("\n");
	}
    }
    if(rank == 0)
	fclose(fp);
	
    
    
    free(sample_per_variable_buffer);
    sample_per_variable_buffer = 0;

    MPI_Finalize();
    return 0;
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:l:t:";
    int one_opt = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        /* postpone error checking for after while loop */
        switch (one_opt) {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
                break;
            case('l'):
                sscanf(optarg, "%dx%dx%d", &count_local[0], &count_local[1], &count_local[2]);
                break;
            case('t'):
                sscanf(optarg, "%d", &time_step);
                break;
            case('?'):
                return (-1);
        }
    }
    /* need positive dimensions */
    if (extents[0] < 1 || extents[1] < 1 || extents[2] < 1 || count_local[0] < 1 || count_local[1] < 1 || count_local[2] < 1) {
        printf("Error: bad dimension specification.\n");
        return (-1);
    }

    /* need global dimension to be larger than the local */
    if (extents[0] < count_local[0] || extents[1] < count_local[1] || extents[2] < count_local[2]) {
        printf("Error: global dimensions and local dimensions aren't evenly divisible\n");
        return (-1);
    }
    return (0);
}

/* prints usage instructions */
static void usage(void) {
    printf("Usage: rst-phase -g 4x4x4 -l 2x2x2 -t 1\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("  -t: number of timesteps\n");
    printf("\n");
    return;
}

static void print_error(char *error_message, char* file, int line) {
    fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(0);
#endif
}
