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
#include "PIDX.h"
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

/* output IDX file Name Template*/
static char output_file_template[512] = {0};
static char **output_file_name;

int main(int argc, char **argv) 
{
    int restructuring_ACTIVE = 0;
    int write_ACTIVE = 1;
    
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

    MPI_File fh;
    char filename[512];
    MPI_Status status;
    int dimension = 3;
    int send_c = 0, send_o = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*Rank 0 parses the command Line Arguments*/
    if (rank == 0) 
    {
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

    MPI_Bcast(extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

    /*Creating the filename*/
    output_file_name = (char**) malloc(sizeof (char*) * time_step);
    for (i = 0; i < time_step; i++) {
        output_file_name[i] = (char*) malloc(sizeof (char) * 512);
        sprintf(output_file_name[i], "%s%04d", output_file_template, i);
    }

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
    
    double *file_create_start, *file_create_end, *file_open_start, *file_open_end, **var_start, **var_end, *file_close_start, *file_close_end;
    double *rst_init_start, *rst_init_end, **rst_start, **rst_end, **internal_rst_start, **internal_rst_end, *rst_cleanup_start, *rst_cleanup_end;

    rst_init_start = (double*) malloc(sizeof (double) * time_step);
    rst_init_end = (double*) malloc(sizeof (double) * time_step);
    rst_cleanup_start = (double*) malloc(sizeof (double) * time_step);
    rst_cleanup_end = (double*) malloc(sizeof (double) * time_step);
    memset(rst_init_start, 0, (sizeof (double) * time_step));
    memset(rst_init_end, 0, (sizeof (double) * time_step));
    memset(rst_cleanup_start, 0, (sizeof (double) * time_step));
    memset(rst_cleanup_end, 0, (sizeof (double) * time_step));

    file_create_start = (double*) malloc(sizeof (double) * time_step);
    file_create_end = (double*) malloc(sizeof (double) * time_step);
    file_open_start = (double*) malloc(sizeof (double) * time_step);
    file_open_end = (double*) malloc(sizeof (double) * time_step);
    file_close_start = (double*) malloc(sizeof (double) * time_step);
    file_close_end = (double*) malloc(sizeof (double) * time_step);

    memset(file_create_start, 0, (sizeof (double) * time_step));
    memset(file_create_end, 0, (sizeof (double) * time_step));
    memset(file_open_start, 0, (sizeof (double) * time_step));
    memset(file_open_end, 0, (sizeof (double) * time_step));
    memset(file_close_start, 0, (sizeof (double) * time_step));
    memset(file_close_end, 0, (sizeof (double) * time_step));

    var_start = (double**) malloc(sizeof (double*) * time_step);
    var_end = (double**) malloc(sizeof (double*) * time_step);
    rst_start = (double**) malloc(sizeof (double*) * time_step);
    rst_end = (double**) malloc(sizeof (double*) * time_step);
    internal_rst_start = (double**) malloc(sizeof (double*) * time_step);
    internal_rst_end = (double**) malloc(sizeof (double*) * time_step);
    memset(var_start, 0, (sizeof (double*) * time_step));
    memset(var_end, 0, (sizeof (double*) * time_step));
    memset(rst_start, 0, (sizeof (double*) * time_step));
    memset(rst_end, 0, (sizeof (double*) * time_step));
    memset(internal_rst_start, 0, (sizeof (double*) * time_step));
    memset(internal_rst_end, 0, (sizeof (double*) * time_step));

    for (i = 0; i < time_step; i++) 
    {
        var_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
        var_end[i] = (double*) malloc(sizeof (double) * number_of_variables);
        rst_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
        rst_end[i] = (double*) malloc(sizeof (double) * number_of_variables);
        internal_rst_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
        internal_rst_end[i] = (double*) malloc(sizeof (double) * number_of_variables);

        memset(var_start[i], 0, (sizeof (double) * number_of_variables));
        memset(var_end[i], 0, (sizeof (double) * number_of_variables));
        memset(rst_start[i], 0, (sizeof (double) * number_of_variables));
        memset(rst_end[i], 0, (sizeof (double) * number_of_variables));
        memset(internal_rst_start[i], 0, (sizeof (double) * number_of_variables));
        memset(internal_rst_end[i], 0, (sizeof (double) * number_of_variables));
    }
    
    if (write_ACTIVE == 1) 
    {
        if (restructuring_ACTIVE == 0) 
	{
            for (ts = 0; ts < time_step; ts++) 
	    {
                //WRITE BUFFER
                write_data = (double**) malloc(sizeof (double*) * number_of_variables);
                if (!write_data)
                    print_error("Error Allocating Buffer", __FILE__, __LINE__);
                memset(write_data, 0, sizeof (double*) * number_of_variables);

                for (vc = 0; vc < number_of_variables; vc++) 
		{
                    write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);
                    if (!write_data[vc])
                        print_error("Error Allocating Buffer", __FILE__, __LINE__);
                    for (k = 0; k < count[2]; k++)
			for (j = 0; j < count[1]; j++)
			    for (i = 0; i < count[0]; i++) 
			    {
                                long long index = (long long) (count[0] * count[1] * k) + (count[0] * j) + i;
                                for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
                                    write_data[vc][index * sample_per_variable_buffer[vc] + spv] = 100 + (0) + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i);
                            }
                }
                
                file_create_start[ts] = MPI_Wtime();
                if (rank == 0) 
		    mkdir(output_file_name[ts], 0770);    
		
                MPI_Barrier(MPI_COMM_WORLD);
                file_create_end[ts] = MPI_Wtime();

                file_open_start[ts] = MPI_Wtime();
                sprintf(filename, "%s%srank%04d", output_file_name[ts], "/", rank);
                MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
                file_open_end[ts] = MPI_Wtime();

                off_t data_offset = 0;
                for (vc = 0; vc < number_of_variables; vc++) 
		{
		    var_start[ts][vc] = MPI_Wtime();
		    MPI_File_write_at(fh, data_offset, write_data[vc], (count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]), MPI_DOUBLE, &status);
		    data_offset = data_offset + (count[0] * count[1] * count[2]) * sample_per_variable_buffer[vc] * sizeof (double);
		    var_end[ts][vc] = MPI_Wtime();
                }
                file_close_start[ts] = MPI_Wtime();
                MPI_File_close(&fh);
                file_close_end[ts] = MPI_Wtime();

                for (vc = 0; vc < number_of_variables; vc++)
		{
                    free(write_data[vc]);
                    write_data[vc] = 0;
	        }
                free(write_data);
                write_data = 0;
            }
        }
        else 
	{
            int *expanded_box_dimension;
            double* expanded_box;
            PIDX_rst_id rst_id;
            off_t data_offset = 0;
            PIDX_Ndim_buffer*** rst_output_buffer;
            PIDX_Ndim_buffer** in_buf;
            int buffer_size;
            int rst_output_buffer_count = 0, var = 0;

            expanded_box_dimension = (int*)malloc(sizeof(int) * 5);
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
                }file_create_start[ts] = MPI_Wtime();
		if (rank == 0)
		    ret = mkdir(output_file_name[ts], 0770);
		MPI_Barrier(MPI_COMM_WORLD);
		file_create_end[ts] = MPI_Wtime();

		rst_init_start[ts] = MPI_Wtime();
		in_buf = (PIDX_Ndim_buffer**) malloc(number_of_variables * sizeof (PIDX_Ndim_buffer*));
		for (i = 0; i < number_of_variables; i++)
		    in_buf[i] = (PIDX_Ndim_buffer*) malloc(sizeof (PIDX_Ndim_buffer));

		rst_output_buffer = (PIDX_Ndim_buffer***) malloc(sizeof (*rst_output_buffer) * number_of_variables);
		memset(rst_output_buffer, 0, sizeof (*rst_output_buffer) * number_of_variables);
		rst_output_buffer_count = 0;

		rst_id = PIDX_rst_init(MPI_COMM_WORLD, dimension, (int*) gextent, count, offset, 0, NULL, &rst_output_buffer_count);
		rst_init_end[ts] = MPI_Wtime();

		file_create_start[ts] = MPI_Wtime();
		if (rst_output_buffer_count != 0) 
		{
                    sprintf(filename, "%s%srank%04d", output_file_name[ts], "/", rank);
                    MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		}
		file_create_end[ts] = MPI_Wtime();

		data_offset = 0;
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
                    if (rst_output_buffer_count != 0) {
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
                        for (r = 0; r < rst_output_buffer_count; r++) {
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
                    
                    /*
                    if (rst_output_buffer_count != 0)
                    {
                        printf("[%d] Data offset = %d and Data Count %d\n", rank, data_offset, (buffer_size * sample_per_variable_buffer[var]));
                        int uu = 0;
                        
                        for (uu = 0; uu < (buffer_size); uu++)
                            for (spv = 0; spv < sample_per_variable_buffer[var]; spv++)
                                printf("[%d : %d] [%d]  Values[%d] = [%d]\n", var, spv, rank, (uu * sample_per_variable_buffer[var] + spv), (int) expanded_box[uu * sample_per_variable_buffer[var] + spv]);
                        
                        //for (uu = 0; uu < (buffer_size * sample_per_variable_buffer[var]); uu++)
                          //  printf("[%d : %d] [%d %d] values[%d] = [%d]\n", rank, var, data_offset, (buffer_size * sample_per_variable_buffer[var]), uu, (int) expanded_box[uu]);
                    }
                    */ 
                    

                    var_start[ts][var] = MPI_Wtime();
                    if (rst_output_buffer_count != 0) {
                        ret = MPI_File_write_at(fh, data_offset, expanded_box, (buffer_size * sample_per_variable_buffer[var]), MPI_DOUBLE, &status);
                        if (ret != MPI_SUCCESS) {
                            return (-1);
                        }
                        data_offset = data_offset + (buffer_size) * sample_per_variable_buffer[var] * sizeof (double);
                        free(expanded_box);
                        expanded_box = 0;
                    }
                    var_end[ts][var] = MPI_Wtime();

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
                    rst_cleanup_end[ts] = MPI_Wtime();
                }
                free(rst_output_buffer);
                rst_output_buffer = 0;

                file_close_start[ts] = MPI_Wtime();
                if (rst_output_buffer_count != 0)
                    MPI_File_close(&fh);
                file_close_end[ts] = MPI_Wtime();


                for (vc = 0; vc < number_of_variables; vc++) {
                    free(write_data[vc]);
                    write_data[vc] = 0;
                }
                free(write_data);
                write_data = 0;
            }
        }
        double max_time = 0, total_time = 0;
        long long total_data = 0;
        int var = 0, sample_sum = 0;
        for (var = 0; var < number_of_variables; var++) {
            sample_sum = sample_sum + sample_per_variable_buffer[var];
        }
        total_data = (gextent[0] * gextent[1] * gextent[2] * sample_sum * sizeof (double));
        for (ts = 0; ts < time_step; ts++) {
            max_time = 0, total_time = 0;
            total_time = file_close_end[ts] - file_create_start[ts];

            MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            if (total_time == max_time) {
                printf("[Rank %d Nproc %d] [Time Step %d] [Total Time %f] [Total Volume %d x %d x %d x %d] [Throughput %f MiB/sec]\n", rank, nprocs, ts, max_time, gextent[0], gextent[1], gextent[2], sample_sum, (float) total_data / (1024 * 1024 * max_time));
                printf("Directory Creation Time: %f\n", (file_create_end[ts] - file_create_start[ts]));
                if (restructuring_ACTIVE == 1)
                    printf("Restructuring Init Time: %f\n", (rst_init_end[ts] - rst_init_start[ts]));

                printf("File Open Time: %f\n", (file_open_end[ts] - file_open_start[ts]));
                for (var = 0; var < number_of_variables; var++) {
                    printf("[%d] Variable IO Time %f\n", var, (var_end[ts][var] - var_start[ts][var]));
                    if (restructuring_ACTIVE == 1) {
                        printf("[%d] Retructuring Time :: Internal Restructuring Time [%f %f]\n", var, (rst_end[ts][var] - rst_start[ts][var]), (internal_rst_end[ts][var] - internal_rst_end[ts][var]));
                    }
                }
                printf("File Close Time: %f\n", (file_close_end[ts] - file_close_start[ts]));
                if (restructuring_ACTIVE == 1) {
                    printf("Cleanup Time: %f\n", (rst_cleanup_end[ts] - rst_cleanup_start[ts]));
                }
                printf("\n");
            }
        }
    } else {
        for (ts = 0; ts < time_step; ts++) {
            
            //WRITE BUFFER
            write_data = (double**) malloc(sizeof (double*) * number_of_variables);
            if (!write_data)
                print_error("Error Allocating Buffer", __FILE__, __LINE__);
            memset(write_data, 0, sizeof (double*) * number_of_variables);

            for (vc = 0; vc < number_of_variables; vc++) {
                write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * count[3] * count[4] * sample_per_variable_buffer[vc]);
                if (!write_data[vc])
                    print_error("Error Allocating Buffer", __FILE__, __LINE__);
            }

            sprintf(filename, "%s%srank%04d", output_file_name[ts], "/", rank);
            //printf("Filename : %s : [%d %d %d %d %d]\n", filename, count[0], count[1], count[2], count[3], count[4]);
            ret = MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            if (ret != MPI_SUCCESS) {

            }
            off_t data_offset = 0;
            for (vc = 0; vc < number_of_variables; vc++) {
                ret = MPI_File_read_at(fh, data_offset, write_data[vc], (count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]), MPI_DOUBLE, &status);
                if (ret != MPI_SUCCESS) {
                    return (-1);
                }
                /*
                {
                    printf("Data offset = %d and Data Count %d\n", data_offset, (count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]));
                    int uu = 0;
                    for (uu = 0; uu < (count[0] * count[1] * count[2]); uu++)
                        for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
                            printf("[%d : %d] [%d]  ValueS[%d] = [%d]\n", vc, spv, rank, uu, (int) write_data[vc][uu * sample_per_variable_buffer[vc] + spv]);
                }
                 */

                data_offset = data_offset + (count[0] * count[1] * count[2]) * sample_per_variable_buffer[vc] * sizeof (double);
                int lost_count = 0;
                for (k = 0; k < count[2]; k++)
                    for (j = 0; j < count[1]; j++)
                        for (i = 0; i < count[0]; i++) {
                            long long index = (long long) (count[0] * count[1] * k) + (count[0] * j) + i;
                            for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++) {
                                if ((int) write_data[vc][index * sample_per_variable_buffer[vc] + spv] != 100 + (ts) + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i)) {
                                    //printf("[BAD %d] [%d %d] [%d] Values : %d = %d\n", lost_count, vc, spv, index, (int) write_data[vc][index * sample_per_variable_buffer[vc] + spv], (100 + (ts) + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i)));
                                    lost_count++;
                                }
                                else {
                                    //printf("[%d] [GOOD %d] Values : %d = %d\n", rank, lost_count, (int) write_data[vc][index * sample_per_variable_buffer[vc] + spv], (100 + (ts) + spv + (extents[0] * extents[1] * extents[2] * extents[3] * (offset[4] + v)) + (extents[0] * extents[1] * extents[2] * (offset[3] + u)) + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i)));
                                }
                            }
                        }
                assert(lost_count == 0);
            }
            MPI_File_close(&fh);

            for (vc = 0; vc < number_of_variables; vc++) {
                free(write_data[vc]);
                write_data[vc] = 0;
            }
            free(write_data);
            write_data = 0;
            if(rank == 0)
                printf("Verifying Time Step %d\n", ts);
        }
    }
    free(sample_per_variable_buffer);
    sample_per_variable_buffer = 0;

    MPI_Finalize();
    return 0;
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:l:f:t:";
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
            case('f'):
                sprintf(output_file_template, "%s", optarg);
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
    printf("Usage: test-PIDX -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("  -f: IDX Filename\n");
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
