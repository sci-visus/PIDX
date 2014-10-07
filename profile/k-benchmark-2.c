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

#ifdef MPI
#include "mpi.h"
#endif

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

static int bpf = 0;
static int aggf = 1;

/* output IDX file Name Template*/
static char output_file_template[512] = {0};
static char *output_file_name;

int main(int argc, char **argv) {
#if 0 /* Sidharth fix */    
    int i = 0, j = 0, k = 0, u = 0, v = 0;
    int spv = 0; /*samples per variable*/
    int ts, vc; /*time step counter and variable counter*/
    int nprocs, rank; /*process count and rank*/
    int slice, ret;
    int sub_div[5], offset_local[5];

    const int total_dimensions = 3; /*Total dimension of the dataset (3)*/
    const char *output_file; /*Output File Name*/
    const int *gextent; /*Global Extensions of the dataset (64 64 64 0 0)*/
    const int *count; /*Local extents of each process*/
    const int *offset; /*Local counts of each process*/
    const int bits_per_block = 15; /*Total number of samples in each block*/

    int number_of_variables = 1;
    char var_name[512];
    int* sample_per_variable_buffer;
    
    PIDX_file idx_ptr;
    PIDX_variable *variable_ptr;

    double** write_data;

#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
    nprocs = 1;
#endif


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
    MPI_Bcast(&bpf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&aggf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

    /*Creating the filename*/
    output_file_name = (char*) malloc(sizeof (char) * 512);
    sprintf(output_file_name, "%s%s", output_file_template, ".idx");
    //blocks_per_file = bpf;

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
    output_file = output_file_name;
    assert(offset[0] < extents[0]);
    assert(offset[1] < extents[1]);
    assert(offset[2] < extents[2]);
    assert(offset[3] < extents[3]);
    assert(offset[4] < extents[4]);

    /*
    Initializing All Variable Buffer for the Simulation
    In this sample run we have 4 variables (following the pattern of S3D simulation)
    First Variable has one sample per data point (Pressure)
    Second Variable has one sample per data point (Temperature)
    Third Variable has three sample per data point (Velocity)
    Fourth Variable has eleven sample per data point (Species)
     
    All time step writes the same data as input.
     */
    sample_per_variable_buffer = (int*) malloc(sizeof (int) * number_of_variables);
    if (!sample_per_variable_buffer)
        print_error("Error Allocating {sample_per_variable_buffer} buffer", __FILE__, __LINE__);

    sample_per_variable_buffer[0] = 1; //Variable 1
    //sample_per_variable_buffer[1] = 1; //Variable 2
    //sample_per_variable_buffer[2] = 1; //Variable 3
    //sample_per_variable_buffer[3] = 1; //Variable 4


    // Allocating variable pointer for all variables.
    variable_ptr = (PIDX_variable*) malloc(sizeof (PIDX_variable) * number_of_variables);
    if (!variable_ptr)
        print_error("Error Allocating Variable pointer", __FILE__, __LINE__);
    memset(variable_ptr, 0, sizeof (PIDX_variable) * number_of_variables);


    PIDX_time_step_define(0, time_step, "time%04d/");
    for (ts = 0; ts < time_step; ts++) {
        /***********************************WRITE IDX FILE****************************************/ 

        //WRITE BUFFER
        write_data = (double**) malloc(sizeof (double*) * number_of_variables);
        if (!write_data)
            print_error("Error Allocating Buffer", __FILE__, __LINE__);
        memset(write_data, 0, sizeof (double*) * number_of_variables);

        for (vc = 0; vc < number_of_variables; vc++) {
            write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * count[3] * count[4] * sample_per_variable_buffer[vc]);
            if (!write_data[vc])
                print_error("Error Allocating Buffer", __FILE__, __LINE__);

            
            //if (ts == 0 || ts == (time_step - 1)) {
                for (v = 0; v < count[4]; v++)
                    for (u = 0; u < count[3]; u++)
                        for (k = 0; k < count[2]; k++)
                            for (j = 0; j < count[1]; j++)
                                for (i = 0; i < count[0]; i++) {
                                    long long index = (long long) (count[0] * count[1] * count[2] * count[3] * v) + (count[0] * count[1] * count[2] * u) + (count[0] * count[1] * k) + (count[0] * j) + i;
                                    for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
                                        write_data[vc][index * sample_per_variable_buffer[vc] + spv] = 100 + spv + ((extents[0] * extents[1] * extents[2] * extents[3] * (offset[4] + v)) + (extents[0] * extents[1] * extents[2] * (offset[3] + u)) + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i));
                                }
            //}
            /*
            else{
              
                for (v = 0; v < count[4]; v++)
                    for (u = 0; u < count[3]; u++)
                        for (k = 0; k < count[2]; k++)
                            for (j = 0; j < count[1]; j++)
                                for (i = 0; i < count[0]; i++) {
                                    long long index = (long long) (count[0] * count[1] * count[2] * count[3] * v) + (count[0] * count[1] * count[2] * u) + (count[0] * count[1] * k) + (count[0] * j) + i;
                                    for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
                                        write_data[vc][index * sample_per_variable_buffer[vc] + spv] = (double) (ts+1) * rank;//) * 256)/nprocs;
                                }
            }
            */
        }

        //MPI_Barrier(MPI_COMM_WORLD);

        idx_ptr = PIDX_create(MPI_COMM_WORLD, ts, output_file, bits_per_block, bpf, total_dimensions, aggf, gextent, 1);
        for (vc = 0; vc < number_of_variables; vc++) {
            sprintf(var_name, "var_%d", vc);
            variable_ptr[vc] = PIDX_variable_global_define(idx_ptr, var_name, sample_per_variable_buffer[vc], MPI_DOUBLE);
            ret = PIDX_variable_local_add(idx_ptr, variable_ptr[vc], (int*) offset, (int*) count);
            ret = PIDX_variable_local_layout(idx_ptr, variable_ptr[vc], write_data[vc], MPI_DOUBLE);
        }
        ret = PIDX_write(idx_ptr);
        ret = PIDX_close(idx_ptr);

        for (vc = 0; vc < number_of_variables; vc++) {
            free(write_data[vc]);
            write_data[vc] = 0;
        }
        free(write_data);
        write_data = 0;
    }
    free(variable_ptr);
    variable_ptr = 0;
    free(sample_per_variable_buffer);
    sample_per_variable_buffer = 0;

    MPI_Finalize();
#endif
    
    return 0;
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:l:f:t:b:";
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
            case('b'):
                sscanf(optarg, "%d", &bpf);
                break;
            //case('a'):
            //    sscanf(optarg, "%d", &aggf);
            //    break;
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
    printf("Usage: test-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4 -b 128\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("  -f: IDX Filename\n");
    printf("  -t: number of timesteps\n");
    printf("  -b: number of blocks per file\n");
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

