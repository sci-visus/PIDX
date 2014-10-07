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
#define MPI 1

#ifdef MPI
#include "mpi.h"
#endif


static int parse_args(int argc, char **argv);
static void usage(void);
static int print_error(char *error_message, char* file, int line);


/* Number of time-steps */
static int time_step = 0;



static MPI_Win win;

int main(int argc, char **argv) {
    
    int i = 0, j = 0, k = 0;
    int spv = 0; /*samples per variable*/
    int ts, var; /*time step counter and variable counter*/
    int nprocs, rank; /*process count and rank*/
    int slice, ret;
  
    int retval = 0;
    

    int number_of_variables = 1;
    int* sample_per_variable_buffer;
    double *start_time, *end_time;
    double** write_data;
    write_data = (double**) malloc(sizeof (double*) * number_of_variables);
    if (!write_data)
        print_error("Error Allocating Buffer", __FILE__, __LINE__);
    memset(write_data, 0, sizeof (double*) * number_of_variables);

    
    
   
    
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
    }
    
#ifdef MPI
    /*   The command line arguments are shared by all processes  */
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    start_time = (double*)malloc(time_step * sizeof(double));
    memset(start_time, 0, time_step * sizeof(double));
    end_time = (double*)malloc(time_step * sizeof(double));
    memset(end_time, 0, time_step * sizeof(double));   
    

    sample_per_variable_buffer = (int*) malloc(sizeof (int) * number_of_variables);
    if (!sample_per_variable_buffer)
        print_error("Error Allocating {sample_per_variable_buffer} buffer", __FILE__, __LINE__);

    sample_per_variable_buffer[0] = 16;
    //sample_per_variable_buffer[1] = 4;
    //sample_per_variable_buffer[2] = 4;
    //sample_per_variable_buffer[3] = 4;

    long long local_process_count = 32 * 32 * 32;
    for (var = 0; var < number_of_variables; var++) {
        write_data[var] = (double*) malloc(sizeof (double) * local_process_count * sample_per_variable_buffer[var]);
        if (!write_data[var])
            print_error("Error Allocating Buffer", __FILE__, __LINE__);
        //memset(write_data[var], rank, sizeof (double) * local_process_count * sample_per_variable_buffer[var]);
        for(i = 0 ; i < local_process_count * sample_per_variable_buffer[var] ; i++)
            write_data[var][i] = rank;
    }
    //printf("[%d] VALUE For : %d\n", rank, (int)write_data[0][0]);
    double* buffer;
    if (rank == nprocs - 2)
    {
        buffer = (double*) malloc(sizeof (double) * local_process_count * sample_per_variable_buffer[0] * nprocs);
        memset(buffer, 0, sizeof (double) * local_process_count * sample_per_variable_buffer[0] * nprocs);
    }
    for(ts = 0 ; ts < time_step ; ts++)
    {
            start_time[ts] = MPI_Wtime();
            if (rank == nprocs - 2) 
            {
                retval = MPI_Win_create(buffer, (sizeof(double) * nprocs * local_process_count), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &(win));
                if (MPI_SUCCESS != retval) 
                {
                    fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                    MPI_Abort(MPI_COMM_WORLD, retval);
                }
            }
            else 
            {
                retval = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &(win));
                if (MPI_SUCCESS != retval) 
                {
                    fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] No agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                    MPI_Abort(MPI_COMM_WORLD, retval);
                }
            }
            retval = MPI_Win_fence(0, win);
            if (MPI_SUCCESS != retval) 
            {
                fprintf(stderr, "Error Creating MPI_Win_fence\n");
                MPI_Abort(MPI_COMM_WORLD, retval);
            }

            for (var = 0; var < number_of_variables; var++) 
            {
                retval = MPI_Put(write_data[var], local_process_count * sample_per_variable_buffer[var] , MPI_DOUBLE, (nprocs - 2), (rank  * local_process_count * sample_per_variable_buffer[var]), (local_process_count * sample_per_variable_buffer[var]), MPI_DOUBLE, win);
                if (MPI_SUCCESS != retval) 
                {
                    fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] MPI_Put\n", rank, __FILE__, __LINE__);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }

            retval = MPI_Win_fence(0, win);
            if (MPI_SUCCESS != retval) 
            {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] MPI_Put\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }

            retval = MPI_Win_free(&(win));
            if (MPI_SUCCESS != retval) 
            {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] MPI_Put\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }
            end_time[ts] = MPI_Wtime();
            
            
            if (rank == nprocs - 2)
            {
                for(i = 0 ; i < nprocs  ; i++)
                {
                    for(j = 0 ; j < local_process_count * sample_per_variable_buffer[0] ; j++)
                    {
                        if((int)buffer[i * local_process_count * sample_per_variable_buffer[0] + j] != i)
                            printf("[%d] [%d : %d] SCREAM!!!!!!!! %d : %d\n", i * local_process_count * sample_per_variable_buffer[0] + j, i, j, (int)buffer[i * local_process_count * sample_per_variable_buffer[0] + j],  j);
                    }
                }
                memset(buffer, 0, sizeof (double) * local_process_count * sample_per_variable_buffer[0] * nprocs);
            }
            MPI_Barrier(MPI_COMM_WORLD);
    }
    
    double max_time = 0, total_time = 0;
    for(ts = 0 ; ts < time_step ; ts++)
    {
        total_time = end_time[ts] - start_time[ts];
        MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (max_time == total_time) {
            printf("[%d] Aggregation Time for Time Step %d = %f\n", rank, ts, total_time);
        }
        
    }
    
    if(rank == nprocs - 2)
    {
        free(buffer);
        buffer = 0;
    }

    for (var = 0; var < number_of_variables; var++) 
    {
        free(write_data[var]);
        write_data[var] = 0;
    }
    free(write_data);
    write_data = 0;
    
    free(sample_per_variable_buffer);
    sample_per_variable_buffer = 0;

    MPI_Finalize();
}

static int parse_args(int argc, char **argv) {
    char flags[] = "t:";
    int one_opt = 0, i = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        /* postpone error checking for after while loop */
        switch (one_opt) {
            case('t'):
                sscanf(optarg, "%d", &time_step);
                break;
            case('?'):
                return (-1);
        }
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

static int print_error(char *error_message, char* file, int line) {
    fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(0);
#endif
}
