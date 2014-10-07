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
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define PIDX_MAX_DIMENSIONS 5

static int agg_count;
static int agg_participate_count;
static int time_step;

static int parse_args(int argc, char **argv);
static void usage(void);

/* global dimensions of 3D volume */
static int* extents = 0;
static char** dir_path;
static char dir_path_template[1024];

static int regular_box_dim[5];
static MPI_Win win;

static int number_of_variables;
static int FACTOR;
static int sub_FACTOR;

//Function to find the power of 2 of an integer value (example 5->8)

int getPowerOftwo(int x) {
    int n = 1;
    while (n < x)
        n <<= 1;
    return n;
}

int main(int argc, char **argv) {
    int i, j, k, ts;
    int retval = 0, ret = 0;
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_File fh;
    MPI_Status status;
    int V_number, F_number;
    int nfiles = 0;
    int agg_slot, proc_slot;
    int agg_count = 0, proc_count = 0;
    int r_is_agg = 0, r_is_proc = 0;
    long long global_volume;
    long long sample_per_file = 512 * 32768;

    int **rank_holder;
    int var_counter = 0, file_counter = 0;
    double* proc_buffer;
    double** proc_buffer_with_chunk;
    double* agg_buffer;
    long long proc_buffer_count;
    long long agg_buffer_count;
    double *agg_start, *agg_end, *init_time_start, *init_time_end, *file_create_start, *file_create_end, *io_start_time, *io_end_time;
    off_t data_offset = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /*Rank 0 parses the command Line Arguments*/
    extents = (int*) malloc(3 * sizeof (int));
    assert(extents);
    memset(extents, 0, 3 * sizeof (int));


    if (rank == 0) {
        ret = parse_args(argc, argv);
        if (ret < 0) {
            usage();
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }


    MPI_Bcast(&extents[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&extents[1], 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&extents[2], 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_variables, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&FACTOR, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sub_FACTOR, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dir_path_template, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);

    agg_start = (double*) malloc(sizeof (double) * time_step);
    agg_end = (double*) malloc(sizeof (double) * time_step);
    init_time_start = (double*) malloc(sizeof (double) * time_step);
    init_time_end = (double*) malloc(sizeof (double) * time_step);
    file_create_start = (double*) malloc(sizeof (double) * time_step);
    file_create_end = (double*) malloc(sizeof (double) * time_step);
    io_start_time = (double*) malloc(sizeof (double) * time_step);
    io_end_time = (double*) malloc(sizeof (double) * time_step);


    dir_path = (char**) malloc(time_step * sizeof (char*));
    for (ts = 0; ts < time_step; ts++) {
        dir_path[ts] = (char*) malloc(1024 * sizeof (char));
        sprintf(dir_path[ts], "%s_%d", dir_path_template, ts);
        //printf("The Paths are %s\n", dir_path[ts]);
    }


    //number_of_variables = 2;
    //int FACTOR = 2;
    //int sub_FACTOR = 8;

    global_volume = (long long) (extents[0] * extents[1] * extents[2]) * number_of_variables;
    long long local_volume = (long long) ((long long) FACTOR * ((long long) extents[0] * extents[1] * extents[2]) * number_of_variables) / nprocs;

    //if(rank == 0)
    //printf("Global and Local :: %lld : %lld\n", local_volume, global_volume);


    proc_count = nprocs / FACTOR;
    //printf("nprocs : FACTOR : PROC_COUNT :: %d : %d : %d\n", nprocs, FACTOR, proc_count);
    assert(nprocs % proc_count == 0);
    proc_slot = nprocs / proc_count;


    nfiles = (getPowerOftwo(extents[0]) * getPowerOftwo(extents[1]) * getPowerOftwo(extents[2])) / ((unsigned long long) sample_per_file);
    if ((getPowerOftwo(extents[0]) * getPowerOftwo(extents[1]) * getPowerOftwo(extents[2])) % ((unsigned long long) sample_per_file))
        nfiles++;

    agg_count = nfiles * number_of_variables;
    assert(nprocs % agg_count == 0);
    agg_slot = nprocs / agg_count;

    /*
    if(rank == 0)
    {
        printf("Process Count : Process Slot :: %d : %d\n", proc_count, proc_slot);
        printf("[Files %d] : Aggregator Count : Aggregator Slot :: %d : %d\n", nfiles, agg_count, agg_slot);
    }
     */

    int afloat = 0, pfloat = 0;
    for (afloat = 0; afloat < nprocs; afloat = afloat + agg_slot) {
        if (rank == (int) afloat)
            r_is_agg = 1;
    }
    for (pfloat = 0; pfloat < nprocs; pfloat = pfloat + proc_slot) {
        if (rank == (int) pfloat)
            r_is_proc = 1;
    }

    rank_holder = (int**) malloc(nfiles * sizeof (int*));
    assert(rank_holder);
    for (i = 0; i < nfiles; i++) {
        rank_holder[i] = (int*) malloc(number_of_variables * sizeof (int));
        assert(rank_holder[i]);
    }

    for (i = 0; i < nprocs; i = i + agg_slot) {
        if (file_counter == (nfiles - 1) && var_counter == (number_of_variables))
            break;

        if (var_counter == number_of_variables) {
            file_counter = file_counter + 1;
            var_counter = 0;
        }
        rank_holder[file_counter][var_counter] = i;

        if (rank == rank_holder[file_counter][var_counter]) {
            V_number = var_counter;
            F_number = file_counter;
            //printf("For Rank [%d] : VC : FC :: %d : %d\n", rank, V_number, F_number);
        }

        var_counter = var_counter + 1;
    }

    /*
    if(rank == 0)
      for(i = 0 ; i < nfiles ; i++)
          for(j = 0 ; j < number_of_variables ; j++)
              printf("[%d] : RANK HOLDER   :   File %d Variable %d : %d\n", rank, i, j, rank_holder[i][j]);
     */
    //creating buffer

    for (ts = 0; ts < time_step; ts++) {
        init_time_start[ts] = MPI_Wtime();
        int MODE = 0;
        if (r_is_proc == 1) {
            proc_buffer_count = local_volume;
            if (MODE == 0) {
                proc_buffer = (double*) malloc(sizeof (double) * (long long) proc_buffer_count);
                assert(proc_buffer);

                for (j = 0; j < proc_buffer_count; j++)
                    proc_buffer[j] = (double) 17.85;
            } else {
                proc_buffer_with_chunk = (double**) malloc(sizeof (double*) * (FACTOR * sub_FACTOR));
                assert(proc_buffer_with_chunk);
                for (i = 0; i < (FACTOR * sub_FACTOR); i++) {
                    proc_buffer_with_chunk[i] = (double*) malloc(sizeof (double) * proc_buffer_count / (FACTOR * sub_FACTOR));
                    assert(proc_buffer_with_chunk[i]);
                }

                for (i = 0; i < (FACTOR * sub_FACTOR); i++) {
                    for (j = 0; j < proc_buffer_count / ((FACTOR * sub_FACTOR)); j++)
                        proc_buffer_with_chunk[i][j] = (double) 17.85;
                }
            }

        }
        if (r_is_agg == 1) {
            agg_buffer_count = sample_per_file;
            //printf("[%d] Aggregating Process : %d\n", rank, agg_buffer_count);

            agg_buffer = (double*) malloc(sample_per_file * sizeof (double));
            assert(agg_buffer);
            memset(agg_buffer, 0, sample_per_file * sizeof (double));

        }


        agg_start[ts] = MPI_Wtime();

        //initiating communication
        if (r_is_agg == 1) {
            //printf("[%d]Aggregation Process : Buffer Size :: %d\n", rank, sample_per_file * sizeof(double));

            retval = MPI_Win_create(agg_buffer, sample_per_file * sizeof (double), sizeof (double), MPI_INFO_NULL, MPI_COMM_WORLD, &(win));
            if (MPI_SUCCESS != retval) {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }
        } else {
            retval = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &(win));
            if (MPI_SUCCESS != retval) {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] No agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }
        }

        retval = MPI_Win_fence(0, win);
        if (MPI_SUCCESS != retval) {
            fprintf(stderr, "Error Creating MPI_Win_fence\n");
            MPI_Abort(MPI_COMM_WORLD, retval);
        }


        assert(proc_count % nfiles == 0);
        if (r_is_proc == 1) {
            for (i = 0; i < number_of_variables; i++) {
                assert(proc_buffer_count % number_of_variables == 0);

                //if(rank == 1)
                //printf("[%d] Local Count %d : [%d %d] Target Rank %d : Target Offset %d\n", rank, (proc_buffer_count/number_of_variables), rank/((FACTOR * proc_count)/nfiles), i, rank_holder[rank/((FACTOR * proc_count)/nfiles)][i], (rank%((FACTOR * proc_count)/nfiles)) / FACTOR);

                if (MODE == 0) {

                    retval = MPI_Put(proc_buffer, (proc_buffer_count / number_of_variables), MPI_DOUBLE, rank_holder[rank / ((FACTOR * proc_count) / nfiles)][i], ((rank % ((FACTOR * proc_count) / nfiles)) / FACTOR) * (proc_buffer_count / number_of_variables), (proc_buffer_count / number_of_variables), MPI_DOUBLE, win);
                    if (MPI_SUCCESS != retval) {
                        fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] No agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                        MPI_Abort(MPI_COMM_WORLD, retval);
                    }

                } else {
                    for (k = 0; k < (FACTOR * sub_FACTOR); k++) {
                        retval = MPI_Put(proc_buffer_with_chunk[k], (proc_buffer_count / (number_of_variables * FACTOR * sub_FACTOR)), MPI_DOUBLE, rank_holder[rank / ((FACTOR * proc_count) / nfiles)][i], (((rank % ((FACTOR * proc_count) / nfiles)) / FACTOR) * (proc_buffer_count / number_of_variables)) + (k * (proc_buffer_count / (number_of_variables * FACTOR * sub_FACTOR))), (proc_buffer_count / (number_of_variables * FACTOR * sub_FACTOR)), MPI_DOUBLE, win);
                        if (MPI_SUCCESS != retval) {
                            fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] No agg : MPI_Win_create\n", rank, __FILE__, __LINE__);
                            MPI_Abort(MPI_COMM_WORLD, retval);
                        }

                    }
                }
            }
        }


        retval = MPI_Win_fence(0, win);
        if (MPI_SUCCESS != retval) {
            fprintf(stderr, "Error Creating MPI_Win_fence\n");
            MPI_Abort(MPI_COMM_WORLD, retval);
        }

        retval = MPI_Win_free(&win);
        if (MPI_SUCCESS != retval) {
            fprintf(stderr, "Error Creating MPI_Win_fence\n");
            MPI_Abort(MPI_COMM_WORLD, retval);
        }
        if (rank == 0)
            printf("Done aggregating for Time Step %d\n", ts);

        file_create_start[ts] = MPI_Wtime();
        if (rank == 0) {
            ret = mkdir(dir_path[ts], 0770);
            if (ret != 0 && errno != EEXIST) {
                fprintf(stderr, "Error: failed to mkdir %s\n", dir_path[ts]);
                MPI_Abort(MPI_COMM_WORLD, errno);
            }
        }

        for (i = 0; i < nfiles; i++) {
            if (rank == 0) {
                char bin_file[1024];
                sprintf(bin_file, "%s/file_%d.bin", dir_path[ts], i);
                //printf("File Name is %s\n", bin_file);
                ret = MPI_File_open(MPI_COMM_SELF, bin_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
                if (ret != MPI_SUCCESS) {
                    fprintf(stderr, "Error: failed to mkdir %d\n", i);
                    MPI_Abort(MPI_COMM_WORLD, errno);
                }

                ret = MPI_File_close(&fh);
                if (ret != MPI_SUCCESS) {
                    fprintf(stderr, "Error: failed to mkdir %d\n", i);
                    MPI_Abort(MPI_COMM_WORLD, errno);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        io_start_time[ts] = MPI_Wtime();


        if (r_is_agg == 1) {
            MPI_File filep;
            char file_name[1024];
            sprintf(file_name, "%s/file_%d.bin", dir_path[ts], F_number);
            //printf(" Rank [%d] : Filename :%s\n", rank, file_name);
            ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &filep);
            if (ret != MPI_SUCCESS) {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] Opening File\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }
            //printf("[Rank %d] : V_number %d\n", rank, V_number);
            for (k = 0; k <= V_number; k++) {
                data_offset = k * agg_buffer_count * sizeof (double);
                //printf("[%d] Data offset %lld [%d : %d x %d] : Length %d\n", rank, data_offset, V_number, k, agg_buffer_count, agg_buffer_count);
            }
            //printf("[%d] : Data offset %lld [%d : %d x %d] : Length %d\n", rank, data_offset, V_number, k, agg_buffer_count, agg_buffer_count);

            ret = MPI_File_write_at(filep, data_offset, agg_buffer, agg_buffer_count, MPI_DOUBLE, &status);
            if (ret != MPI_SUCCESS) {
                fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] Opening File\n", rank, __FILE__, __LINE__);
                MPI_Abort(MPI_COMM_WORLD, retval);
            }

            MPI_File_close(&filep);
        }
        io_end_time[ts] = MPI_Wtime();

        if (r_is_agg == 1) {
            free(agg_buffer);
            agg_buffer = 0;
        }
        if (r_is_proc == 1) {
            if (MODE == 0) {
                free(proc_buffer);
                proc_buffer = 0;
            } else {
                for (i = 0; i < (FACTOR * sub_FACTOR); i++) {
                    free(proc_buffer_with_chunk[i]);
                    proc_buffer_with_chunk[i] = 0;
                }
                free(proc_buffer_with_chunk);
                proc_buffer_with_chunk = 0;
            }
        }
    }

    double avg_throughput = 0, sdev = 0, throughput_sum = 0, throughput_max = 0, sdev_sum = 0, sdev_sum1 = 0, sdev_sum2 = 0;
    double sdev_agg, sdev_file, sdev_io;

    double *total_time, *max_time, *throughput;
    double agg_ts_sum = 0, file_ts_sum = 0, io_ts_sum = 0, agg_ts_max, file_ts_max, io_ts_max;
    double *agg_ts, *io_ts, *file_ts;

    agg_ts = (double*) malloc(sizeof (double) * time_step);
    memset(agg_ts, 0, sizeof (double) * time_step);
    io_ts = (double*) malloc(sizeof (double) * time_step);
    memset(io_ts, 0, sizeof (double) * time_step);
    file_ts = (double*) malloc(sizeof (double) * time_step);
    memset(file_ts, 0, sizeof (double) * time_step);

    total_time = (double*) malloc(sizeof (double) * time_step);
    memset(total_time, 0, sizeof (double) * time_step);

    max_time = (double*) malloc(sizeof (double) * time_step);
    memset(max_time, 0, sizeof (double) * time_step);

    throughput = (double*) malloc(sizeof (double) * time_step);
    memset(throughput, 0, sizeof (double) * time_step);

    for (ts = 0; ts < time_step; ts++) {
        double total_agg = 0, max_agg = 0, total_file = 0, max_file = 0, total_io = 0, max_io = 0;
        double avg_agg, avg_file, avg_io;
        total_time[ts] = io_end_time[ts] - agg_start[ts];
        MPI_Allreduce(&total_time[ts], &max_time[ts], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        total_agg = (file_create_start[ts] - agg_start[ts]);
        total_file = (io_start_time[ts] - file_create_start[ts]);
        total_io = (io_end_time[ts] - io_start_time[ts]);
        MPI_Allreduce(&total_agg, &max_agg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&total_agg, &avg_agg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        avg_agg = avg_agg / (nprocs / FACTOR);

        MPI_Allreduce(&total_file, &max_file, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&total_file, &avg_file, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        avg_file = avg_file / (nprocs / FACTOR);

        MPI_Allreduce(&total_io, &max_io, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&total_io, &avg_io, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        avg_io = avg_io / (nprocs / FACTOR);

        if (max_time[ts] == total_time[ts]) {
            long long total_data = (long long) extents[0] * extents[1] * extents[2] * number_of_variables * sizeof (double);

            throughput[ts] = (double) total_data / (1024 * 1024 * max_time[ts]);
            throughput_sum = throughput_sum + (double) total_data / (1024 * 1024 * max_time[ts]); //throughput;

            agg_ts_sum = agg_ts_sum + (file_create_start[ts] - agg_start[ts]);
            agg_ts[ts] = (file_create_start[ts] - agg_start[ts]);

            file_ts_sum = file_ts_sum + (io_start_time[ts] - file_create_start[ts]);
            file_ts[ts] = (io_start_time[ts] - file_create_start[ts]);

            io_ts_sum = io_ts_sum + (io_end_time[ts] - io_start_time[ts]);
            io_ts[ts] = (io_end_time[ts] - io_start_time[ts]);


            printf("\n[R%d T%d]Time Taken [%d %d %d] : [Proc %d] : [Factor %d] : [sub_FACTOR %d] : [Variables %d]: %f Seconds Throughput %f MiB/sec\n", rank, ts, extents[0], extents[1], extents[2], nprocs, FACTOR, sub_FACTOR, number_of_variables, max_time[ts], (float) total_data / (1024 * 1024 * max_time[ts]));
            printf("Agg Time : [%f] : [%f : %f]\n", (file_create_start[ts] - agg_start[ts]), max_agg, avg_agg);
            printf("File Time : [%f] : [%f : %f]\n", (io_start_time[ts] - file_create_start[ts]), max_file, avg_file);
            printf("I/O Time : [%f] : [%f : %f]\n", (io_end_time[ts] - io_start_time[ts]), max_io, avg_io);
        }
    }


    MPI_Allreduce(&throughput_sum, &throughput_max, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    throughput_max = throughput_max / time_step;

    MPI_Allreduce(&agg_ts_sum, &agg_ts_max, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    agg_ts_max = agg_ts_max / time_step;
    MPI_Allreduce(&file_ts_sum, &file_ts_max, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    file_ts_max = file_ts_max / time_step;
    MPI_Allreduce(&io_ts_sum, &io_ts_max, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    io_ts_max = io_ts_max / time_step;

    for (ts = 0; ts < time_step; ts++) {
        if (max_time[ts] == total_time[ts]) {
            sdev = sdev + (throughput[ts] - throughput_max) * (throughput[ts] - throughput_max);
            sdev_agg = sdev_agg + (agg_ts[ts] - agg_ts_max) * (agg_ts[ts] - agg_ts_max);
            sdev_file = sdev_file + (file_ts[ts] - file_ts_max) * (file_ts[ts] - file_ts_max);
            sdev_io = sdev_io + (io_ts[ts] - io_ts_max) * (io_ts[ts] - io_ts_max);
        }
    }
    MPI_Allreduce(&sdev, &sdev_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdev_sum1 = sdev_sum / time_step;
    //       sdev_sum1 = sqrt(sdev_sum1);

    sdev_sum2 = sdev_sum / (time_step - 1);
    //       sdev_sum2 = sqrt(sdev_sum2);

    double sdev_agg_sum = 0, sdev_file_sum = 0, sdev_io_sum = 0;
    double sdev_agg_sum1 = 0, sdev_file_sum1 = 0, sdev_io_sum1 = 0;
    double sdev_agg_sum2 = 0, sdev_file_sum2 = 0, sdev_io_sum2 = 0;

    MPI_Allreduce(&sdev_agg, &sdev_agg_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdev_agg_sum1 = sdev_agg_sum / time_step;
    //       sdev_agg_sum1 = sqrt(sdev_agg_sum1);

    sdev_agg_sum2 = sdev_agg_sum / (time_step - 1);
    //       sdev_agg_sum2 = sqrt(sdev_agg_sum2);

    MPI_Allreduce(&sdev_file, &sdev_file_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdev_file_sum1 = sdev_file_sum / time_step;
    //       sdev_file_sum1 = sqrt(sdev_file_sum1);

    sdev_file_sum2 = sdev_file_sum / (time_step - 1);
    //       sdev_file_sum2 = sqrt(sdev_file_sum2);

    MPI_Allreduce(&sdev_io, &sdev_io_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdev_io_sum1 = sdev_io_sum / time_step;
    //       sdev_io_sum1 = sqrt(sdev_io_sum1);

    sdev_io_sum2 = sdev_io_sum / (time_step - 1);
    //       sdev_io_sum2 = sqrt(sdev_io_sum2);


    if (rank == 0) {
        printf("Average Throughput acroos %d simulation runs = %f with a standard deviation of %f [%f] : %f [%f]\n", time_step, throughput_max, sdev_sum1, (sdev_sum1 * sdev_sum1), sdev_sum2, (sdev_sum2 * sdev_sum2));
        printf("Average Aggregation Time acroos %d simulation runs = %f with a standard deviation of %f [%f] : %f [%f]\n", time_step, agg_ts_max, sdev_agg_sum1, (sdev_agg_sum1 * sdev_agg_sum1), sdev_agg_sum2, (sdev_agg_sum2 * sdev_agg_sum2));
        printf("Average File Creation Time acroos %d simulation runs = %f with a standard deviation of %f [%f] : %f [%f]\n", time_step, file_ts_max, sdev_file_sum1, (sdev_file_sum1 * sdev_file_sum1), sdev_file_sum2, (sdev_file_sum2 * sdev_file_sum2));
        printf("Average IO Time acroos %d simulation runs = %f with a standard deviation of %f [%f] : %f [%f]\n", time_step, io_ts_max, sdev_io_sum1, (sdev_io_sum1 * sdev_io_sum1), sdev_io_sum2, (sdev_io_sum2 * sdev_io_sum2));
    }

    /*
    char file_name[512] = {0};
    FILE* fp;
    double *buffer;
    int counter = 0;
    size_t result;
    if(rank == 0)
    {
        for(ts = 0 ; ts <  1 ; ts++)
        {
            for(i = 0 ; i < nfiles ; i++)
            {
                sprintf(file_name,"%s/file_%d.bin", dir_path[ts], i);
                //printf("File Name %s\n", file_name);
                fp = fopen(file_name, "rb");
                if(!fp)
                {
                    //
                }
                buffer = (double*)malloc(sizeof(double) * 512 * 32768 * number_of_variables);
                assert(buffer);
                memset(buffer, 0 , (sizeof(double) * 512 * 32768 * number_of_variables));
		  
                result = fread(buffer, sizeof(double), 512 * 32768 * number_of_variables, fp);
                assert(result == 512 * 32768 * number_of_variables);
		  
                for( j = 0 ; j < 512 * 32768 * number_of_variables ; j++)
                {
                      if(buffer[j] != 17.85)
                      {
                        counter++;
                        printf("j = %d : %f\n", j, buffer[j]);
                      }
                }
                printf("Number of unequal elements in file %s [%d] are %d\n", file_name, i, counter);
                fclose(fp);
            }
        }
    }
     */
    MPI_Finalize();
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:f:t:v:e:r:";
    int one_opt = 0, i = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        switch (one_opt) {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
                break;
            case('f'):
                sprintf(dir_path_template, "%s", optarg);
                break;
            case('t'):
                sscanf(optarg, "%d", &time_step);
                break;
            case('v'):
                sscanf(optarg, "%d", &number_of_variables);
            case('e'):
                sscanf(optarg, "%d", &FACTOR);
            case('r'):
                sscanf(optarg, "%d", &sub_FACTOR);
                break;


            case('?'):
                return (-1);
        }
    }
    return (0);
}

static void usage(void) {
    printf("Usage: test-PIDX -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("  -f: IDX Filename\n");
    printf("  -t: number of timesteps\n");
    printf("\n");
    return;
} 
