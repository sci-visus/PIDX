/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */

#include <mpi.h>
#include <pnetcdf.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "PIDX.h"

static int parse_args(int argc, char **argv);
static void usage(void);
static int print_error(char *error_message, char* file, int line);

/* global dimensions of 3D volume */
static int extents[3] = {0, 0, 0};
static int extents1[4] = {0, 0, 0, 0};
static int extents2[4] = {0, 0, 0, 0};

/* per-process dimensions of each sub-block */
static int count[5] = {0, 0, 0, 0, 0};
static int offset[5] = {0, 0, 0, 0, 0};

static MPI_Offset count1[3] = {0, 0, 0};
static MPI_Offset count2[3] = {0, 0, 0};
static MPI_Offset count3[4] = {0, 0, 0, 0};
static MPI_Offset count4[4] = {0, 0, 0, 0};

static MPI_Offset offset1[3] = {0, 0, 0};
static MPI_Offset offset2[3] = {0, 0, 0};
static MPI_Offset offset3[4] = {0, 0, 0, 0};
static MPI_Offset offset4[4] = {0, 0, 0, 0};

/* Number of time-steps */
static int time_step = 0;

/* output IDX file Name Template*/
static char output_file_template[512] = {0};

static void handle_error(int status, int lineno) {
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {

    int restructuring_ACTIVE = 0;
    int write_ACTIVE = 1;

    int ret, ncfile, nprocs, rank, varid0 = 0, varid1 = 1, varid2 = 2, varid3 = 3;
    int dimid[4];

    //lfs setstripe --size 4M --index -1 --count 156 $SCRATCH2/READ/parallel_idx/source/PIDX-lib/PNET_R/

    int i = 0, j = 0, k = 0;
    int spv = 0; /*samples per variable*/
    int ts, vc; /*time step counter and variable counter*/
    int slice;
    int sub_div[4];
    char **output_file; /*Output File Name*/
    int gextent[5]; /*Global Extensions of the dataset (64 64 64 0 0)*/

    int number_of_variables = 4;
    int* sample_per_variable_buffer;
    double** write_data;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /*Rank 0 parses the command Line Arguments*/
    if (rank == 0) {
        ret = parse_args(argc, argv);
        if (ret < 0) {
            usage();
            print_error("ret error", __FILE__, __LINE__);
        }
        if (count[0] == 0 || count[1] == 0 || count[2] == 0) {
            usage();
            print_error("Local Dimension cannot be 0!!!!!!!!!\n", __FILE__, __LINE__);
        } else {
            if ((extents[0] / count[0]) * (extents[1] / count[1]) * (extents[2] / count[2]) != nprocs) {
                usage();
                print_error("Wrong Number of Processes\n", __FILE__, __LINE__);
            }
        }
    }
    count[3] = 1;
    count[4] = 1;
    offset[3] = 0;
    offset[4] = 0;

    /*   The command line arguments are shared by all processes  */
    MPI_Bcast(extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(count, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

    /*Creating the filename*/
    output_file = (char**) malloc(sizeof (char*) * time_step);
    for (i = 0; i < time_step; i++)
        output_file[i] = (char*) malloc(sizeof (char) * 512);

    /* Calculating every process's offset and count */
    gextent[0] = extents[0];
    gextent[1] = extents[1];
    gextent[2] = extents[2];
    gextent[3] = 1;
    gextent[4] = 1;

    sub_div[0] = (extents[0] / count[0]);
    sub_div[1] = (extents[1] / count[1]);
    sub_div[2] = (extents[2] / count[2]);
    offset[2] = (rank / (sub_div[0] * sub_div[1])) * count[2];
    slice = rank % (sub_div[0] * sub_div[1]);
    offset[1] = (slice / sub_div[0]) * count[1];
    offset[0] = (slice % sub_div[0]) * count[0];

    for(i = 0 ; i < 3 ; i++)
    {
        count1[i] = count[i];
        offset1[i] = offset[i];
        count2[i] = count[i];
        offset2[i] = offset[i];
        count3[i] = count[i];
        offset3[i] = offset[i];
        count4[i] = count[i];
        offset4[i] = offset[i];
        
        extents1[i] = extents[i];
        extents2[i] = extents[i];
    }
    offset3[3] = 0;
    count3[3] = 3;

    offset4[3] = 0;
    count4[3] = 11;

    extents1[3] = 3;
    extents2[3] = 11;
    
    assert(offset[0] < extents[0]);
    assert(offset[1] < extents[1]);
    assert(offset[2] < extents[2]);
    
    sample_per_variable_buffer = (int*) malloc(sizeof (int) * number_of_variables);
    if (!sample_per_variable_buffer)
        print_error("Error Allocating {sample_per_variable_buffer} buffer", __FILE__, __LINE__);

    sample_per_variable_buffer[0] = 1; //Variable 1
    sample_per_variable_buffer[1] = 1; //Variable 2
    sample_per_variable_buffer[2] = 3; //Variable 3
    sample_per_variable_buffer[3] = 11; //Variable 4

    write_data = (double**) malloc(sizeof (double*) * number_of_variables);
    if (!write_data)
        print_error("Error Allocating Buffer", __FILE__, __LINE__);
    memset(write_data, 0, sizeof (double*) * number_of_variables);


    double *start_time = 0, *end_time = 0, *var1_rst_start = 0, *var1_rst_end, *var2_rst_start, *var2_rst_end, *var3_rst_start, *var3_rst_end, *var4_rst_start, *var4_rst_end;
    double *create_start = 0, *create_end = 0, *var1_io_start, *var1_io_end, *var2_io_start, *var2_io_end, *var3_io_start, *var3_io_end, *var4_io_start, *var4_io_end;
    double *close_start = 0, *close_end = 0;
    start_time = (double*) malloc(sizeof (double) * time_step);
    end_time = (double*) malloc(sizeof (double) * time_step);


    if (write_ACTIVE == 1) 
    {
        
        create_start = (double*) malloc(sizeof (double) * time_step);
        create_end = (double*) malloc(sizeof (double) * time_step);
        close_start = (double*) malloc(sizeof (double) * time_step);
        close_end = (double*) malloc(sizeof (double) * time_step);
        var1_io_start = (double*) malloc(sizeof (double) * time_step);
        var1_io_end = (double*) malloc(sizeof (double) * time_step);
        var2_io_start = (double*) malloc(sizeof (double) * time_step);
        var2_io_end = (double*) malloc(sizeof (double) * time_step);
        var3_io_start = (double*) malloc(sizeof (double) * time_step);
        var3_io_end = (double*) malloc(sizeof (double) * time_step);
        var4_io_start = (double*) malloc(sizeof (double) * time_step);
        var4_io_end = (double*) malloc(sizeof (double) * time_step);
        memset(var1_io_start, 0, sizeof (double) * time_step);
        memset(var1_io_end, 0, sizeof (double) * time_step);
        memset(var2_io_start, 0, sizeof (double) * time_step);
        memset(var2_io_end, 0, sizeof (double) * time_step);
        memset(var3_io_start, 0, sizeof (double) * time_step);
        memset(var3_io_end, 0, sizeof (double) * time_step);
        memset(var4_io_start, 0, sizeof (double) * time_step);
        memset(var4_io_end, 0, sizeof (double) * time_step);
        
        for (vc = 0; vc < number_of_variables; vc++) 
        {
            write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);
            if (!write_data[vc])
                print_error("Error Allocating Buffer", __FILE__, __LINE__);

            for (k = 0; k < count[2]; k++)
                for (j = 0; j < count[1]; j++)
                    for (i = 0; i < count[0]; i++) 
                    {
                        int64_t index = (int64_t) (count[0] * count[1] * k) + (count[0] * j) + i;
                        for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++) 
                                write_data[vc][index * sample_per_variable_buffer[vc] + spv] = 100 + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i);
                    }
        }
        
        if (restructuring_ACTIVE == 0) 
        {
// 	    if(offset1[0]+count1[0] == gextent[0] || offset1[1]+count1[1] == gextent[1] || offset1[2]+count1[2] == gextent[2])
// 	    //if(rank < 32)
// 	    {
// 		printf("Edge Processes (Rank) = %d\n", rank);
// 		offset1[0] = 0;
// 		offset1[1] = 0;
// 		offset1[2] = 0;
// 		
// 		count1[0] = 1;
// 		count1[1] = 1;
// 		count1[2] = 1;
// 	
// 	    }
            for (ts = 0; ts < time_step; ts++) 
            {
                start_time[ts] = MPI_Wtime();
                sprintf(output_file[ts], "%s%s%d%s", output_file_template, "_", ts, ".nc");
                
                create_start[ts] = MPI_Wtime();
                ret = ncmpi_create(MPI_COMM_WORLD, output_file[ts], NC_CLOBBER | NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

		
                ret = ncmpi_def_dim(ncfile, "nx", extents[0], &dimid[0]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "ny", extents[1], &dimid[1]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "nz", extents[2], &dimid[2]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		
                ret = ncmpi_def_var(ncfile, "temp", NC_DOUBLE, 3, dimid, &varid0);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
		

		/*
                ret = ncmpi_def_var(ncfile, "pressure", NC_DOUBLE, 3, dimid, &varid1);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "number_of_velocity_components", extents1[3], &dimid[3]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "velocity", NC_DOUBLE, 4, dimid, &varid2);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "number_of_species", extents2[3], &dimid[3]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "yspecies", NC_DOUBLE, 4, dimid, &varid3);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                */

                ret = ncmpi_enddef(ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                create_end[ts] = MPI_Wtime();

                var1_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid0, offset1, count1, write_data[0]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		var1_io_end[ts] = MPI_Wtime();
                
		/*
                var2_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid1, offset2, count2, write_data[1]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                var2_io_end[ts] = MPI_Wtime();

                var3_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid2, offset3, count3, write_data[2]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                var3_io_end[ts] = MPI_Wtime();

                var4_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid3, offset4, count4, write_data[3]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                var4_io_end[ts] = MPI_Wtime();
		 */
                
                close_start[ts] = MPI_Wtime();
                ret = ncmpi_close(ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                close_end[ts] = MPI_Wtime();
                end_time[ts] = MPI_Wtime();
            }
        } 
        else 
        {

            int *expanded_box_dimension;
            double* expanded_box = NULL;
            PIDX_rst_id rst_id;

            PIDX_Ndim_buffer*** rst_output_buffer;
            PIDX_Ndim_buffer** in_buf;
            int buffer_size, send_c = 0, send_o = 0;
            int rst_output_buffer_count = 0, var = 0;
            int dimension = 3;
            int k1, j1, i1, r, index = 0, recv_o = 0;
            
            expanded_box_dimension = (int*)malloc(sizeof(int) * 5);
            var1_rst_start = (double*) malloc(sizeof (double) * time_step);
            var1_rst_end = (double*) malloc(sizeof (double) * time_step);
            var2_rst_start = (double*) malloc(sizeof (double) * time_step);
            var2_rst_end = (double*) malloc(sizeof (double) * time_step);
            var3_rst_start = (double*) malloc(sizeof (double) * time_step);
            var3_rst_end = (double*) malloc(sizeof (double) * time_step);
            var4_rst_start = (double*) malloc(sizeof (double) * time_step);
            var4_rst_end = (double*) malloc(sizeof (double) * time_step);
            
            memset(var1_rst_start, 0, sizeof (double) * time_step);
            memset(var1_rst_end, 0, sizeof (double) * time_step);
            memset(var2_rst_start, 0, sizeof (double) * time_step);
            memset(var2_rst_end, 0, sizeof (double) * time_step);
            memset(var3_rst_start, 0, sizeof (double) * time_step);
            memset(var3_rst_end, 0, sizeof (double) * time_step);
            memset(var4_rst_start, 0, sizeof (double) * time_step);
            memset(var4_rst_end, 0, sizeof (double) * time_step);
            
            for (ts = 0; ts < time_step; ts++) 
            {    
                start_time[ts] = MPI_Wtime();
                sprintf(output_file[ts], "%s%s%d%s", output_file_template, "_", ts, ".nc");
              
                create_start[ts] = MPI_Wtime();
                ret = ncmpi_create(MPI_COMM_WORLD, output_file[ts], NC_CLOBBER | NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);


                ret = ncmpi_def_dim(ncfile, "nx", extents[0], &dimid[0]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "ny", extents[1], &dimid[1]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "nz", extents[2], &dimid[2]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "temp", NC_DOUBLE, 3, dimid, &varid0);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "pressure", NC_DOUBLE, 3, dimid, &varid1);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "number_of_velocity_components", extents1[3], &dimid[3]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "velocity", NC_DOUBLE, 4, dimid, &varid2);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_dim(ncfile, "number_of_species", extents2[3], &dimid[3]);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_def_var(ncfile, "yspecies", NC_DOUBLE, 4, dimid, &varid3);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                ret = ncmpi_enddef(ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);


                in_buf = (PIDX_Ndim_buffer**) malloc(number_of_variables * sizeof (PIDX_Ndim_buffer*));
                for (i = 0; i < number_of_variables; i++)
                    in_buf[i] = (PIDX_Ndim_buffer*) malloc(sizeof (PIDX_Ndim_buffer));

                rst_output_buffer = (PIDX_Ndim_buffer***) malloc(sizeof (*rst_output_buffer) * number_of_variables);
                memset(rst_output_buffer, 0, sizeof (*rst_output_buffer) * number_of_variables);
                rst_output_buffer_count = 0;


                rst_id = PIDX_rst_init(MPI_COMM_WORLD, dimension, (int*) gextent, count, offset, 0, NULL, &rst_output_buffer_count);
                expanded_box_dimension = PIDX_rst_get_box_dimension(rst_id);
                
                buffer_size = 0;
                buffer_size = expanded_box_dimension[0] * expanded_box_dimension[1] * expanded_box_dimension[2];
                
                if(rank == 0)
                    printf("[RST] : Time Step %d Restructured Box Size : %d %d %d\n", ts, expanded_box_dimension[0], expanded_box_dimension[1], expanded_box_dimension[2]);
                create_end[ts] = MPI_Wtime();

                var1_rst_start[ts] = MPI_Wtime();
                var = 0;
                PIDX_rst_buf_init(in_buf[var], 3, offset, count, write_data[var], sample_per_variable_buffer[var], MPI_DOUBLE, "var_name", var);
                rst_output_buffer[var] = (PIDX_Ndim_buffer**) malloc((rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                memset(rst_output_buffer[var], 0, (rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                PIDX_rst_restructure(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count));
                PIDX_rst_restructure_IO(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count), 1);
                
                if (rst_output_buffer_count != 0) 
                {
                    index = 0, recv_o = 0;
                    expanded_box = (double*) malloc(sizeof (double) * buffer_size * sample_per_variable_buffer[var]);
                    memset(expanded_box, 0, (sizeof (double) * buffer_size * sample_per_variable_buffer[var]));

                    for (r = 0; r < rst_output_buffer_count; r++) 
                        for (k1 = rst_output_buffer[var][r]->lower_bounds[2]; k1 < rst_output_buffer[var][r]->upper_bounds[2]; k1++)
                            for (j1 = rst_output_buffer[var][r]->lower_bounds[1]; j1 < rst_output_buffer[var][r]->upper_bounds[1]; j1++)
                                for (i1 = rst_output_buffer[var][r]->lower_bounds[0]; i1 < rst_output_buffer[var][r]->upper_bounds[0]; i1 = i1 + rst_output_buffer[var][r]->upper_bounds[0] - rst_output_buffer[var][r]->lower_bounds[0]) 
                                {
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
                if (rst_output_buffer_count != 0) 
                {
                    for(i = 0 ; i < 3 ; i++)
                    {
                        offset1[i] = rst_output_buffer[0][0]->regular_lower_bounds[i];
                        count1[i] = rst_output_buffer[0][0]->regular_upper_bounds[i] - rst_output_buffer[0][0]->regular_lower_bounds[i];
                    }
                } 
                else 
                {
                    for(i = 0 ; i < 3 ; i++)
                    {
                        offset1[i] = 0;
                        count1[i] = 0;
                    }
                }
                var1_rst_end[ts] = MPI_Wtime();

                var1_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid0, offset1, count1, expanded_box);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);


                ret = PIDX_rst_buf_destroy(in_buf[var]);
                for (j = 0; j < rst_output_buffer_count; j++) 
                {
                    ret = PIDX_rst_buf_destroy(rst_output_buffer[var][j]);
                    if (ret == -1) 
                    {
                        fprintf(stderr, "PIDX : [%d] Error in PIDX_rst_buf_destroy\n", rank);
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                }
                if (rst_output_buffer_count != 0) 
                {
                    free(rst_output_buffer[var]);
                    rst_output_buffer[var] = 0;
                }
                if (rst_output_buffer_count != 0)
                    free(expanded_box);
                var1_io_end[ts] = MPI_Wtime();

                var2_rst_start[ts] = MPI_Wtime();
                var = 1;
                PIDX_rst_buf_init(in_buf[var], 3, offset, count, write_data[var], sample_per_variable_buffer[var], MPI_DOUBLE, "var_name", var);
                rst_output_buffer[var] = (PIDX_Ndim_buffer**) malloc((rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                memset(rst_output_buffer[var], 0, (rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                PIDX_rst_restructure(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count));
                PIDX_rst_restructure_IO(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count), 1);
                //HELPER_rst(rst_output_buffer[var], rst_id, rst_output_buffer_count, sample_per_variable_buffer[var]);



                if (rst_output_buffer_count != 0) {
                    
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


                if (rst_output_buffer_count != 0) {
                    offset2[0] = rst_output_buffer[var][0]->regular_lower_bounds[0];
                    offset2[1] = rst_output_buffer[var][0]->regular_lower_bounds[1];
                    offset2[2] = rst_output_buffer[var][0]->regular_lower_bounds[2];

                    count2[0] = rst_output_buffer[var][0]->regular_upper_bounds[0] - rst_output_buffer[var][0]->regular_lower_bounds[0];
                    count2[1] = rst_output_buffer[var][0]->regular_upper_bounds[1] - rst_output_buffer[var][0]->regular_lower_bounds[1];
                    count2[2] = rst_output_buffer[var][0]->regular_upper_bounds[2] - rst_output_buffer[var][0]->regular_lower_bounds[2];
                } else {
                    offset2[0] = 0;
                    offset2[1] = 0;
                    offset2[2] = 0;

                    count2[0] = 0;
                    count2[1] = 0;
                    count2[2] = 0;
                }
                //    if (rst_output_buffer_count != 0)
                //	printf("[%d] [%d] [E%d] OFFSET :: COUNT %d %d %d : %d %d %d\n", rank, rst_output_buffer_count, expanded_box_dimension, offset2[0], offset2[1], offset2[2], count2[0], count2[1], count2[2]);

                var2_rst_end[ts] = MPI_Wtime();

                var2_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid1, offset2, count2, expanded_box);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

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
                if (rst_output_buffer_count != 0)
                    free(expanded_box);
                var2_io_end[ts] = MPI_Wtime();

                var3_rst_start[ts] = MPI_Wtime();
                var = 2;
                PIDX_rst_buf_init(in_buf[var], 3, offset, count, write_data[var], sample_per_variable_buffer[var], MPI_DOUBLE, "var_name", var);
                rst_output_buffer[var] = (PIDX_Ndim_buffer**) malloc((rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                memset(rst_output_buffer[var], 0, (rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                PIDX_rst_restructure(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count));
                PIDX_rst_restructure_IO(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count), 1);
                //HELPER_rst(rst_output_buffer[var], rst_id, rst_output_buffer_count, sample_per_variable_buffer[var]);



                if (rst_output_buffer_count != 0) {
                    

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
                if (rst_output_buffer_count != 0) {
                    offset3[0] = rst_output_buffer[var][0]->regular_lower_bounds[0];
                    offset3[1] = rst_output_buffer[var][0]->regular_lower_bounds[1];
                    offset3[2] = rst_output_buffer[var][0]->regular_lower_bounds[2];
                    offset3[3] = 0;

                    count3[0] = rst_output_buffer[var][0]->regular_upper_bounds[0] - rst_output_buffer[var][0]->regular_lower_bounds[0];
                    count3[1] = rst_output_buffer[var][0]->regular_upper_bounds[1] - rst_output_buffer[var][0]->regular_lower_bounds[1];
                    count3[2] = rst_output_buffer[var][0]->regular_upper_bounds[2] - rst_output_buffer[var][0]->regular_lower_bounds[2];
                    count3[3] = 3;
                } else {
                    offset3[0] = 0;
                    offset3[1] = 0;
                    offset3[2] = 0;
                    offset3[3] = 0;

                    count3[0] = 0;
                    count3[1] = 0;
                    count3[2] = 0;
                    count3[3] = 0;
                }
                var3_rst_end[ts] = MPI_Wtime();
                //    if (rst_output_buffer_count != 0)
                //	printf("[%d] [%d] [E%d] OFFSET :: COUNT %d %d %d : %d %d %d\n", rank, rst_output_buffer_count, expanded_box_dimension, offset3[0], offset3[1], offset3[2], count3[0], count3[1], count3[2]);

                var3_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid2, offset3, count3, expanded_box);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

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
                if (rst_output_buffer_count != 0)
                    free(expanded_box);
                var3_io_end[ts] = MPI_Wtime();

                var4_rst_start[ts] = MPI_Wtime();
                var = 3;
                PIDX_rst_buf_init(in_buf[var], 3, offset, count, write_data[var], sample_per_variable_buffer[var], MPI_DOUBLE, "var_name", var);
                rst_output_buffer[var] = (PIDX_Ndim_buffer**) malloc((rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                memset(rst_output_buffer[var], 0, (rst_output_buffer_count) * sizeof (PIDX_Ndim_buffer*));
                PIDX_rst_restructure(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count));
                PIDX_rst_restructure_IO(rst_id, in_buf[var], rst_output_buffer[var], (rst_output_buffer_count), 1);
                //HELPER_rst(rst_output_buffer[var], rst_id, rst_output_buffer_count, sample_per_variable_buffer[var]);



                if (rst_output_buffer_count != 0) {
                    

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

                if (rst_output_buffer_count != 0) {
                    offset4[0] = rst_output_buffer[var][0]->regular_lower_bounds[0];
                    offset4[1] = rst_output_buffer[var][0]->regular_lower_bounds[1];
                    offset4[2] = rst_output_buffer[var][0]->regular_lower_bounds[2];
                    offset4[3] = 0;

                    count4[0] = rst_output_buffer[var][0]->regular_upper_bounds[0] - rst_output_buffer[var][0]->regular_lower_bounds[0];
                    count4[1] = rst_output_buffer[var][0]->regular_upper_bounds[1] - rst_output_buffer[var][0]->regular_lower_bounds[1];
                    count4[2] = rst_output_buffer[var][0]->regular_upper_bounds[2] - rst_output_buffer[var][0]->regular_lower_bounds[2];
                    count4[3] = 11;
                } else {
                    offset4[0] = 0;
                    offset4[1] = 0;
                    offset4[2] = 0;
                    offset4[3] = 0;

                    count4[0] = 0;
                    count4[1] = 0;
                    count4[2] = 0;
                    count4[3] = 0;
                }
                var4_rst_end[ts] = MPI_Wtime();

                var4_io_start[ts] = MPI_Wtime();
                ret = ncmpi_put_vara_double_all(ncfile, varid3, offset4, count4, expanded_box);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);

                // if (rst_output_buffer_count != 0)
                //		printf("[%d] [%d] [E%d] OFFSET :: COUNT %d %d %d : %d %d %d\n", rank, rst_output_buffer_count, expanded_box_dimension, offset4[0], offset4[1], offset4[2], count4[0], count4[1], count4[2]);      

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
                
                if (rst_output_buffer_count != 0)
                    free(expanded_box);
                var4_io_end[ts] = MPI_Wtime();
                
                close_start[ts] = MPI_Wtime();

                free(rst_output_buffer);
                rst_output_buffer = 0;

                ret = ncmpi_close(ncfile);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                close_end[ts] = MPI_Wtime();

                end_time[ts] = MPI_Wtime();
            }

        }
        for (vc = 0; vc < number_of_variables; vc++) {
            free(write_data[vc]);
            write_data[vc] = 0;
        }
        free(write_data);
        write_data = 0;


    } else {
        write_data = (double**) malloc(sizeof (double*) * number_of_variables);
        if (!write_data)
            print_error("Error Allocating Buffer", __FILE__, __LINE__);
        memset(write_data, 0, sizeof (double*) * number_of_variables);

        for (vc = 0; vc < number_of_variables; vc++) {
            write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);
            if (!write_data[vc])
                print_error("Error Allocating Buffer", __FILE__, __LINE__);
            memset(write_data[vc], 0, sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);
        }
        //printf("[%d] : %d\n", rank, time_step);
        for (ts = 0; ts < time_step; ts++) {

            start_time[ts] = MPI_Wtime();
            for (vc = 0; vc < number_of_variables; vc++)
                memset(write_data[vc], 0, sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);

            sprintf(output_file[ts], "%s%s%d%s", output_file_template, "_", ts, ".nc");
            //printf("[%d] Filename for Time-step %d : O : C :: %d %d %d :: %d %d %d = %s\n", rank,  ts, offset[0], offset[1], offset[2], count[0], count[1], count[2], output_file[ts]);
            ret = ncmpi_open(MPI_COMM_WORLD, output_file[ts], NC_NOWRITE, MPI_INFO_NULL, &ncfile);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);

            ret = ncmpi_get_vara_double_all(ncfile, varid0, offset1, count1, write_data[0]);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);

            ret = ncmpi_get_vara_double_all(ncfile, varid1, offset2, count2, write_data[1]);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);


            //count[3] = 3;
            ret = ncmpi_get_vara_double_all(ncfile, varid2, offset3, count3, write_data[2]);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);

            //count[3] = 11;
            ret = ncmpi_get_vara_double_all(ncfile, varid3, offset4, count4, write_data[3]);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);

            ret = ncmpi_close(ncfile);
            if (ret != NC_NOERR) handle_error(ret, __LINE__);

            int lost_element = 0;
            for (vc = 0; vc < number_of_variables; vc++) {

                /*
                for (k = 0; k < count[2]; k++)
            for (j = 0; j < count[1]; j++)
                for (i = 0; i < count[0]; i++) {
                    int64_t index = (int64_t) (count[0] * count[1] * k) + (count[0] * j) + i;
                    for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++){
                        //write_data[vc][index * sample_per_variable_buffer[vc] + spv] = 100 + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i);
                        //if(vc == 0)
                        //printf("[%d] [%d] [%d] Data at %d %d %d = %d\n", rank, vc, ts, i + offset[0], j + offset[1], k + offset[2], (int)write_data[vc][index * sample_per_variable_buffer[vc] + spv]);
                    }
                }
                 */

                //if(rank == 1)
                for (k = 0; k < count[2]; k++)
                    for (j = 0; j < count[1]; j++)
                        for (i = 0; i < count[0]; i++) {
                            int64_t index = (int64_t) (count[0] * count[1] * k) + (count[0] * j) + i;
                            for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++) {
                                if ((int) write_data[vc][index * sample_per_variable_buffer[vc] + spv] != 100 + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i)) {
                                    printf("[%d] SCREAM!!!!!!!! : %d : %d\n", vc, (int) write_data[vc][index * sample_per_variable_buffer[vc] + spv], (100 + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i)));
                                    lost_element++;
                                } else {
                                    //printf("[%d] VALUE : %d\n", vc, (int)write_data[vc][index * sample_per_variable_buffer[vc] + spv]);
                                }
                            }
                        }


                assert(lost_element == 0);
            }

            if (rank == 0)
                printf("Done verifying for Time-step %d\n", ts);
            end_time[ts] = MPI_Wtime();
        }
        for (vc = 0; vc < number_of_variables; vc++) {
            free(write_data[vc]);
            write_data[vc] = 0;
        }
        free(write_data);
        write_data = 0;
    }
    double max_time = 0, total_time = 0;
    int64_t total_data = 0;
    int var = 0, sample_sum = 0;
    for (var = 0; var < number_of_variables; var++) {
        sample_sum = sample_sum + sample_per_variable_buffer[var];
    }
    total_data = (gextent[0] * gextent[1] * gextent[2] * sample_sum * sizeof (double));
    for (ts = 0; ts < time_step; ts++) 
    {
        max_time = 0, total_time = 0;
        total_time = end_time[ts] - start_time[ts];
        MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (total_time == max_time) 
        {
            if(restructuring_ACTIVE == 0)
                printf("[Rank %d Nproc %d] [Time Step %d] [Total Time : %2.4f [%2.4f = %2.4f] = %2.4f + [%2.4f]  + [%2.4f]  + [%2.4f]  + [%2.4f] + %2.4f] [Total Volume %d x %d x %d x %d] [Throughput %2.4f MiB/sec]\n", rank, nprocs, ts, max_time, 
                        ((create_end[ts] - create_start[ts]) + (var1_io_end[ts] - var1_io_start[ts]) + (var2_io_end[ts] - var2_io_start[ts]) + (var3_io_end[ts] - var3_io_start[ts]) + (var4_io_end[ts] - var4_io_start[ts]) + (close_end[ts] - close_start[ts])),
                        ((var1_io_end[ts] - var1_io_start[ts]) + (var2_io_end[ts] - var2_io_start[ts]) + (var3_io_end[ts] - var3_io_start[ts]) + (var4_io_end[ts] - var4_io_start[ts])),
                        (create_end[ts] - create_start[ts]), 
                        (var1_io_end[ts] - var1_io_start[ts]), 
                        (var2_io_end[ts] - var2_io_start[ts]), 
                        (var3_io_end[ts] - var3_io_start[ts]), 
                        (var4_io_end[ts] - var4_io_start[ts]), 
                        (close_end[ts] - close_start[ts]), 
                        gextent[0], gextent[1], gextent[2], sample_sum, (float) ((int64_t)total_data / (1024 * 1024 * max_time)));
            else
                printf("[R] [Rank %d Nproc %d] [Time Step %d] [Total Time : %2.4f [%2.4f = %2.4f + %2.4f] = %2.4f + [%2.4f + %2.4f]  + [%2.4f  + %2.4f]  + [%2.4f + %2.4f]  + [%2.4f  + %2.4f] + %2.4f] [Total Volume %d x %d x %d x %d] [Throughput %2.4f MiB/sec]\n", 
                        rank, nprocs, ts, max_time, 
                        ((create_end[ts] - create_start[ts]) + (var1_rst_end[ts] - var1_rst_start[ts]) + (var1_io_end[ts] - var1_io_start[ts]) + (var2_rst_end[ts] - var2_rst_start[ts]) + (var2_io_end[ts] - var2_io_start[ts]) + (var3_rst_end[ts] - var3_rst_start[ts]) + (var3_io_end[ts] - var3_io_start[ts]) + (var4_rst_end[ts] - var4_rst_start[ts]) + (var4_io_end[ts] - var4_io_start[ts]) + (close_end[ts] - close_start[ts])),
                        ((var1_rst_end[ts] - var1_rst_start[ts]) + (var2_rst_end[ts] - var2_rst_start[ts]) + (var3_rst_end[ts] - var3_rst_start[ts]) + (var4_rst_end[ts] - var4_rst_start[ts])),
                        ((var1_io_end[ts] - var1_io_start[ts]) + (var2_io_end[ts] - var2_io_start[ts]) + (var3_io_end[ts] - var3_io_start[ts]) + (var4_io_end[ts] - var4_io_start[ts])),
                        (create_end[ts] - create_start[ts]), 
                        (var1_rst_end[ts] - var1_rst_start[ts]), (var1_io_end[ts] - var1_io_start[ts]), 
                        (var2_rst_end[ts] - var2_rst_start[ts]), (var2_io_end[ts] - var2_io_start[ts]), 
                        (var3_rst_end[ts] - var3_rst_start[ts]), (var3_io_end[ts] - var3_io_start[ts]), 
                        (var4_rst_end[ts] - var4_rst_start[ts]), (var4_io_end[ts] - var4_io_start[ts]), 
                        (close_end[ts] - close_start[ts]), 
                        gextent[0], gextent[1], gextent[2], sample_sum, (float) ((int64_t)total_data / (1024 * 1024 * max_time)));
        }
    }



    MPI_Finalize();
//0.0161 + 0.0001 + 0.0083  + 0.0001  + 0.0071  + 0.0011 + 0.0142  + 0.0022  + 0.0471
    return 0;
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:l:f:t:";
    int one_opt = 0, i = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        /* postpone error checking for after while loop */
        switch (one_opt) {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
                break;
            case('l'):
                sscanf(optarg, "%dx%dx%d", &count[0], &count[1], &count[2]);
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
    if (extents[0] < 1 || extents[1] < 1 || extents[2] < 1 || count[0] < 1 || count[1] < 1 || count[2] < 1) {
        printf("Error: bad dimension specification.\n");
        return (-1);
    }

    /* need global dimension to be larger than the local */
    if (extents[0] < count[0] || extents[1] < count[1] || extents[2] < count[2]) {
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

static int print_error(char *error_message, char* file, int line) {
    fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(0);
#endif
}
