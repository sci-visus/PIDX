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

static int parse_args(int argc, char **argv);
static void usage(void);
static int print_error(char *error_message, char* file, int line);

/* global dimensions of 3D volume */
static int extents[3] = {0, 0, 0};


/* per-process dimensions of each sub-block */
static int count[3] = {0, 0, 0};
static int offset[3] = {0, 0, 0};

/* Number of time-steps */
static int time_step = 0;

/* output IDX file Name Template*/
static char output_file_template[512] = {0};


static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {

    int ret, ncfile, nprocs, rank, varid0 = 0, varid1 = 1, varid2 = 2, varid3 = 3;
    int dimid[4];
    
    MPI_Offset offset1[3] = {0, 0, 0};
    MPI_Offset offset2[3] = {0, 0, 0};
    MPI_Offset offset3[4] = {0, 0, 0, 0};
    MPI_Offset offset4[4] = {0, 0, 0, 0};
       
    MPI_Offset count1[3] = {0, 0, 0};
    MPI_Offset count2[3] = {0, 0, 0};
    MPI_Offset count3[4] = {0, 0, 0, 0};
    MPI_Offset count4[4] = {0, 0, 0, 0};
    
    int i = 0, j = 0, k = 0;
    int spv = 0; /*samples per variable*/
    int ts, vc; /*time step counter and variable counter*/
    int slice;
    int sub_div[4];
    char **output_file; /*Output File Name*/
    
    
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
        }
        else {
            if ((extents[0] / count[0]) * (extents[1] / count[1]) * (extents[2] / count[2]) != nprocs) {
                usage();
                print_error("Wrong Number of Processes\n", __FILE__, __LINE__);
            }
        }
    }
    

    /*   The command line arguments are shared by all processes  */
    MPI_Bcast(extents, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(count, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);


    /*Creating the filename*/
    output_file = (char**) malloc(sizeof (char*) * time_step);
    for(i = 0 ; i < time_step ; i++)
        output_file[i] = (char*) malloc(sizeof (char) * 512);
    
    



    
    sub_div[0] = (extents[0] / count[0]);
    sub_div[1] = (extents[1] / count[1]);
    sub_div[2] = (extents[2] / count[2]);
    offset[2] = (rank / (sub_div[0] * sub_div[1])) * count[2];
    slice = rank % (sub_div[0] * sub_div[1]);
    offset[1] = (slice / sub_div[0]) * count[1];
    offset[0] = (slice % sub_div[0]) * count[0];

    
    offset1[0] = offset[0];
    offset1[1] = offset[1];
    offset1[2] = offset[2];
    count1[0] = count[0];
    count1[1] = count[1];
    count1[2] = count[2];
    
    offset2[0] = offset[0];
    offset2[1] = offset[1];
    offset2[2] = offset[2];
    count2[0] = count[0];
    count2[1] = count[1];
    count2[2] = count[2];

    offset3[0] = offset[0];
    offset3[1] = offset[1];
    offset3[2] = offset[2];
    offset3[3] = 0;
    count3[0] = count[0];
    count3[1] = count[1];
    count3[2] = count[2];
    count3[3] = 3;
    
    offset4[0] = offset[0];
    offset4[1] = offset[1];
    offset4[2] = offset[2];
    offset4[3] = 0;
    count4[0] = count[0];
    count4[1] = count[1];
    count4[2] = count[2];
    count4[3] = 11;
    
    //offset = offset_local;
    //count = count_local;
    
    assert(offset[0] < extents[0]);
    assert(offset[1] < extents[1]);
    assert(offset[2] < extents[2]);
    //assert(offset[3] < extents[3]);
    

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
    sample_per_variable_buffer[1] = 1; //Variable 2
    sample_per_variable_buffer[2] = 3; //Variable 3
    sample_per_variable_buffer[3] = 11; //Variable 4

    //printf("Number of time steps %d\n", time_step);

    write_data = (double**) malloc(sizeof (double*) * number_of_variables);
    if (!write_data)
        print_error("Error Allocating Buffer", __FILE__, __LINE__);
    memset(write_data, 0, sizeof (double*) * number_of_variables);

    for (vc = 0; vc < number_of_variables; vc++) {
        write_data[vc] = (double*) malloc(sizeof (double) * count[0] * count[1] * count[2] * sample_per_variable_buffer[vc]);
        if (!write_data[vc])
            print_error("Error Allocating Buffer", __FILE__, __LINE__);


    }
    
    for (ts = 0; ts < time_step; ts++) {
        
        
        sprintf(output_file[ts], "%s%s%d%s", output_file_template, "_", ts, ".nc");
        //printf("Filename for Time-step %d = %s\n", ts, output_file[ts]);
        ret = ncmpi_open(MPI_COMM_WORLD, output_file[ts], NC_NOWRITE, MPI_INFO_NULL, &ncfile);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        /*
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
        
        extents[3] = 3;
        ret = ncmpi_def_dim(ncfile, "number_of_velocity_components", extents[3], &dimid[3]);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        
        ret = ncmpi_def_var(ncfile, "velocity", NC_DOUBLE, 4, dimid, &varid2);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        
        extents[3] = 11;
        ret = ncmpi_def_dim(ncfile, "number_of_species", extents[3], &dimid[3]);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        
        ret = ncmpi_def_var(ncfile, "yspecies", NC_DOUBLE, 4, dimid, &varid3);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
      
        ret = ncmpi_enddef(ncfile);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        */ 
        
        
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
                            
        for (k = 0; k < count[2]; k++)
            for (j = 0; j < count[1]; j++)
                for (i = 0; i < count[0]; i++) {
                    long long index = (long long) (count[0] * count[1] * k) + (count[0] * j) + i;
                    for (spv = 0; spv < sample_per_variable_buffer[vc]; spv++)
                    {
                        if((int)write_data[vc][index * sample_per_variable_buffer[vc] + spv] != 100 + spv + (extents[0] * extents[1]*(offset[2] + k))+(extents[0]*(offset[1] + j)) + (offset[0] + i))
                        {
                            printf("[%d] SCREAM!!!!!!!! : %d\n", vc, (int)write_data[vc][index * sample_per_variable_buffer[vc] + spv]);
                            lost_element++;
                        }
                        else
                        {
                            //printf("[%d] [%d] VALUE : %d\n", ts, vc, (int)write_data[vc][index * sample_per_variable_buffer[vc] + spv]);
                        }
                            
                    }
                }
        }       
        assert(lost_element == 0);
        if(rank == 0)
                printf("Done verifying for Time-step %d\n", ts);
        
        
         
    }
    for (vc = 0; vc < number_of_variables; vc++) {
        free(write_data[vc]);
        write_data[vc] = 0;
    }
    free(write_data);
    write_data = 0;
    MPI_Finalize();

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
