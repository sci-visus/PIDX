/*
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
//#include <unistd.h>
//#include <stdarg.h>
//#include <stdint.h>
#include <string>
#include <PIDX.h>

#define VIS_DUMP // << ENABLE Visualization dump with downcasting


static int parse_args(int argc, char **argv);
static void usage(void);
static void report_error(std::string func_name, std::string file_name, int line_no);

static const int bits_per_block = 16;
static const int blocks_per_file = 256;
static int global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static int local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static PIDX_point global_size, local_offset, local_size;
static int local_box_offset[3];
static int restart_time = -1;
static int nprocs = 1, rank = 0;

int init(int argc, char **argv){

    int ret;

    // MPI initialization
#if PIDX_HAVE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0)
    {
        ret = parse_args(argc, argv);
        if (ret < 0)
        {
            usage();
#if PIDX_HAVE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#else
            exit(-1);
#endif
        }

        // check if the num procs is appropriate
        int num_bricks = (global_box_size[0] / local_box_size[0]) * (global_box_size[1] / local_box_size[1]) * (global_box_size[2] / local_box_size[2]);
        if(num_bricks != nprocs)
        {
            fprintf(stderr, "Error: number of sub-blocks (%d) doesn't match number of procs (%d)\n", num_bricks, nprocs);
            fprintf(stderr, "Incorrect distribution of data across processes i.e.\n(global_x / local_x) X (global_x / local_x) X (global_x / local_x) != nprocs\n(%d/%d) X (%d/%d) X (%d/%d) != %d\n", global_box_size[0], local_box_size[0], global_box_size[1], local_box_size[1], global_box_size[2], local_box_size[2], nprocs);

#if PIDX_HAVE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#else
            exit(-1);
#endif
        }
    }

#if PIDX_HAVE_MPI
    MPI_Bcast(global_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(local_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&restart_time, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    // Calculating every process data's offset and size
    int sub_div[3];
    sub_div[0] = (global_box_size[0] / local_box_size[0]);
    sub_div[1] = (global_box_size[1] / local_box_size[1]);
    sub_div[2] = (global_box_size[2] / local_box_size[2]);
    local_box_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
    int slice = rank % (sub_div[0] * sub_div[1]);
    local_box_offset[1] = (slice / sub_div[0]) * local_box_size[1];
    local_box_offset[0] = (slice % sub_div[0]) * local_box_size[0];

    PIDX_set_point_5D(global_size, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);
    PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);
    PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);

    return 0;

}

// Create a 3D sin function as example of synthetic simulation data (change phase over time)
template<typename T>
int create_data(T** data, unsigned variable_count, int* vps1, int ts){
    double fsize = 16.0;
    double fratiox = (fsize)/global_box_size[0];
    double fratioy = (fsize)/global_box_size[1];
    double fratioz = (fsize)/global_box_size[2];

    for(unsigned var = 0; var < variable_count; var++)
    {
        fsize = fsize * (var+1); // change freq for different variables

        data[var] = (double*)malloc(sizeof (T) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * vps1[var]);

        for (int k = 0; k < local_box_size[2]; k++)
        {
            for (int j = 0; j < local_box_size[1]; j++)
            {
                for (int i = 0; i < local_box_size[0]; i++)
                {
                    unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
                    for (int vps = 0; vps < vps1[var]; vps++){
                        double xc = (i + local_box_offset[0])*fratiox - fsize/2;
                        double yc = (j + local_box_offset[1])*fratioy - fsize/2;
                        double zc = (k + local_box_offset[2])*fratioz;
                        double zv = sin((xc*xc + yc*yc + ts%180)*3.14159265359/180)* fsize/2 + fsize/2;

                        if(zc > zv-0.7 && zc < zv+0.7)
                            data[var][index * vps1[var] + vps] = zv*global_box_size[2];
                        else
                            data[var][index * vps1[var] + vps] = 0;
                    }
                }
            }
        }
    }

    return 0;

}

// Read the data in the current IDX file at a defined time step,
// useful to check the information about the last checkpoint and then restart the simulation
int check_last_state(std::string file_name, int ts){
    //  Creating access
    PIDX_access access;
    PIDX_create_access(&access);
#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif
    int variable_count, ret;

    PIDX_file file;            // IDX file descriptor
    PIDX_variable* variable;   // variable descriptor
    PIDX_point dims;

    ret = PIDX_file_open(file_name.c_str(), PIDX_MODE_RDONLY, access, dims, &file);
    if (ret != PIDX_success)  report_error("PIDX_file_open", __FILE__, __LINE__);

    ret = PIDX_get_variable_count(file, &variable_count);
    if (ret != PIDX_success)  report_error("PIDX_set_variable_count", __FILE__, __LINE__);

    if(rank == 0){
        printf("Restart reading last state:\n\tvariable count: %d\n", variable_count);
    }

    ret = PIDX_set_current_time_step(file, ts);
    if (ret != PIDX_success)  report_error("PIDX_get_current_time_step", __FILE__, __LINE__);

    variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);

    int* vps = (int*)malloc(sizeof(*vps) * variable_count);
    memset(vps, 0, sizeof(*vps) * variable_count);

    double** data = (double**)malloc(sizeof(*data) * variable_count);
    memset(data, 0, sizeof(*data) * variable_count);

    for (int var = 0; var < variable_count; var++)
    {
        ret = PIDX_get_next_variable(file, &variable[var]);
        if (ret != PIDX_success) report_error("PIDX_get_next_variable", __FILE__, __LINE__);

        vps[var] = variable[var]->vps;

        int bits_per_sample = 0;
        ret = PIDX_default_bits_per_datatype(variable[var]->type_name, &bits_per_sample);
        if (ret != PIDX_success)  report_error("PIDX_default_bytes_per_datatype", __FILE__, __LINE__);

        if(rank == 0){
            printf("\t variable %d values per sample %d bits per sample %d \n", var, vps[var], bits_per_sample);
        }

        data[var] = (double*)malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable[var]->vps);
        memset(data[var], 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2] * variable[var]->vps);
    }

    ret = PIDX_reset_variable_counter(file);
    if (ret != PIDX_success)  report_error("PIDX_reset_variable_counter", __FILE__, __LINE__);

    for (int var = 0; var < variable_count; var++)
    {
        ret = PIDX_get_next_variable(file, &variable[var]);
        if (ret != PIDX_success)  report_error("PIDX_get_next_variable", __FILE__, __LINE__);

        ret = PIDX_variable_read_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
        if (ret != PIDX_success)  report_error("PIDX_variable_read_data_layout", __FILE__, __LINE__);

        ret = PIDX_read_next_variable(file, variable[var]);
        if (ret != PIDX_success)  report_error("PIDX_read_next_variable", __FILE__, __LINE__);
    }

    ret = PIDX_close(file);
    if (ret != PIDX_success)  report_error("PIDX_close", __FILE__, __LINE__);

    ret = PIDX_close_access(access);
    if (ret != PIDX_success)  report_error("PIDX_close_access", __FILE__, __LINE__);

    for(int var = 0; var < variable_count; var++)
    {
        free(data[var]);
        data[var] = 0;
    }
    free(data);

    return 0;

}

int open_file(PIDX_file& file, PIDX_access& access, std::string output_filename){
    int ret = 1;

    if(restart_time > 0)
    {
        ret = PIDX_file_open(output_filename.c_str(), PIDX_MODE_RDWR, access, global_size, &file);
        if (ret != PIDX_success)  report_error("PIDX_file_open", __FILE__, __LINE__);
    }
    else
    {
        ret = PIDX_file_create(output_filename.c_str(), PIDX_MODE_CREATE, access, global_size, &file);
        if (ret != PIDX_success)  report_error("PIDX_file_create", __FILE__, __LINE__);
        PIDX_set_block_size(file, bits_per_block);
        PIDX_set_block_count(file, blocks_per_file);

    }

    return ret;

}

// Dump the data with PIDX
template<typename T>
int write_data(T** data, PIDX_file& file, PIDX_access& access, int ts, int variable_count, int* vps, std::string output_file_name, std::string var_prefix, PIDX_data_type data_type){

    int ret = 0;

    PIDX_variable* variable;   // variable descriptor
    variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);

    open_file(file, access, output_file_name);

    ret = PIDX_set_current_time_step(file, ts);
    if (ret != PIDX_success)  report_error("PIDX_set_current_time_step", __FILE__, __LINE__);
    ret = PIDX_set_variable_count(file, variable_count);
    if (ret != PIDX_success)  report_error("PIDX_set_variable_count", __FILE__, __LINE__);

    char var_name[512];
    for (int var = 0; var < variable_count; var++)
    {
        sprintf(var_name, var_prefix.c_str(), var);

        ret = PIDX_variable_create(var_name, vps[var] * sizeof(T) * 8, data_type, &variable[var]);
        if (ret != PIDX_success)  report_error("PIDX_variable_create", __FILE__, __LINE__);
        ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
        if (ret != PIDX_success)  report_error("PIDX_variable_data_layout", __FILE__, __LINE__);
        ret = PIDX_append_and_write_variable(file, variable[var]);
        if (ret != PIDX_success)  report_error("PIDX_append_and_write_variable", __FILE__, __LINE__);

    }

    ret = PIDX_close(file);
    if (ret != PIDX_success)  report_error("PIDX_close", __FILE__, __LINE__);

    delete [] variable;

    return ret;

}


// Simple downcast from the checkpoint_data (double) to a downcasted version (float)
void downcast_checkpoint(double** checkpoint_data, float** vis_data,
                            int cp_var_idx, int vis_var_idx, int len){

    vis_data[vis_var_idx] = (float*)malloc(sizeof(float) * len);

    for(int i=0; i < len; i++){
        vis_data[vis_var_idx][i] = (float)(checkpoint_data[cp_var_idx][i]);
    }
}

int main(int argc, char** argv){
    // Initialize PIDX
    init(argc, argv);

    const unsigned checkpoint_t_period = 10;    // checkpoint period
    const unsigned vis_t_period = 5;            // visualization (downcasted) dump period
    const unsigned max_t = 500;                 // total number of timesteps of the simulation

    std::string checkpoint_filename = "checkpoint.idx";
    std::string vis_filename = "vis.idx";

    unsigned t = 0;

    // if a restart timestep has been specified
    // check the checkpoint data and restart the simulation from this time
    if(restart_time > 0){
        t = restart_time; // set current simulation time to the restarting time
        check_last_state(checkpoint_filename, t/checkpoint_t_period);

    }else{ // No restart, so the simulation will start from t=0
        t = 0;
    }

    if(rank == 0){
        printf("Start from t = %d\n", t);
    }

#ifdef VIS_DUMP

    // Define file descriptor and access for visualization dump
    PIDX_file vis_file;
    PIDX_access vis_access;
    PIDX_create_access(&vis_access);

 #if PIDX_HAVE_MPI
    PIDX_set_mpi_access(vis_access, MPI_COMM_WORLD);
 #endif

#endif

    // Define file descriptor and access for checkpoint dump
    PIDX_file checkpoint_file;
    PIDX_access checkpoint_access;
    PIDX_create_access(&checkpoint_access);

#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(checkpoint_access, MPI_COMM_WORLD);
#endif

    // Main simulation loop
    for(; t < max_t; t++){

        // Generate data for Checkpoint dump
        int checkpoint_variable_count = 3;                           //< How many variables
        int checkpoint_vps[checkpoint_variable_count];      //< How many values per sample
        checkpoint_vps[0] = 1;
        checkpoint_vps[1] = 1;
        checkpoint_vps[2] = 1;

        double** checkpoint_data = (double**)malloc(sizeof(*checkpoint_data) * checkpoint_variable_count);
        memset(checkpoint_data, 0, sizeof(*checkpoint_data) * checkpoint_variable_count);

        // Create syntethic data for the checkpoint
        double t1, t2;
        t1 = MPI_Wtime();
        create_data<double>(checkpoint_data, checkpoint_variable_count, checkpoint_vps, t);
        t2 = MPI_Wtime();

        if(rank == 0)
            printf( "Generating state data time is %f\n", t2 - t1 );

        if(t % checkpoint_t_period == 0){
            write_data(checkpoint_data, checkpoint_file, checkpoint_access, t/checkpoint_t_period, checkpoint_variable_count, checkpoint_vps, checkpoint_filename, "state_var_%d", FLOAT64);
        }

#ifdef VIS_DUMP
        // Vis DUMP with downcasting

        int vis_variable_count = 1;                             //< How many variables
        int vis_vps[vis_variable_count];          //< How many values per sample
        vis_vps[0] = 1;

        float** vis_data = (float**)malloc(sizeof(*vis_data) * vis_variable_count);
        memset(vis_data, 0, sizeof(*vis_data) * vis_variable_count);

        // downcast check point data
        downcast_checkpoint(checkpoint_data, vis_data, 0, 0, local_box_size[0]*local_box_size[1]*local_box_size[2]);

        // dump visualization data
        if(t % vis_t_period == 0){
            write_data(vis_data, vis_file, vis_access, t/vis_t_period, vis_variable_count, vis_vps, vis_filename, "data_var_%d", FLOAT32);
        }

        delete [] vis_data[0];
        delete [] vis_data;
#endif

        if(rank == 0)
        {
            printf("time %d\n", t);
        }
        //sleep(1);

        for(int varidx=0; varidx < checkpoint_variable_count; varidx++)
            delete [] checkpoint_data[varidx];

        delete [] checkpoint_data;

    }

    int ret = PIDX_close_access(checkpoint_access);
    if (ret != PIDX_success)  report_error("PIDX_close_access state", __FILE__, __LINE__);

#ifdef VIS_DUMP
    ret = PIDX_close_access(vis_access);
    if (ret != PIDX_success)  report_error("PIDX_close_access data", __FILE__, __LINE__);
#endif

    //  MPI finalize
#if PIDX_HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}


///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
    char flags[] = "g:l:r:";
    int one_opt = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF)
    {
        /* postpone error checking for after while loop */
        switch (one_opt)
        {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &global_box_size[0], &global_box_size[1], &global_box_size[2]);
                break;
            case('l'):
                sscanf(optarg, "%dx%dx%d", &local_box_size[0], &local_box_size[1], &local_box_size[2]);
                break;
            case('r'):
                sscanf(optarg, "%dx", &restart_time);
                break;
            case('?'):
                return (-1);
        }
    }
    /* need positive dimensions */
    if (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1 || local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1)
    {
        fprintf(stderr, "Error: bad dimension specification.\n");
        return (-1);
    }

    /* need global dimension to be larger than the local */
    if (global_box_size[0] < local_box_size[0] || global_box_size[1] < local_box_size[1] || global_box_size[2] < local_box_size[2])
    {
        fprintf(stderr, "Error: Per-process local box size cannot be greater than the global box\n");
        return (-1);
    }

    if (local_box_size[0] == 0 || local_box_size[1] == 0 || local_box_size[2] == 0)
    {
        fprintf(stderr, "Local Dimension cannot be 0!!!!!!!!!\n");
        return (-1);
    }

    if (global_box_size[0] == 0 || global_box_size[1] == 0 || global_box_size[2] == 0)
    {
        fprintf(stderr, "Global Dimension cannot be 0!!!!!!!!!\n");
        return (-1);
    }

    return (0);
}

///< How to use this progam
static void usage(void)
{
    printf("Serial Usage: ./checkpoint-restart -g 4x4x4 -l 4x4x4 [-r restart_timestep]\n");
    printf("Parallel Usage: mpirun -n 8 ./checkpoint-restart -g 4x4x4 -l 2x2x2 [-r restart_timestep]\n");
    printf("  -g: global dimensions\n");
    printf("  -l: local (per-process) dimensions\n");
    printf("\n");

    return;
}

///< Print error and exit program
static void report_error(std::string func_name, std::string file_name, int line_no)
{
    fprintf(stderr, "Error in function %s Program %s Line %d\n", func_name.c_str(), file_name.c_str(), line_no);
#if PIDX_HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(-1);
#endif
}
