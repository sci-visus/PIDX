/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

/*

             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        IDX Format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*

*/

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static unsigned long long global_box_size[3] = {162, 162, 42};
static unsigned long long local_box_offset[3];
static unsigned long long local_box_size[3] = {20, 20, 20};
int sub_div[NUM_DIMS] = {1,1,1};
static int time_step_count = 1;
static int variable_count = 1;
static char output_file_template[512] = "test";
static double **data;
static int *int_data;
static char output_file_name[512] = "test.idx";
static int vps = 1;
static char *usage = "Serial Usage: ./checkpoint -g 32x32x32 -l 32x32x32 -v 3 -t 16 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./checkpoint -g 32x32x32 -l 16x16x16 -f output_idx_file_name -v 3 -t 16\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables\n";


//----------------------------------------------------------------
static void terminate()
{
#if PIDX_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
}

//----------------------------------------------------------------
static void terminate_with_error_msg(const char *format, ...)
{
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
  terminate();
}

//----------------------------------------------------------------
static void rank_0_print(const char *format, ...)
{
  if (rank != 0) return;
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
}

//----------------------------------------------------------------
static void init_mpi(int argc, char **argv)
{
#if PIDX_HAVE_MPI
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
#endif
}

//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}

static void calculate_per_process_offsets()
{
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
}

static void create_synthetic_simulation_data()
{
  int var = 0;
  unsigned long long i, j, k, vps = 0;
  data = malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  // Synthetic simulation data

  for(var = 0; var < variable_count; var++)
  {
    if (var == 0)
    {
      vps = 1;
      int_data = malloc(sizeof (*(int_data)) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
      for (k = 0; k < local_box_size[2]; k++)
        for (j = 0; j < local_box_size[1]; j++)
          for (i = 0; i < local_box_size[0]; i++)
          {
            unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
            for (vps = 0; vps < vps; vps++)
              int_data[index * vps + vps] = var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i));
          }
    }

    else
    {
      if (var == 2)
        vps = 3;
      else
        vps = 1;

      data[var] = malloc(sizeof (*(data[var])) * local_box_size[0] * local_box_size[1] * local_box_size[2] * vps);
      for (k = 0; k < local_box_size[2]; k++)
        for (j = 0; j < local_box_size[1]; j++)
          for (i = 0; i < local_box_size[0]; i++)
          {
            unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
            for (vps = 0; vps < vps; vps++)
              data[var][index * vps + vps] = var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i));
          }
    }
  }
}

static void destroy_synthetic_simulation_data()
{
  int var = 0;
  for(var = 0; var < variable_count; var++)
  {
    if (var == 0)
      continue;

    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;
}

///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) || (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||(local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
        terminate_with_error_msg("Invalid local dimension\n%s", usage);
      break;

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // number of timesteps
      if (sscanf(optarg, "%d", &time_step_count) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('v'): // number of variables
      if (sscanf(optarg, "%d", &variable_count) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

static void check_args()
{
  if (global_box_size[X] < local_box_size[X] || global_box_size[Y] < local_box_size[Y] || global_box_size[Z] < local_box_size[Z])
    terminate_with_error_msg("ERROR: Global box is smaller than local box in one of the dimensions\n");

  // check if the number of processes given by the user is consistent with the actual number of processes needed
  int brick_count = (int)((global_box_size[X] + local_box_size[X] - 1) / local_box_size[X]) *
                    (int)((global_box_size[Y] + local_box_size[Y] - 1) / local_box_size[Y]) *
                    (int)((global_box_size[Z] + local_box_size[Z] - 1) / local_box_size[Z]);
  if(brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);
}

int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();
  calculate_per_process_offsets();
  create_synthetic_simulation_data();

  rank_0_print("Simulation Data Created\n");

  int ret;
  int var;
  int ts;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

  variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_point global_size, local_offset, local_size;
  PIDX_set_point_5D(global_size, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);
  PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);
  PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);

  //  Creating access
  PIDX_access access;
  PIDX_create_access(&access);
#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  for (ts = 0; ts < time_step_count; ts++)
  {
    //  PIDX mandatory calls
    ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, global_size, &file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

    ret = PIDX_set_current_time_step(file, ts);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");
    ret = PIDX_set_variable_count(file, variable_count);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

    PIDX_set_block_count(file, 256);
    PIDX_set_block_size(file, 15);

    //PIDX_set_compression_type(file, PIDX_CHUNKING_ONLY);
    //PIDX_set_compression_type(file, PIDX_CHUNKING_ZFP);
    //PIDX_set_lossy_compression_bit_rate(file, 16);

    char var_name[512];
    for (var = 0; var < variable_count; var++)
    {
      sprintf(var_name, "variable_%d", var);

      if (var == 0)
      {
        vps = 1;
        ret = PIDX_variable_create(var_name,  vps * sizeof(int) * 8, INT32 , &variable[var]);
        if (ret != PIDX_success)  terminate_with_error_msg("A PIDX_variable_create");

        ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, int_data, PIDX_row_major);
        if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_data_layout");
      }
      else
      {
        if (var == 2)
        {
          vps = 3;
          ret = PIDX_variable_create(var_name,  vps * sizeof(double) * 8, FLOAT64_RGB , &variable[var]);
          if (ret != PIDX_success)  terminate_with_error_msg("B[%d] [%d] PIDX_variable_create", var, ret);
        }
        else
        {
          vps = 1;
          ret = PIDX_variable_create(var_name,  vps * sizeof(double) * 8, FLOAT64 , &variable[var]);
          if (ret != PIDX_success)  terminate_with_error_msg("B[%d] [%d] PIDX_variable_create", var, ret);
        }

        ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
        if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_data_layout");
      }

      ret = PIDX_append_and_write_variable(file, variable[var]);
      if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

      //ret = PIDX_flush(file);
      //if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");
    }

    ret = PIDX_close(file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");
  }

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

  free(variable);
  variable = 0;

  destroy_synthetic_simulation_data();
  shutdown_mpi();

  return 0;
}
