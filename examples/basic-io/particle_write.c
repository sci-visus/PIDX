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

/*
  PIDX write example

  In this example we show how to write data using the PIDX library.

  We consider a global 3D regular grid domain that we will call
  global domain (g).
  This global domain represents the grid space where all the data are stored.

  In a parallel environment each core (e.g. MPI rank) owns a portion of the data
  that has to be written on the disk. We refer to this portion of the domain as
  local domain (l).

  In this example we well see how to execute parallel write with PIDX of a
  syntethic dataset.

  In the following picture is represented a sample domain decomposition
  of the global domain (l) in per-core local domains (l), sometimes referred
  as patches.
  In this example all the local domains have same dimesions for simplicity.
  PIDX supports different number and sizes of patches per core.

             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        IDX Data format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*

*/
#if !defined _MSC_VER
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#endif
#include <stdarg.h>
#include <stdint.h>
#include <ctype.h>
#include <PIDX.h>

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

#define TYPE_COUNT 5
#define COLOR_COUNT 2
#define MAX_VAR_COUNT 256
enum { X, Y, Z, NUM_DIMS };

static int process_count = 1, rank = 0;
static unsigned long long logical_global_box_size[NUM_DIMS];
static unsigned long long logical_local_box_offset[NUM_DIMS];
static unsigned long long logical_local_box_size[NUM_DIMS];

static double physical_global_box_size[NUM_DIMS] = {1.0, 1.0, 1.0};
static double physical_local_box_offset[NUM_DIMS];
static double physical_local_box_size[NUM_DIMS];

/*
 * position double vector (3)
 * color double scalar
 * density double scalar
 * volume double scalar
 * type int scalar
 * stress double tensor (9)
 */
static int variable_count = 6;



static int time_step_count = 1;
static size_t particle_count = 32;
static char output_file_template[512];
static unsigned char **data;
static char output_file_name[512];
static char var_name[MAX_VAR_COUNT][512];
static int bpv[MAX_VAR_COUNT];
static char type_name[MAX_VAR_COUNT][512];
static int vps[MAX_VAR_COUNT];

static PIDX_point logical_global_size;
static PIDX_physical_point physical_global_size, physical_local_offset, physical_local_size;
static PIDX_access p_access;
static PIDX_file file;
static PIDX_variable* variable;

static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static int generate_vars();
static void check_args();
static void calculate_per_process_offsets();
static void create_synthetic_simulation_data();
static void terminate_with_error_msg(const char *format, ...);
static void terminate();
static void set_pidx_file(int ts);
static void set_pidx_variable();
static void create_pidx_var_point_and_access();
static void destroy_pidx_var_point_and_access();
static void destroy_synthetic_simulation_data();
static void shutdown_mpi();
static char *usage = "Serial Usage: ./idx_write -g 32x32x32 -l 32x32x32 -v 2 -t 4 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./idx_write -g 64x64x64 -l 32x32x32 -v 2 -t 4 -f output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -r: restructured box dimension\n"
                     "  -f: file name template (without .idx)\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables (or file containing a list of variables)\n";

int main(int argc, char **argv)
{
  int ts = 0, var = 0;

  // Init MPI and MPI vars (e.g. rank and process_count)
  init_mpi(argc, argv);

  // Parse input arguments and initialize
  // corresponing variables
  parse_args(argc, argv);

  // Verify that the domain decomposition is valid
  // for the given number of cores
  check_args();

  // Initialize per-process local domain
  calculate_per_process_offsets();


  // Generate synthetic data
  create_synthetic_simulation_data();


  // Create variables
  create_pidx_var_point_and_access();

  for (ts = 0; ts < time_step_count; ts++)
  {
    // Set PIDX_file for this timestep
    set_pidx_file(ts);

    // Set all the PIDX_variable that we want to write
    for (var = 0; var < variable_count; var++)
      set_pidx_variable(var);

    // PIDX_close triggers the actual write on the disk
    // of the variables that we just set
    PIDX_close(file);
  }

  // Clean up our mess
  destroy_pidx_var_point_and_access();

  destroy_synthetic_simulation_data();

  shutdown_mpi();

  return 0;
}

//----------------------------------------------------------------
static void init_mpi(int argc, char **argv)
{
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
}


//----------------------------------------------------------------
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:p:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &logical_global_box_size[X], &logical_global_box_size[Y], &logical_global_box_size[Z]) == EOF) || (logical_global_box_size[X] < 1 || logical_global_box_size[Y] < 1 || logical_global_box_size[Z] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &logical_local_box_size[X], &logical_local_box_size[Y], &logical_local_box_size[Z]) == EOF) ||(logical_local_box_size[X] < 1 || logical_local_box_size[Y] < 1 || logical_local_box_size[Z] < 1))
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

    case('p'): // number of particles per patch
      if (sscanf(optarg, "%lu", &particle_count) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }

  generate_vars();

}

static int generate_vars(){

 /*
  * position double vector (3)
  * color double scalar
  * density double scalar
  * volume double scalar
  * type int scalar
  * stress double tensor (9)
  */

  bpv[0] = sizeof(double) * CHAR_BIT;
  vps[0] = 3;
  strcpy(type_name[0], "3*float64");
  strcpy(var_name[0], "position");

  bpv[1] = sizeof(double) * CHAR_BIT;
  vps[1] = 1;
  strcpy(type_name[1], "1*float64");
  strcpy(var_name[1], "color");

  bpv[2] = sizeof(double) * CHAR_BIT;
  vps[2] = 1;
  strcpy(type_name[2], "1*float64");
  strcpy(var_name[2], "density");

  bpv[3] = sizeof(double) * CHAR_BIT;
  vps[3] = 1;
  strcpy(type_name[3], "1*float64");
  strcpy(var_name[3], "volume");

  bpv[4] = sizeof(int) * CHAR_BIT;
  vps[4] = 1;
  strcpy(type_name[4], "1*int32");
  strcpy(var_name[4], "type");

  bpv[5] = sizeof(double) * CHAR_BIT;
  vps[5] = 9;
  strcpy(type_name[5], "9*float64");
  strcpy(var_name[5], "stress");

  return 0;
}



//----------------------------------------------------------------
static void check_args()
{
  if (logical_global_box_size[X] < logical_local_box_size[X] || logical_global_box_size[Y] < logical_local_box_size[Y] || logical_global_box_size[Z] < logical_local_box_size[Z])
    terminate_with_error_msg("ERROR: Global box is smaller than local box in one of the dimensions\n");

  // check if the number of processes given by the user is consistent with the actual number of processes needed
  int brick_count = (int)((logical_global_box_size[X] + logical_local_box_size[X] - 1) / logical_local_box_size[X]) *
                    (int)((logical_global_box_size[Y] + logical_local_box_size[Y] - 1) / logical_local_box_size[Y]) *
                    (int)((logical_global_box_size[Z] + logical_local_box_size[Z] - 1) / logical_local_box_size[Z]);
  if(brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);
}

//----------------------------------------------------------------
static void calculate_per_process_offsets()
{
  int sub_div[NUM_DIMS];
  sub_div[X] = (logical_global_box_size[X] / logical_local_box_size[X]);
  sub_div[Y] = (logical_global_box_size[Y] / logical_local_box_size[Y]);
  sub_div[Z] = (logical_global_box_size[Z] / logical_local_box_size[Z]);
  logical_local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * logical_local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  logical_local_box_offset[Y] = (slice / sub_div[X]) * logical_local_box_size[Y];
  logical_local_box_offset[X] = (slice % sub_div[X]) * logical_local_box_size[X];

  physical_local_box_size[X] = physical_global_box_size[X] / sub_div[X];
  physical_local_box_size[Y] = physical_global_box_size[Y] / sub_div[Y];
  physical_local_box_size[Z] = physical_global_box_size[Z] / sub_div[Z];

  physical_local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * physical_local_box_size[Z];
  physical_local_box_offset[Y] = (slice / sub_div[X]) * physical_local_box_size[Y];
  physical_local_box_offset[X] = (slice % sub_div[X]) * physical_local_box_size[X];

  printf("[%d] [%lld %lld %lld : %f %f %f] - [%lld %lld %lld - %lld %lld %lld] [%f %f %f - %f %f %f]\n", rank,
         logical_global_box_size[X], logical_global_box_size[Y], logical_global_box_size[Z], physical_global_box_size[X], physical_global_box_size[Y], physical_global_box_size[Z],
         logical_local_box_offset[X], logical_local_box_offset[Y], logical_local_box_offset[Z], logical_local_box_size[X], logical_local_box_size[Y], logical_local_box_size[Z],
         physical_local_box_offset[X], physical_local_box_offset[Y], physical_local_box_offset[Z], physical_local_box_size[X], physical_local_box_size[Y], physical_local_box_size[Z]);
}

//----------------------------------------------------------------
static void create_synthetic_simulation_data()
{
  /*
   * position double vector (3)
   * color double scalar
   * density double scalar
   * volume double scalar
   * type int scalar
   * stress double vector (9)
   */

  data = malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  // Synthetic simulation data

  // data[0] corresponds to the location of the particle
  data[0] = malloc (particle_count * 3 * sizeof(double));

  // data[1] corresponds to colors doubles (1)
  data[1] = malloc (particle_count * sizeof(double));

  // data[2] corresponds to density doubles (1)
  data[2] = malloc (particle_count * sizeof(double));

  // data[3] corresponds to volume doubles (1)
  data[3] = malloc (particle_count * sizeof(double));

  // data[4] corresponds to type int (1)
  data[4] = malloc (particle_count * sizeof(int));

  // data[5] corresponds to tensor volume doubles (9)
  data[5] = malloc (particle_count * sizeof(double) * 9);

  double scale;
  double color;
  int type;
  int particle_type[TYPE_COUNT] = {1,2,3,4,5};
  double particle_color[COLOR_COUNT] = {0.25, 0.75};

  srand((unsigned int)time(NULL));
  for (int k = 0; k < particle_count; k++)
  {
    for (int j = 0; j < 3; j++)
    {
      scale = physical_local_box_offset[j] + (rand() / (float) RAND_MAX) * physical_local_box_size[j];
      memcpy(data[0] + (k * 3 + j) * sizeof(double), &scale, sizeof(double));
    }

    color = particle_color[rand() % COLOR_COUNT];
    memcpy(data[1] + (k) * sizeof(double), &color, sizeof(double));

    scale = rand() / (float) RAND_MAX;
    memcpy(data[2] + (k) * sizeof(double), &scale, sizeof(double));

    scale = rand() / (float) RAND_MAX;
    memcpy(data[3] + (k) * sizeof(double), &scale, sizeof(double));

    type = particle_type[rand() % TYPE_COUNT];
    memcpy(data[4] + (k) * sizeof(int), &type, sizeof(int));

    for (int j = 0; j < 9; j++)
    {
      scale = rand() / (float) RAND_MAX;
      memcpy(data[5] + (k * 9 + j) * sizeof(double), &scale, sizeof(double));
    }
    //printf("[%d] -> %f %f %f\n", k, dvalue_X, dvalue_Y, dvalue_Z);
  }
}

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
static void create_pidx_var_point_and_access()
{
  // Allocate a PIDX_variable array where we store the information
  // of all the variables
  variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  // Set variables that define the global and local domain information
  PIDX_set_point(logical_global_size, logical_global_box_size[X], logical_global_box_size[Y], logical_global_box_size[Z]);

  PIDX_set_physical_point(physical_global_size, physical_global_box_size[X], physical_global_box_size[Y], physical_global_box_size[Z]);
  PIDX_set_physical_point(physical_local_offset, physical_local_box_offset[X], physical_local_box_offset[Y], physical_local_box_offset[Z]);
  PIDX_set_physical_point(physical_local_size, physical_local_box_size[X], physical_local_box_size[Y], physical_local_box_size[Z]);

  // Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  return;
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  // Create IDX file
  ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, p_access, logical_global_size, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create\n");

  PIDX_set_physical_dims(file, physical_global_size);

  // Set the current timestep
  // WILL: Always set it to 0, to reduce disk consumption
  PIDX_set_current_time_step(file, 0);
  // Set the number of variables
  PIDX_set_variable_count(file, variable_count);

  // Select I/O mode (PIDX_IDX_IO for the multires, PIDX_RAW_IO for non-multires)
  PIDX_set_io_mode(file, PIDX_IDX_IO); // TODO:

  return;
}

//----------------------------------------------------------------
static void set_pidx_variable(int var)
{
  PIDX_return_code ret = 0;

  // Set variable name, number of bits, typename
  ret = PIDX_variable_create(var_name[var], bpv[var] * vps[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  // Set the variable offset and size of the local domain,
  // where the data is in memory (data) and what is its layout in memory (row major)
  ret = PIDX_variable_write_particle_data_physical_layout(variable[var], physical_local_offset, physical_local_size, data[var], particle_count, PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

  // Tell PIDX that we want to write this variable
  ret = PIDX_append_and_write_variable(file, variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  return;
}

//----------------------------------------------------------------
static void destroy_pidx_var_point_and_access()
{
  if (PIDX_close_access(p_access) != PIDX_success)
    terminate_with_error_msg("PIDX_close_access");

  free(variable);
  variable = 0;

  return;
}

//----------------------------------------------------------------
static void destroy_synthetic_simulation_data()
{
  int var = 0;
  for(var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;
}

//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}

