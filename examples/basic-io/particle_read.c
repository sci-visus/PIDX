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
  PIDX read example

  In this example we show how to read data using the PIDX library.

  We consider a global 3D regular grid domain that we will call
  global domain (g).
  This global domain represents the grid space where all the data are stored.

  In a parallel environment each core (e.g. MPI rank) owns a portion of the data
  that has to be written on the disk. We refer to this portion of the domain as
  local domain (l).

  In this example we well see how to execute parallel read with PIDX of a
  syntethic dataset.

  In the following picture is represented a sample domain decomposition
  of the global domain (l) in per-core local domains (l), sometimes referred
  as patches.
  In this example all the local domains have same dimesions for simplicity.
  PIDX supports different number and sizes of patches per core.
  This also means that you can actually read the same data from a different
  core configurations.

                                         *---------*--------*
                                       /         /         /| P7
                                      *---------*---------* |
                                     /         /         /| |
                                    *---------*---------* | *
                                    |         |         | |/|
                                    |         |         | * |
                                    | P4      | P5      |/| | P3
        IDX Format      ------>     *---------*---------* | *
                                    |         |         | |/
                                    |         |         | *
                                    | P0      | P1      |/
                                    *---------*---------*

*/

#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define DEBUG_PRINT_OUTPUT 1

#define PARTICLES_POSITION_VAR 0
#define PARTICLES_COLOR_VAR 1

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static unsigned long long logical_local_box_offset[3];
static unsigned long long logical_global_box_size[3] = {0, 0, 0};
static unsigned long long logical_local_box_size[3] = {0, 0, 0};

static double physical_global_box_size[NUM_DIMS] = {0};
static double physical_local_box_offset[NUM_DIMS] = {0};
static double physical_local_box_size[NUM_DIMS] = {0};

static int current_ts = 1;
static int variable_index = 0;
static char output_file_template[512] = "test_idx";
static unsigned char *data = NULL;
static uint64_t particle_count = 0;
static int checkpoint_restart = 1;

static PIDX_point local_offset, local_size;
static PIDX_physical_point physical_local_offset, physical_local_size;

static PIDX_access p_access;
static int variable_count;
static PIDX_file file;
static PIDX_variable variable;
static int bits_per_sample = 0;
static int values_per_sample = 0;
static char type_name[512];
static char output_file_name[512] = "test.idx";
static char *usage = "Serial Usage: ./particle_read -g 32x32x32 -l 32x32x32 -v 0 -f input_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./particle_read -g 32x32x32 -l 16x16x16 -v 0 -f input_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX input filename\n"
                     "  -t: time step index to read\n"
                     "  -v: variable index to read\n"
                     "  -c: read all vars to simulate checkpoint/restart";


static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static void check_args();
static void calculate_per_process_offsets();
static void terminate_with_error_msg(const char *format, ...);
static void terminate();
static void create_pidx_var_point_and_access();
static void set_pidx_file(int ts);
static void set_pidx_variable_and_create_buffer();
static void shutdown_mpi();

int main(int argc, char **argv)
{
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

  // Create variables
  create_pidx_var_point_and_access();

  // Set PIDX_file for this timestep
  set_pidx_file(current_ts);

  PIDX_variable *checkpoint_vars = NULL;
  uint64_t *checkpoint_particle_counts = NULL;
  void **checkpoint_data = NULL;

  if (!checkpoint_restart)
  {
    // Get all the information about the variable that we want to read
    set_pidx_variable_and_create_buffer();

    // Read the data into a local buffer (data) in row major order
    // TODO WILL: the data should be a pointer to a pointer, since PIDX internally
    // will be in charge of allocating the data buffer. So this data param will take
    // a double pointer.
    PIDX_variable_read_particle_data_layout(variable, physical_local_offset, physical_local_size,
        (void**)&data, &particle_count, PIDX_row_major);
  }
  else
  {
    if (rank == 0) {
      printf("Reading checkpoint-restart style\n");
    }

    checkpoint_vars = calloc(variable_count, sizeof(PIDX_variable));
    checkpoint_particle_counts = calloc(variable_count, sizeof(uint64_t));
    checkpoint_data = calloc(variable_count, sizeof(void*));
    PIDX_return_code ret = PIDX_set_current_variable_index(file, 0);
    if (ret != PIDX_success) terminate_with_error_msg("PIDX_set_current_variable_index");

    for (int i = 0; i < variable_count; ++i)
    {
      ret = PIDX_get_next_variable(file, &checkpoint_vars[i]);
      if (ret != PIDX_success) terminate_with_error_msg("PIDX_get_next_variable");

      ret = PIDX_variable_read_particle_data_layout(checkpoint_vars[i], physical_local_offset, physical_local_size,
          &checkpoint_data[i], &checkpoint_particle_counts[i], PIDX_row_major);
      if (ret != PIDX_success) terminate_with_error_msg("PIDX_variable_read_particle_data_layout");
      if (i + 1 < variable_count)
      {
        ret = PIDX_read_next_variable(file, checkpoint_vars[i]);
        if (ret != PIDX_success) terminate_with_error_msg("PIDX_read_next_variable");
      }
    }
  }

  // PIDX_close triggers the actual write on the disk
  // of the variables that we just set
  PIDX_return_code ret = PIDX_close(file);
  if (ret != PIDX_success) terminate_with_error_msg("PIDX_close");


  /*
  printf("Rank %d:\nPhysical box {[%f, %f, %f], [%f, %f, %f]}\n",
      rank,
      physical_local_box_offset[0], physical_local_box_offset[1], physical_local_box_offset[2],
      physical_local_box_offset[0] + physical_local_box_size[0],
      physical_local_box_offset[1] + physical_local_box_size[1],
      physical_local_box_offset[2] + physical_local_box_size[2]);

  if (checkpoint_restart)
  {
    for (int i = 0; i < variable_count; ++i) {
      printf("Read %lu attribs for variable %d\n", checkpoint_particle_counts[i], i);
    }
    particle_count = checkpoint_particle_counts[0];
  }
  */

  /*
  printf("Logical box {[%lld, %lld, %lld], [%lld, %lld, %lld]}\nLoaded %lu particles\n",
      logical_local_box_offset[0], logical_local_box_offset[1], logical_local_box_offset[2],
      logical_local_box_offset[0] + logical_local_box_size[0],
      logical_local_box_offset[1] + logical_local_box_size[1],
      logical_local_box_offset[2] + logical_local_box_size[2],
      particle_count);
      */

  // Close PIDX_access
  PIDX_close_access(p_access);

#if DEBUG_PRINT_OUTPUT
  char rank_filename[PATH_MAX];
  sprintf(rank_filename, "%s_r_%d", output_file_template, rank);
  FILE *fp = fopen(rank_filename, "w");

  int ch1, ch2;
  
  for (uint32_t p = 0; p < checkpoint_particle_counts[0]; p++)
  {
    double px, py, pz;
    memcpy(&px, (char*)checkpoint_data[PARTICLES_POSITION_VAR] + (p * 3 + 0) * sizeof(double), sizeof(double));
    memcpy(&py, (char*)checkpoint_data[PARTICLES_POSITION_VAR] + (p * 3 + 1) * sizeof(double), sizeof(double));
    memcpy(&pz, (char*)checkpoint_data[PARTICLES_POSITION_VAR] + (p * 3 + 2) * sizeof(double), sizeof(double));

    double d1;
    memcpy(&d1, (char*)checkpoint_data[PARTICLES_COLOR_VAR] + (p) * sizeof(double), sizeof(double));

    double d2;
    memcpy(&d2, (char*)checkpoint_data[2] + (p) * sizeof(double), sizeof(double));

    double d3;
    memcpy(&d3, (char*)checkpoint_data[3] + (p) * sizeof(double), sizeof(double));

    int i1;
    memcpy(&i1, (char*)checkpoint_data[4] + (p) * sizeof(int), sizeof(int));

    double e1, e2, e3, e4, e5, e6, e7, e8, e9;
    memcpy(&e1, (char*)checkpoint_data[5] + (p * 9 + 0) * sizeof(double), sizeof(double));
    memcpy(&e2, (char*)checkpoint_data[5] + (p * 9 + 1) * sizeof(double), sizeof(double));
    memcpy(&e3, (char*)checkpoint_data[5] + (p * 9 + 2) * sizeof(double), sizeof(double));
    memcpy(&e4, (char*)checkpoint_data[5] + (p * 9 + 3) * sizeof(double), sizeof(double));
    memcpy(&e5, (char*)checkpoint_data[5] + (p * 9 + 4) * sizeof(double), sizeof(double));
    memcpy(&e6, (char*)checkpoint_data[5] + (p * 9 + 5) * sizeof(double), sizeof(double));
    memcpy(&e7, (char*)checkpoint_data[5] + (p * 9 + 6) * sizeof(double), sizeof(double));
    memcpy(&e8, (char*)checkpoint_data[5] + (p * 9 + 7) * sizeof(double), sizeof(double));
    memcpy(&e9, (char*)checkpoint_data[5] + (p * 9 + 8) * sizeof(double), sizeof(double));

    int64_t particle_id;
    memcpy(&particle_id, (char*)checkpoint_data[6] + (p) * sizeof(int64_t), sizeof(int64_t));

    fprintf(fp, "%f %f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %lld\n", px, py, pz, d1, d2, d3, i1, e1, e2, e3, e4, e5, e6, e7, e8, e9, particle_id);
  }
  
  fclose(fp);
  
  // Verify output of read and write (TODO check for RST mode)
  char rank_filename_verify[PATH_MAX];
  sprintf(rank_filename_verify, "%s_w_%d", output_file_template, rank);
  FILE *fpv = fopen(rank_filename_verify, "r");
  fp = fopen(rank_filename, "r");
  
  int error_count = 0;
  if (fpv == NULL || fp == NULL) {
    printf("Cannot verify rank %d looking for file %s and %s for reading ", rank, rank_filename, rank_filename_verify);
    exit(1);
  } else {
    ch1 = getc(fp);
    ch2 = getc(fpv);
    
    while ((ch1 != EOF) && (ch2 != EOF) && (ch1 == ch2)) {
      ch1 = getc(fp);
      ch2 = getc(fpv);
      
      if(ch1 != ch2)
        error_count++;
    }
    
    fclose(fpv);
    fclose(fp);
  }
  
#endif
  
  // TODO fix error count
  int total_errors = 0;
//  MPI_Allreduce(&error_count, &total_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//  
//  if(rank == 0)
    printf("Error Count: %d this %d\n", total_errors, error_count);
  
  free(data);
  if (checkpoint_restart)
  {
    for (int i = 0; i < variable_count; ++i)
    {
      free(checkpoint_data[i]);
    }
    free(checkpoint_data);
    free(checkpoint_vars);
    free(checkpoint_particle_counts);
  }
  shutdown_mpi();
  
  return error_count;

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



static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:c";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &logical_global_box_size[0], &logical_global_box_size[1], &logical_global_box_size[2]) == EOF) ||
          (logical_global_box_size[0] < 1 || logical_global_box_size[1] < 1 || logical_global_box_size[2] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &logical_local_box_size[0], &logical_local_box_size[1], &logical_local_box_size[2]) == EOF) ||
          (logical_local_box_size[0] < 1 || logical_local_box_size[1] < 1 || logical_local_box_size[2] < 1))
        terminate_with_error_msg("Invalid local dimension\n%s", usage);
      break;

    case('f'): // input file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // number of timesteps
      if (sscanf(optarg, "%d", &current_ts) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('v'): // number of variables
      if (sscanf(optarg, "%d", &variable_index) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('c'): // checkpoint restart style read
      checkpoint_restart = 1;
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }


}

static void check_args()
{
  if (logical_global_box_size[X] < logical_local_box_size[X] || logical_global_box_size[Y] < logical_local_box_size[Y] || logical_global_box_size[Z] < logical_local_box_size[Z])
    terminate_with_error_msg("ERROR: Global box is smaller than local box in one of the dimensions\n");

  // check if the number of processes given by the user is consistent with the actual number of processes needed
  int brick_count = (int)((logical_global_box_size[X] + logical_local_box_size[X] - 1) / logical_local_box_size[X]) *
                    (int)((logical_global_box_size[Y] + logical_local_box_size[Y] - 1) / logical_local_box_size[Y]) *
                    (int)((logical_global_box_size[Z] + logical_local_box_size[Z] - 1) / logical_local_box_size[Z]);
  if (brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);

}

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

  physical_global_box_size[X] = logical_global_box_size[X];
  physical_global_box_size[Y] = logical_global_box_size[Y];
  physical_global_box_size[Z] = logical_global_box_size[Z];

  physical_local_box_size[X] = ((double)logical_local_box_size[X] / logical_global_box_size[X]) * physical_global_box_size[X];
  physical_local_box_size[Y] = ((double)logical_local_box_size[Y] / logical_global_box_size[Y]) * physical_global_box_size[Y];
  physical_local_box_size[Z] = ((double)logical_local_box_size[Z] / logical_global_box_size[Z]) * physical_global_box_size[Z];

  physical_local_box_offset[X] = ((double)logical_local_box_offset[X] / logical_global_box_size[X]) * physical_global_box_size[X];
  physical_local_box_offset[Y] = ((double)logical_local_box_offset[Y] / logical_global_box_size[Y]) * physical_global_box_size[Y];
  physical_local_box_offset[Z] = ((double)logical_local_box_offset[Z] / logical_global_box_size[Z]) * physical_global_box_size[Z];

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

  PIDX_set_point(local_offset, logical_local_box_offset[X], logical_local_box_offset[Y], logical_local_box_offset[Z]);
  PIDX_set_point(local_size, logical_local_box_size[X], logical_local_box_size[Y], logical_local_box_size[Z]);

  PIDX_set_physical_point(physical_local_offset, physical_local_box_offset[X], physical_local_box_offset[Y], physical_local_box_offset[Z]);
  PIDX_set_physical_point(physical_local_size, physical_local_box_size[X], physical_local_box_size[Y], physical_local_box_size[Z]);

  //  Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  return;
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  // Open IDX file
  ret = PIDX_file_open(output_file_name, PIDX_MODE_RDONLY, p_access, NULL, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_open\n");

  // Set the current timestep
  PIDX_set_current_time_step(file, ts);
  // Get the total number of variables
  PIDX_get_variable_count(file, &variable_count);

  //PIDX_disable_agg(file);
  return;
}

//----------------------------------------------------------------
static void set_pidx_variable_and_create_buffer()
{
  PIDX_return_code ret;

  // Check if the index variable requested is valid (< num variables in the dataset)
  if (variable_index >= variable_count) terminate_with_error_msg("Variable index more than variable count\n");

  // Set the variable index that we want to read
  ret = PIDX_set_current_variable_index(file, variable_index);
  //ret = PIDX_set_current_variable_by_name(file, "var_1");
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_variable_index\n");

  // Get corresponding PIDX_variable
  PIDX_get_current_variable(file, &variable);

  // Get some information about this variable (typename, number of values per sample,
  // number of bits per sample)
  PIDX_values_per_datatype(variable->type_name, &values_per_sample, &bits_per_sample);
  strcpy(type_name, variable->type_name);

  // TODO WILL: The user will not provide us a buffer, we'll have to allocate it on the fly
  //data = malloc((uint64_t)(bits_per_sample/8) * logical_local_box_size[0] * logical_local_box_size[1] * logical_local_box_size[2]  * values_per_sample);
  //memset(data, 0, (uint64_t)(bits_per_sample/8) * logical_local_box_size[0] * logical_local_box_size[1] * logical_local_box_size[2]  * values_per_sample);
}



//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}
