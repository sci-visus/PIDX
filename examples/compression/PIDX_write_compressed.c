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

#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

#include "pidx_examples_utils.h"

static int parse_args(int argc, char **argv);
static void usage(void);
static void report_error(char* func_name, char* file_name, int line_no);

//static int global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
//static int local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
//static int time_step_count = 1;                       ///< Number of time-steps
//static int variable_count = 1;                        ///< Number of fields
static char output_file_template[512];   ///< output IDX file Name Template
static int compression_bit_rate = 64;

double **data;             // variable data buffer
int *values_per_sample;    // Example: 1 for scalar 3 for vector

int main(int argc, char **argv)
{
  int ret;
  int i, j, k;
  int var, vps;
  char output_file_name[512];

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

  values_per_sample = malloc(sizeof(*values_per_sample) * variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);

  // Creating the filename
  sprintf(output_file_name, "%s%s", output_file_template,".idx");

  calculate_per_process_offsets();

  data = (double**)malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  // Synthetic simulation data
  for(var = 0; var < variable_count; var++)
  {
    if (var % 2  == 1)
      values_per_sample[var] =  3;
    else
      values_per_sample[var] =  1;
    data[var] = malloc(sizeof (uint64_t) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample[var]);

    for (k = 0; k < local_box_size[2]; k++)
      for (j = 0; j < local_box_size[1]; j++)
        for (i = 0; i < local_box_size[0]; i++)
        {
          int64_t index = (int64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
          for (vps = 0; vps < values_per_sample[var]; vps++)
            data[var][index * values_per_sample[var] + vps] = 100 + var + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i));
        }
  }

  
  //  Creating access
  create_pidx_point_and_access();
  
  int ts;
  PIDX_variable* variable;   // variable descriptor

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);
  
  ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, p_access, global_size, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_open\n");
  
  for (ts = 0; ts < time_step_count; ts++)
  {
    if (ret != PIDX_success)  report_error("PIDX_set_dims", __FILE__, __LINE__);
    ret = PIDX_set_current_time_step(file, ts);
    if (ret != PIDX_success)  report_error("PIDX_set_current_time_step", __FILE__, __LINE__);
    ret = PIDX_set_variable_count(file, variable_count);
    if (ret != PIDX_success)  report_error("PIDX_set_variable_count", __FILE__, __LINE__);

    /* PIDX compression related calls */
    PIDX_set_compression_type(file, PIDX_CHUNKING_ZFP);
    PIDX_set_lossy_compression_bit_rate(file, compression_bit_rate);

    char var_name[512];
    char data_type[512];
    for (var = 0; var < variable_count; var++)
    {
      sprintf(var_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);

      ret = PIDX_variable_create(var_name, values_per_sample[var] * sizeof(uint64_t) * 8, data_type, &variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_variable_create", __FILE__, __LINE__);

      ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
      if (ret != PIDX_success)  report_error("PIDX_variable_data_layout", __FILE__, __LINE__);

      ret = PIDX_append_and_write_variable(file, variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_append_and_write_variable", __FILE__, __LINE__);
    }

    ret = PIDX_close(file);
    if (ret != PIDX_success)  report_error("PIDX_close", __FILE__, __LINE__);
  }

  // Close access and free memory
  if (PIDX_close_access(p_access) != PIDX_success)
    terminate_with_error_msg("PIDX_close_access");

  for(var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;

  shutdown_mpi();

  return 0;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:b:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
      case('g'): // global dimension
        if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) ||
            (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
          terminate_with_error_msg("Invalid global dimensions\n%s", usage);
        break;
        
      case('l'): // local dimension
        if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||
            (local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
          terminate_with_error_msg("Invalid local dimension\n%s", usage);
        break;
        
      case('f'):
          sprintf(output_file_template, "%s", optarg);
          break;
      case('t'):
          sscanf(optarg, "%d", &time_step_count);
          break;
      case('b'):
          sscanf(optarg, "%d", &compression_bit_rate);
          break;
      case('v'):
          sscanf(optarg, "%d", &variable_count);
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
  printf("Serial Usage: ./pidx-s3d-checkpoint -g 4x4x4 -l 4x4x4 -f Filename -t 1 -v 1\n");
  printf("Parallel Usage: mpirun -n 8 ./pidx-s3d-checkpoint -g 4x4x4 -l 2x2x2 -f Filename_ -t 1 -v 1\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("  -v: number of fields\n");

  printf("pidx-s3d-checkpoint generates a 3D volume of size g_x g_y g_z specified by -g command line parameter\nEach process writes a sub-volume of size l_x l_y l_z specified by -l command line parameter. \nData is written in the idx format, with the filename Filename specified by -f command line parameter. \nThe number of time-steps and the number of fields can be optionally provided by -t and -v command line parameters.\n");

  printf("\n");

  return;
}

///< Print error and exit program
static void report_error(char* func_name, char* file_name, int line_no)
{
  fprintf(stderr, "Error in function %s Program %s Line %d\n", func_name, file_name, line_no);
#if PIDX_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
}
