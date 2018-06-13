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

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static unsigned long long local_box_offset[3];
static unsigned long long global_box_size[3] = {0, 0, 0};
static unsigned long long local_box_size[3] = {0, 0, 0};
static int current_ts = 0;
static int variable_index = 0;
static char output_file_template[512] = "test_idx";
static unsigned char *data;
static PIDX_point global_size, local_offset, local_size, global_bounds;
static PIDX_access p_access;
static int variable_count;
static PIDX_file file;
static PIDX_variable variable;
static int bits_per_sample = 0;
static int values_per_sample = 0;
static char type_name[512];
static char output_file_name[512] = "test.idx";
static char *usage = "Serial Usage: ./minmax file_name variable_index\n";

                     //"Parallel Usage: mpirun -n 8 ./restart -g 32x32x32 -l 16x16x16 -f output_idx_file_name\n"
                     //"  -g: global dimensions\n"
                     //"  -l: local (per-process) dimensions\n"
                     //"  -f: IDX filename\n"
                     //"  -t: time step index that needs to be read\n"
                     //"  -v: Variable Index";


static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static void check_args();
static void calculate_per_process_offsets();
static void terminate_with_error_msg(const char *format, ...);
static void terminate();
static void create_pidx_var_point_and_access();
static void set_pidx_file(int ts);
static void set_pidx_variable_and_create_buffer();
static void verify_read_results();
static void shutdown_mpi();

int main(int argc, char **argv)
{
  init_mpi(argc, argv);

  if (argc == 3)
  {
    sprintf(output_file_name, "%s%s", argv[1], ".idx");
    variable_index = atoi(argv[2]);
  }
  else if (argc > 3)
  {
    parse_args(argc, argv);
    check_args();
  }
  else
      terminate_with_error_msg("Wrong Usage\n%s", usage);

  create_pidx_var_point_and_access();

  set_pidx_file(current_ts);

  calculate_per_process_offsets();

  set_pidx_variable_and_create_buffer();
  PIDX_variable_read_data_layout(variable, local_offset, local_size, data, PIDX_row_major);
  PIDX_close(file);

  PIDX_close_access(p_access);
  verify_read_results();

  free(data);
  shutdown_mpi();

  return 0;
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
      if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) ||
          (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||
          (local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
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
  if (brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);

}

static void calculate_per_process_offsets()
{
  if (local_box_size[X] == 0)
      local_box_size[X] = global_bounds[X];
  if (local_box_size[Y] == 0)
      local_box_size[Y] = global_bounds[Y];
  if (local_box_size[Z] == 0)
      local_box_size[Z] = global_bounds[Z];

  int sub_div[NUM_DIMS];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];

  PIDX_set_point(local_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z]);
  PIDX_set_point(local_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);
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
  //  Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  return;
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  ret = PIDX_file_open(output_file_name, PIDX_MODE_RDONLY, p_access, global_bounds, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);

  int last_ts;
  PIDX_get_last_time_step(file, &last_ts);

  if (ts != 0)
    PIDX_set_current_time_step(file, ts);
  else
    PIDX_set_current_time_step(file, last_ts);

  PIDX_get_variable_count(file, &variable_count);

  if (global_box_size[X] == 0 && global_box_size[Y] == 0 && global_box_size[Z] == 0)
  {
    global_box_size[X] = global_bounds[X];
    global_box_size[Y] = global_bounds[Y];
    global_box_size[Z] = global_bounds[Z];
    PIDX_query_box(file, global_bounds);
  }
  else
    PIDX_query_box(file, global_size);

  return;
}

//----------------------------------------------------------------
static void set_pidx_variable_and_create_buffer()
{
  PIDX_return_code ret;

  if (variable_index >= variable_count) terminate_with_error_msg("Variable index more than variable count\n");
  ret = PIDX_set_current_variable_index(file, variable_index);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_variable_index");

  PIDX_get_current_variable(file, &variable);

  PIDX_values_per_datatype(variable->type_name, &values_per_sample, &bits_per_sample);
  strcpy(type_name, variable->type_name);

  data = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
  memset(data, 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
}

//----------------------------------------------------------------
static void verify_read_results()
{
  int i, j, k, vps;
  int int_val = 0;

  double double_val = 0, double_min = 0, double_max = 0;
  float float_val = 0, float_min = 0, float_max = 0;
  int int_min = 0, int_max = 0;
  double *min_temp, *max_temp;
  double min_magnitude = 7000000;
  double max_magnitude = 0;

  bits_per_sample = bits_per_sample / 8;

  if (strcmp(type_name, INT32) == 0)
  {
    memcpy(&int_min, data, bits_per_sample);
    memcpy(&int_max, data, bits_per_sample);
  }
  else if (strcmp(type_name, FLOAT32) == 0)
  {
    memcpy(&float_min, data, bits_per_sample);
    memcpy(&float_max, data, bits_per_sample);
  }
  else if (strcmp(type_name, FLOAT64) == 0)
  {
    memcpy(&double_min, data, bits_per_sample);
    memcpy(&double_max, data, bits_per_sample);
  }
  else if (strcmp(type_name, FLOAT64_RGB) == 0)
  {

    min_temp = malloc(values_per_sample * sizeof (*min_temp));
    max_temp = malloc(values_per_sample * sizeof (*max_temp));

    memcpy(min_temp, data, values_per_sample * sizeof (*min_temp));
    memcpy(max_temp, data, values_per_sample * sizeof (*min_temp));
  }

  //fprintf(stderr, "values_per_sample = %d %d [ %f %f %f ] [ %f %f %f ]\n", values_per_sample, bits_per_sample, min_temp[X], min_temp[Y], min_temp[Z], max_temp[X], max_temp[Y], max_temp[Z]);


  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        uint64_t index = (uint64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;

        if (strcmp(type_name, INT32) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&int_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (int_val < int_min)
              int_min = int_val;

            if (int_val > int_max)
              int_max = int_val;
          }
        }

        else if (strcmp(type_name, FLOAT64) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&double_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (double_val < double_min)
              double_min = double_val;

            if (double_val > double_max)
              double_max = double_val;
          }
        }

        else if (strcmp(type_name, FLOAT64_RGB) == 0)
        {
          double *temp = malloc(values_per_sample * sizeof (*temp));


          double c_magnitude = 0;
          for (vps = 0; vps < values_per_sample; vps++)
            memcpy(&temp[vps], data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);

          c_magnitude = (temp[0]*temp[0]) + (temp[1]*temp[1]) + (temp[2]*temp[2]);
          if (c_magnitude < min_magnitude)
          {
            min_magnitude = c_magnitude;
            min_temp[0] = temp[0];
            min_temp[1] = temp[1];
            min_temp[2] = temp[2];
          }

          if (c_magnitude > max_magnitude)
          {
              max_magnitude = c_magnitude;
              max_temp[0] = temp[0];
              max_temp[1] = temp[1];
              max_temp[2] = temp[2];
          }

          free(temp);

        }

        else if (strcmp(type_name, FLOAT32) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&float_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (float_val < float_min)
              float_min = float_val;

            if (float_val > float_max)
              float_max = float_val;
          }
        }

      }


  if (strcmp(type_name, INT32) == 0)
  {
    int int_final_min = 0, int_final_max = 0;
#if PIDX_HAVE_MPI
    MPI_Allreduce(&int_min, &int_final_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&int_max, &int_final_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    int_final_min = int_min;
    int_final_max = int_max;
#endif
    if (rank == 0)
      fprintf(stderr, "[INT] Min %d Max %d\n", int_final_min, int_final_max);
  }
  else if (strcmp(type_name, FLOAT32) == 0)
  {
    float float_final_min = 0, float_final_max = 0;
#if PIDX_HAVE_MPI
    MPI_Allreduce(&float_min, &float_final_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&float_max, &float_final_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#else
    float_final_min = float_min;
    float_final_max = float_max;
#endif
    if (rank == 0)
      fprintf(stderr, "[FLOAT32] Min %.13f Max %f\n", float_final_min, float_final_max);
  }
  else if (strcmp(type_name, FLOAT64) == 0)
  {
    double double_final_min = 0, double_final_max = 0;
#if PIDX_HAVE_MPI
    MPI_Allreduce(&double_min, &double_final_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&double_max, &double_final_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    double_final_max = final_max;
    double_final_min = final_min;
#endif
    if (rank == 0)
      fprintf(stderr, "Min %f Max %f\n", double_final_min, double_final_max);
  }

  else if (strcmp(type_name, FLOAT64_RGB) == 0)
  {
    if (rank == 0)
      fprintf(stderr, "Min %f %f %f Max %f %f %f\n", min_temp[0], min_temp[1], min_temp[2], max_temp[0], max_temp[1], max_temp[2]);

    free(min_temp);
    free(max_temp);
  }


}


//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}
