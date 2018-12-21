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

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

#include "pidx_examples_utils.h"

static char input_file_name[512];
int current_ts = 0;
int variable_index = 0;
int values_per_sample = 0;
int bits_per_sample = 0;
PIDX_variable variable;
char type_name[512];
static PIDX_point global_bounds;
unsigned char *data;

static char *usage = "Serial Usage: ./idx_read -g 32x32x32 -l 32x32x32 -v 0 -f input_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./idx_read -g 32x32x32 -l 16x16x16 -f -v 0 input_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX input filename\n"
                     "  -t: time step index to read\n"
                     "  -v: variable index to read";

static void parse_args(int argc, char **argv);
static void set_pidx_variable_and_create_buffer();
static int verify_read_results();
static void set_pidx_file(int ts);


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
  create_pidx_point_and_access();

  // Set PIDX_file for this timestep
  set_pidx_file(current_ts);

  // Get all the information about the variable that we want to read
  set_pidx_variable_and_create_buffer();

  // Read the data into a local buffer (data) in row major order
  PIDX_variable_read_data_layout(variable, local_offset, local_size, data, PIDX_row_major);

  // PIDX_close triggers the actual write on the disk
  // of the variables that we just set
  PIDX_close(file);

  // Close PIDX_access
  PIDX_close_access(p_access);

  // Compare the data that we just against the syntethic data
  int ret = verify_read_results();
  
  free(data);
  shutdown_mpi();
  
  return ret;
}

static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:";
  int one_opt = 0;
  char input_file_template[512];

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
      if (sprintf(input_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(input_file_name, "%s%s", input_file_template, ".idx");
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

#if 0
  if (rank == 0)
  {
    local_box_size[X] = 1025;
    local_box_size[Y] = 1025;
    local_box_size[Z] = 2612;

    local_box_offset[X] = 0;
    local_box_offset[Y] = 0;
    local_box_offset[Z] = 0;
  }

  else if (rank == 1)
  {
    local_box_size[X] = 1025;
    local_box_size[Y] = 1025;
    local_box_size[Z] = 2612;

    local_box_offset[X] = 1023;
    local_box_offset[Y] = 0;
    local_box_offset[Z] = 0;
  }

  else if (rank == 2)
  {
    local_box_size[X] = 1025;
    local_box_size[Y] = 1025;
    local_box_size[Z] = 2612;

    local_box_offset[X] = 0;
    local_box_offset[Y] = 1023;
    local_box_offset[Z] = 0;
  }

  else if (rank == 3)
  {
    local_box_size[X] = 1025;
    local_box_size[Y] = 1025;
    local_box_size[Z] = 2612;

    local_box_offset[X] = 1023;
    local_box_offset[Y] = 1023;
    local_box_offset[Z] = 0;
  }
#endif

}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  // Open IDX file
  ret = PIDX_file_open(input_file_name, PIDX_MODE_RDONLY, p_access, global_bounds, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_open\n");

  PIDX_query_box(file, global_size);

  // Set the current timestep
  PIDX_set_current_time_step(file, ts);
  // Get the total number of variables
  PIDX_get_variable_count(file, &variable_count);
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

  data = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
  memset(data, 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
}

//----------------------------------------------------------------
static int verify_read_results()
{
  int read_error_count = 0, read_count = 0;
  int total_read_error_count = 0, total_read_count = 0;
  int i, j, k, vps;
  int int_val = 0;
  double double_val = 0;
  int64_t int64_val = 0;
  uint64_t uint64_val = 0;
  float float_val = 0;
  int var = variable_index;

  bits_per_sample = bits_per_sample / 8;
  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        uint64_t index = (uint64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;

        if (strcmp(type_name, PIDX_DType.INT32) == 0 || strcmp(type_name, PIDX_DType.INT32_GA) == 0 || strcmp(type_name, PIDX_DType.INT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&int_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (int_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %d\n", i,j ,k, vps, int_val);//, var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

        else if (strcmp(type_name, PIDX_DType.FLOAT64) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_RGB) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_GA) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&double_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (double_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val,100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val, var /*+ vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
            }
          }
        }

        else if (strcmp(type_name, PIDX_DType.FLOAT32) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_GA) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&float_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (float_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
//              if (rank == 1)
//                printf("W [%d] [%d %d %d] [%d] Read error %f %d Diff %f\n", rank, i,j ,k, vps, float_val, var + vps + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)), (float) (float_val - (var + vps + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, float_val, 100 + var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }
        
        else if (strcmp(type_name, PIDX_DType.INT64) == 0 || strcmp(type_name, PIDX_DType.INT64_RGB) == 0 || strcmp(type_name, PIDX_DType.INT64_GA) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&int64_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (int64_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              if (rank == 0)
                printf("W[%d %d %d] [%d] Read error %lld %lld\n", i,j ,k, vps, int64_val,100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %lld %lld\n", i,j ,k, vps, int64_val, var /*+ vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
            }
          }
        }
        
        else if (strcmp(type_name, PIDX_DType.UINT64) == 0 || strcmp(type_name, PIDX_DType.UINT64_RGB) == 0 || strcmp(type_name, PIDX_DType.UINT64_GA) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&uint64_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (uint64_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //printf("W[%d %d %d] [%d] Read error %lld %lld\n", i,j ,k, vps, int64_val,100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %lld %lld\n", i,j ,k, vps, int64_val, var /*+ vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
            }
          }
        }

      }

  MPI_Allreduce(&read_count, &total_read_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&read_error_count, &total_read_error_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0)
    printf("Correct Sample Count %d Incorrect Sample Count %d\n", total_read_count, total_read_error_count);
  
  if(total_read_error_count == 0 && (total_read_count == global_bounds[2]*global_bounds[1]*global_bounds[0]*values_per_sample)){
    return 0;
  }
  else{
    return 1;
  }

}
