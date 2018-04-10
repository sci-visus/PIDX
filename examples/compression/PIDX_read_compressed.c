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

/// main
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
  verify_read_results();
  
  free(data);
  shutdown_mpi();

  return 0;
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
  
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;
  
  // Open IDX file
  ret = PIDX_file_open(input_file_name, PIDX_MODE_RDONLY, p_access, global_bounds, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_open\n");
  
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
  float float_val = 0;
  int var = variable_index;
  
  bits_per_sample = bits_per_sample / 8;
  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
        
        if (strcmp(type_name, INT32) == 0 || strcmp(type_name, INT32_GA) == 0 || strcmp(type_name, INT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&int_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (int_val != var + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
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
        
        else if (strcmp(type_name, FLOAT64) == 0 || strcmp(type_name, FLOAT64_RGB) == 0 || strcmp(type_name, FLOAT64_GA) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&double_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (double_val != var + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val,100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val, var /*+ vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
            }
          }
        }
        
        else if (strcmp(type_name, FLOAT32) == 0 || strcmp(type_name, FLOAT32_GA) == 0 || strcmp(type_name, FLOAT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&float_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (float_val != var + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 1)
              //printf("W [%d] [%d %d %d] [%d] Read error %f %d Diff %f\n", rank, i,j ,k, vps, float_val, var + vps + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)), (float) (float_val - (var + vps + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }
        
      }
  
  MPI_Allreduce(&read_count, &total_read_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&read_error_count, &total_read_error_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  if (rank == 0)
    printf("Correct Sample Count %d Incorrect Sample Count %d\n", total_read_count, total_read_error_count);
  
  if(total_read_error_count == 0 && (total_read_count == global_bounds[2]*global_bounds[1]*global_bounds[0]*values_per_sample))
    return 0;
  else
    return 1;
  
}

