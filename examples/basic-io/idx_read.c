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
#include <unistd.h>
#include <getopt.h>

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static unsigned long long local_box_offset[3];
static unsigned long long global_box_size[3] = {0, 0, 0};
static unsigned long long local_box_size[3] = {0, 0, 0};
static int current_ts = 1;
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
static char *usage = "Serial Usage: ./idx_read -g 32x32x32 -l 32x32x32 -v 0 -f input_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./idx_read -g 32x32x32 -l 16x16x16 -f -v 0 input_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX input filename\n"
                     "  -t: time step index to read\n"
                     "  -v: variable index to read";


static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static void check_args();
static void calculate_per_process_offsets();
static void terminate_with_error_msg(const char *format, ...);
static void terminate();
static void create_pidx_var_point_and_access();
static void set_pidx_file(int ts);
static void set_pidx_variable_and_create_buffer();
static int verify_read_results();
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

#if 0
  if (rank == 0)
  {
    local_box_size[X] = 514;
    local_box_size[Y] = 1024;
    local_box_size[Z] = 512;

    local_box_offset[X] = 0;
    local_box_offset[Y] = 0;
    local_box_offset[Z] = 0;
  }
  else if (rank == 1)
  {
    local_box_size[X] = 512;
    local_box_size[Y] = 1024;
    local_box_size[Z] = 512;

    local_box_offset[X] = 512;
    local_box_offset[Y] = 0;
    local_box_offset[Z] = 0;
  }
#endif

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

static void calculate_per_process_offsets()
{
  int sub_div[NUM_DIMS];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
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
  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);
  PIDX_set_point(local_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z]);
  PIDX_set_point(local_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);

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
  ret = PIDX_file_open(output_file_name, PIDX_MODE_RDONLY, p_access, global_bounds, NULL, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_open\n");

  // Set the current timestep
  PIDX_set_current_time_step(file, ts);
  // Get the total number of variables
  PIDX_get_variable_count(file, &variable_count);

  PIDX_query_box(file, global_size);

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

  data = malloc((uint64_t)(bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
  memset(data, 0, (uint64_t)(bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample);
}

//----------------------------------------------------------------
static int verify_read_results()
{
  int read_error_count = 0, read_count = 0;
  int total_read_error_count = 0, total_read_count = 0;
  int i, j, k, vps;
  unsigned char uchar_val = 0;
  short short_val = 0;
  int int_val = 0;
  double double_val = 0;
  float float_val = 0;
  int var = variable_index;

  //if (rank == 100)
  //    fprintf(stderr, "LBO %d %d %d LBS %d %d %d\n", local_box_offset[0], local_box_offset[1], local_box_offset[2], local_box_size[0], local_box_size[1], local_box_size[2]);

  bits_per_sample = bits_per_sample / 8;
  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;

        if (strcmp(type_name, UINT8) == 0 || strcmp(type_name, UINT8_GA) == 0 || strcmp(type_name, UINT8_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&uchar_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (uchar_val != var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //fprintf(stderr, "W[%d %d %d] [%d] Read error %d\n", i,j ,k, vps, int_val);//, var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  fprintf(stderr, "C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

        else if (strcmp(type_name, INT16) == 0 || strcmp(type_name, INT16_GA) == 0 || strcmp(type_name, INT16_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&short_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (short_val != var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //fprintf(stderr, "W[%d %d %d] [%d] Read error %d\n", i,j ,k, vps, int_val);//, var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  fprintf(stderr, "C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

        else if (strcmp(type_name, INT32) == 0 || strcmp(type_name, INT32_GA) == 0 || strcmp(type_name, INT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&int_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (int_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //fprintf(stderr, "W[%d %d %d] [%d] Read error %d\n", i,j ,k, vps, int_val);//, var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  fprintf(stderr, "C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

        else if (strcmp(type_name, FLOAT64) == 0 || strcmp(type_name, FLOAT64_RGB) == 0 || strcmp(type_name, FLOAT64_GA) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&double_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (double_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 0)
              //  fprintf(stderr, "W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val,100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 0)
              //  fprintf(stderr, "C[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val,100+ var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

        else if (strcmp(type_name, FLOAT32) == 0 || strcmp(type_name, FLOAT32_GA) == 0 || strcmp(type_name, FLOAT32_RGB) == 0)
        {
          for (vps = 0; vps < values_per_sample; vps++)
          {
            memcpy(&float_val, data + (index * values_per_sample + vps) * bits_per_sample, bits_per_sample);
            if (float_val != var + 100 + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //if (rank == 100)
              //fprintf(stderr, "W [%d] [%d %d %d] [%d] Read error %f %d \n", rank, i,j ,k, vps, float_val, var + vps + 100 + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
            else
            {
              read_count++;
              //if (rank == 100)
              //  fprintf(stderr, "C [%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, float_val, 100 + var + vps + ((global_bounds[0] * global_bounds[1]*(local_box_offset[2] + k))+(global_bounds[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
          }
        }

      }

  //fprintf(stderr, "RANK %d CC %d WC %d\n", rank, read_count, read_error_count);
  MPI_Allreduce(&read_count, &total_read_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&read_error_count, &total_read_error_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0)
    fprintf(stderr, "Correct Sample Count %d Incorrect Sample Count %d\n", total_read_count, total_read_error_count);

  if(total_read_error_count == 0 && (total_read_count == global_bounds[2]*global_bounds[1]*global_bounds[0]*values_per_sample))
    return 0;
  else
    return 1;
}


//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}
