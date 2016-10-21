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
static unsigned long long global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static unsigned long long local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_index = 0;
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static unsigned char *data;
static char output_file_name[512] = "test.idx";
static char *usage = "Serial Usage: ./vis_read -g 32x32x32 -l 32x32x32 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./restart -g 32x32x32 -l 16x16x16 -f output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -t: time step index that needs to be read\n"
                     "  -v: Variable Index";

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
  int sub_div[NUM_DIMS];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
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
      if (sscanf(optarg, "%d", &time_step_count) < 0)
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
  if(brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);

}

int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();
  calculate_per_process_offsets();

  int ret;
  int variable_count;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable variable;   // variable descriptor

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

  //  PIDX mandatory calls
  ret = PIDX_file_open(output_file_name, PIDX_MODE_RDONLY, access, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  ret = PIDX_get_dims(file, global_size);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_dims");

  ret = PIDX_get_variable_count(file, &variable_count);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

  ret = PIDX_set_current_time_step(file, time_step_count);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");

  if (variable_index >= variable_count) terminate_with_error_msg("Variable index more than variable count\n");
  ret = PIDX_set_current_variable_index(file, variable_index);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_variable_index");

  ret = PIDX_get_current_variable(file, &variable);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_current_variable");

  int bits_per_sample = 0;
  ret = PIDX_default_bits_per_datatype(variable->type_name, &bits_per_sample);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_default_bytes_per_datatype");

  data = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable->vps);
  memset(data, 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable->vps);

  int v_per_sample = 0;
  char type_name[512];
  PIDX_values_per_datatype(variable->type_name, &v_per_sample, &bits_per_sample);
  strcpy(type_name, variable->type_name);

  ret = PIDX_variable_read_data_layout(variable, local_offset, local_size, data, PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_read_data_layout");

  ret = PIDX_close(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

  int read_error_count = 0, read_count = 0;
  /*
  int i = 0, j = 0, k = 0, vps = 0;
  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
        for (vps = 0; vps < v_per_sample; vps++)
        {
          if (data[index * v_per_sample + vps] != variable_index + vps + ((global_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            read_error_count++;
          else
            read_count++;
        }
      }
  */


  int i, j, k, vps;
  int int_val = 0;
  double double_val = 0;
  float float_val = 0;
  int var = variable_index;
  {
    bits_per_sample = bits_per_sample / 8;
    for (k = 0; k < local_box_size[2]; k++)
      for (j = 0; j < local_box_size[1]; j++)
        for (i = 0; i < local_box_size[0]; i++)
        {
          unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;

          if (strcmp(type_name, INT32) == 0)
          {
            for (vps = 0; vps < v_per_sample; vps++)
            {
              memcpy(&int_val, data + (index * v_per_sample + vps) * bits_per_sample, bits_per_sample);
              if (int_val != var + vps + 100 + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                if (rank == 0)
                  printf("W[%d %d %d] [%d] Read error %d\n", i,j ,k, vps, int_val);//, var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
            }
          }

          else if (strcmp(type_name, FLOAT64) == 0 || strcmp(type_name, FLOAT64_RGB) == 0)
          {

            for (vps = 0; vps < v_per_sample; vps++)
            {
              memcpy(&double_val, data + (index * v_per_sample + vps) * bits_per_sample, bits_per_sample);

              if (double_val != vps + 100 + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val,100 + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %f %d\n", i,j ,k, vps, double_val, var /*+ vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
              }
            }
          }

          else if (strcmp(type_name, FLOAT32) == 0)
          {
            for (vps = 0; vps < v_per_sample; vps++)
            {
              memcpy(&float_val, data + (index * v_per_sample + vps) * bits_per_sample, bits_per_sample);
              if (float_val != var + vps + 100 + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %d %f\n", i,j ,k, vps, var, float_val /*+ vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i))*/);
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[index * vps[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
            }
          }

        }
  }


  printf("Correct Sample Count %d Incorrect Sample Count %d\n", read_count, read_error_count);
  free(data);
  shutdown_mpi();

  return 0;
}
