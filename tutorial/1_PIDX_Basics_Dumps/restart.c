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
static int variable_count = 0;                        ///< Number of fields
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static unsigned char **data;
static char output_file_name[512] = "test.idx";
static char *usage = "Serial Usage: ./restart -g 32x32x32 -l 32x32x32 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./restart -g 32x32x32 -l 16x16x16 -f output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -t: the timestep to read\n";

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


static void destroy_synthetic_simulation_data()
{
  int var;
  for(var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }

  free(data);
  data = 0;

}

///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) ||
          (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
        terminate_with_error_msg("[g] Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||
          (local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
        terminate_with_error_msg("[l] Invalid local dimension\n%s", usage);
      break;

    case('f'): // input file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("[f] Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // number of timesteps
      if (sscanf(optarg, "%d", &time_step_count) < 0)
        terminate_with_error_msg("[t] Invalid variable file\n%s", usage);
      break;

    default:
      terminate_with_error_msg("[D] Wrong arguments\n%s", usage);
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

  int var;
  int ret;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

  int *v_per_sample;
  int *bits_per_sample;
  char **type_name;

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

  PIDX_debug_output(file);

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  data = malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  int v = 0;
  bits_per_sample = malloc(sizeof(*bits_per_sample) * variable_count);
  v_per_sample = malloc(sizeof(*v_per_sample) * variable_count);
  type_name = malloc(sizeof(*type_name) * variable_count);
  for (v = 0; v < variable_count; v++)
    type_name[v] = malloc(sizeof(*type_name[v]) * 512);

  /*
  PIDX_enable_raw_io(file);
  PIDX_point reg_patch_size;
  PIDX_set_point_5D(reg_patch_size, 64, 64, 64, 1, 1);
  PIDX_set_restructuring_box(file, reg_patch_size);
  */

  for (var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

    int bits_per_sample = 0;
    ret = PIDX_default_bits_per_datatype(variable[var]->type_name, &bits_per_sample);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_default_bytes_per_datatype");

    data[var] = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable[var]->values_per_sample);
    memset(data[var], 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable[var]->values_per_sample);

    ret = PIDX_read_next_variable(file, variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_read_next_variable");
  }

  ret = PIDX_reset_variable_counter(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_reset_variable_counter");

  for (var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

    ret = PIDX_variable_read_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_read_data_layout");

    ret = PIDX_read_next_variable(file, variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_read_next_variable");

    PIDX_values_per_datatype(variable[var]->type_name, &v_per_sample[var], &bits_per_sample[var]);
    strcpy(type_name[var], variable[var]->type_name);
  }

  ret = PIDX_close(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

  int i, j, k, vps;
  int read_error_count = 0, read_count = 0;
  int int_val = 0;
  double double_val = 0;
  float float_val = 0;


  for(var = 0; var < variable_count; var++)
  {
    bits_per_sample[var] = bits_per_sample[var] / 8;
    for (k = 0; k < local_box_size[2]; k++)
      for (j = 0; j < local_box_size[1]; j++)
        for (i = 0; i < local_box_size[0]; i++)
        {
          unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;

          if (strcmp(type_name[var], INT32) == 0)
          {
            for (vps = 0; vps < v_per_sample[var]; vps++)
            {
              memcpy(&int_val, data[var] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
              if (int_val != var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %d %lld\n", i,j ,k, vps, int_val, var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[var][index * values_per_sample[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
            }
          }

          else if (strcmp(type_name[var], FLOAT64) == 0 || strcmp(type_name[var], FLOAT64_RGB) == 0)
          {
            for (vps = 0; vps < v_per_sample[var]; vps++)
            {
              memcpy(&double_val, data[var] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
              if (double_val != var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %f %lld\n", i,j ,k, vps, double_val, var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[var][index * values_per_sample[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
            }
          }

          else if (strcmp(type_name[var], FLOAT32) == 0)
          {
            for (vps = 0; vps < v_per_sample[var]; vps++)
            {
              memcpy(&float_val, data[var] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
              if (float_val != var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
              {
                read_error_count++;
                //if (rank == 0)
                //  printf("W[%d %d %d] [%d] Read error %d %lld\n", i,j ,k, vps, int_val, var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
              else
              {
                read_count++;
                //if (rank == 0)
                //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[var][index * values_per_sample[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
              }
            }
          }

        }
  }
  //

  int var_sample_count = 0;
  for(var = 0; var < variable_count; var++)
    var_sample_count = var_sample_count + v_per_sample[var];

  printf("[%d] Read Error Count + Right Count : [%d + %d] Total Count = %lld\n", rank, read_error_count, read_count, global_box_size[0]*global_box_size[1]*global_box_size[2] * var_sample_count);

  free(bits_per_sample);
  free(v_per_sample);
  for (v = 0; v < variable_count; v++)
    free(type_name[v]);
  free(type_name);

  free(variable);
  variable = 0;

  destroy_synthetic_simulation_data();

  shutdown_mpi();

  return 0;
}
