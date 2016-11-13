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
static char raw_file_template[512] = "test_idx";   ///< output IDX file Name Template
int partition_size[3] = {1, 1, 1};
static unsigned char *data;
static char raw_file_name[512] = "test.idx";
static char idx_file_name[512] = "test.idx";
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
  char flags[] = "g:l:p:f:t:v:";
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

    case('p'): // local dimension
      if ((sscanf(optarg, "%dx%dx%d", &partition_size[0], &partition_size[1], &partition_size[2]) == EOF) ||(partition_size[0] < 1 || partition_size[1] < 1 || partition_size[2] < 1))
        terminate_with_error_msg("Invalid partition dimension\n%s", usage);
      break;

    case('f'): // input file name
      if (sprintf(raw_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(raw_file_name, "%s%s", raw_file_template, ".idx");
      sprintf(idx_file_name, "%s_NEW%s", raw_file_template, ".idx");
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
  char var_name[1024];
  int ret;
  int variable_count;
  PIDX_file file_r;
  PIDX_variable variable_r;

  PIDX_point global_size, local_offset, local_size, global_bounds;
  PIDX_set_point_5D(global_size, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);
  PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);
  PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);

  PIDX_access access;
  PIDX_create_access(&access);
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);

  PIDX_file_open(raw_file_name, PIDX_MODE_RDONLY, access, global_bounds, &file_r);
  PIDX_query_box(file_r, global_size);
  PIDX_get_variable_count(file_r, &variable_count);
  PIDX_set_current_time_step(file_r, time_step_count);

  if (variable_index >= variable_count)
    terminate_with_error_msg("Variable index more than variable count\n");

  int bits_per_sample = 0;
  int v_per_sample = 0;
  PIDX_set_current_variable_index(file_r, variable_index);
  PIDX_get_current_variable(file_r, &variable_r);
  PIDX_default_bits_per_datatype(variable_r->type_name, &bits_per_sample);
  strcpy(var_name, variable_r->var_name);

  data = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable_r->vps);
  memset(data, 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable_r->vps);

  PIDX_values_per_datatype(variable_r->type_name, &v_per_sample, &bits_per_sample);
  PIDX_variable_read_data_layout(variable_r, local_offset, local_size, data, PIDX_row_major);
  PIDX_close(file_r);


  PIDX_file file_w;
  PIDX_variable variable_w;
  PIDX_point reg_patch_size;
  PIDX_set_point_5D(reg_patch_size, 512, 256, 256, 1, 1);

  int ts = 0;
  time_step_count = 1;
  for (ts = 0; ts < time_step_count; ts++)
  {
    PIDX_file_create(idx_file_name, PIDX_MODE_CREATE, access, global_size, &file_w);

    PIDX_set_current_time_step(file_w, ts);
    PIDX_set_variable_count(file_w, 1);
    PIDX_set_partition_size(file_w, partition_size[0], partition_size[1], partition_size[2]);
    PIDX_set_block_count(file_w, 512);
    PIDX_set_block_size(file_w, 16);
    PIDX_set_restructuring_box(file_w, reg_patch_size);
    PIDX_disable_agg(file_w);
    PIDX_save_little_endian(file_w);

    PIDX_variable_create(var_name,  bits_per_sample, FLOAT32 , &variable_w);
    PIDX_variable_write_data_layout(variable_w, local_offset, local_size, data, PIDX_row_major);
    ret = PIDX_append_and_write_variable(file_w, variable_w);

    PIDX_close(file_w);
  }
  PIDX_close_access(access);


  free(data);
  shutdown_mpi();

  return 0;
}
