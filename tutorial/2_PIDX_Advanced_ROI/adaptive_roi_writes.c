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

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static unsigned long long global_box_size[3] = {0, 0, 0};
static unsigned long long local_box_offset[3];
static unsigned long long local_box_size[3] = {0, 0, 0};
static int time_step_count = 1;
static char output_file_template[512] = "test";
static double *data;
static char output_file_name[512] = "test.idx";
static PIDX_point global_size, local_offset, local_size;
static int roi_type = 0;
static int current_time_step = 0;
static int reduced_resolution = 2;
static char *usage = "Serial Usage: ./adaptive_roi_writes -g 32x32x128 -l 32x32x128 -f file_name -r 4\n"
                     "Parallel Usage: mpirun -n 32 ./adaptive_roi_writes -g 32x32x128 -l 32x32x4 -f file_name -r 4\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -r: ROI type (1,2,4)\n";


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
static void rank_0_print(const char *format, ...)
{
  if (rank != 0) return;
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
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

static void calculate_per_process_offset()
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



///< Parse the input argumencurrent_time_step
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:r:";
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

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('r'): // ROI type
      if (sscanf(optarg, "%d", &roi_type) < 0)
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

static void set_ROI_extents(PIDX_file file)
{
  switch (roi_type)
  {
    case 0:
      break;

    case 1:
      PIDX_set_resolution(file, 0, reduced_resolution);
      break;

    case 2:
      if (rank % 3 == 0)
      {
        local_box_size[0] = 0;
        local_box_size[1] = 0;
        local_box_size[2] = 0;
      }
      break;

    case 3:
      if (rank % 4 == 0)
      {
        local_box_size[0] = local_box_size[0];
        local_box_size[1] = local_box_size[1];
        local_box_size[2] = local_box_size[2];
      }
      else if (rank % 4 == 1)
      {
        local_box_offset[1] = (3 * local_box_offset[1]) / 4;
        local_box_size[1] = (3 * local_box_size[1]) / 4;
      }
      else if (rank % 4 == 2)
      {
        local_box_offset[1] = (2 * local_box_offset[1]) / 4;
        local_box_size[1] = (2 * local_box_size[1]) / 4;
      }
      else if (rank % 4 == 3)
      {
        local_box_offset[1] = (1 * local_box_offset[1]) / 4;
        local_box_size[1] = (1 * local_box_size[1]) / 4;
      }
      //printf ("%d: %d %d %d\n", rank, local_box_offset[0], local_box_offset[1], local_box_offset[2]);
      break;

    case 4:
      if (rank >= 0 && rank < 16)
      {
        local_box_size[0] = local_box_size[0];
        local_box_size[1] = local_box_size[1];
        local_box_size[2] = local_box_size[2];
      }
      else if (rank >= 16 && rank < 32)
      {
        local_box_offset[1] = local_box_offset[1] / 2;
        local_box_size[1] = local_box_size[1] / 2;
      }
      else if (rank >= 32 && rank < 48)
      {
        local_box_offset[1] = local_box_offset[1] / 4;
        local_box_size[1] = local_box_size[1] / 4;
      }
      else if (rank >= 48 && rank < 64)
      {
        local_box_offset[1] = local_box_offset[1] / 8;
        local_box_size[1] = local_box_size[1] / 8;
      }
      break;

    case 5:
      if (rank != current_time_step * (global_box_size[X] / local_box_size[X]) * (global_box_size[Y] / local_box_size[Y]))
      {
        PIDX_set_point_5D(local_size, 0, 0, 0, 1, 1);
      }
      else
        PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
      break;

    default:
      break;
  }
}

static void create_synthetic_simulation_data()
{
  int i = 0, j = 0, k = 0;
  const double pi = acos(-1.0);

  data = malloc(sizeof (unsigned long long) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
  for (k = 0; k < local_box_size[2]; k++)
    for (j = 0; j < local_box_size[1]; j++)
      for (i = 0; i < local_box_size[0]; i++)
      {
        unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
        data[index] = 100;
        if (roi_type == 0)
          data[index] = cos(2 * pi * i / local_box_size[0]) * cos(2 * pi * j / local_box_size[1]) * cos(2 * pi * k / local_box_size[2]);
      }
}

int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();
  calculate_per_process_offset();

  rank_0_print("Simulation Data Created\n");

  int ret;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable variable;   // variable descriptor

  //  Creating access
  PIDX_access access;
  PIDX_create_access(&access);
#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  PIDX_set_point_5D(global_size, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);

  if (roi_type == 5)
    time_step_count = global_box_size[Z]/local_box_size[Z];

  for (current_time_step = 0; current_time_step < time_step_count; current_time_step++)
  {
    ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, global_size, &file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

    //PIDX_point reg_patch_size;
    //PIDX_set_point_5D(reg_patch_size, 32, 32, 32, 1, 1);
    //PIDX_set_restructuring_box(file, reg_patch_size);

    set_ROI_extents(file);
    create_synthetic_simulation_data();
    PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);
    if (roi_type != 5)
      PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);

    PIDX_set_block_size(file, 12);

    ret = PIDX_set_current_time_step(file, current_time_step);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");

    ret = PIDX_set_variable_count(file, 1);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

    ret = PIDX_set_ROI_type(file, 1);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_ROI_writes");

    ret = PIDX_variable_create("ROI_Var", sizeof(unsigned long long) * 8, FLOAT64, &variable);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

    ret = PIDX_variable_write_data_layout(variable, local_offset, local_size, data, PIDX_row_major);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_data_layout");

    ret = PIDX_append_and_write_variable(file, variable);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

    ret = PIDX_close(file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");
  }

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

  //free(data);
  shutdown_mpi();

  return 0;
}
