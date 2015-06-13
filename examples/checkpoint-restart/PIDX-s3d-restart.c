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

#include <PIDX.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

static int parse_args(int argc, char **argv);
static void usage(void);
static void report_error(char* func_name, char* file_name, int line_no);

static int global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static int local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 0;                        ///< Number of fields
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template

int main(int argc, char **argv)
{
  int ret;
  int i, j, k;
  int var, vps;
  int slice = 0;
  int nprocs = 1, rank = 0;
  char output_file_name[512];

  // MPI initialization
#if PIDX_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  double **data;             // variable data buffer
  int *values_per_sample;    // Example: 1 for scalar 3 for vector

  int local_box_offset[3];

  if (rank == 0)
  {
    ret = parse_args(argc, argv);
    if (ret < 0)
    {
      usage();
#if PIDX_HAVE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }

    // check if the num procs is appropriate
    int num_bricks = (global_box_size[0] / local_box_size[0]) * (global_box_size[1] / local_box_size[1]) * (global_box_size[2] / local_box_size[2]);
    if(num_bricks != nprocs)
    {
      fprintf(stderr, "Error: number of sub-blocks (%d) doesn't match number of procs (%d)\n", num_bricks, nprocs);
      fprintf(stderr, "Incorrect distribution of data across processes i.e.\n(global_x / local_x) X (global_x / local_x) X (global_x / local_x) != nprocs\n(%d/%d) X (%d/%d) X (%d/%d) != %d\n", global_box_size[0], local_box_size[0], global_box_size[1], local_box_size[1], global_box_size[2], local_box_size[2], nprocs);

#if PIDX_HAVE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }
  }

  //  The command line arguments are shared by all processes
#if PIDX_HAVE_MPI
  MPI_Bcast(global_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

  values_per_sample = malloc(sizeof(*values_per_sample) * variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);

  // Creating the filename
  sprintf(output_file_name, "%s%s", output_file_template,".idx");

  // Calculating every process data's offset and size
  int sub_div[3];
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  local_box_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_box_offset[1] = (slice / sub_div[0]) * local_box_size[1];
  local_box_offset[0] = (slice % sub_div[0]) * local_box_size[0];

  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

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
  ret = PIDX_file_open(output_file_name, PIDX_file_trunc, access, &file);
  if (ret != PIDX_success)  report_error("PIDX_file_create", __FILE__, __LINE__);

  ret = PIDX_get_dims(file, global_size);
  if (ret != PIDX_success)  report_error("PIDX_set_dims", __FILE__, __LINE__);

  ret = PIDX_get_variable_count(file, &variable_count);
  if (ret != PIDX_success)  report_error("PIDX_set_variable_count", __FILE__, __LINE__);

  ret = PIDX_set_current_time_step(file, 0);
  if (ret != PIDX_success)  report_error("PIDX_set_current_time_step", __FILE__, __LINE__);

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  data = (double**)malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  for (var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  report_error("PIDX_get_next_variable", __FILE__, __LINE__);

    values_per_sample[var] = variable[var]->values_per_sample;

    int bits_per_sample = 0;
    ret = PIDX_default_bits_per_datatype(variable[var]->type_name, &bits_per_sample);
    if (ret != PIDX_success)  report_error("PIDX_default_bytes_per_datatype", __FILE__, __LINE__);

    data[var] = malloc((bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable[var]->values_per_sample);
    memset(data[var], 0, (bits_per_sample/8) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * variable[var]->values_per_sample);
  }

  ret = PIDX_reset_variable_counter(file);
  if (ret != PIDX_success)  report_error("PIDX_reset_variable_counter", __FILE__, __LINE__);

  for (var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  report_error("PIDX_get_next_variable", __FILE__, __LINE__);

    ret = PIDX_variable_read_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
    if (ret != PIDX_success)  report_error("PIDX_variable_read_data_layout", __FILE__, __LINE__);

    ret = PIDX_read_next_variable(file, variable[var]);
    if (ret != PIDX_success)  report_error("PIDX_read_next_variable", __FILE__, __LINE__);
  }

  ret = PIDX_close(file);
  if (ret != PIDX_success)  report_error("PIDX_close", __FILE__, __LINE__);

   ret = PIDX_close_access(access);
  if (ret != PIDX_success)  report_error("PIDX_close_access", __FILE__, __LINE__);

  int read_error_count = 0;
  for(var = 0; var < variable_count; var++)
  {
    for (k = 0; k < local_box_size[2]; k++)
      for (j = 0; j < local_box_size[1]; j++)
        for (i = 0; i < local_box_size[0]; i++)
        {
          int64_t index = (int64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
          for (vps = 0; vps < values_per_sample[var]; vps++)
            if (data[var][index * values_per_sample[var] + vps] != 100 + var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)))
            {
              read_error_count++;
              //printf("[%d %d %d] Read error %f %d\n", i,j ,k, data[var][index * values_per_sample[var] + vps], 100 + var + vps + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i)));
            }
        }
  }
  printf("Read Error Count = %d Total Error Count = %d\n", read_error_count, global_box_size[0]*global_box_size[1]*global_box_size[2]);

  for(var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;

  //  MPI finalize
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
      case('g'):
          sscanf(optarg, "%dx%dx%d", &global_box_size[0], &global_box_size[1], &global_box_size[2]);
          break;
      case('l'):
          sscanf(optarg, "%dx%dx%d", &local_box_size[0], &local_box_size[1], &local_box_size[2]);
          break;
      case('f'):
          sprintf(output_file_template, "%s", optarg);
          break;
      case('t'):
          sscanf(optarg, "%d", &time_step_count);
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
