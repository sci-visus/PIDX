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

#include <sys/stat.h> 
#include <fcntl.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>

#include <PIDX.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

#define PIDX_IO 1

static int parse_args(int argc, char **argv);
static void usage(void);
#if PIDX_IO
static void report_error(char* func_name, char* file_name, int line_no);
#endif
static void delete_buffers();

static unsigned long long global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static unsigned long long local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 1;                        ///< Number of fields
static char **var_name;
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static char var_file[512];
static char input_file[512];
static char **file_name;
static double **buffer;
static int *values_per_sample;    // Example: 1 for scalar 3 for vector
static uint64_t restructured_box_size[5] = {32, 32, 32, 1, 1};

int main(int argc, char **argv)
{
  int ret;
  int t;
  int var;
  int slice = 0;
  int nprocs = 1, rank = 0;
  char output_file_name[1024];
  unsigned long long local_box_offset[3];

  // MPI initialization
#if PIDX_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

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
    fprintf(stderr, "Incorrect distribution of data across processes i.e.\n(global_x / local_x) X (global_x / local_x) X (global_x / local_x) != nprocs\n(%d/%d) X (%d/%d) X (%d/%d) != %d\n", (int)global_box_size[0], (int)local_box_size[0], (int)global_box_size[1], (int)local_box_size[1], (int)global_box_size[2], (int)local_box_size[2], nprocs);

#if PIDX_HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(-1);
#endif
  }

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

#if PIDX_IO
  PIDX_point global_bounding_box, local_offset, local_size;
  PIDX_set_point(global_bounding_box, global_box_size[0], global_box_size[1], global_box_size[2]);
  PIDX_set_point(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2]);
  PIDX_set_point(local_size, local_box_size[0], local_box_size[1], local_box_size[2]);

  PIDX_file file;
  PIDX_access access;
  PIDX_variable *variable;
#endif

#if PIDX_IO
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif
#endif


  for (t = 0; t < time_step_count; t++)
  {
    int fp = open(file_name[t], O_RDONLY);
    for(var = 0; var < variable_count; var++)
    {
      values_per_sample[var] =  1;
      int64_t variable_offset = var * global_box_size[0] * global_box_size[1] * global_box_size[2] * sizeof(double);
      buffer[var] = malloc(sizeof (uint64_t) * local_box_size[0] * local_box_size[1] * local_box_size[2]  * values_per_sample[var]);
      pread(fp, buffer[var], local_box_size[0] * local_box_size[1] * local_box_size[2] * sizeof(double), variable_offset + (rank * local_box_size[0] * local_box_size[1] * local_box_size[2] * sizeof(double)));

      /*
      const double pi = acos(-1.0);
      for (k = 0; k < local_box_size[2]; k++)
        for (j = 0; j < local_box_size[1]; j++)
          for (i = 0; i < local_box_size[0]; i++)
          {
            int64_t index = (int64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
            for (spv = 0; spv < values_per_sample[var]; spv++)
              double_data[var][index * values_per_sample[var] + spv] = 100 + var + ((global_box_size[0] * global_box_size[1]*(local_offset[2] + k))+(global_box_size[0]*(local_offset[1] + j)) + (local_offset[0] + i));
          }
      */
    }
    close(fp);

#if PIDX_IO
    variable = malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);
    PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, global_bounding_box, &file);
    PIDX_set_variable_count(file, variable_count);
    PIDX_set_current_time_step(file, t);

	PIDX_point rst_box_size;
	PIDX_set_point(rst_box_size, restructured_box_size[0], restructured_box_size[1],
			restructured_box_size[2]);
    PIDX_set_restructuring_box(file, rst_box_size);

	PIDX_set_io_mode(file, PIDX_RAW_IO);
    for(var = 0; var < variable_count; var++)
    {
      ret = PIDX_variable_create(var_name[var], sizeof(double) * 8, "1*float64", &variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_variable_data_layout", __FILE__, __LINE__);

      ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, buffer[var], PIDX_column_major);
      if (ret != PIDX_success)  report_error("PIDX_variable_data_layout", __FILE__, __LINE__);

      ret = PIDX_append_and_write_variable(file, variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_append_and_write_variable", __FILE__, __LINE__);
    }
    PIDX_close(file);
    free(variable);
#endif

  }

#if PIDX_IO
  PIDX_close_access(access);
#endif

  delete_buffers();

#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

static void delete_buffers()
{
  int i;

  for (i = 0; i < variable_count; i++)
  {
    free(buffer[i]);
    buffer[i] = 0;
  }
  free(buffer);
  buffer = 0;

  free(values_per_sample);
  values_per_sample = 0;

  for (i = 0; i < variable_count; i++)
  {
    free(var_name[i]);
    var_name[i] = 0;
  }
  free(var_name);
  var_name = 0;

  for (i = 0; i < time_step_count; i++)
  {
    free(file_name[i]);
    file_name[i] = 0;
  }
  free(file_name);
  file_name = 0;

  return;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  int ret, i;
  char flags[] = "g:l:f:i:v:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
      case('g'):
          sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]);
          break;
      case('l'):
          sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]);
          break;
      case('f'):
          sprintf(output_file_template, "%s", optarg);
          break;
      case('i'):
          sprintf(input_file, "%s", optarg);
          break;
      case('v'):
          sprintf(var_file, "%s", optarg);
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

  FILE *fp = fopen(var_file, "r");
  ret = fscanf(fp, "%d", &variable_count);
  if (ret != EOF && ret != 1)
    return (-1);

  var_name = malloc(sizeof(char*) * variable_count);
  memset(var_name, 0, sizeof(char*) * variable_count);

  buffer = malloc(sizeof(double*) * variable_count);
  memset(buffer, 0, sizeof(double*) * variable_count);
  for (i = 0; i < variable_count; i++)
  {
    buffer[i] = malloc(sizeof(double) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
    memset(buffer[i], 0, sizeof(double) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
  }

  values_per_sample = malloc(sizeof(*values_per_sample) * variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);

  for (i = 0; i < variable_count; i++)
  {
    char temp_var_name[1024];
    ret = fscanf(fp, "%s %d", temp_var_name, &values_per_sample[i]);
    if (ret != 2 || ret == EOF)
      return (-1);
    var_name[i] = strdup(temp_var_name);
  }
  fclose(fp);


  fp = fopen(input_file, "r");
  ret = fscanf(fp, "%d", &time_step_count);
  if (ret != EOF && ret != 1)
    return (-1);
  file_name = malloc(sizeof(char*) * time_step_count);
  memset(file_name, 0, sizeof(char*) * time_step_count);

  for (i = 0; i < time_step_count; i++)
  {
    char temp_file_name[1024];
    ret = fscanf(fp, "%s", temp_file_name);
    if (ret != 1 || ret == EOF)
      return (-1);
    file_name[i] = strdup(temp_file_name);
  }
  fclose(fp);

  return (0);
}


///< How to use this progam
static void usage(void)
{
  printf("Serial Usage: ./hdf5-to-idx -g 4x4x4 -l 4x4x4 -v var_list -i input_file_list -f output_idx_file_name\n");
  printf("Parallel Usage: mpirun -n 8 ./hdf5-to-idx -g 4x4x4 -l 2x2x2 -f Filename_ -t 1 -v 1\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -i: list of input\n");
  printf("  -v: list of input fields\n");

  printf("\n");

  return;
}

#if PIDX_IO
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
#endif
