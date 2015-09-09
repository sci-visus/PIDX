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

#if PIDX_HAVE_MPI
 #include <mpi.h>
#endif

static int parse_args(int argc, char **argv);
static void usage(void);
static void report_error(char* func_name, char* file_name, int line_no);

static int global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static int local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 1;                        ///< Number of fields
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static int patch_count = 1;


/// main
int main(int argc, char **argv)
{
#if 1
  int ret;
  int i, j, k;
  int d, p, var, vps, ts;
  int slice = 0;
  int nprocs = 1, rank = 0;
  char output_file_name[512];

  // MPI initialization
#if PIDX_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

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
  MPI_Bcast(&variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&patch_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

  values_per_sample = (int*) malloc(sizeof (int) * variable_count);
  for (var = 0; var < variable_count; var++)
    values_per_sample[var] = 1;

  // Creating the filename
  sprintf(output_file_name, "%s%s", output_file_template,".idx");

  //   Calculating every process's offset and count
  int sub_div[3];
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  local_box_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_box_offset[1] = (slice / sub_div[0]) * local_box_size[1];
  local_box_offset[0] = (slice % sub_div[0]) * local_box_size[0];

  int ***var_count;
  int ***var_offset;
  var_count = malloc(sizeof(int**) * variable_count);
  var_offset = malloc(sizeof(int**) * variable_count);

  for(var = 0; var < variable_count; var++)
  {
    var_count[var] = malloc(sizeof(int*) * patch_count);
    var_offset[var] = malloc(sizeof(int*) * patch_count);
    for(i = 0; i < patch_count ; i++)
    {
      var_count[var][i] = malloc(sizeof(int) * 3);
      var_offset[var][i] = malloc(sizeof(int) * 3);
    }
  }

  for(var = 0; var < variable_count; var++)
  {
    // One patch for this variable
    if (patch_count == 1)
    {
      for(d = 0; d < 3; d++)
      {
        var_count[var][0][d] = local_box_size[d];
        var_offset[var][0][d] = local_box_offset[d];
      }
    }

    // two patches for this variable
    else if (patch_count == 2)
    {
      for(d = 0; d < 3; d++)
      {
        var_count[var][0][d] = local_box_size[d];
        var_offset[var][0][d] = local_box_offset[d];
        var_count[var][1][d] = local_box_size[d];
        var_offset[var][1][d] = local_box_offset[d];
      }

      var_count[var][0][0] = local_box_size[0]/2;
      if(local_box_size[0] % 2 == 0)
        var_count[var][1][0] = local_box_size[0]/2;
      else
        var_count[var][1][0] = local_box_size[0]/2 + 1;

      var_offset[var][1][0] = local_box_offset[0] + local_box_size[0]/2;
    }

    // four patches for this variable
    else if (patch_count == 4)
    {
      for(i = 0; i < patch_count; i++)
      {
        for(d = 0; d < 3; d++)
        {
          var_count[var][i][d] = local_box_size[d];
          var_offset[var][i][d] = local_box_offset[d];
        }
      }
      var_count[var][0][0] = local_box_size[0]/2;
      var_count[var][0][1] = local_box_size[1]/2;

      var_count[var][1][1] = local_box_size[1]/2;
      if(local_box_size[0] % 2 == 0)
      {
        var_count[var][1][0] = local_box_size[0]/2;
        var_count[var][3][0] = local_box_size[0]/2;
        var_offset[var][1][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][3][0] = var_offset[var][0][0] + local_box_size[0]/2;
      }
      else
      {
        var_count[var][1][0] = local_box_size[0]/2 + 1;
        var_count[var][3][0] = local_box_size[0]/2 + 1;
        var_offset[var][1][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][3][0] = var_offset[var][0][0] + local_box_size[0]/2;
      }

      var_count[var][2][0] = local_box_size[0]/2;
      if(local_box_size[1] % 2 == 0)
      {
        var_count[var][2][1] = local_box_size[1]/2;
        var_count[var][3][1] = local_box_size[1]/2;
        var_offset[var][2][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][3][1] = var_offset[var][0][1] + local_box_size[1]/2;
      }
      else
      {
        var_count[var][2][1] = local_box_size[1]/2 + 1;
        var_count[var][3][1] = local_box_size[1]/2 + 1;
        var_offset[var][2][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][3][1] = var_offset[var][0][1] + local_box_size[1]/2;
      }
    }

    // eight patches for this variable
    else if (patch_count == 8)
    {
      for(i = 0; i < patch_count; i++)
      {
        for(d = 0; d < 3; d++)
        {
          var_count[var][i][d] = local_box_size[d];
          var_offset[var][i][d] = local_box_offset[d];
        }
      }
      var_count[var][0][0] = local_box_size[0]/2;
      var_count[var][0][1] = local_box_size[1]/2;

      var_count[var][4][0] = local_box_size[0]/2;
      var_count[var][4][1] = local_box_size[1]/2;

      var_count[var][1][1] = local_box_size[1]/2;
      var_count[var][5][1] = local_box_size[1]/2;

      if(local_box_size[0] % 2 == 0)
      {
        var_count[var][1][0] = local_box_size[0]/2;
        var_count[var][3][0] = local_box_size[0]/2;
        var_offset[var][1][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][3][0] = var_offset[var][0][0] + local_box_size[0]/2;

        var_count[var][5][0] = local_box_size[0]/2;
        var_count[var][7][0] = local_box_size[0]/2;
        var_offset[var][5][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][7][0] = var_offset[var][0][0] + local_box_size[0]/2;
      }
      else
      {
        var_count[var][1][0] = local_box_size[0]/2 + 1;
        var_count[var][3][0] = local_box_size[0]/2 + 1;
        var_offset[var][1][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][3][0] = var_offset[var][0][0] + local_box_size[0]/2;

        var_count[var][5][0] = local_box_size[0]/2 + 1;
        var_count[var][7][0] = local_box_size[0]/2 + 1;
        var_offset[var][5][0] = var_offset[var][0][0] + local_box_size[0]/2;
        var_offset[var][7][0] = var_offset[var][0][0] + local_box_size[0]/2;
      }

      var_count[var][2][0] = local_box_size[0]/2;
      var_count[var][6][0] = local_box_size[0]/2;

      if(local_box_size[1] % 2 == 0)
      {
        var_count[var][2][1] = local_box_size[1]/2;
        var_count[var][3][1] = local_box_size[1]/2;
        var_offset[var][2][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][3][1] = var_offset[var][0][1] + local_box_size[1]/2;

        var_count[var][6][1] = local_box_size[1]/2;
        var_count[var][7][1] = local_box_size[1]/2;
        var_offset[var][6][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][7][1] = var_offset[var][0][1] + local_box_size[1]/2;
      }
      else
      {
        var_count[var][2][1] = local_box_size[1]/2 + 1;
        var_count[var][3][1] = local_box_size[1]/2 + 1;
        var_offset[var][2][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][3][1] = var_offset[var][0][1] + local_box_size[1]/2;

        var_count[var][6][1] = local_box_size[1]/2 + 1;
        var_count[var][7][1] = local_box_size[1]/2 + 1;
        var_offset[var][6][1] = var_offset[var][0][1] + local_box_size[1]/2;
        var_offset[var][7][1] = var_offset[var][0][1] + local_box_size[1]/2;
      }

      var_count[var][0][2] = local_box_size[2]/2;
      var_count[var][1][2] = local_box_size[2]/2;
      var_count[var][2][2] = local_box_size[2]/2;
      var_count[var][3][2] = local_box_size[2]/2;
      if(local_box_size[1] % 2 == 0)
      {
        var_count[var][4][2] = local_box_size[2]/2;
        var_count[var][5][2] = local_box_size[2]/2;
        var_count[var][6][2] = local_box_size[2]/2;
        var_count[var][7][2] = local_box_size[2]/2;

        var_offset[var][4][2] = var_offset[var][0][2] + local_box_size[2]/2;
        var_offset[var][5][2] = var_offset[var][1][2] + local_box_size[2]/2;
        var_offset[var][6][2] = var_offset[var][2][2] + local_box_size[2]/2;
        var_offset[var][7][2] = var_offset[var][3][2] + local_box_size[2]/2;
      }
      else
      {
        var_count[var][4][2] = local_box_size[2]/2 + 1;
        var_count[var][5][2] = local_box_size[2]/2 + 1;
        var_count[var][6][2] = local_box_size[2]/2 + 1;
        var_count[var][7][2] = local_box_size[2]/2 + 1;

        var_offset[var][4][2] = var_offset[var][0][2] + local_box_size[2]/2;
        var_offset[var][5][2] = var_offset[var][1][2] + local_box_size[2]/2;
        var_offset[var][6][2] = var_offset[var][2][2] + local_box_size[2]/2;
        var_offset[var][7][2] = var_offset[var][3][2] + local_box_size[2]/2;
      }
    }
    else
      printf("This patch count not supported !!!!\n");
  }

  double   ***double_data;
  double_data = malloc(sizeof(double**) * variable_count);
  for (var = 0; var < variable_count; var++)
  {
    double_data[var] = malloc(sizeof(double*) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      double_data[var][p] = malloc(sizeof (double) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * values_per_sample[var]);
      for (k = 0; k < var_count[var][p][2]; k++)
        for (j = 0; j < var_count[var][p][1]; j++)
          for (i = 0; i < var_count[var][p][0]; i++)
          {
            int64_t index = (int64_t) (var_count[var][p][0] * var_count[var][p][1] * k) + (var_count[var][p][0] * j) + i;
            for (vps = 0; vps < values_per_sample[var]; vps++)
              double_data[var][p][index * values_per_sample[var] + vps] = 100 +                 (global_box_size[0] * global_box_size[1]*(var_offset[var][p][2] + k))+(global_box_size[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i);
          }
    }
  }

  PIDX_point global_bounding_box, **local_offset_point, **local_box_count_point;

  local_offset_point = malloc(sizeof(PIDX_point*) * variable_count);
  local_box_count_point = malloc(sizeof(PIDX_point*) * variable_count);
  for(var = 0; var < variable_count; var++)
  {
    local_offset_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    local_box_count_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      PIDX_set_point_5D(local_offset_point[var][p], (int64_t)var_offset[var][p][0], (int64_t)var_offset[var][p][1], (int64_t)var_offset[var][p][2], 0, 0);
      PIDX_set_point_5D(local_box_count_point[var][p], (int64_t)var_count[var][p][0], (int64_t)var_count[var][p][1], (int64_t)var_count[var][p][2], 1, 1);
    }
  }
  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_set_point_5D(global_bounding_box, (int64_t)global_box_size[0], (int64_t)global_box_size[1], (int64_t)global_box_size[2], 1, 1);

  PIDX_access access;
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  for (ts = 0; ts < time_step_count; ts++)
  {
    PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_variable_count(file, variable_count);

    for (var = 0; var < variable_count; var++)
    {
      char variable_name[512];
      char data_type[512];
      sprintf(variable_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);
      PIDX_variable_create(variable_name, values_per_sample[var] * sizeof(double) * 8, data_type, &variable[var]);

      for (p = 0 ; p < patch_count ; p++)
        PIDX_variable_write_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], double_data[var][p], PIDX_row_major);


      PIDX_append_and_write_variable(file, variable[var]/*, local_offset_point[var][p], local_box_count_point[var][p], double_data[var][p], PIDX_row_major*/);
    }

    PIDX_close(file);
  }

  PIDX_close_access(access);
  for (var = 0; var < variable_count; var++)
  {
    for(p = 0 ; p < patch_count ; p++)
      free(double_data[var][p]);
    free(double_data[var]);
  }
  free(double_data);  double_data = 0;

  for(var = 0; var < variable_count; var++)
  {
    free(local_offset_point[var]);
    free(local_box_count_point[var]);
  }
  free(local_offset_point);
  free(local_box_count_point);

  free(variable);
  variable = 0;
  free(values_per_sample);
  values_per_sample = 0;
#endif

  return 0;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:p:";
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
      case('p'):
          sscanf(optarg, "%d", &patch_count);
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
  printf("Serial Usage: ./pidx-uintah-checkpoint -g 4x4x4 -l 4x4x4 -f Filename -t 1 -v 1\n");
  printf("Parallel Usage: mpirun -n 8 ./pidx-s3d-checkpoint -g 4x4x4 -l 2x2x2 -f Filename_ -t 1 -v 1\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("  -v: number of fields\n");
  printf("  -p: number of patches\n");

  printf("pidx-s3d-uintah generates a 3D volume of size g_x g_y g_z specified by -g command line parameter\nEach process writes a sub-volume of size l_x l_y l_z specified by -l command line parameter. \nData is written in the idx format, with the filename Filename specified by -f command line parameter. \nThe number of time-steps and the number of fields can be optionally provided by -t and -v command line parameters.\n");

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
