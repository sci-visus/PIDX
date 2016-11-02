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
static int process_count = 1, pidx_rank = 0;
static unsigned long long global_box_size[5] = {0, 0, 0, 0, 0};         ///< global dimensions of 3D volume
static unsigned long long local_box_size[5] = {0, 0, 0, 0, 0};          ///< local dimensions of the per-process block
static unsigned long long local_box_offset[5] = {0, 0, 0, 0, 0};
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 0;                         ///< Number of fields
static int patch_count = 1;
static int ***var_count;
static int ***var_offset;
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static unsigned char ***data;
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
  if (MPI_Comm_rank(MPI_COMM_WORLD, &pidx_rank) != MPI_SUCCESS)
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
  char flags[] = "g:l:f:t:p:";
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

    case('p'): // patch count
      if (sscanf(optarg, "%d", &patch_count) < 0)
        terminate_with_error_msg("[p] Invalid patch count\n%s", usage);
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

static void calculate_per_process_offsets()
{
  //   Calculating every process's offset and count
  int var = 0, d = 0, i = 0;
  int sub_div[3];
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  local_box_offset[2] = (pidx_rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  int slice = pidx_rank % (sub_div[0] * sub_div[1]);
  local_box_offset[1] = (slice / sub_div[0]) * local_box_size[1];
  local_box_offset[0] = (slice % sub_div[0]) * local_box_size[0];

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

    // One patch for this variable
    if (patch_count == 1)
    {
      for(d = 0; d < 3; d++)
      {
        var_count[var][0][d] = local_box_size[d];
        var_offset[var][0][d] = local_box_offset[d];
      }
    }
#if 1
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

#else
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

      if(var_offset[var][0][0] == 0){
     //   printf("moving a from %d to %lld\n", var_offset[var][0][0], var_offset[var][0][0]+ local_box_size[0]);
        var_offset[var][0][0] = local_box_offset[0] + local_box_size[0];
      }
      else if(var_offset[var][0][0] == local_box_size[0]){
      //  printf("moving b from %d to %lld\n", var_offset[var][0][0], var_offset[var][0][0] - local_box_size[0]);
        var_offset[var][0][0] = local_box_offset[0] - local_box_size[0];
      }
    }
#endif

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
}


int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();

  int var, p;
  int ret;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor
  int *v_per_sample;
  int *bits_per_sample;
  char **type_name;

  PIDX_point **local_offset_point, **local_box_count_point;
  //  Creating access
  PIDX_access access;
  PIDX_create_access(&access);
#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  //  PIDX mandatory calls
  ret = PIDX_file_open(output_file_name, PIDX_MODE_RDONLY, access, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  ret = PIDX_get_dims(file, global_box_size);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_dims");

  ret = PIDX_get_variable_count(file, &variable_count);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

  ret = PIDX_set_current_time_step(file, time_step_count);
  if (ret == PIDX_err_file_exists)  terminate_with_error_msg("PIDX_file_open: file does not exist");
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");

  calculate_per_process_offsets();

  local_offset_point = malloc(sizeof(PIDX_point*) * variable_count);
  local_box_count_point = malloc(sizeof(PIDX_point*) * variable_count);

  for(var = 0; var < variable_count; var++)
  {
    local_offset_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    local_box_count_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      PIDX_set_point_5D(local_offset_point[var][p], (unsigned long long)var_offset[var][p][0], (unsigned long long)var_offset[var][p][1], (unsigned long long)var_offset[var][p][2], 0, 0);
      PIDX_set_point_5D(local_box_count_point[var][p], (unsigned long long)var_count[var][p][0], (unsigned long long)var_count[var][p][1], (unsigned long long)var_count[var][p][2], 1, 1);
    }
  }

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  int v = 0;
  bits_per_sample = malloc(sizeof(*bits_per_sample) * variable_count);
  v_per_sample = malloc(sizeof(*v_per_sample) * variable_count);
  type_name = malloc(sizeof(*type_name) * variable_count);
  for (v = 0; v < variable_count; v++)
    type_name[v] = malloc(sizeof(*type_name[v]) * 512);

  data = malloc(sizeof(**data) * variable_count);
  for (var = 0; var < variable_count; var++)
  {

    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

    int bits_per_sample = 0;
    ret = PIDX_default_bits_per_datatype(variable[var]->type_name, &bits_per_sample);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_default_bytes_per_datatype");

    data[var] = malloc(sizeof(*data) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      data[var][p] = malloc((bits_per_sample/8) * local_box_count_point[var][p][0] * local_box_count_point[var][p][1] * local_box_count_point[var][p][2]  * variable[var]->vps);
      memset(data[var][p], 0, (bits_per_sample/8) * local_box_count_point[var][p][0] * local_box_count_point[var][p][1] * local_box_count_point[var][p][2]  * variable[var]->vps);
    }

    ret = PIDX_read_next_variable(file, variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_read_next_variable");
  }

  ret = PIDX_reset_variable_counter(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_reset_variable_counter");

  for (var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

    for (p = 0 ; p < patch_count ; p++){
      ret = PIDX_variable_read_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], data[var][p], PIDX_row_major);
      if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_read_data_layout");
    }

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

#if 1
  for(var = 0; var < variable_count; var++)
  {
    bits_per_sample[var] = bits_per_sample[var] / 8;
    for (p = 0 ; p < patch_count ; p++)
    {
      for (k = 0; k < local_box_count_point[var][p][2]; k++)
        for (j = 0; j < local_box_count_point[var][p][1]; j++)
          for (i = 0; i < local_box_count_point[var][p][0]; i++)
          {
            unsigned long long index = (unsigned long long) (local_box_count_point[var][p][0] * local_box_count_point[var][p][1] * k) + (local_box_count_point[var][p][0] * j) + i;

            if (strcmp(type_name[var], INT32) == 0)
            {
              for (vps = 0; vps < v_per_sample[var]; vps++)
              {
                memcpy(&int_val, data[var][p] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
                if (int_val != var + vps + 100 + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)))
                {
                  read_error_count++;
                  //if (pidx_rank == 0)
                  //  printf("W[%d %d %d] [%d] Read error %d %lld\n", i,j ,k, vps, int_val, var + vps + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)));
                }
                else
                {
                  read_count++;
                  //if (pidx_rank == 0)
                  //  printf("C[%d %d %d] [%d] Read %f %lld\n", i,j ,k, vps, data[var][p][index * vps[var] + vps], var + vps + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)));
                }
              }
            }

            else if (strcmp(type_name[var], FLOAT64) == 0 || strcmp(type_name[var], FLOAT64_RGB) == 0)
            {
              for (vps = 0; vps < v_per_sample[var]; vps++)
              {
                memcpy(&double_val, data[var][p] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
                if (double_val != (double)(var + 100 + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i))))
                {
                  read_error_count++;
                  //if (pidx_rank == 0)
                  //  printf("[%d] W[%d %d %d] [%d] Read error %f %d index %d [%d]\n", p, i,j ,k, vps, double_val, var + 100+ ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)), index, bits_per_sample[var]);
                }
                else
                {
                  read_count++;
                  //if (pidx_rank == 0)
                  //  printf("C[%d %d %d] [%d] Read error %f %d index %d indexv %d\n", i,j ,k, vps, double_val, var + 100+ ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)), index, indexv);
                }
              }
            }

            else if (strcmp(type_name[var], FLOAT32) == 0)
            {
              for (vps = 0; vps < v_per_sample[var]; vps++)
              {
                memcpy(&float_val, data[var][p] + (index * v_per_sample[var] + vps) * bits_per_sample[var], bits_per_sample[var]);
                if (float_val != var + vps + 100 + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)))
                {
                  read_error_count++;
                  //if (pidx_rank == 0)
                  //  printf("%d ", var/* + vps + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i))*/);
                    //printf("W %d [%d %d %d] [%d] Read error %f %d\n", p, i,j ,k, vps, float_val, 100 + var + vps + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)));
                }
                else
                {
                  read_count++;
                  //printf("C %d [%d %d %d] [%d] Read error %f %d\n", p, i,j ,k, vps, float_val, 100 + var + vps + ((global_box_size[0] * global_box_size[1]*(local_offset_point[var][p][2] + k))+(global_box_size[0]*(local_offset_point[var][p][1] + j)) + (local_offset_point[var][p][0] + i)));
                }
              }
            }

          }
    }
  }
  //

  int var_sample_count = 0;
  for(var = 0; var < variable_count; var++)
    var_sample_count = var_sample_count + v_per_sample[var];

  printf("[%d] Read Error Count + Right Count : [%d + %d] Total Count = %lld\n", pidx_rank, read_error_count, read_count, global_box_size[0]*global_box_size[1]*global_box_size[2] * var_sample_count);

  free(bits_per_sample);
  free(v_per_sample);
  for (v = 0; v < variable_count; v++)
    free(type_name[v]);
  free(type_name);

  free(variable);
  variable = 0;

  destroy_synthetic_simulation_data();
#endif
  shutdown_mpi();

  return 0;
}
