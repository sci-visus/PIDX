/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

#define MAX_VAR_COUNT 128
enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;

static unsigned long long global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static unsigned long long local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 1;                        ///< Number of fields
static char output_file_template[512] = "test_idx";   ///< output IDX file Name Template
static int patch_count = 1;
static int ***var_count;
static int ***var_offset;
static float   ***float_data;
static int partition_size[3] = {1, 1, 1};
static char output_file_name[512] = "test.idx";

static char var_list[512] = "var_list";
static char var_name[MAX_VAR_COUNT][512];
static char type_name[MAX_VAR_COUNT][512];
static int bpv[MAX_VAR_COUNT];
static int vps[MAX_VAR_COUNT];

static char *usage = "Serial Usage: ./hdf5-to-idx -g 4x4x4 -l 4x4x4 -v var_list -i hdf5_file_names_list -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./hdf5-to-idx -g 4x4x4 -l 2x2x2 -f Filename_ -v var_list -i hdf5_file_names_list\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -i: file containing list of input hdf5 files\n"
                     "  -v: file containing list of input fields\n";

static int parse_var_list();
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

// //----------------------------------------------------------------
// static void rank_0_print(const char *format, ...)
// {
//   if (rank != 0) return;
//   va_list arg_ptr;
//   va_start(arg_ptr, format);
//   vfprintf(stderr, format, arg_ptr);
//   va_end(arg_ptr);
// }

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
  //   Calculating every process's offset and count
  int var = 0, d = 0, i = 0;
  int sub_div[3];
  int local_box_offset[3];
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  local_box_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  int slice = rank % (sub_div[0] * sub_div[1]);
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

static void create_synthetic_simulation_data()
{
  int var = 0, p = 0;
  unsigned long long i, j, k;
  //printf("creating data for rank %d\n", rank);

  float_data = malloc(sizeof(float**) * variable_count);
  for (var = 0; var < variable_count; var++)
  {
    float_data[var] = malloc(sizeof(float*) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      float_data[var][p] = malloc(sizeof (float) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2]);
      for (k = 0; k < var_count[var][p][2]; k++)
        for (j = 0; j < var_count[var][p][1]; j++)
          for (i = 0; i < var_count[var][p][0]; i++)
          {
            unsigned long long index = (unsigned long long) (var_count[var][p][0] * var_count[var][p][1] * k) + (var_count[var][p][0] * j) + i;
            float_data[var][p][index] = var + 100 + (global_box_size[0] * global_box_size[1]*(var_offset[var][p][2] + k))+(global_box_size[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i);
          }
    }
  }
}

static void destroy_synthetic_simulation_data()
{
  int var, p;
  for (var = 0; var < variable_count; var++)
  {
    for(p = 0 ; p < patch_count ; p++)
      free(float_data[var][p]);
    free(float_data[var]);
  }
  free(float_data);  float_data = 0;
}

///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:p:f:t:v:r:";
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

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // a file with a list of variables
      if (sscanf(optarg, "%d", &time_step_count) < 0)
        terminate_with_error_msg("Invalid time step count\n%s", usage);
      break;

    case('v'): // number of variables
      if (sprintf(var_list, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      parse_var_list();
      break;

    case('r'): // a file with a list of variables
      if (sscanf(optarg, "%d", &patch_count) < 0)
        terminate_with_error_msg("Invalid patch count\n%s", usage);
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

/// main
int main(int argc, char **argv)
{
  int i = 0;
  init_mpi(argc, argv);
  parse_args(argc, argv);

#if 1
  if (rank % 4 == 0)
    patch_count = 1;
  else if (rank % 4 == 1)
    patch_count = 2;
  else if (rank % 4 == 2)
    patch_count = 4;
  else
    patch_count = 8;
#endif

  check_args();
  calculate_per_process_offsets();
  create_synthetic_simulation_data();

#if 0
  for (i = 0; i < patch_count; i++)
  {
    printf("[%d] : [%d %d] --> %d %d %d : %d %d %d\n", rank, patch_count, i, var_offset[0][i][0], var_offset[0][i][1], var_offset[0][i][2], var_count[0][i][0], var_count[0][i][1], var_count[0][i][2]);
  }
#endif

  int var, p;
  PIDX_point global_bounding_box, **local_offset_point, **local_box_count_point;

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
  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_set_point_5D(global_bounding_box, (unsigned long long)global_box_size[0], (unsigned long long)global_box_size[1], (unsigned long long)global_box_size[2], 1, 1);

  PIDX_access access;
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  int ts;
  for (ts = 0; ts < time_step_count; ts++)
  {
    PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, global_bounding_box, &file);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_variable_count(file, variable_count);
    PIDX_set_block_count(file, 128);
    PIDX_set_io_mode(file, PIDX_RAW_IO);

    PIDX_point rst_box;
    PIDX_set_point_5D(rst_box, 64,64,64,1,1);
    PIDX_set_restructuring_box(file, rst_box);

    PIDX_set_partition_size(file, partition_size[0], partition_size[1], partition_size[2]);


    for (var = 0; var < variable_count; var++)
    {
      char variable_name[512];
      sprintf(variable_name, "variable_%d", var);
      PIDX_variable_create(variable_name, sizeof(float) * 8, FLOAT32, &variable[var]);

      for (p = 0 ; p < patch_count ; p++)
        PIDX_variable_write_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], float_data[var][p], PIDX_row_major);

      PIDX_append_and_write_variable(file, variable[var]/*, local_offset_point[var][p], local_box_count_point[var][p], float_data[var][p], PIDX_row_major*/);
    }

    PIDX_close(file);
  }
  PIDX_close_access(access);

  destroy_synthetic_simulation_data();

  for(var = 0; var < variable_count; var++)
  {
    free(local_offset_point[var]);
    free(local_box_count_point[var]);
  }
  free(local_offset_point);
  free(local_box_count_point);


  for(var = 0; var < variable_count; var++)
  {
    for(i = 0; i < patch_count ; i++)
    {
      free(var_count[var][i]);
      free(var_offset[var][i]);
    }
    free(var_count[var]);
    free(var_offset[var]);
  }
  free(var_count);
  free(var_offset);

  free(variable);
  variable = 0;

  shutdown_mpi();

  return 0;
}


static int parse_var_list()
{
  FILE *fp = fopen(var_list, "r");
  if (fp == NULL)
  {
    fprintf(stdout, "Error Opening %s\n", var_list);
    return PIDX_err_file;
  }

  int variable_counter = 0, count = 0, len = 0;
  char *pch1;
  char line [ 512 ];

  while (fgets(line, sizeof (line), fp) != NULL)
  {
    //printf("%s", line);
    line[strcspn(line, "\r\n")] = 0;

    if (strcmp(line, "(fields)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      count = 0;
      variable_counter = 0;

      while (line[0] != '(')
      {
        pch1 = strtok(line, " +");
        while (pch1 != NULL)
        {
          if (count == 0)
          {
            char* temp_name = strdup(pch1);
            strcpy(var_name[variable_counter], temp_name);
            free(temp_name);
          }

          if (count == 1)
          {
            len = strlen(pch1) - 1;
            if (pch1[len] == '\n')
              pch1[len] = 0;

            strcpy(type_name[variable_counter], pch1);
            int ret;
            int bits_per_sample = 0;
            ret = PIDX_default_bits_per_datatype(type_name[variable_counter], &bits_per_sample);
            if (ret != PIDX_success)  return PIDX_err_file;

            bpv[variable_counter] = bits_per_sample;
            vps[variable_counter] = 1;
          }
          count++;
          pch1 = strtok(NULL, " +");
        }
        count = 0;

        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        variable_counter++;
      }
      variable_count = variable_counter;
    }
  }
  fclose(fp);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    int v = 0;
    for(v = 0; v < variable_count; v++)
      printf("[%d] -> %s %d %d\n", v, var_name[v], bpv[v], vps[v]);
  }

  return PIDX_success;
}
