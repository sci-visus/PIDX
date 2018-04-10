/*
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
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

/*
             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        RAW Data format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*
*/

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

#define MAX_VAR_COUNT 128
enum { X, Y, Z, NUM_DIMS };

static int process_count = 1, rank = 0;
static unsigned long long rst_box_size[NUM_DIMS];
static unsigned long long global_box_size[NUM_DIMS];
static unsigned long long local_box_offset[NUM_DIMS];
static unsigned long long local_box_size[NUM_DIMS];
static int patch_count = 1;
static unsigned long long ***var_count;
static unsigned long long ***var_offset;
int sub_div[NUM_DIMS];
static int time_step_count = 1;
static int variable_count = 1;
static char output_file_template[512];
static char var_list[512];
static unsigned char ***data;
static char output_file_name[512];
static char var_name[MAX_VAR_COUNT][512];
static int bpv[MAX_VAR_COUNT];
static char type_name[MAX_VAR_COUNT][512];
static int vps[MAX_VAR_COUNT];

static PIDX_point **local_offset_point, **local_box_count_point;
static PIDX_point global_size, local_offset, local_size, reg_size;
static PIDX_access p_access;
static PIDX_file file;
static PIDX_variable* variable;

static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static int parse_var_list();
static void check_args();
static void calculate_per_process_offsets();
static void create_synthetic_simulation_data();
static void terminate_with_error_msg(const char *format, ...);
static void terminate();
static void rank_0_print(const char *format, ...);
static void set_pidx_file(int ts);
static void set_pidx_variable();
static void create_pidx_var_point_and_access();
static void destroy_pidx_var_point_and_access();
static void destroy_synthetic_simulation_data();
static void shutdown_mpi();

static char *usage = "Serial Usage: ./multi_buffer_raw_write -g 32x32x32 -l 32x32x32 -v VL -t 4 -f -p 4 output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./multi_buffer_raw_write -g 64x64x64 -l 32x32x32 -v VL -t 4 -f -p 4 output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: Raw file name template\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables\n"
                     "  -p: number of patches per process\n";

int main(int argc, char **argv)
{
  int ts = 0, var = 0;
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();
  calculate_per_process_offsets();
  create_synthetic_simulation_data();

  rank_0_print("Simulation Data Created\n");

  create_pidx_var_point_and_access();
  for (ts = 0; ts < time_step_count; ts++)
  {
    set_pidx_file(ts + 5);
    for (var = 0; var < variable_count; var++)
      set_pidx_variable(var);
    PIDX_close(file);
  }
  destroy_pidx_var_point_and_access();

  destroy_synthetic_simulation_data();
  shutdown_mpi();

  return 0;
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
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:r:f:t:v:p:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[X], &global_box_size[Y], &global_box_size[Z]) == EOF) || (global_box_size[X] < 1 || global_box_size[Y] < 1 || global_box_size[Z] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", usage);
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[X], &local_box_size[Y], &local_box_size[Z]) == EOF) ||(local_box_size[X] < 1 || local_box_size[Y] < 1 || local_box_size[Z] < 1))
        terminate_with_error_msg("Invalid local dimension\n%s", usage);
      break;

    case('r'): // local dimension
      if ((sscanf(optarg, "%lldx%lldx%lld", &rst_box_size[X], &rst_box_size[Y], &rst_box_size[Z]) == EOF) ||(rst_box_size[X] < 1 || rst_box_size[Y] < 1 || rst_box_size[Z] < 1))
        terminate_with_error_msg("Invalid partition dimension\n%s", usage);
      break;

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // number of timesteps
      if (sscanf(optarg, "%d", &time_step_count) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('v'): // number of variables
      if (sprintf(var_list, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      parse_var_list();
      break;

    case('p'): // number of timesteps
      if (sscanf(optarg, "%d", &patch_count) < 0)
        terminate_with_error_msg("Invalid patch count\n%s", usage);
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

//----------------------------------------------------------------
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
    line[strcspn(line, "\r\n")] = 0;

    if (strcmp(line, "(fields)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      count = 0;
      variable_counter = 0;

      while (line[X] != '(')
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

//----------------------------------------------------------------
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

//----------------------------------------------------------------
static void calculate_per_process_offsets()
{
  int var = 0, i = 0, d = 0;

  //   Calculating every process's offset and count
  var_count = malloc(sizeof(unsigned long long**) * variable_count);
  var_offset = malloc(sizeof(unsigned long long**) * variable_count);

  for(var = 0; var < variable_count; var++)
  {
    var_count[var] = malloc(sizeof(unsigned long long*) * patch_count);
    var_offset[var] = malloc(sizeof(unsigned long long*) * patch_count);
    for(i = 0; i < patch_count ; i++)
    {
      var_count[var][i] = malloc(sizeof(unsigned long long) * 3);
      var_offset[var][i] = malloc(sizeof(unsigned long long) * 3);
    }
  }

#if 0
  int vi = 0;
  int pi = 0;
  int vo_x = 0, vo_y = 0, vo_z = 0, vc_x = 0, vc_y = 0, vc_z = 0;
  char fn[512];
  sprintf(fn, "/home/sid/uintah/PIDX/PIDX/build/CCVars_1_process_state_dump/rank_%d", rank);
  //printf("FN %s\n", fn);
  FILE *fp = fopen(fn, "r");
  if (fp == NULL)
    fprintf(stdout, "Error Opening %s\n", var_list);

  for (i = 0; i < 16; i++)
  {
    assert(8 == fscanf(fp, "[%d] [%d] %d %d %d %d %d %d\n", &vi, &pi, &vo_x, &vo_y, &vo_z, &vc_x, &vc_y, &vc_z));
    assert(pi == i);
    assert(vi == 0);
    for(var = 0; var < variable_count; var++)
    {
      var_offset[var][i][X] = vo_x;
      var_offset[var][i][Y] = vo_y;
      var_offset[var][i][Z] = vo_z;

      var_count[var][i][X] = vc_x;
      var_count[var][i][Y] = vc_y;
      var_count[var][i][Z] = vc_z;
      //printf("%d : %d %d %d %d %d %d\n", i, (int)var_offset[var][i][X], (int)var_offset[var][i][Y], (int)var_offset[var][i][Z], (int)var_count[var][i][X], (int)var_count[var][i][Y], (int)var_count[var][i][Z]);
    }
  }
  fclose(fp);
#else
  int sub_div[3];
  int local_box_offset[3];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];

  for(var = 0; var < variable_count; var++)
  {
    // One patch for this variable
    if (patch_count == 1)
    {
      for(d = 0; d < 3; d++)
      {
        var_count[var][X][d] = local_box_size[d];
        var_offset[var][X][d] = local_box_offset[d];
      }
    }
#if 1
    // two patches for this variable
    else if (patch_count == 2)
    {
      for(d = 0; d < 3; d++)
      {
        var_count[var][X][d] = local_box_size[d];
        var_offset[var][X][d] = local_box_offset[d];
        var_count[var][Y][d] = local_box_size[d];
        var_offset[var][Y][d] = local_box_offset[d];
      }
      var_count[var][X][X] = local_box_size[X]/2;
      if(local_box_size[X] % 2 == 0)
        var_count[var][Y][X] = local_box_size[X]/2;
      else
        var_count[var][Y][X] = local_box_size[X]/2 + 1;
      var_offset[var][Y][X] = local_box_offset[X] + local_box_size[X]/2;
    }

#else
    else if (patch_count == 2)
    {

      for(d = 0; d < 3; d++)
      {
        var_count[var][X][d] = local_box_size[d];
        var_offset[var][X][d] = local_box_offset[d];
        var_count[var][Y][d] = local_box_size[d];
        var_offset[var][Y][d] = local_box_offset[d];
      }

      var_count[var][X][X] = local_box_size[X]/2;

      if(local_box_size[X] % 2 == 0)
        var_count[var][Y][X] = local_box_size[X]/2;
      else
        var_count[var][Y][X] = local_box_size[X]/2 + 1;

      var_offset[var][Y][X] = local_box_offset[X] + local_box_size[X]/2;

      if(var_offset[var][X][X] == 0){
        var_offset[var][X][X] = local_box_offset[X] + local_box_size[X];
      }
      else if(var_offset[var][X][X] == local_box_size[X]){
        var_offset[var][X][X] = local_box_offset[X] - local_box_size[X];
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
      var_count[var][X][X] = local_box_size[X]/2;
      var_count[var][X][Y] = local_box_size[Y]/2;

      var_count[var][Y][Y] = local_box_size[Y]/2;
      if(local_box_size[X] % 2 == 0)
      {
        var_count[var][Y][X] = local_box_size[X]/2;
        var_count[var][3][X] = local_box_size[X]/2;
        var_offset[var][Y][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][3][X] = var_offset[var][X][X] + local_box_size[X]/2;
      }
      else
      {
        var_count[var][Y][X] = local_box_size[X]/2 + 1;
        var_count[var][3][X] = local_box_size[X]/2 + 1;
        var_offset[var][Y][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][3][X] = var_offset[var][X][X] + local_box_size[X]/2;
      }

      var_count[var][Z][X] = local_box_size[X]/2;
      if(local_box_size[Y] % 2 == 0)
      {
        var_count[var][Z][Y] = local_box_size[Y]/2;
        var_count[var][3][Y] = local_box_size[Y]/2;
        var_offset[var][Z][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][3][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
      }
      else
      {
        var_count[var][Z][Y] = local_box_size[Y]/2 + 1;
        var_count[var][3][Y] = local_box_size[Y]/2 + 1;
        var_offset[var][Z][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][3][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
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
      var_count[var][X][X] = local_box_size[X]/2;
      var_count[var][X][Y] = local_box_size[Y]/2;

      var_count[var][4][X] = local_box_size[X]/2;
      var_count[var][4][Y] = local_box_size[Y]/2;

      var_count[var][Y][Y] = local_box_size[Y]/2;
      var_count[var][5][Y] = local_box_size[Y]/2;

      if(local_box_size[X] % 2 == 0)
      {
        var_count[var][Y][X] = local_box_size[X]/2;
        var_count[var][3][X] = local_box_size[X]/2;
        var_offset[var][Y][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][3][X] = var_offset[var][X][X] + local_box_size[X]/2;

        var_count[var][5][X] = local_box_size[X]/2;
        var_count[var][7][X] = local_box_size[X]/2;
        var_offset[var][5][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][7][X] = var_offset[var][X][X] + local_box_size[X]/2;
      }
      else
      {
        var_count[var][Y][X] = local_box_size[X]/2 + 1;
        var_count[var][3][X] = local_box_size[X]/2 + 1;
        var_offset[var][Y][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][3][X] = var_offset[var][X][X] + local_box_size[X]/2;

        var_count[var][5][X] = local_box_size[X]/2 + 1;
        var_count[var][7][X] = local_box_size[X]/2 + 1;
        var_offset[var][5][X] = var_offset[var][X][X] + local_box_size[X]/2;
        var_offset[var][7][X] = var_offset[var][X][X] + local_box_size[X]/2;
      }

      var_count[var][Z][X] = local_box_size[X]/2;
      var_count[var][6][X] = local_box_size[X]/2;

      if(local_box_size[Y] % 2 == 0)
      {
        var_count[var][Z][Y] = local_box_size[Y]/2;
        var_count[var][3][Y] = local_box_size[Y]/2;
        var_offset[var][Z][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][3][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;

        var_count[var][6][Y] = local_box_size[Y]/2;
        var_count[var][7][Y] = local_box_size[Y]/2;
        var_offset[var][6][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][7][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
      }
      else
      {
        var_count[var][Z][Y] = local_box_size[Y]/2 + 1;
        var_count[var][3][Y] = local_box_size[Y]/2 + 1;
        var_offset[var][Z][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][3][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;

        var_count[var][6][Y] = local_box_size[Y]/2 + 1;
        var_count[var][7][Y] = local_box_size[Y]/2 + 1;
        var_offset[var][6][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
        var_offset[var][7][Y] = var_offset[var][X][Y] + local_box_size[Y]/2;
      }

      var_count[var][X][Z] = local_box_size[Z]/2;
      var_count[var][Y][Z] = local_box_size[Z]/2;
      var_count[var][Z][Z] = local_box_size[Z]/2;
      var_count[var][3][Z] = local_box_size[Z]/2;
      if(local_box_size[Y] % 2 == 0)
      {
        var_count[var][4][Z] = local_box_size[Z]/2;
        var_count[var][5][Z] = local_box_size[Z]/2;
        var_count[var][6][Z] = local_box_size[Z]/2;
        var_count[var][7][Z] = local_box_size[Z]/2;

        var_offset[var][4][Z] = var_offset[var][X][Z] + local_box_size[Z]/2;
        var_offset[var][5][Z] = var_offset[var][Y][Z] + local_box_size[Z]/2;
        var_offset[var][6][Z] = var_offset[var][Z][Z] + local_box_size[Z]/2;
        var_offset[var][7][Z] = var_offset[var][3][Z] + local_box_size[Z]/2;
      }
      else
      {
        var_count[var][4][Z] = local_box_size[Z]/2 + 1;
        var_count[var][5][Z] = local_box_size[Z]/2 + 1;
        var_count[var][6][Z] = local_box_size[Z]/2 + 1;
        var_count[var][7][Z] = local_box_size[Z]/2 + 1;

        var_offset[var][4][Z] = var_offset[var][X][Z] + local_box_size[Z]/2;
        var_offset[var][5][Z] = var_offset[var][Y][Z] + local_box_size[Z]/2;
        var_offset[var][6][Z] = var_offset[var][Z][Z] + local_box_size[Z]/2;
        var_offset[var][7][Z] = var_offset[var][3][Z] + local_box_size[Z]/2;
      }
    }
    else
      printf("This patch count not supported !!!!\n");
  }
#endif
}

//----------------------------------------------------------------
static void create_synthetic_simulation_data()
{
  int var = 0, p = 0;
  unsigned long long i, j, k, vps = 0;
  float fvalue = 0;

  data = malloc(sizeof(float**) * variable_count);
  for (var = 0; var < variable_count; var++)
  {
    int sample_count = 1;
    if ((bpv[var]) == 32)
      sample_count = 1;
    else if ((bpv[var]) == 192)
      sample_count = 3;
    else if ((bpv[var]) == 64)
      sample_count = 1;

    data[var] = malloc(sizeof(float*) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      data[var][p] = malloc(sizeof (float) * var_count[var][p][X] * var_count[var][p][Y] * var_count[var][p][Z]);
      for (k = 0; k < var_count[var][p][Z]; k++)
        for (j = 0; j < var_count[var][p][Y]; j++)
          for (i = 0; i < var_count[var][p][X]; i++)
          {
            unsigned long long index = (unsigned long long) (var_count[var][p][X] * var_count[var][p][Y] * k) + (var_count[var][p][X] * j) + i;
            for (vps = 0; vps < sample_count; vps++)
            {
              fvalue = var + 100 + (global_box_size[X] * global_box_size[Y]*(var_offset[var][p][Z] + k))+(global_box_size[X]*(var_offset[var][p][Y] + j)) + (var_offset[var][p][X] + i);

              memcpy(data[var][p] + (index * sample_count + vps) * sizeof(float), &fvalue, sizeof(float));
            }
          }
    }
  }
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
static void rank_0_print(const char *format, ...)
{
  if (rank != 0) return;
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
}

//----------------------------------------------------------------
static void create_pidx_var_point_and_access()
{
  int var, p;
  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);
  PIDX_set_point(local_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z]);
  PIDX_set_point(local_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);
  PIDX_set_point(reg_size, rst_box_size[X], rst_box_size[Y], rst_box_size[Z]);

  local_offset_point = malloc(sizeof(PIDX_point*) * variable_count);
  local_box_count_point = malloc(sizeof(PIDX_point*) * variable_count);

  for(var = 0; var < variable_count; var++)
  {
    local_offset_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    local_box_count_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      PIDX_set_point(local_offset_point[var][p], var_offset[var][p][X], var_offset[var][p][Y], var_offset[var][p][Z]);
      PIDX_set_point(local_box_count_point[var][p], var_count[var][p][X], var_count[var][p][Y], var_count[var][p][Z]);
    }
  }

  //  Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  return;
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, p_access, global_size, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  PIDX_set_current_time_step(file, ts);
  PIDX_set_variable_count(file, variable_count);

  // Seting the restructuring box size
  PIDX_set_restructuring_box(file, reg_size);

  // Selecting raw I/O mode
  PIDX_set_io_mode(file, PIDX_RAW_IO);

  PIDX_dump_state(file, PIDX_META_DATA_DUMP_ONLY);

  return;
}

//----------------------------------------------------------------
static void set_pidx_variable(int var)
{
  int p = 0;
  PIDX_return_code ret = 0;

  PIDX_variable_create(var_name[var], bpv[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  for (p = 0 ; p < patch_count ; p++)
  {
    ret = PIDX_variable_write_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], data[var][p], PIDX_row_major);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");
  }

  ret = PIDX_append_and_write_variable(file, variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  return;
}

//----------------------------------------------------------------
static void destroy_pidx_var_point_and_access()
{
  if (PIDX_close_access(p_access) != PIDX_success)
    terminate_with_error_msg("PIDX_close_access");

  free(variable);
  variable = 0;

  return;
}

//----------------------------------------------------------------
static void destroy_synthetic_simulation_data()
{
  int var = 0;
  for(var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;
}

//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}
