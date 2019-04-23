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

/*
             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        Global Partition IDX Data format
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
static unsigned long long global_box_size[NUM_DIMS];
static unsigned long long local_box_offset[NUM_DIMS];
static unsigned long long local_box_size[NUM_DIMS];
static int partition_box_size[NUM_DIMS];
int sub_div[NUM_DIMS];
static int time_step_count = 1;
static int variable_count = 1;
static int bit_string_type = 0;
static int blocks_per_file = 256;
static char output_file_template[512];
static char raw_file_name[512];
static char var_list[512];
static unsigned char **data;
static unsigned char* sub_sampled_buffer;
static char output_file_name[512];
static char output_file_name2[512];
static char var_name[MAX_VAR_COUNT][512];
static int bpv[MAX_VAR_COUNT];
static char type_name[MAX_VAR_COUNT][512];
static int vps[MAX_VAR_COUNT];

static int sf = 1;
static float bit_rate1 = 0;
static float bit_rate2 = 0;
static int wavelet_type = 0;
static int wavelet_level = WAVELET_STENCIL;

static PIDX_point global_size, local_offset, local_size;

static PIDX_access p_access;
static PIDX_file file;
static PIDX_variable* variable;

static PIDX_access p_access2;
static PIDX_file file2;
static PIDX_variable* variable2;

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
static void set_pidx_file2(int ts);
static void set_pidx_variable();
static void set_pidx_variable2();
static void create_pidx_var_point_and_access();
static void create_pidx_var_point_and_access2();
static void destroy_pidx_var_point_and_access();
static void destroy_pidx_var_point_and_access2();
static void destroy_synthetic_simulation_data();
static void shutdown_mpi();

static char *usage = "Serial Usage: ./single_buffer_global_partitioned_idx_write -g 32x32x32 -l 32x32x32 -v VL -t 4 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./single_buffer_global_partitioned_idx_write -g 64x64x64 -l 32x32x32 -v VL -t 4 -f output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX file name template\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables\n";

int main(int argc, char **argv)
{
  int ts = 0, var = 0;
  init_mpi(argc, argv);
  parse_args(argc, argv);
  check_args();
  calculate_per_process_offsets();
  create_synthetic_simulation_data();
  double s_time = 0, e_time = 0;

  MPI_Barrier(MPI_COMM_WORLD);

#if 1
  double sampling_start = MPI_Wtime();
  int **tpatch;
  int **allign_offset;
  int **allign_count;
  int **nsamples_per_level;
  tpatch = (int**) malloc(2 * sizeof (int*));
  memset(tpatch, 0, 2 * sizeof (int*));
  tpatch[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  tpatch[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[0], 0, PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[1], 0, PIDX_MAX_DIMENSIONS * sizeof (int));

  int maxh = 28;
  allign_offset = malloc(sizeof (int*) * maxh);
  allign_count = malloc(sizeof (int*) * maxh);
  memset(allign_offset, 0, sizeof (int*) * maxh);
  memset(allign_count, 0, sizeof (int*) * maxh);

  nsamples_per_level = malloc(sizeof (int*) * maxh);
  memset(nsamples_per_level, 0, sizeof (int*) * maxh);

  int j = 0, p = 0;
  for (j = 0; j < maxh; j++)
  {
    allign_offset[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    memset(allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

    allign_count[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    memset(allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

    nsamples_per_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    memset(nsamples_per_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
  }

  tpatch[0][0] = local_box_offset[0];
  tpatch[0][1] = local_box_offset[1];
  tpatch[0][2] = local_box_offset[2];

  tpatch[1][0] = local_box_offset[0] + local_box_size[0] - 1;
  tpatch[1][1] = local_box_offset[1] + local_box_size[1] - 1;
  tpatch[1][2] = local_box_offset[2] + local_box_size[2] - 1;

  int i1 = 0, j1 = 0, k1 = 0;
  int count = 0;

  if (sf == 1)
  {
    sub_sampled_buffer = malloc(sizeof(float) * local_box_size[X] * local_box_size[Y] * local_box_size[Z]);
    memcpy(sub_sampled_buffer, data[0], sizeof(float) * local_box_size[X] * local_box_size[Y] * local_box_size[Z]);
  }
  else
  {
    char bitSequence[512] = "V221100210210210210210210210";
    char bitPattern[512];
    int level = maxh - (int)log2(sf) * 3;
    int i = 0;

    for (i = 0; i <= maxh; i++)
      bitPattern[i] = RegExBitmaskBit(bitSequence, i);
    Align((maxh - 1), level, bitPattern, tpatch, allign_offset, allign_count, nsamples_per_level);
    //if (rank == 0)
    //  fprintf(stderr, "[%s %d] [%d %d %d] XXXXX %d %d %d\n", bitSequence, maxh, allign_offset[level][0], allign_offset[level][1], allign_offset[level][2], nsamples_per_level[level][0], nsamples_per_level[level][1], nsamples_per_level[level][2]);

    sub_sampled_buffer = malloc(sizeof(float) * nsamples_per_level[level][0] * nsamples_per_level[level][1] * nsamples_per_level[level][2]);
    for (k1 = allign_offset[level][2] - local_box_offset[2]; k1 < local_box_size[2]; k1 = k1 + sf)
    {
      for (j1 = allign_offset[level][1] -  local_box_offset[1]; j1 < local_box_size[1]; j1 = j1 + sf)
      {
        for (i1 = allign_offset[level][0] - local_box_offset[0]; i1 < local_box_size[0]; i1 = i1 + sf)
        {
          // (r_p->size[0] * r_p->size[1] * k1) + (r_p->size[0] * j1) + i1
          ((float*)sub_sampled_buffer)[count] = ((float*)data[0])[(local_box_size[0] * local_box_size[1] * k1) + (local_box_size[0] * j1) + i1];
          assert (count < nsamples_per_level[level][0] * nsamples_per_level[level][1] * nsamples_per_level[level][2]);
          count++;
        }
      }
    }

    global_box_size[0] = global_box_size[0] / sf;
    global_box_size[1] = global_box_size[1] / sf;
    global_box_size[2] = global_box_size[2] / sf;

    //partition_box_size[0] = global_box_size[0];
    //partition_box_size[1] = global_box_size[1];
    //partition_box_size[2] = global_box_size[2];

    local_box_offset[0] = local_box_offset[0] / sf;
    local_box_offset[1] = local_box_offset[1] / sf;
    local_box_offset[2] = local_box_offset[2] / sf;

    local_box_size[0] = nsamples_per_level[level][0];
    local_box_size[1] = nsamples_per_level[level][1];
    local_box_size[2] = nsamples_per_level[level][2];
  }

  free(tpatch[0]);
  free(tpatch[1]);
  free(tpatch);

  for (j = 0; j < maxh; j++)
  {
    free(allign_offset[j]);
    free(allign_count[j]);
    free(nsamples_per_level[j]);
  }
  free(allign_offset);
  free(allign_count);
  free(nsamples_per_level);

  double sampling_end = MPI_Wtime();

#endif


#if 1
  //rank_0_print("Simulation Data Created\n");


  create_pidx_var_point_and_access();
  //create_pidx_var_point_and_access2();
  for (ts = 0; ts < time_step_count; ts++)
  {
    s_time = MPI_Wtime();
    set_pidx_file(ts);
    for (var = 0; var < variable_count; var++)
      set_pidx_variable(var);
    PIDX_close(file);

    //MPI_Barrier(MPI_COMM_WORLD);

    //set_pidx_file2(ts);
    //for (var = 0; var < variable_count; var++)
    //  set_pidx_variable2(var);
    //PIDX_close(file2);

    e_time = MPI_Wtime();
    double max_time;
    double t_time = e_time - s_time;
    MPI_Allreduce(&t_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //if (max_time == t_time)
    //  fprintf(stderr, "[%d] [%d] [%d] [%d %d %d | %d %d %d] SAMPLING %f BOX %f SETUP %f RST %f INSITU %f CLEANUP %f TOTAL TIME %f [%f] = [%f]\n", ts, rank, sf, global_box_size[0], global_box_size[1], global_box_size[2], partition_box_size[0] * local_box_size[0], partition_box_size[1] * local_box_size[1], partition_box_size[2] * local_box_size[2], (sampling_end - sampling_start), (pa2 - pa1), (pa3 - pa2), (pa4 - pa3), (pa5 - pa4), (pa6 - pa5), max_time, ((pa2 - pa1) + (pa3 - pa2) + (pa4 - pa3) + (pa5 - pa4) + (pa6 - pa5)), max_time + (sampling_end - sampling_start));
  }
  destroy_pidx_var_point_and_access();
  //destroy_pidx_var_point_and_access2();

  free(sub_sampled_buffer);
#endif
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
  char flags[] = "g:l:c:f:t:v:b:a:w:x:p:q:r:s:";
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

    case('c'): // partition box dimension
      if ((sscanf(optarg, "%dx%dx%d", &partition_box_size[X], &partition_box_size[Y], &partition_box_size[Z]) == EOF) ||(partition_box_size[X] < 1 || partition_box_size[Y] < 1 || partition_box_size[Z] < 1))
        terminate_with_error_msg("Invalid local dimension\n%s", usage);
      break;

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s%s", output_file_template, "_63.idx");
      sprintf(output_file_name2, "%s%s", output_file_template, "_avg.idx");
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

    case('b'): // bit string type
      if (sscanf(optarg, "%d", &bit_string_type) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('a'): // blocks per file
      if (sscanf(optarg, "%d", &blocks_per_file) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('w'): // wavelet type
      if (sscanf(optarg, "%d", &wavelet_type) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('x'): // wavelet level
      if (sscanf(optarg, "%d", &wavelet_level) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('p'): // zfp bit rate
      if (sscanf(optarg, "%f", &bit_rate1) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('q'): // zfp bit rate
      if (sscanf(optarg, "%f", &bit_rate2) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    case('r'): // output file name
      if (sprintf(raw_file_name, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      break;

    case('s'): // output file name
      if (sscanf(optarg, "%d", &sf) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
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
    fprintf(stderr, "Error Opening %s\n", var_list);
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

  /*
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    int v = 0;
    for(v = 0; v < variable_count; v++)
      fprintf(stderr, "[%d] -> %s %d %d\n", v, var_name[v], bpv[v], vps[v]);
  }
  */

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
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
}

//----------------------------------------------------------------
static void create_synthetic_simulation_data()
{
  int i1, j1, k1;
  int index = 0;
  int var = 0;
  data = malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  // Synthetic simulation data
  for(var = 0; var < variable_count; var++)
  {
    int sample_count = 1;
    unsigned long long i, j, k, vps = 0;
    if ((bpv[var]) == 32)
      sample_count = 1;
    else if ((bpv[var]) == 192)
      sample_count = 3;
    else if ((bpv[var]) == 64)
      sample_count = 1;

    data[var] = malloc(sizeof (*(data[var])) * local_box_size[X] * local_box_size[Y] * local_box_size[Z] * (bpv[var]/8));

#if 0
    float fvalue = 0;
    double dvalue = 0;
    for (k = 0; k < local_box_size[Z]; k++)
      for (j = 0; j < local_box_size[Y]; j++)
        for (i = 0; i < local_box_size[X]; i++)
        {
          unsigned long long index = (unsigned long long) (local_box_size[X] * local_box_size[Y] * k) + (local_box_size[X] * j) + i;

          for (vps = 0; vps < sample_count; vps++)
          {
            if ((bpv[var]) == 32)
            {
              fvalue = 100 + var + vps + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              //fvalue = (local_box_offset[X] + i) * (local_box_offset[X] + i) + (local_box_offset[Y] + j) * (local_box_offset[Y] + j) + (local_box_offset[Z] + k) * (local_box_offset[Z] + k);
              //fvalue = (local_box_offset[X] + i) * (local_box_offset[X] + i) + (local_box_offset[Y] + j) * (local_box_offset[Y] + j) + (local_box_offset[Z] + k) * (local_box_offset[Z] + k);
              memcpy(data[var] + (index * sample_count + vps) * sizeof(float), &fvalue, sizeof(float));
            }

            else if ((bpv[var]) == 64)
            {
              dvalue = 100 + var + vps + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              memcpy(data[var] + (index * sample_count + vps) * sizeof(double), &dvalue, sizeof(double));
            }

            else if ((bpv[var]) == 192)
            {
              dvalue = 100 + var + vps + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              memcpy(data[var] + (index * sample_count + vps) * sizeof(double), &dvalue, sizeof(double));
            }
          }
        }
#else
#if 1
    //float* temp_buffer = malloc(local_box_size[0] * sizeof(*temp_buffer));
    //memset(temp_buffer, 0, local_box_size[0] * sizeof(*temp_buffer));

    //
    //int fp = open("magnetic-512-volume.raw", O_RDONLY);
    if (strcmp(raw_file_name, "sine") == 0)
    {
      for (k1 = 0; k1 < local_box_size[Z]; k1++)
      {
        for (j1 = 0; j1 < local_box_size[Y]; j1++)
        {
          for (i1 = 0; i1 < local_box_size[X]; i1 ++)
          {
            index = (local_box_size[0]* local_box_size[1] * k1) +
                    (local_box_size[0] * j1) +
                     i1;
            float i = (float) i1;
            float j = (float) j1;
            float k = (float) k1;
            const float PI = 3.1415926535897932384;
            float PI_2 = 2 * PI;
            float xx = i == -1 ? 1 : sin(PI_2 * i / 256);
            float yy = j == -1 ? 1 : sin(PI_2 * j / 256);
            float zz = k == -1 ? 1 : sin(PI_2 * k / 256);
            //return xx * yy * zz;

            float x = xx * yy * zz;
            memcpy(data[0] + index * sizeof(float), &x, sizeof (float));
          }
        }
      }
    }
    else
    {
      int fp = open(raw_file_name, O_RDONLY | O_BINARY);
      int send_o = 0;
      int send_c = 0;
      int recv_o = 0;

      for (k1 = local_box_offset[Z]; k1 < local_box_offset[Z] + local_box_size[Z]; k1++)
      {
        for (j1 = local_box_offset[Y]; j1 < local_box_offset[Y] + local_box_size[Y]; j1++)
        {
          for (i1 = local_box_offset[X]; i1 < local_box_offset[X] + local_box_size[X]; i1 = i1 + local_box_size[X])
          {
            index = (local_box_size[0]* local_box_size[1] * (k1 - local_box_offset[2])) +
                    (local_box_size[0] * (j1 - local_box_offset[1])) +
                    (i1 - local_box_offset[0]);
            send_o = index;
            send_c = local_box_size[0];
            recv_o = (global_box_size[X] * global_box_size[Y] * k1) + (global_box_size[X] * j1) + i1;

            pread(fp, data[0] + send_o * sizeof(float), send_c * sizeof(float), recv_o * sizeof(float));
          }
        }
      }
      close(fp);
    }
#endif
#endif
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
  variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);
  PIDX_set_point(local_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z]);
  PIDX_set_point(local_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);

  //  Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  return;
}

static void create_pidx_var_point_and_access2()
{
  variable2 = (PIDX_variable*)malloc(sizeof(*variable2) * variable_count);
  memset(variable2, 0, sizeof(*variable2) * variable_count);

  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);
  PIDX_set_point(local_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z]);
  PIDX_set_point(local_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);

  //  Creating access
  PIDX_create_access(&p_access2);
  PIDX_set_mpi_access(p_access2, MPI_COMM_WORLD);

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

  //PIDX_debug_disable_agg(file);

  PIDX_point reg_size;
  PIDX_set_point(reg_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);
  PIDX_set_restructuring_box(file, reg_size);

  PIDX_set_process_decomposition(file, global_box_size[X]/local_box_size[X], global_box_size[Y]/local_box_size[Y], global_box_size[Z]/local_box_size[Z]);
  //fprintf(stderr, "XYZ %lld / %lld   %lld / %lld    %lld / %lld\n", global_box_size[X], local_box_size[X], global_box_size[Y], local_box_size[Y], global_box_size[Z], local_box_size[Z]);
  //PIDX_set_partition_size(file, partition_box_size[0], partition_box_size[1], partition_box_size[2]);
  PIDX_set_partition_size(file, global_box_size[0], global_box_size[1], global_box_size[2]);

  PIDX_set_block_count(file, blocks_per_file);
  PIDX_set_block_size(file, 15);
  PIDX_set_bit_string_type(file, bit_string_type);

  // Selecting idx I/O mode
  //PIDX_set_io_mode(file, PIDX_WAVELET_IO);
  PIDX_set_io_mode(file, PIDX_MERGE_TREE_ANALYSIS);

  //PIDX_set_wavelet_implementation_type(file, wavelet_type);
  //PIDX_set_wavelet_level(file, wavelet_level);

  //PIDX_set_compression_type(file, PIDX_CHUNKING_ZFP_63_COEFFICIENT);
  //PIDX_set_lossy_compression_bit_rate(file, bit_rate1);
  //PIDX_set_zfp_precisison(file, precisison);


  return;
}

static void set_pidx_file2(int ts)
{
  PIDX_return_code ret;

  ret = PIDX_file_create(output_file_name2, PIDX_MODE_CREATE, p_access2, global_size, &file2);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  PIDX_set_current_time_step(file2, ts);
  PIDX_set_variable_count(file2, variable_count);

  //PIDX_debug_disable_agg(file);

  PIDX_point reg_size;
  //PIDX_set_point(reg_size, local_box_size[X], local_box_size[Y], local_box_size[Z]);
  PIDX_set_point(reg_size, partition_box_size[X] * local_box_size[0], partition_box_size[Y] * local_box_size[1], partition_box_size[Z] * local_box_size[2]);
  PIDX_set_restructuring_box(file2, reg_size);

  PIDX_set_process_decomposition(file2, global_box_size[X]/local_box_size[X], global_box_size[Y]/local_box_size[Y], global_box_size[Z]/local_box_size[Z]);
  //fprintf(stderr, "XYZ %lld / %lld   %lld / %lld    %lld / %lld\n", global_box_size[X], local_box_size[X], global_box_size[Y], local_box_size[Y], global_box_size[Z], local_box_size[Z]);
  PIDX_set_partition_size(file2, partition_box_size[0], partition_box_size[1], partition_box_size[2]);

  PIDX_set_block_count(file2, blocks_per_file);
  PIDX_set_block_size(file2, 15);
  PIDX_set_bit_string_type(file2, bit_string_type);

  // Selecting idx I/O mode
  PIDX_set_io_mode(file2, PIDX_WAVELET_IO);
  //PIDX_set_io_mode(file, PIDX_WAVELET_IO);
  PIDX_set_wavelet_implementation_type(file2, wavelet_type);
  PIDX_set_wavelet_level(file2, wavelet_level);

  PIDX_set_compression_type(file2, PIDX_CHUNKING_AVERAGE);
  PIDX_set_average_compression_factor(file2, 64, bit_rate2);
  //PIDX_set_zfp_precisison(file, precisison);

  return;
}


//----------------------------------------------------------------
static void set_pidx_variable(int var)
{
  PIDX_return_code ret = 0;

  ret = PIDX_variable_create(var_name[var],  bpv[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, /*data[var]*/sub_sampled_buffer, PIDX_row_major);
  //ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

  ret = PIDX_append_and_write_variable(file, variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  return;
}

static void set_pidx_variable2(int var)
{
  PIDX_return_code ret = 0;

  ret = PIDX_variable_create(var_name[var],  bpv[var], type_name[var], &variable2[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  ret = PIDX_variable_write_data_layout(variable2[var], local_offset, local_size, data[var], PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

  ret = PIDX_append_and_write_variable(file2, variable2[var]);
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

static void destroy_pidx_var_point_and_access2()
{
  if (PIDX_close_access(p_access2) != PIDX_success)
    terminate_with_error_msg("PIDX_close_access");

  free(variable2);
  variable2 = 0;

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
