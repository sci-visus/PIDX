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
  PIDX write example

  In this example we show how to write data using the PIDX library.

  We consider a global 3D regular grid domain that we will call
  global domain (g).
  This global domain represents the grid space where all the data are stored.

  In a parallel environment each core (e.g. MPI rank) owns a portion of the data
  that has to be written on the disk. We refer to this portion of the domain as
  local domain (l).

  In this example we well see how to execute parallel write with PIDX of a
  syntethic dataset.

  In the following picture is represented a sample domain decomposition
  of the global domain (l) in per-core local domains (l), sometimes referred
  as patches.
  In this example all the local domains have same dimesions for simplicity.
  PIDX supports different number and sizes of patches per core.

             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        IDX Data format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*

*/
#if !defined _MSC_VER
#include <unistd.h>
#endif
#include <stdarg.h>
#include <stdint.h>
#include <ctype.h>
#include <PIDX.h>

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

#include "pidx_examples_utils.h"

char var_name[MAX_VAR_COUNT][512];
int bpv[MAX_VAR_COUNT];
char type_name[MAX_VAR_COUNT][512];
int vps[MAX_VAR_COUNT];
PIDX_variable* variable;
char output_file_template[512];
char var_list[512];
char output_file_name[512];
unsigned char **data;

char *usage = "Serial Usage: ./idx_write -g 32x32x32 -l 32x32x32 -v 2 -t 4 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./idx_write -g 64x64x64 -l 32x32x32 -v 2 -t 4 -f output_idx_file_name\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -r: restructured box dimension\n"
                     "  -f: file name template (without .idx)\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables (or file containing a list of variables)\n";

static int generate_vars();
static void parse_args(int argc, char **argv);
static int parse_var_list();
static void create_synthetic_simulation_data();
static void set_pidx_variable(int var);
static void set_pidx_file(int ts);
static void destroy_synthetic_simulation_data();

int main(int argc, char **argv)
{
  int ts = 0, var = 0;

  // Init MPI and MPI vars (e.g. rank and process_count)
  init_mpi(argc, argv);

  // Parse input arguments and initialize
  // corresponing variables
  parse_args(argc, argv);

  // Verify that the domain decomposition is valid
  // for the given number of cores
  check_args();

  // Initialize per-process local domain
  calculate_per_process_offsets();

  // Generate synthetic data
  create_synthetic_simulation_data();

  // Create variables
  create_pidx_point_and_access();
  variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  for (ts = 0; ts < time_step_count; ts++)
  {
    // Set PIDX_file for this timestep
    set_pidx_file(ts);

    // Set all the PIDX_variable that we want to write
    for (var = 0; var < variable_count; var++)
      set_pidx_variable(var);

    // PIDX_close triggers the actual write on the disk
    // of the variables that we just set
    PIDX_close(file);
  }

  // Close access and free memory
  if (PIDX_close_access(p_access) != PIDX_success)
    terminate_with_error_msg("PIDX_close_access");

  if (PIDX_free_metadata_cache(cache) != PIDX_success)
    terminate_with_error_msg("PIDX_free_meta_data_cache");

  free(variable);
  variable = 0;

  destroy_synthetic_simulation_data();

  shutdown_mpi();

  return 0;
}


//----------------------------------------------------------------
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:t:v:";
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
      if (!isNumber(optarg)){ // the param is a file with the list of variables
        if (sprintf(var_list, "%s", optarg) > 0)
          parse_var_list();
        else
          terminate_with_error_msg("Invalid variable list file\n%s", usage);
      }else { // the param is a number of variables (default: 1*float32)
        if (sscanf(optarg, "%d", &variable_count) > 0)
          generate_vars();
        else
          terminate_with_error_msg("Invalid number of variables\n%s", usage);
      }
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

static void set_pidx_variable(int var)
{
  PIDX_return_code ret = 0;

  // Set variable name, number of bits, typename
  ret = PIDX_variable_create(var_name[var], bpv[var] * vps[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  // Set the variable offset and size of the local domain,
  // where the data is in memory (data) and what is its layout in memory (row major)
  ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

  // Tell PIDX that we want to write this variable
  ret = PIDX_append_and_write_variable(file, variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  return;
}

static int generate_vars(){

  int variable_counter = 0;

  for (variable_counter = 0; variable_counter < variable_count; variable_counter++){
    int ret;
    int bits_per_sample = 0;
    int sample_count = 0;
    char temp_name[512];
    char* temp_type_name = "1*float32";
    //char* temp_type_name = "1*int8";
    sprintf(temp_name, "var_%d", variable_counter);
    strcpy(var_name[variable_counter], temp_name);
    strcpy(type_name[variable_counter], temp_type_name);

    ret = PIDX_values_per_datatype(temp_type_name, &sample_count, &bits_per_sample);
    if (ret != PIDX_success)  return PIDX_err_file;

    bpv[variable_counter] = bits_per_sample;
    vps[variable_counter] = sample_count;
  }

  return 0;
}

//----------------------------------------------------------------
static void create_synthetic_simulation_data()
{
  int var = 0;
  data = malloc(sizeof(*data) * variable_count);
  memset(data, 0, sizeof(*data) * variable_count);

  // Synthetic simulation data
  for (var = 0; var < variable_count; var++)
  {
    uint64_t i, j, k, val_per_sample = 0;
    data[var] = malloc(sizeof (*(data[var])) * local_box_size[X] * local_box_size[Y] * local_box_size[Z] * (bpv[var]/8) * vps[var]);

    unsigned char cvalue = 0;
    short svalue = 0;
    float fvalue = 0;
    double dvalue = 0;
    int ivalue = 0;
    uint64_t u64ivalue = 0;
    int64_t i64value = 0;
#if 0
    FILE* fp = fopen("OF_data.txt", "r");
    int o = 0;
    /*
     for (o = 0; o < 1400; o++)
     {
     fscanf(fp, "%f\n", &fvalue);
     printf("[%d] val = %.16f\n", o, fvalue);
     memcpy(data[var] + o * sizeof(float), &fvalue, sizeof(float));
     }
     */
    fclose(fp);
#else
    for (k = 0; k < local_box_size[Z]; k++)
      for (j = 0; j < local_box_size[Y]; j++)
        for (i = 0; i < local_box_size[X]; i++)
        {
          uint64_t index = (uint64_t) (local_box_size[X] * local_box_size[Y] * k) + (local_box_size[X] * j) + i;

          for (val_per_sample = 0; val_per_sample < vps[var]; val_per_sample++)
          {
            if (strcmp(type_name[var], PIDX_DType.UINT8) == 0 || strcmp(type_name[var], PIDX_DType.UINT8_GA) == 0 || strcmp(type_name[var], PIDX_DType.UINT8_RGB) == 0)
            {
              cvalue = (int)(var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i)));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(unsigned char), &cvalue, sizeof(unsigned char));
            }
            if (strcmp(type_name[var], PIDX_DType.INT16) == 0 || strcmp(type_name[var], PIDX_DType.INT16_GA) == 0 || strcmp(type_name[var], PIDX_DType.INT16_RGB) == 0)
            {
              svalue = (int)(var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i)));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(short), &svalue, sizeof(short));
            }
            if (strcmp(type_name[var], PIDX_DType.INT32) == 0 || strcmp(type_name[var], PIDX_DType.INT32_GA) == 0 || strcmp(type_name[var], PIDX_DType.INT32_RGB) == 0)
            {
              ivalue = (int)( 100 + var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i)));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(int), &ivalue, sizeof(int));
            }
            else if (strcmp(type_name[var], PIDX_DType.FLOAT32) == 0 || strcmp(type_name[var], PIDX_DType.FLOAT32_GA) == 0 || strcmp(type_name[var], PIDX_DType.FLOAT32_RGB) == 0)
            {
              fvalue = (float)( 100 + var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i)));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(float), &fvalue, sizeof(float));
            }
            else if (strcmp(type_name[var], PIDX_DType.FLOAT64) == 0 || strcmp(type_name[var], PIDX_DType.FLOAT64_GA) == 0 || strcmp(type_name[var], PIDX_DType.FLOAT64_RGB) == 0)
            {
              dvalue = (double) 100 + var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(double), &dvalue, sizeof(double));
            }
            else if (strcmp(type_name[var], PIDX_DType.INT64) == 0 || strcmp(type_name[var], PIDX_DType.INT64_GA) == 0 || strcmp(type_name[var], PIDX_DType.INT64_RGB) == 0)
            {
              i64value = (int64_t) 100 + var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(int64_t), &i64value, sizeof(int64_t));
            }
            else if (strcmp(type_name[var], PIDX_DType.UINT64) == 0 || strcmp(type_name[var], PIDX_DType.UINT64_GA) == 0 || strcmp(type_name[var], PIDX_DType.UINT64_RGB) == 0)
            {
              u64ivalue = (uint64_t) 100 + var + val_per_sample + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              memcpy(data[var] + (index * vps[var] + val_per_sample) * sizeof(uint64_t), &u64ivalue, sizeof(uint64_t));
            }
          }
        }
#endif
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
      if ( fgets(line, sizeof line, fp) == NULL)
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
            int sample_count = 0;
            ret = PIDX_values_per_datatype(type_name[variable_counter], &sample_count, &bits_per_sample);
            if (ret != PIDX_success)  return PIDX_err_file;

            bpv[variable_counter] = bits_per_sample;
            vps[variable_counter] = sample_count;
          }
          count++;
          pch1 = strtok(NULL, " +");
        }
        count = 0;

        if ( fgets(line, sizeof line, fp) == NULL)
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
   for (v = 0; v < variable_count; v++)
   fprintf(stderr, "[%d] -> %s %d %d\n", v, var_name[v], bpv[v], vps[v]);
   }
   */

  return PIDX_success;
}

//----------------------------------------------------------------
static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  // Create IDX file
  ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, p_access, global_size, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create\n");

  // Set the current timestep
  PIDX_set_current_time_step(file, ts);
  // Set the number of variables
  PIDX_set_variable_count(file, variable_count);

  // Advanced settings
  PIDX_set_meta_data_cache(file, cache);

  // Select I/O mode (PIDX_IDX_IO for the multires, PIDX_RAW_IO for non-multires)
  PIDX_set_io_mode(file, PIDX_IDX_IO);

  // Set how many blocks we want to write in a single file
  PIDX_set_block_count(file, 128);

  // Set the size of a block: how many 2^N samples we want to put in a single block
  PIDX_set_block_size(file, 15);

  // If the domain decomposition and the cores configuration do not change over time
  // we can instruct PIDX to cache and reuse these information for the next timesteps
  PIDX_set_cache_time_step(file, 0);

  return;
}

static void destroy_synthetic_simulation_data()
{
  int var = 0;
  for (var = 0; var < variable_count; var++)
  {
    free(data[var]);
    data[var] = 0;
  }
  free(data);
  data = 0;
}
