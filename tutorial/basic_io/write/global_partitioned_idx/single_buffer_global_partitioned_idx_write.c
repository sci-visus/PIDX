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
static char var_list[512];
static unsigned char **data;
static char output_file_name[512];
static char var_name[MAX_VAR_COUNT][512];
static int bpv[MAX_VAR_COUNT];
static char type_name[MAX_VAR_COUNT][512];
static int vps[MAX_VAR_COUNT];

static int wavelet_type = 0;
static int wavelet_level = WAVELET_STENCIL;

static PIDX_point global_size, local_offset, local_size;
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
#if 1
  //rank_0_print("Simulation Data Created\n");

  create_pidx_var_point_and_access();
  for (ts = 0; ts < time_step_count; ts++)
  {
    set_pidx_file(ts);
    for (var = 0; var < variable_count; var++)
      set_pidx_variable(var);
    PIDX_close(file);
  }
  destroy_pidx_var_point_and_access();
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
  char flags[] = "g:l:c:f:t:v:b:a:w:x:";
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

  /*
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    int v = 0;
    for(v = 0; v < variable_count; v++)
      printf("[%d] -> %s %d %d\n", v, var_name[v], bpv[v], vps[v]);
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
              //fvalue = 100 + var + vps + ((global_box_size[X] * global_box_size[Y]*(local_box_offset[Z] + k))+(global_box_size[X]*(local_box_offset[Y] + j)) + (local_box_offset[X] + i));
              //fvalue = (local_box_offset[X] + i) * (local_box_offset[X] + i) + (local_box_offset[Y] + j) * (local_box_offset[Y] + j) + (local_box_offset[Z] + k) * (local_box_offset[Z] + k);
              fvalue = (local_box_offset[X] + i) * (local_box_offset[X] + i) + (local_box_offset[Y] + j) * (local_box_offset[Y] + j) + (local_box_offset[Z] + k) * (local_box_offset[Z] + k);
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
    int fp = open("magnetic-512-volume.raw", O_RDONLY);
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

          //printf("[%d] : read offset %d write offset %d Count %d\n", rank, send_o, recv_o, send_c);
          pread(fp, data[0] + send_o * sizeof(float), send_c * sizeof(float), recv_o * sizeof(float));
        }
      }
    }
    close(fp);
#endif
#if 0
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
  //printf("XYZ %lld / %lld   %lld / %lld    %lld / %lld\n", global_box_size[X], local_box_size[X], global_box_size[Y], local_box_size[Y], global_box_size[Z], local_box_size[Z]);
  PIDX_set_partition_size(file, partition_box_size[0], partition_box_size[1], partition_box_size[2]);

  PIDX_set_block_count(file, blocks_per_file);
  PIDX_set_block_size(file, 15);
  PIDX_set_bit_string_type(file, bit_string_type);

  // Selecting idx I/O mode
  //PIDX_set_io_mode(file, PIDX_GLOBAL_PARTITION_IDX_IO);
  PIDX_set_io_mode(file, PIDX_WAVELET_IO);
  PIDX_set_wavelet_implementation_type(file, wavelet_type);
  PIDX_set_wavelet_level(file, wavelet_level);

  PIDX_set_compression_type(file, PIDX_ZFP_COMPRESSION);

  return;
}

//----------------------------------------------------------------
static void set_pidx_variable(int var)
{
  PIDX_return_code ret = 0;

  ret = PIDX_variable_create(var_name[var],  bpv[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

  ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

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
