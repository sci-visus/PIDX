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

static int create_patches();
static void destroy_patches();
static int parse_args(int argc, char **argv);
static void usage(void);
static int verify_parameters();
static void report_error(char* func_name, char* file_name, int line_no);

static int ***local_patch_size;
static int ***local_patch_offset;
static double ***double_data;
static uint64_t ***ulong_data;
static int *values_per_sample;    // Example: 1 for scalar 3 for vector
static int local_box_size[3] = {0, 0, 0};             ///< local dimensions of the per-process block

static int global_box_size[3] = {0, 0, 0};            ///< global dimensions of 3D volume
static int local_box_offset[3];
static int time_step_count = 1;                       ///< Number of time-steps
static int variable_count = 1;                        ///< Number of fields
static char variable_type[512];   ///< output IDX file Name Template
static char output_file_template[512];   ///< output IDX file Name Template
static int patch_count;
static int bits_per_block;
static int blocks_per_file;
static int64_t restructured_box_size[5] = {1,1,1,1,1};
static int compression_type;
static int compression_bit_rate;
static int perform_phases[6];
static int debug_rst_hz[2];
static int dump_agg_io[2];
static int idx_count[3];
static int aggregation_factor;
static int is_rank_z_ordering;
static int hz_from_to[2];
static int rst_agg[2];

int main(int argc, char **argv)
{
  int ret;
  int i, j, k, p, ts;
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

#if 1
  //  The command line arguments are shared by all processes
#if PIDX_HAVE_MPI
  MPI_Bcast(global_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_box_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&patch_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&variable_type, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&blocks_per_file, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bits_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(restructured_box_size, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&compression_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&compression_bit_rate, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(perform_phases, 6, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(debug_rst_hz, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dump_agg_io, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(idx_count, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&aggregation_factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_rank_z_ordering, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(hz_from_to, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rst_agg, 2, MPI_INT, 0, MPI_COMM_WORLD);
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

  unsigned int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice;
  if (is_rank_z_ordering == 0)
  {
    rank_z = rank / (sub_div[0] * sub_div[1]);
    rank_slice = rank % (sub_div[0] * sub_div[1]);
    rank_y = (rank_slice / sub_div[0]);
    rank_x = (rank_slice % sub_div[0]);
  }

  create_patches();

  if (strcmp(variable_type, "double") == 0)
  {
    double_data = malloc(sizeof(double**) * variable_count);
    for (var = 0; var < variable_count; var++)
    {
      values_per_sample[var] = 1;

      double_data[var] = malloc(sizeof(double*) * patch_count);
      for(p = 0 ; p < patch_count ; p++)
      {
        double_data[var][p] = malloc(sizeof (double) * local_patch_size[var][p][0] * local_patch_size[var][p][1] * local_patch_size[var][p][2] * local_patch_size[var][p][3]);

        for (k = 0; k < local_patch_size[var][p][2]; k++)
          for (j = 0; j < local_patch_size[var][p][1]; j++)
            for (i = 0; i < local_patch_size[var][p][0]; i++)
            {
              int64_t index = (int64_t) (local_patch_size[var][p][0] * local_patch_size[var][p][1] * k) + (local_patch_size[var][p][0] * j) + i;
              for (vps = 0; vps < values_per_sample[var]; vps++)
                double_data[var][p][index * values_per_sample[var] + vps] = (100 + var + ((global_box_size[0] * global_box_size[1]*(local_patch_offset[var][p][2] + k))+(global_box_size[0]*(local_patch_offset[var][p][1] + j)) + (local_patch_offset[var][p][0] + i)));
            }
      }
    }
  }
  else if (strcmp(variable_type, "unsigned long long") == 0)
  {
    ulong_data = malloc(sizeof(double**) * variable_count);
    for (var = 0; var < variable_count; var++)
    {
      values_per_sample[var] = 1;

      ulong_data[var] = malloc(sizeof(double*) * patch_count);
      for(p = 0 ; p < patch_count ; p++)
      {
        ulong_data[var][p] = malloc(sizeof (double) * local_patch_size[var][p][0] * local_patch_size[var][p][1] * local_patch_size[var][p][2] * local_patch_size[var][p][3]);

        for (k = 0; k < local_patch_size[var][p][2]; k++)
          for (j = 0; j < local_patch_size[var][p][1]; j++)
            for (i = 0; i < local_patch_size[var][p][0]; i++)
            {
              int64_t index = (int64_t) (local_patch_size[var][p][0] * local_patch_size[var][p][1] * k) + (local_patch_size[var][p][0] * j) + i;
              for (vps = 0; vps < values_per_sample[var]; vps++)
              {
                ulong_data[var][p][index * values_per_sample[var] + vps] = (100 + var + ((global_box_size[0] * global_box_size[1]*(local_patch_offset[var][p][2] + k))+(global_box_size[0]*(local_patch_offset[var][p][1] + j)) + (local_patch_offset[var][p][0] + i)));
                //printf("Value at %lld %lld %lld = %lld\n", (unsigned long long)i, (unsigned long long)j, (unsigned long long)k, (unsigned long long)ulong_data[var][p][index * values_per_sample[var] + vps]);
              }
            }
      }
    }
  }

  PIDX_file file;            // IDX file descriptor
  PIDX_variable* variable;   // variable descriptor

  variable = malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);

  PIDX_point global_bounding_box, **local_offset_point, **local_box_count_point;

  local_offset_point = malloc(sizeof(PIDX_point*) * variable_count);
  local_box_count_point = malloc(sizeof(PIDX_point*) * variable_count);
  for(var = 0; var < variable_count; var++)
  {
    local_offset_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    local_box_count_point[var] = malloc(sizeof(PIDX_point) * patch_count);
    for(p = 0 ; p < patch_count ; p++)
    {
      PIDX_set_point_5D(local_offset_point[var][p], (int64_t)local_patch_offset[var][p][0], (int64_t)local_patch_offset[var][p][1], (int64_t)local_patch_offset[var][p][2], 0, 0);
      PIDX_set_point_5D(local_box_count_point[var][p], (int64_t)local_patch_size[var][p][0], (int64_t)local_patch_size[var][p][1], (int64_t)local_patch_size[var][p][2], 1, 1);
    }
  }
  PIDX_set_point_5D(global_bounding_box, (int64_t)global_box_size[0], (int64_t)global_box_size[1], (int64_t)global_box_size[2], 1, 1);

  PIDX_access access;
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
  PIDX_set_idx_count(access, idx_count[0], idx_count[1], idx_count[2]);
  PIDX_set_process_extent(access, sub_div[0], sub_div[1], sub_div[2]);
  PIDX_set_process_rank_decomposition(access, rank_x, rank_y, rank_z);
#endif

  for (ts = 0; ts < time_step_count; ts++)
  {
    PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_variable_count(file, variable_count);

    PIDX_set_resolution(file, hz_from_to[0], hz_from_to[1]);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);

    // PIDX set restructuring box size
    PIDX_set_restructuring_box(file, restructured_box_size);

    // PIDX compression related calls
    if (compression_type == 0)
      PIDX_set_compression_type(file, PIDX_NO_COMPRESSION);
    if (compression_type == 1)
      PIDX_set_compression_type(file, PIDX_CHUNKING_ONLY);
    if (compression_type == 2)
    {
      PIDX_set_compression_type(file, PIDX_CHUNKING_ZFP);
      PIDX_set_lossy_compression_bit_rate(file, compression_bit_rate);
    }

    PIDX_debug_rst(file, debug_rst_hz[0]);
    PIDX_debug_hz(file, debug_rst_hz[1]);
    PIDX_dump_agg_info(file, dump_agg_io[0]);
    PIDX_dump_io_info(file, dump_agg_io[1]);

    if (perform_phases[0] == 0)
      PIDX_debug_disable_restructuring(file);
    if (perform_phases[1] == 0)
      PIDX_debug_disable_chunking(file);
    if (perform_phases[2] == 0)
      PIDX_debug_disable_hz(file);
    if (perform_phases[3] == 0)
      PIDX_debug_disable_compression(file);
    if (perform_phases[4] == 0)
      PIDX_debug_disable_agg(file);
    if (perform_phases[5] == 0)
      PIDX_debug_disable_io(file);

    if (rst_agg[0] == 0)
      PIDX_disable_rst(file);
    if (rst_agg[1] == 0)
      PIDX_disable_agg(file);

    for (var = 0; var < variable_count; var++)
    {
      char variable_name[512];
      sprintf(variable_name, "variable_%d", var);

      PIDX_data_type type;
      if (strcmp(variable_type, "double") == 0)
        strcpy(type, FLOAT64);
      else if (strcmp(variable_type, "unsigned long long") == 0)
        strcpy(type, UINT64);

      ret = PIDX_variable_create(variable_name, values_per_sample[var] * sizeof(uint64_t) * 8, type, &variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_variable_create", __FILE__, __LINE__);

      for (p = 0 ; p < patch_count ; p++)
      {
        if (strcmp(variable_type, "double") == 0)
          ret = PIDX_variable_write_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], double_data[var][p], PIDX_row_major);
        else if (strcmp(variable_type, "unsigned long long") == 0)
          ret = PIDX_variable_write_data_layout(variable[var], local_offset_point[var][p], local_box_count_point[var][p], ulong_data[var][p], PIDX_row_major);

        if (ret != PIDX_success)  report_error("PIDX_variable_data_layout", __FILE__, __LINE__);
      }

      ret = PIDX_append_and_write_variable(file, variable[var]);
      if (ret != PIDX_success)  report_error("PIDX_append_and_write_variable", __FILE__, __LINE__);
    }
    ret = PIDX_close(file);
    if (ret != PIDX_success)  report_error("PIDX_close", __FILE__, __LINE__);
  }
  PIDX_close_access(access);

  if (strcmp(variable_type, "double") == 0)
  {
    for (var = 0; var < variable_count; var++)
    {
      for(p = 0 ; p < patch_count ; p++)
        free(double_data[var][p]);
      free(double_data[var]);
    }
    free(double_data);  double_data = 0;
  }
  if (strcmp(variable_type, "unsigned long long") == 0)
  {
    for (var = 0; var < variable_count; var++)
    {
      for(p = 0 ; p < patch_count ; p++)
        free(ulong_data[var][p]);
      free(ulong_data[var]);
    }
    free(ulong_data);  ulong_data = 0;
  }

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

  destroy_patches();
#endif

#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  FILE* fp;
  char *pch;
  int count = 0, len = 0;
  fp = fopen(argv[1], "r");
  if (!fp)
  {
    fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
    return (-1);
  }

  printf("Parsing the Input File\n");
  char line [ 512 ];
  while (fgets(line, sizeof (line), fp) != NULL)
  {
    line[strcspn(line, "\r\n")] = 0;

    //printf("%s: ", line);
    if (strcmp(line, "(global box)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        global_box_size[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
      //printf("%%d %d %d\n", global_box_size[0], global_box_size[1], global_box_size[2]);
    }

    if (strcmp(line, "(rst agg)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        rst_agg[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(local box)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        local_box_size[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(restructured box size)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        restructured_box_size[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(file name)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;

      len = strlen(line);
      strncpy(output_file_template, line, len);
    }

    if (strcmp(line, "(time steps)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      time_step_count = atoi(line);
    }

    if (strcmp(line, "(fields)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      variable_count= atoi(line);
    }

    if (strcmp(line, "(field type)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      len = strlen(line);
      strncpy(variable_type, line, len);
    }

    if (strcmp(line, "(patch count)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      patch_count= atoi(line);
    }

    if (strcmp(line, "(idx count)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        idx_count[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(debug rst:hz)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        debug_rst_hz[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(dump agg:io)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        dump_agg_io[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }
    if (strcmp(line, "(perform rst:brst:hz:comp:agg:io)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        perform_phases[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(compression type)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      compression_type= atoi(line);
    }

    if (strcmp(line, "(compression bit rate)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      compression_bit_rate= atoi(line);
    }

    if (strcmp(line, "(blocks per file)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      blocks_per_file= atoi(line);
    }

    if (strcmp(line, "(bits per block)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      bits_per_block = atoi(line);
    }

    if (strcmp(line, "(aggregator multiplier)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      aggregation_factor = atoi(line);
    }

    if (strcmp(line, "(rank z ordering)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      is_rank_z_ordering = atoi(line);
    }

    if (strcmp(line, "(hz from:to)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return (-1);
      line[strcspn(line, "\r\n")] = 0;
      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        hz_from_to[count] = atoi(pch) ;
        count++;
        pch = strtok(NULL, " ");
      }
    }
  }
  fclose(fp);

  if (verify_parameters() < 0)
    return (-1);

  return (0);
}

static int verify_parameters()
{
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

  if (time_step_count <= 0)
  {
    fprintf(stderr, "Time step count cannot be 0!!!!!!!!!\n");
    return (-1);
  }

  if (variable_count <= 0)
  {
    fprintf(stderr, "Variable count cannot be 0!!!!!!!!!\n");
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

static void destroy_patches()
{
  int i = 0, var = 0;
  for(var = 0; var < variable_count; var++)
  {
    for(i = 0; i < patch_count ; i++)
    {
      free(local_patch_size[var][i]);
      free(local_patch_offset[var][i]);
    }
    free(local_patch_size[var]);
    free(local_patch_offset[var]);
  }
  free (local_patch_size);
  free (local_patch_offset);
}

static int create_patches()
{
  int var = 0, d = 0, i = 0;
  local_patch_size = malloc(sizeof(int**) * variable_count);
  local_patch_offset = malloc(sizeof(int**) * variable_count);

  for(var = 0; var < variable_count; var++)
  {
    local_patch_size[var] = malloc(sizeof(int*) * patch_count);
    local_patch_offset[var] = malloc(sizeof(int*) * patch_count);
    for(i = 0; i < patch_count ; i++)
    {
      local_patch_size[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
      local_patch_offset[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
    }
  }

  for(var = 0; var < variable_count; var++)
  {
    // One patch for this variable
    if (patch_count == 1)
    {
      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        local_patch_size[var][0][d] = local_box_size[d];
        local_patch_offset[var][0][d] = local_box_offset[d];
      }
    }

    // two patches for this variable
    else if (patch_count == 2)
    {
      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        local_patch_size[var][0][d] = local_box_size[d];
        local_patch_offset[var][0][d] = local_box_offset[d];
        local_patch_size[var][1][d] = local_box_size[d];
        local_patch_offset[var][1][d] = local_box_offset[d];
      }

      local_patch_size[var][0][0] = local_box_size[0]/2;
      if(local_box_size[0] % 2 == 0)
        local_patch_size[var][1][0] = local_box_size[0]/2;
      else
        local_patch_size[var][1][0] = local_box_size[0]/2 + 1;

      local_patch_offset[var][1][0] = local_box_offset[0] + local_box_size[0]/2;
    }

    // four patches for this variable
    else if (patch_count == 4)
    {
      for(i = 0; i < patch_count; i++)
      {
        for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          local_patch_size[var][i][d] = local_box_size[d];
          local_patch_offset[var][i][d] = local_box_offset[d];
        }
      }
      local_patch_size[var][0][0] = local_box_size[0]/2;
      local_patch_size[var][0][1] = local_box_size[1]/2;

      local_patch_size[var][1][1] = local_box_size[1]/2;
      if(local_box_size[0] % 2 == 0)
      {
        local_patch_size[var][1][0] = local_box_size[0]/2;
        local_patch_size[var][3][0] = local_box_size[0]/2;
        local_patch_offset[var][1][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][3][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
      }
      else
      {
        local_patch_size[var][1][0] = local_box_size[0]/2 + 1;
        local_patch_size[var][3][0] = local_box_size[0]/2 + 1;
        local_patch_offset[var][1][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][3][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
      }

      local_patch_size[var][2][0] = local_box_size[0]/2;
      if(local_box_size[1] % 2 == 0)
      {
        local_patch_size[var][2][1] = local_box_size[1]/2;
        local_patch_size[var][3][1] = local_box_size[1]/2;
        local_patch_offset[var][2][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][3][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
      }
      else
      {
        local_patch_size[var][2][1] = local_box_size[1]/2 + 1;
        local_patch_size[var][3][1] = local_box_size[1]/2 + 1;
        local_patch_offset[var][2][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][3][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
      }
    }

    // eight patches for this variable
    else if (patch_count == 8)
    {
      for(i = 0; i < patch_count; i++)
      {
        for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          local_patch_size[var][i][d] = local_box_size[d];
          local_patch_offset[var][i][d] = local_box_offset[d];
        }
      }
      local_patch_size[var][0][0] = local_box_size[0]/2;
      local_patch_size[var][0][1] = local_box_size[1]/2;

      local_patch_size[var][4][0] = local_box_size[0]/2;
      local_patch_size[var][4][1] = local_box_size[1]/2;

      local_patch_size[var][1][1] = local_box_size[1]/2;
      local_patch_size[var][5][1] = local_box_size[1]/2;
      if(local_box_size[0] % 2 == 0)
      {
        local_patch_size[var][1][0] = local_box_size[0]/2;
        local_patch_size[var][3][0] = local_box_size[0]/2;
        local_patch_offset[var][1][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][3][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;

        local_patch_size[var][5][0] = local_box_size[0]/2;
        local_patch_size[var][7][0] = local_box_size[0]/2;
        local_patch_offset[var][5][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][7][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
      }
      else
      {
        local_patch_size[var][1][0] = local_box_size[0]/2 + 1;
        local_patch_size[var][3][0] = local_box_size[0]/2 + 1;
        local_patch_offset[var][1][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][3][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;

        local_patch_size[var][5][0] = local_box_size[0]/2 + 1;
        local_patch_size[var][7][0] = local_box_size[0]/2 + 1;
        local_patch_offset[var][5][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
        local_patch_offset[var][7][0] = local_patch_offset[var][0][0] + local_box_size[0]/2;
      }

      local_patch_size[var][2][0] = local_box_size[0]/2;
      local_patch_size[var][6][0] = local_box_size[0]/2;

      if(local_box_size[1] % 2 == 0)
      {
        local_patch_size[var][2][1] = local_box_size[1]/2;
        local_patch_size[var][3][1] = local_box_size[1]/2;
        local_patch_offset[var][2][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][3][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;

        local_patch_size[var][6][1] = local_box_size[1]/2;
        local_patch_size[var][7][1] = local_box_size[1]/2;
        local_patch_offset[var][6][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][7][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
      }
      else
      {
        local_patch_size[var][2][1] = local_box_size[1]/2 + 1;
        local_patch_size[var][3][1] = local_box_size[1]/2 + 1;
        local_patch_offset[var][2][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][3][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;

        local_patch_size[var][6][1] = local_box_size[1]/2 + 1;
        local_patch_size[var][7][1] = local_box_size[1]/2 + 1;
        local_patch_offset[var][6][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
        local_patch_offset[var][7][1] = local_patch_offset[var][0][1] + local_box_size[1]/2;
      }

      local_patch_size[var][0][2] = local_box_size[2]/2;
      local_patch_size[var][1][2] = local_box_size[2]/2;
      local_patch_size[var][2][2] = local_box_size[2]/2;
      local_patch_size[var][3][2] = local_box_size[2]/2;
      if(local_box_size[1] % 2 == 0)
      {
        local_patch_size[var][4][2] = local_box_size[2]/2;
        local_patch_size[var][5][2] = local_box_size[2]/2;
        local_patch_size[var][6][2] = local_box_size[2]/2;
        local_patch_size[var][7][2] = local_box_size[2]/2;

        local_patch_offset[var][4][2] = local_patch_offset[var][0][2] + local_box_size[2]/2;
        local_patch_offset[var][5][2] = local_patch_offset[var][1][2] + local_box_size[2]/2;
        local_patch_offset[var][6][2] = local_patch_offset[var][2][2] + local_box_size[2]/2;
        local_patch_offset[var][7][2] = local_patch_offset[var][3][2] + local_box_size[2]/2;
      }
      else
      {
        local_patch_size[var][4][2] = local_box_size[2]/2 + 1;
        local_patch_size[var][5][2] = local_box_size[2]/2 + 1;
        local_patch_size[var][6][2] = local_box_size[2]/2 + 1;
        local_patch_size[var][7][2] = local_box_size[2]/2 + 1;

        local_patch_offset[var][4][2] = local_patch_offset[var][0][2] + local_box_size[2]/2;
        local_patch_offset[var][5][2] = local_patch_offset[var][1][2] + local_box_size[2]/2;
        local_patch_offset[var][6][2] = local_patch_offset[var][2][2] + local_box_size[2]/2;
        local_patch_offset[var][7][2] = local_patch_offset[var][3][2] + local_box_size[2]/2;
      }
    }
    else
      printf("This patch count not supported !!!!\n");
  }
  return 0;
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
