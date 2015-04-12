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

// This program reads an original IDX file and convert it into a compressed IDX file
#include "pidxtest.h"

#ifdef __cplusplus
extern "C" {
#endif

int64_t mylog2(int64_t num)
{
  if (num <= 0) return -1;
  int64_t log = 0;
  int64_t n = 1;
  while (n < num)
  {
    n *= 2;
    ++log;
  }
  return log;
}

int test_converter(struct Args args, int rank)
{
  /// The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_steps, args.time_step_count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.idx_count, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.compression_block_size, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.restructured_box_size, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.compression_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.compression_bit_rate, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.perform_brst, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.perform_hz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.perform_compression, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.perform_agg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.perform_io, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.debug_rst, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.debug_hz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.dump_agg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.blocks_per_file, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.bits_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.aggregation_factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variables, args.variable_count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.hz_from, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.hz_to, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.is_rank_z_ordering, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.is_global_indexing, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Creating the filenames
  char *input_file_name = (char *)malloc(sizeof(char) * 512);
  sprintf(input_file_name, "%s%s", args.output_file_template, ".idx");
  args.output_file_name = (char *)malloc(sizeof(char) * 512);
  sprintf(args.output_file_name, "%scompressed%s", args.output_file_template, ".idx");

  // Calculating every process's offset and count
  unsigned int sub_div[3], local_offset[3];
  sub_div[0] = (args.extents[0] / args.count_local[0]);
  sub_div[1] = (args.extents[1] / args.count_local[1]);
  sub_div[2] = (args.extents[2] / args.count_local[2]);
  int slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * args.count_local[2];
  local_offset[1] = (slice / sub_div[0]) * args.count_local[1];
  local_offset[0] = (slice % sub_div[0]) * args.count_local[0];

  unsigned int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice = 0;
  rank_z = rank / (sub_div[0] * sub_div[1]);
  rank_slice = rank % (sub_div[0] * sub_div[1]);
  rank_y = (rank_slice / sub_div[0]);
  rank_x = (rank_slice % sub_div[0]);
  

  PIDX_point global_bounding_box, local_offset_point, local_box_count_point, compression_block_size_point, restructured_box_size_point;
  PIDX_set_point_5D(global_bounding_box, args.extents[0], args.extents[1], args.extents[2], 1, 1);
  PIDX_set_point_5D(local_offset_point, local_offset[0], local_offset[1], local_offset[2], 0, 0);
  PIDX_set_point_5D(local_box_count_point, args.count_local[0], args.count_local[1], args.count_local[2], 1, 1);
  PIDX_set_point_5D(compression_block_size_point, args.compression_block_size[0], args.compression_block_size[1], args.compression_block_size[2], 1, 1);
  PIDX_set_point_5D(restructured_box_size_point, args.restructured_box_size[0], args.restructured_box_size[1], args.restructured_box_size[2], 1, 1);

  PIDX_time_step_caching_ON();
  int time_step = 0;
  for (time_step = 0; time_step < args.time_step_count; time_step++)
  {
    PIDX_access access;
    PIDX_create_access(&access);
    PIDX_set_mpi_access(access, args.idx_count[0], args.idx_count[1], args.idx_count[2], MPI_COMM_WORLD);
    PIDX_set_global_indexing_order(access, 0);
    PIDX_set_process_extent(access, sub_div[0], sub_div[1], sub_div[2]);
    PIDX_set_process_rank_decomposition(access, rank_x, rank_y, rank_z);

    int bits_per_block = 15;
    int variable_count = 1;
    int blocks_per_file = 256;
    int var = 0;
    
    double **double_data = (double **)malloc(sizeof(*double_data) * variable_count); // DUONG_HARDCODE?
    memset(double_data, 0, sizeof(*double_data) * variable_count);
    for (var = 0; var < variable_count; var++)
    {
      //double_data[var] = (double*)malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * 1/*read_variable[var]->values_per_sample*/);
    }
    
    
    PIDX_variable *write_variable;
    // get read parameters
    PIDX_file input_file;
    PIDX_file_open(input_file_name, PIDX_file_rdonly, access, &input_file);
    
    PIDX_get_dims(input_file, global_bounding_box);
    PIDX_set_current_time_step(input_file, time_step/*args.time_steps[time_step]*/);
    PIDX_get_block_size(input_file, &bits_per_block);
    PIDX_get_block_count(input_file, &blocks_per_file);
    PIDX_get_variable_count(input_file, &variable_count);

    // enable read phases
    PIDX_enable_block_restructuring(input_file, args.perform_brst);
    PIDX_enable_hz(input_file, args.perform_hz);
    // NOTE: we don't enable compression for input
    PIDX_enable_agg(input_file, args.perform_agg);
    PIDX_enable_io(input_file, args.perform_io);
    

    // Allocate memory to read variables
    PIDX_variable *read_variable = (PIDX_variable *)malloc(sizeof(*read_variable) * variable_count);
    memset(read_variable, 0, sizeof(*read_variable) * variable_count);
    

    
    
    // read
    for (var = 0; var < variable_count; var++)
    {
      PIDX_get_next_variable(input_file, &read_variable[var]);
      double_data[var] = (double*)malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * 1);
      memset(double_data[var], 0, sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * 1);
      PIDX_read_next_variable(read_variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);
    }
    PIDX_close(input_file);
    
    
    for(var = 0; var < variable_count; var++)
    {
      int i, j, k, spv;
      for (k = 0; k < args.count_local[2]; k++)
        for (j = 0; j < args.count_local[1]; j++)
          for (i = 0; i < args.count_local[0]; i++)
          {
            int64_t index = (int64_t) (args.count_local[0] * args.count_local[1] * k) + (args.count_local[0] * j) + i;
            
            if ((int)double_data[var][index] == (int)100 + var + ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i)))
            {
              //printf("Y values = %f %lld\n", double_data[var][index * variable[var]->values_per_sample + spv], (long long)100 + var + ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i)));
            }
            else
            {
              //printf("X values = %f %lld\n", double_data[var][index * variable[var]->values_per_sample + spv], (long long)100 + var + ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i)));
              printf("Shout !!!!!\n");
            }
            
          }
    }
    
    
    
    // set write parameters
    PIDX_file output_file;
    PIDX_file_create(args.output_file_name, PIDX_file_trunc, access, &output_file);
    PIDX_set_dims(output_file, global_bounding_box);
    PIDX_set_current_time_step(output_file, /*args.time_steps[time_step]*/time_step);
    PIDX_set_block_size(output_file, bits_per_block);
    
    PIDX_set_block_count(output_file, blocks_per_file);
    PIDX_set_variable_count(output_file, variable_count);
    PIDX_set_aggregation_factor(output_file, args.aggregation_factor);
    PIDX_set_compression_type(output_file, args.compression_type);
    PIDX_set_compression_block_size(output_file, compression_block_size_point);
    PIDX_set_lossy_compression_bit_rate(output_file, args.compression_bit_rate);
    PIDX_set_restructuring_box(output_file, restructured_box_size_point); // DUONG_TODO: do we have to set this for input?

    // debug write phases
    PIDX_debug_rst(output_file, args.debug_rst);
    PIDX_debug_hz(output_file, args.debug_hz);
    PIDX_dump_agg_info(output_file, args.dump_agg);

    // enable write phases
    PIDX_enable_block_restructuring(output_file, args.perform_brst);
    PIDX_enable_hz(output_file, args.perform_hz);
    PIDX_enable_compression(output_file, args.perform_compression);
    PIDX_enable_agg(output_file, args.perform_agg);
    PIDX_enable_io(output_file, args.perform_io);
    
    write_variable = (PIDX_variable *)malloc(sizeof(*write_variable) * variable_count);
    memset(write_variable, 0, sizeof(*write_variable) * variable_count);
    
    // write
    for (var = 0; var < variable_count; var++)
    {
      PIDX_variable_create(output_file, "var_name", 1 * sizeof(uint64_t) * 8, "1*float64", &write_variable[var]);
      PIDX_append_and_write_variable(write_variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);
    }
    PIDX_close(output_file);
    
    
    PIDX_close_access(access);

    // free stuffs
    for (int i = 0; i < variable_count; ++i)
    {
      free(double_data[i]);
    }
    free(double_data);
    double_data = 0;
    
    free(read_variable);
    read_variable = 0;
    
    free(write_variable);
    write_variable = 0;
  }
  PIDX_time_step_caching_OFF();
  free(args.output_file_name);
  free(input_file_name);

  return 0;
}

#ifdef __cplusplus
}
#endif
