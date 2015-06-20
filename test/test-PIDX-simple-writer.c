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

#include "pidxtest.h"

#ifdef __cplusplus
extern "C" {
#endif

int test_simple_idx_writer(struct Args args, int rank, int nprocs)
{
#if PIDX_HAVE_MPI
  int i = 0, j = 0, k = 0;
  int ts, var, spv;
  int slice;
  unsigned int sub_div[3], local_offset[3];

  /// IDX file descriptor
  PIDX_file file;
  
  /// IDX File Name
  const char *output_file; 
  
  PIDX_variable* variable;                                       // variable descriptor
  //uint64_t     **long_data;
  double     **double_data;
  int* values_per_sample;
  
  PIDX_point restructured_box_size_point;
  PIDX_point chunk_size_point;
  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
    
  /// The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.chunk_size, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
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
  MPI_Bcast(&args.time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.idx_count, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.blocks_per_file, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.bits_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.aggregation_factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.hz_from, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.hz_to, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.is_rank_z_ordering, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.is_global_indexing, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  variable = (PIDX_variable*)malloc(sizeof(*variable) * args.variable_count);
  memset(variable, 0, sizeof(*variable) * args.variable_count);
  
  values_per_sample = (int*)malloc(sizeof(*values_per_sample) * args.variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * args.variable_count);      
  
  /// Creating the filename 
  args.output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(args.output_file_name, "%s%s", args.output_file_template,".idx");
  output_file = args.output_file_name;
  
  /// Calculating every process's offset and count  
  sub_div[0] = (args.extents[0] / args.count_local[0]);
  sub_div[1] = (args.extents[1] / args.count_local[1]);
  sub_div[2] = (args.extents[2] / args.count_local[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * args.count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * args.count_local[1];
  local_offset[0] = (slice % sub_div[0]) * args.count_local[0];
  
  unsigned int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice;

  //printf("Rank Row ordering\n");
  rank_z = rank / (sub_div[0] * sub_div[1]);
  rank_slice = rank % (sub_div[0] * sub_div[1]);
  rank_y = (rank_slice / sub_div[0]);
  rank_x = (rank_slice % sub_div[0]);
  
  PIDX_set_point_5D(global_bounding_box, (int64_t)args.extents[0], (int64_t)args.extents[1], (int64_t)args.extents[2], 1, 1);
  PIDX_set_point_5D(chunk_size_point, (int64_t)args.chunk_size[0], (int64_t)args.chunk_size[1], (int64_t)args.chunk_size[2], 1, 1);
  PIDX_set_point_5D(local_offset_point, (int64_t)local_offset[0], (int64_t)local_offset[1], (int64_t)local_offset[2], 0, 0);
  PIDX_set_point_5D(local_box_count_point, (int64_t)args.count_local[0], (int64_t)args.count_local[1], (int64_t)args.count_local[2], 1, 1);
  PIDX_set_point_5D(restructured_box_size_point, (int64_t)args.restructured_box_size[0], (int64_t)args.restructured_box_size[1], (int64_t)args.restructured_box_size[2], 1, 1);
  
  PIDX_time_step_caching_ON();  
  for (ts = 0; ts < args.time_step_count; ts++) 
  {
#if long_data
    long_data = malloc(sizeof(*long_data) * args.variable_count);
    memset(long_data, 0, sizeof(*long_data) * args.variable_count);
#else
    double_data = (double**)malloc(sizeof(*double_data) * args.variable_count);
    memset(double_data, 0, sizeof(*double_data) * args.variable_count);
#endif
    
    /// Creating access type (parallel here)
    PIDX_access access;
    PIDX_create_access(&access);
#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, MPI_COMM_WORLD);
    PIDX_set_idx_count(access, args.idx_count[0], args.idx_count[1], args.idx_count[2]);
    PIDX_set_process_extent(access, sub_div[0], sub_div[1], sub_div[2]);
    PIDX_set_process_rank_decomposition(access, rank_x, rank_y, rank_z);
    //PIDX_enable_topology_aware_io(access, args.topology_aware);
#else
    PIDX_set_default_access(access);
#endif
    //printf("compression type = %d\n", args.compression_type);  
    /// PIDX mandatory calls
    PIDX_file_create(output_file, PIDX_MODE_CREATE, access, &file);
    
    /// PIDX calls to set different parameters (optional)
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_block_size(file, args.bits_per_block);
    PIDX_set_aggregation_factor(file, args.aggregation_factor);
    PIDX_set_block_count(file, args.blocks_per_file);
    PIDX_set_variable_count(file, args.variable_count);
    PIDX_set_resolution(file, args.hz_from, args.hz_to);
    
    /// PIDX set restructuring box size
    PIDX_set_restructuring_box(file, restructured_box_size_point);
    
    /// PIDX compression related calls
    PIDX_set_compression_type(file, args.compression_type);
    PIDX_set_lossy_compression_bit_rate(file, args.compression_bit_rate);
    
    /// PIDX debuging different phases
    PIDX_debug_rst(file, args.debug_rst);
    PIDX_debug_hz(file, args.debug_hz);
    PIDX_dump_agg_info(file, args.dump_agg);
    
    /// PIDX enabling/disabling different phases
    /*
    PIDX_disable_rst(file);
    PIDX_debug_disable_chunking(file);
    PIDX_debug_disable_hz(file);
    PIDX_debug_disable_compression(file);
    PIDX_debug_disable_agg(file);
    PIDX_debug_disable_io(file);
    */

    char variable_name[512];
    char data_type[512];
    

#if 1
    for(var = 0; var < args.variable_count; var++)
    {
      values_per_sample[var] =  1;
      double_data[var] = (double*)malloc(sizeof (uint64_t) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * values_per_sample[var]);
      
      //const double pi = acos(-1.0);
      for (k = 0; k < args.count_local[2]; k++)
        for (j = 0; j < args.count_local[1]; j++)
          for (i = 0; i < args.count_local[0]; i++)
          {
            int64_t index = (int64_t) (args.count_local[0] * args.count_local[1] * k) + (args.count_local[0] * j) + i;
            for (spv = 0; spv < values_per_sample[var]; spv++)
              //double_data[var][index * values_per_sample[var] + spv] = 100 + var + ((args.count_local[0] * args.count_local[1]*(k))+(args.count_local[0]*(j)) + (i));
              double_data[var][index * values_per_sample[var] + spv] = 100 + var + ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i));
              //double_data[var][index * values_per_sample[var] + spv] = cos(2 * pi * i / args.count_local[0]) * cos(2 * pi * j / args.count_local[1]) * cos(2 * pi * k / args.count_local[2]);
          }
    }
    
    for(var = 0; var < args.variable_count; var++)
    {
      sprintf(variable_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);

      PIDX_variable_create(variable_name, values_per_sample[var] * sizeof(uint64_t) * 8, data_type, &variable[var]);

      PIDX_set_current_variable_index(file, var);
      PIDX_set_current_variable(file, variable[var]);

      PIDX_write_variable(file, variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);

      //PIDX_append_and_write_variable(file, variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);

    }
    
    PIDX_close(file);
    PIDX_close_access(access);
    
#if long_data
    for(var = 0; var < args.variable_count; var++)
    {
      free(long_data[var]);
      long_data[var] = 0;
    }
#else
    for(var = 0; var < args.variable_count; var++)
    {
      free(double_data[var]);
      double_data[var] = 0;
    }
#endif

#endif

#if long_data
    free(long_data);
    long_data = 0;
#else
    free(double_data);
    double_data = 0;
#endif
    
  }
  PIDX_time_step_caching_OFF();
  
  free(variable);
  free(values_per_sample);
  free(args.output_file_name);
  
#endif
  return 0;
}

/*   prints usage instructions   */
void usage_simple_idx_writer(void)
{
  printf("Usage: test-multi-idx-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}

#ifdef __cplusplus
}
#endif
