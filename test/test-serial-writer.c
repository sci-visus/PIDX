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
 **  http://www.sci.utah.edu/pidx                   **
 **  or contact: pascucci@sci.utah.edu              **
 **                                                 **
 *****************************************************/

#include "pidxtest.h"

int serial_writer(struct Args args)
{
#if 1
  int i = 0, j = 0, k = 0;
  int spv = 0, var = 0;
  int ts;

  /* IDX file */
  PIDX_file file = 0;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  
  /* IDX variables */
  PIDX_variable *variable;                                       // variable descriptor
  int *values_per_sample;
  double** double_data;
    
  variable = (PIDX_variable*)malloc(sizeof(*variable) * args.variable_count);
  memset(variable, 0, sizeof(*variable) * args.variable_count);
  
  values_per_sample = (int*) malloc(sizeof(*values_per_sample) * args.variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * args.variable_count);      
  
  //   Creating the filename
  args.output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(args.output_file_name, "%s%s", args.output_file_template, ".idx");

  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  PIDX_set_point_5D(global_bounding_box, args.extents[0], args.extents[1], args.extents[2], 1, 1);
  PIDX_set_point_5D(local_offset_point, 0, 0, 0, 0, 0);
  PIDX_set_point_5D(local_box_count_point, args.extents[0], args.extents[1], args.extents[2], 1, 1);
  
  output_file = args.output_file_name;
   
  PIDX_time_step_caching_ON();
  for (ts = 0; ts < args.time_step_count; ts++)
  {
    double_data = (double**)malloc(sizeof(*double_data) * args.variable_count);
    for (var = 0; var < args.variable_count; var++)
    {
      values_per_sample[var] = 1;
      double_data[var] = (double*)malloc(sizeof (*double_data[var]) * args.extents[0] * args.extents[1] * args.extents[2] * values_per_sample[var]);
      for (k = 0; k < args.extents[2]; k++)
        for (j = 0; j < args.extents[1]; j++)
          for (i = 0; i < args.extents[0]; i++) 
          {
            int64_t index = (int64_t) (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i;
            for (spv = 0; spv < values_per_sample[var]; spv++)
              double_data[var][index * values_per_sample[var] + spv] = (double)100 + var + (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i;
          }
    }
    
    PIDX_access access;
    PIDX_create_access(&access);
    //PIDX_set_default_access(access);
    
    /// PIDX mandatory calls
    PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
    
    /// PIDX calls to set different parameters (optional)
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_block_size(file, args.bits_per_block);
    PIDX_set_aggregation_factor(file, args.aggregation_factor);
    PIDX_set_block_count(file, args.blocks_per_file);
    PIDX_set_variable_count(file, args.variable_count);
    
    /// PIDX debuging different phases
    PIDX_debug_rst(file, args.debug_rst);
    PIDX_debug_hz(file, args.debug_hz);
    PIDX_dump_agg_info(file, args.dump_agg);
    PIDX_dump_io_info(file, args.dump_io);


    /// PIDX disabling different phases
    // PIDX_debug_disable_restructuring(file);
    // PIDX_debug_disable_chunking(file);
    // PIDX_debug_disable_hz(file);
    // PIDX_debug_disable_compression(file);
    // PIDX_debug_disable_agg(file);
    // PIDX_debug_disable_io(file);


    for(var = 0; var < args.variable_count; var++)
    {
      char variable_name[1024];
      char data_type[1024];
      sprintf(variable_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);
      PIDX_variable_create(variable_name, values_per_sample[var] * sizeof(uint64_t) * 8, data_type, &variable[var]);
      PIDX_append_and_write_variable(file, variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);
    }
    PIDX_close(file);

    for(var = 0; var < args.variable_count; var++)
    {
      free(double_data[var]);
      double_data[var] = 0;
    }
    free(double_data);
    double_data = 0;

  }
  PIDX_time_step_caching_OFF();
  
  free(args.output_file_name);
#endif
  return 0;
}

/*   prints usage instructions   */
void usage_serial(void) 
{
  printf("Usage: test-PIDX-one-var-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}
