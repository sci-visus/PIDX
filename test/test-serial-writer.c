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
#include <PIDX.h>

int serial_writer(struct Args args)
{
  int i = 0, j = 0, k = 0;
  int spv = 0;
  int ts;

  /* IDX file */
  PIDX_file file = 0;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 256;                                               // Total number of blocks per file
  
  /* IDX variables */
  PIDX_variable variable = 0;                                       // variable descriptor
  double     *var1_double_scalar_data;
  int sample_count = 1;
  //int one_opt = 0;
  //int variable_count = 1;

  //   Creating the filename
  args.output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(args.output_file_name, "%s%s", args.output_file_template, ".idx");

  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  //PIDX_create_point(&global_bounding_box);
  //PIDX_create_point(&local_offset_point);
  //PIDX_create_point(&local_box_count_point);
  
  PIDX_set_point_5D(args.extents[0], args.extents[1], args.extents[2], 1, 1, global_bounding_box);
  PIDX_set_point_5D(0, 0, 0, 0, 0, local_offset_point);
  PIDX_set_point_5D(args.extents[0], args.extents[1], args.extents[2], 1, 1, local_box_count_point);
  
  output_file = args.output_file_name;
    
  for (ts = 0; ts < args.time_step; ts++) 
  {
    var1_double_scalar_data = malloc(sizeof (double) * args.extents[0] * args.extents[1] * args.extents[2] * sample_count);
    for (k = 0; k < args.extents[2]; k++)
      for (j = 0; j < args.extents[1]; j++)
	for (i = 0; i < args.extents[0]; i++) 
	{
	  long long index = (long long) (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i;
	  for (spv = 0; spv < sample_count; spv++)
	    var1_double_scalar_data[index * sample_count + spv] = (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i;
	}
    
    PIDX_file_create(output_file, PIDX_file_trunc, NULL, &file);
    PIDX_set_dims(file, global_bounding_box);
    
    PIDX_set_current_time_step(file, ts);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    
    PIDX_variable_create(file, "var1_double_scalar_data", sample_count * sizeof(double) * 8, "1*float64", &variable);
    PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, var1_double_scalar_data, PIDX_row_major);
    //PIDX_variable_set_box_metadata_on(variable);
    
    PIDX_close(file);
    free(var1_double_scalar_data);
  }
  
  //PIDX_delete_point(&global_bounding_box);
  //PIDX_delete_point(&local_offset_point);
  //PIDX_delete_point(&local_box_count_point);
  
  free(args.output_file_name);
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
