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

int test_one_var_writer(struct Args args, int rank, int nprocs) 
{
#if 0
#if PIDX_HAVE_MPI
  int i = 0, j = 0, k = 0;
  int spv = 0;
  int ts;
  int slice;
  int sub_div[3], local_offset[3];

  /* IDX file */
  PIDX_file file = 0;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 256;                                               // Total number of blocks per file
  
  /* IDX variables */
  PIDX_variable variable = 0;                                       // variable descriptor
  double     *var1_double_scalar_data;
  int sample_count = 1;

  //The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  //   Creating the filename 
  args.output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(args.output_file_name, "%s%s", args.output_file_template, ".idx");

  //   Calculating every process's offset and count  
  sub_div[0] = (args.extents[0] / args.count_local[0]);
  sub_div[1] = (args.extents[1] / args.count_local[1]);
  sub_div[2] = (args.extents[2] / args.count_local[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * args.count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * args.count_local[1];
  local_offset[0] = (slice % sub_div[0]) * args.count_local[0];

  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  //PIDX_create_point(&global_bounding_box);
  //PIDX_create_point(&local_offset_point);
  //PIDX_create_point(&local_box_count_point);
  
  PIDX_set_point_5D(global_bounding_box, args.extents[0], args.extents[1], args.extents[2], 1, 1);
  PIDX_set_point_5D(local_offset_point, local_offset[0], local_offset[1], local_offset[2], 0, 0);
  PIDX_set_point_5D(local_box_count_point, args.count_local[0], args.count_local[1], args.count_local[2], 1, 1);
  
  output_file = args.output_file_name;
    
  for (ts = 0; ts < args.time_step_count; ts++) 
  {
    var1_double_scalar_data = malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2] * sample_count);
    for (k = 0; k < args.count_local[2]; k++)
      for (j = 0; j < args.count_local[1]; j++)
	for (i = 0; i < args.count_local[0]; i++) 
	{
	  int64_t index = (int64_t) (args.count_local[0] * args.count_local[1] * k) + (args.count_local[0] * j) + i;
	  for (spv = 0; spv < sample_count; spv++)
	    var1_double_scalar_data[index * sample_count + spv] = (args.extents[0] * args.extents[1] * (local_offset[2] + k)) + (args.extents[0] * (local_offset[1] + j)) + (local_offset[0] + i);
	}
    
    PIDX_access access;
    PIDX_create_access(&access);

#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, 1, 1, 1, MPI_COMM_WORLD);
#endif
    
    PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    
    PIDX_variable_create(file, "var1_double_scalar_data", sample_count * sizeof(double) * 8, "1*float64", &variable);
    PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, var1_double_scalar_data, PIDX_row_major);
    //PIDX_variable_set_box_metadata_on(variable);
    
    PIDX_close(file);
    PIDX_close_access(access);
    free(var1_double_scalar_data);
  }
  
  //PIDX_delete_point(&global_bounding_box);
  //PIDX_delete_point(&local_offset_point);
  //PIDX_delete_point(&local_box_count_point);
  
  free(args.output_file_name);
#endif
#endif
  return 0;
}

/*   prints usage instructions   */
void usage_one_var_writer(void) 
{
  printf("Usage: test-PIDX-one-var-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}
