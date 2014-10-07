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
#include <PIDX.h>
#include <Generic_data_structs.h>

int test_multi_var_writer(struct Args args, int rank, int nprocs) 
{
  /*
  int i = 0, j = 0, k = 0, u = 0, v = 0;
  int ts, var, spv;
  int slice;
  int variable_count;
  int sub_div[5], offset_local[5];

  PIDX_file_descriptor pidx_ptr;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  const int *gextent;                                                           // Global Extensions of the dataset (64 64 64 0 0)
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 256;                                               // Total number of blocks per file
  
  PIDX_variable_descriptor* variable_ptr;                                       // variable descriptor
  double     **double_data;
  int* values_per_sample;
  
  //The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  variable_count = args.variable_count;
  variable_ptr = malloc(sizeof(*variable_ptr) * variable_count);
  memset(variable_ptr, 0, sizeof(*variable_ptr) * variable_count);
  
  values_per_sample = malloc(sizeof(*values_per_sample) * variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);
  
  //   Creating the filename 
  args.output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(args.output_file_name, "%s%s", args.output_file_template, ".idx");

  //   Calculating every process's offset and count  
  gextent = args.extents;
  sub_div[0] = (args.extents[0] / args.count_local[0]);
  sub_div[1] = (args.extents[1] / args.count_local[1]);
  sub_div[2] = (args.extents[2] / args.count_local[2]);
  offset_local[2] = (rank / (sub_div[0] * sub_div[1])) * args.count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  offset_local[1] = (slice / sub_div[0]) * args.count_local[1];
  offset_local[0] = (slice % sub_div[0]) * args.count_local[0];

  offset_local[3] = 0;
  offset_local[4] = 0;
  args.count_local[3] = 1;
  args.count_local[4] = 1;

  output_file = args.output_file_name;    
  for (ts = 0; ts < args.time_step; ts++) 
  {
    double_data = malloc(sizeof(*double_data) * variable_count);
    memset(double_data, 0, sizeof(*double_data) * variable_count);
    
    for(var = 0; var < variable_count; var++)
    {
      values_per_sample[var] = var + 1;
      double_data[var] = malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2] * args.count_local[3] * args.count_local[4] * values_per_sample[var]);
      for (v = 0; v < args.count_local[4]; v++)
	for (u = 0; u < args.count_local[3]; u++)
	  for (k = 0; k < args.count_local[2]; k++)
	    for (j = 0; j < args.count_local[1]; j++)
	      for (i = 0; i < args.count_local[0]; i++) 
	      {
		long long index = (long long) (args.count_local[0] * args.count_local[1] * args.count_local[2] * args.count_local[3] * v) + 
				  (args.count_local[0] * args.count_local[1] * args.count_local[2] * u) + (args.count_local[0] * args.count_local[1] * k) + 
				  (args.count_local[0] * j) + i;
		for (spv = 0; spv < values_per_sample[var]; spv++)
		  double_data[var][index * values_per_sample[var] + spv] = (100 + 
		    ((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (offset_local[4] + v)) + 
		    (args.extents[0] * args.extents[1] * args.extents[2] * (offset_local[3] + u)) + 
		    (args.extents[0] * args.extents[1]*(offset_local[2] + k))+(args.extents[0]*(offset_local[1] + j)) + (offset_local[0] + i)));
	      }
    }
    
    pidx_ptr = PIDX_file_create(output_file, gextent, 0);
    PIDX_set_communicator(pidx_ptr, MPI_COMM_WORLD);
    PIDX_set_current_time_step(pidx_ptr, ts);
    PIDX_set_block_size(pidx_ptr, bits_per_block);
    PIDX_set_block_count(pidx_ptr, blocks_per_file);
    
    for(var = 0; var < variable_count; var++)
    {
      variable_ptr[var] = PIDX_variable_create(pidx_ptr, "var1_double_scalar_data", values_per_sample[var], sizeof(double) * 8, "float64");
      PIDX_write_variable(variable_ptr[var], offset_local, args.count_local, double_data[var], "row");
    }
    
    PIDX_close(pidx_ptr);
    
    for(var = 0; var < variable_count; var++)
    {
      free(double_data[var]);
      double_data[var] = 0;
    }
    
    free(double_data);
    double_data = 0;
  }
  
  free(variable_ptr);
  free(values_per_sample);
  free(args.output_file_name);
  */
  return 0;
}

/*   prints usage instructions   */
void usage_multi_var_writer(void) 
{
  printf("Usage: test-multi-var-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}