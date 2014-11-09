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

int test_multi_patch_writer(struct Args args, int rank, int nprocs) 
{
  int i = 0, j = 0, k = 0, u = 0, v = 0, p = 0, var, d;
  int spv = 0;
  int ts;
  int slice;
  int sub_div[5], offset_local[5];

  PIDX_file file;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  const int *gextent;                                                           // Global Extensions of the dataset (64 64 64 0 0)
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 32;                                               // Total number of blocks per file
  
  PIDX_variable *variable;                                       // variable descriptor
  int *values_per_sample;                                                       // values per variable (example, scalar=1, vector=3)
  int *var_patch_count;                                                         // Number of patches a variable has (data blocks per variable that needs to be written to).
  
  int      ***var_count;                                                              // Local extents of the variables in each process
  int      ***var_offset;                                                              // Local counts of the variables in each process
  
  double   ***var_double_scalar_data;
  
  //The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  variable = malloc(sizeof (PIDX_variable) * args.variable_count);        // variable pointers
  memset(variable, 0, sizeof (PIDX_variable) * args.variable_count);
  
  values_per_sample = (int*) malloc(sizeof (int) * args.variable_count);                  
  var_patch_count = (int*) malloc(sizeof (int) * args.variable_count);                  
  for (var = 0; var < args.variable_count; var++)
  {
    values_per_sample[var] = 1;
    
    if(var % 4 == 3)
      var_patch_count[var] = 1;
    else if(var % 4 == 2)
      var_patch_count[var] = 2;
    else if(var % 4 == 1)
      var_patch_count[var] = 4;
    else if(var % 4 == 0)
      var_patch_count[var] = 8;
  }
  
  
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

  var_count = malloc(sizeof(int**) * args.variable_count);
  var_offset = malloc(sizeof(int**) * args.variable_count);
  
  for(var = 0; var < args.variable_count; var++)
  {
    var_count[var] = malloc(sizeof(int*) * var_patch_count[var]);
    var_offset[var] = malloc(sizeof(int*) * var_patch_count[var]);
    for(i = 0; i < var_patch_count[var] ; i++)
    {
      var_count[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
      var_offset[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);    
    }
  }
  
  for(var = 0; var < args.variable_count; var++)
  {
    if(var % 4 == 3)                                                       // One patch for this variable
    {
      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
	var_count[var][0][d] = args.count_local[d];
	var_offset[var][0][d] = offset_local[d];
      }
    }
    else if(var % 4 == 2)                                                  // two patches for this variable
    {
      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
	var_count[var][0][d] = args.count_local[d];
	var_offset[var][0][d] = offset_local[d];
	var_count[var][1][d] = args.count_local[d];
	var_offset[var][1][d] = offset_local[d];
      }
      
      var_count[var][0][0] = args.count_local[0]/2;
      if(args.count_local[0] % 2 == 0)
	var_count[var][1][0] = args.count_local[0]/2;
      else
	var_count[var][1][0] = args.count_local[0]/2 + 1;
      
      var_offset[var][1][0] = offset_local[0] + args.count_local[0]/2;
    }
    else if(var % 4 == 1)
    {
      for(i = 0; i < var_patch_count[var]; i++)
      {
	for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
	{
	  var_count[var][i][d] = args.count_local[d];
	  var_offset[var][i][d] = offset_local[d];
	}
      }
      var_count[var][0][0] = args.count_local[0]/2;
      var_count[var][0][1] = args.count_local[1]/2;
      
      var_count[var][1][1] = args.count_local[1]/2;
      if(args.count_local[0] % 2 == 0)
      {
	var_count[var][1][0] = args.count_local[0]/2;
	var_count[var][3][0] = args.count_local[0]/2;
	var_offset[var][1][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][3][0] = var_offset[var][0][0] + args.count_local[0]/2;
      }
      else
      {
	var_count[var][1][0] = args.count_local[0]/2 + 1;
	var_count[var][3][0] = args.count_local[0]/2 + 1;
	var_offset[var][1][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][3][0] = var_offset[var][0][0] + args.count_local[0]/2;
      }
      
      var_count[var][2][0] = args.count_local[0]/2;
      if(args.count_local[1] % 2 == 0)
      {
	var_count[var][2][1] = args.count_local[1]/2;
	var_count[var][3][1] = args.count_local[1]/2;
	var_offset[var][2][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][3][1] = var_offset[var][0][1] + args.count_local[1]/2;
      }
      else
      {
	var_count[var][2][1] = args.count_local[1]/2 + 1;
	var_count[var][3][1] = args.count_local[1]/2 + 1;
	var_offset[var][2][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][3][1] = var_offset[var][0][1] + args.count_local[1]/2;
      }
    }
    else if(var % 4 == 0)
    {
      for(i = 0; i < var_patch_count[var]; i++)
      {
	for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
	{
	  var_count[var][i][d] = args.count_local[d];
	  var_offset[var][i][d] = offset_local[d];
	}
      }
      var_count[var][0][0] = args.count_local[0]/2;
      var_count[var][0][1] = args.count_local[1]/2;
      
      var_count[var][4][0] = args.count_local[0]/2;
      var_count[var][4][1] = args.count_local[1]/2;
      
      var_count[var][1][1] = args.count_local[1]/2;
      var_count[var][5][1] = args.count_local[1]/2;
      if(args.count_local[0] % 2 == 0)
      {
	var_count[var][1][0] = args.count_local[0]/2;
	var_count[var][3][0] = args.count_local[0]/2;
	var_offset[var][1][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][3][0] = var_offset[var][0][0] + args.count_local[0]/2;
	
	var_count[var][5][0] = args.count_local[0]/2;
	var_count[var][7][0] = args.count_local[0]/2;
	var_offset[var][5][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][7][0] = var_offset[var][0][0] + args.count_local[0]/2;
      }
      else
      {
	var_count[var][1][0] = args.count_local[0]/2 + 1;
	var_count[var][3][0] = args.count_local[0]/2 + 1;
	var_offset[var][1][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][3][0] = var_offset[var][0][0] + args.count_local[0]/2;
	
	var_count[var][5][0] = args.count_local[0]/2 + 1;
	var_count[var][7][0] = args.count_local[0]/2 + 1;
	var_offset[var][5][0] = var_offset[var][0][0] + args.count_local[0]/2;
	var_offset[var][7][0] = var_offset[var][0][0] + args.count_local[0]/2;
      }
      
      var_count[var][2][0] = args.count_local[0]/2;
      var_count[var][6][0] = args.count_local[0]/2;
      
      if(args.count_local[1] % 2 == 0)
      {
	var_count[var][2][1] = args.count_local[1]/2;
	var_count[var][3][1] = args.count_local[1]/2;
	var_offset[var][2][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][3][1] = var_offset[var][0][1] + args.count_local[1]/2;
	
	var_count[var][6][1] = args.count_local[1]/2;
	var_count[var][7][1] = args.count_local[1]/2;
	var_offset[var][6][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][7][1] = var_offset[var][0][1] + args.count_local[1]/2;
      }
      else
      {
	var_count[var][2][1] = args.count_local[1]/2 + 1;
	var_count[var][3][1] = args.count_local[1]/2 + 1;
	var_offset[var][2][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][3][1] = var_offset[var][0][1] + args.count_local[1]/2;
	
	var_count[var][6][1] = args.count_local[1]/2 + 1;
	var_count[var][7][1] = args.count_local[1]/2 + 1;
	var_offset[var][6][1] = var_offset[var][0][1] + args.count_local[1]/2;
	var_offset[var][7][1] = var_offset[var][0][1] + args.count_local[1]/2;
      }
      
      var_count[var][0][2] = args.count_local[2]/2;
      var_count[var][1][2] = args.count_local[2]/2;
      var_count[var][2][2] = args.count_local[2]/2;
      var_count[var][3][2] = args.count_local[2]/2;
      if(args.count_local[1] % 2 == 0)
      {
	var_count[var][4][2] = args.count_local[2]/2;
	var_count[var][5][2] = args.count_local[2]/2;
	var_count[var][6][2] = args.count_local[2]/2;
	var_count[var][7][2] = args.count_local[2]/2;
	
	var_offset[var][4][2] = var_offset[var][0][2] + args.count_local[2]/2;
	var_offset[var][5][2] = var_offset[var][1][2] + args.count_local[2]/2;
	var_offset[var][6][2] = var_offset[var][2][2] + args.count_local[2]/2;
	var_offset[var][7][2] = var_offset[var][3][2] + args.count_local[2]/2;
      }
      else
      {
	var_count[var][4][2] = args.count_local[2]/2 + 1;
	var_count[var][5][2] = args.count_local[2]/2 + 1;
	var_count[var][6][2] = args.count_local[2]/2 + 1;
	var_count[var][7][2] = args.count_local[2]/2 + 1;
	
	var_offset[var][4][2] = var_offset[var][0][2] + args.count_local[2]/2;
	var_offset[var][5][2] = var_offset[var][1][2] + args.count_local[2]/2;
	var_offset[var][6][2] = var_offset[var][2][2] + args.count_local[2]/2;
	var_offset[var][7][2] = var_offset[var][3][2] + args.count_local[2]/2;
      }
    }
  }
  
  output_file = args.output_file_name;
  
  PIDX_point global_bounding_box, **local_offset_point, **local_box_count_point;
  
  local_offset_point = malloc(sizeof(PIDX_point*) * args.variable_count);
  local_box_count_point = malloc(sizeof(PIDX_point*) * args.variable_count);
  for(var = 0; var < args.variable_count; var++)
  {
    local_offset_point[var] = malloc(sizeof(PIDX_point) * var_patch_count[var]);
    local_box_count_point[var] = malloc(sizeof(PIDX_point) * var_patch_count[var]);
    for(p = 0 ; p < var_patch_count[var] ; p++)
    {
      PIDX_set_point_5D((long long)var_offset[var][p][0], (long long)var_offset[var][p][1], (long long)var_offset[var][p][2], 0, 0, local_offset_point[var][p]);
      PIDX_set_point_5D((long long)var_count[var][p][0], (long long)var_count[var][p][1], (long long)var_count[var][p][2], 1, 1, local_box_count_point[var][p]);       
    }
  }
  
  PIDX_set_point_5D((long long)args.extents[0], (long long)args.extents[1], (long long)args.extents[2], 1, 1, global_bounding_box);
  
  for (ts = 0; ts < args.time_step; ts++) 
  {
    PIDX_access access;
    PIDX_create_access(&access);

#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif
    
    PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    PIDX_set_variable_count(file, args.variable_count);
    
#if 0
    var_double_scalar_data = malloc(sizeof(double**) * args.variable_count);
    for (var = 0; var < args.variable_count; var++)
    {
      var_double_scalar_data[var] = malloc(sizeof(double*) * var_patch_count[var]);
      for(p = 0 ; p < var_patch_count[var] ; p++)
      {
	var_double_scalar_data[var][p] = malloc(sizeof (double) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * values_per_sample[var]);
	for (v = 0; v < var_count[var][p][4]; v++)
	  for (u = 0; u < var_count[var][p][3]; u++)
	    for (k = 0; k < var_count[var][p][2]; k++)
	      for (j = 0; j < var_count[var][p][1]; j++)
		for (i = 0; i < var_count[var][p][0]; i++) 
		{
		  long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				    (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				    (var_count[var][p][0] * j) + i;
		  for (spv = 0; spv < values_per_sample[var]; spv++)
		    var_double_scalar_data[var][p][index * values_per_sample[var] + spv] = (100 + 
		      ((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + 
		      (args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + 
		      (args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i)));
		}
      }
    }
    
    char variable_name[512];
    char data_type[512];
    for(var = 0; var < args.variable_count; var++)
    {
      sprintf(variable_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);
      PIDX_variable_create(file, variable_name, values_per_sample[var] * sizeof(double) * 8, data_type, &variable[var]);
      
      for(p = 0 ; p < var_patch_count[var] ; p++)
      {
	if (rank == 0)
	  printf("Writing patch %d\n", p);
	PIDX_append_and_write_variable(variable[var], local_offset_point[var][p], local_box_count_point[var][p], var_double_scalar_data[var][p], PIDX_row_major);
      }
    }
    
    PIDX_close(file);
    PIDX_close_access(access);
    
    for(var = 0; var < args.variable_count; var++)
    {
      for(p = 0; p < var_patch_count[0]; p++)
      {
	free(var_double_scalar_data[var][p]);  
	var_double_scalar_data[var][p] = 0;
      }
      free(var_double_scalar_data[var]);  
      var_double_scalar_data[var] = 0;
    }
#endif

#if 1
    var_double_scalar_data = malloc(sizeof(double**) * args.variable_count);
    for (var = 0; var < args.variable_count; var++)
    {
      var_double_scalar_data[var] = malloc(sizeof(double*) * var_patch_count[var]);
      for(p = 0 ; p < var_patch_count[var] ; p++)
      {
	var_double_scalar_data[var][p] = malloc(sizeof (double) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * values_per_sample[var]);
	for (v = 0; v < var_count[var][p][4]; v++)
	  for (u = 0; u < var_count[var][p][3]; u++)
	    for (k = 0; k < var_count[var][p][2]; k++)
	      for (j = 0; j < var_count[var][p][1]; j++)
		for (i = 0; i < var_count[var][p][0]; i++) 
		{
		  long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				    (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				    (var_count[var][p][0] * j) + i;
		  for (spv = 0; spv < values_per_sample[var]; spv++)
		    var_double_scalar_data[var][p][index * values_per_sample[var] + spv] = (100 + 
		      ((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + 
		      (args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + 
		      (args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i)));
		}
      }
    
    
      char variable_name[512];
      char data_type[512];
      sprintf(variable_name, "variable_%d", var);
      sprintf(data_type, "%d*float64", values_per_sample[var]);
      PIDX_variable_create(file, variable_name, values_per_sample[var] * sizeof(double) * 8, data_type, &variable[var]);
      
      for(p = 0 ; p < var_patch_count[var] ; p++)
      {
	if (rank == 0)
	  printf("Writing patch %d\n", p);
	PIDX_append_and_write_variable(variable[var], local_offset_point[var][p], local_box_count_point[var][p], var_double_scalar_data[var][p], PIDX_row_major);
      }
      
      PIDX_flush(file);
      for(p = 0; p < var_patch_count[var]; p++)
      {
	free(var_double_scalar_data[var][p]);  
	var_double_scalar_data[var][p] = 0;
      }
      
      free(var_double_scalar_data[var]);  
      var_double_scalar_data[var] = 0;
    }
    
    PIDX_close(file);
    PIDX_close_access(access);
#endif
    
    free(var_double_scalar_data);  var_double_scalar_data = 0;
  }
  
  for(var = 0; var < args.variable_count; var++)
  {
    free(local_offset_point[var]);
    free(local_box_count_point[var]);
  }
  free(local_offset_point);
  free(local_box_count_point);
  
  free(args.output_file_name);
  free(variable);
  variable = 0;
  free(values_per_sample);
  values_per_sample = 0;
  
  return 0;
}

/*   prints usage instructions   */
void usage_multi_patch_writer(void) 
{
  printf("Usage: test-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}