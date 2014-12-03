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

int test_writer(struct Args args, int rank, int nprocs) 
{
  /*
  int i = 0, j = 0, k = 0, u = 0, v = 0, p = 0, var, d;
  int spv = 0;
  int ts;
  int slice;
  int sub_div[5], offset_local[5];

  PIDX_file_descriptor pidx_ptr;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  const int *gextent;                                                           // Global Extensions of the dataset (64 64 64 0 0)
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 32;                                               // Total number of blocks per file
  
  int variable_count = 4;                                                            // Number of variables
  PIDX_variable_descriptor *variable_ptr;                                       // variable descriptor
  int *var_sample_counts;                                                       // values per variable (example, scalar=1, vector=3)
  int *var_patch_count;                                                         // Number of patches a variable has (data blocks per variable that needs to be written to).
  
  int      ***var_count;                                                              // Local extents of the variables in each process
  int      ***var_offset;                                                              // Local counts of the variables in each process
  
  double     **var1_double_scalar_data;
  int      **var2_int_scalar_data;
  float    **var3_float_scalar_data;
  int   **var4_char_scalar_data;
  
  
  variable_ptr = malloc(sizeof (PIDX_variable_descriptor) * variable_count);        // variable pointers
  memset(variable_ptr, 0, sizeof (PIDX_variable_descriptor) * variable_count);
  
  var_sample_counts = (int*) malloc(sizeof (int) * variable_count);                  
  var_patch_count = (int*) malloc(sizeof (int) * variable_count);                  
  for (var = 0; var < variable_count; var++)
  {
    if(var % 4 == 0)
    {
      var_patch_count[var] = 1;
      var_sample_counts[var] = 1;
    }
    else if(var % 4 == 1)
    {
      var_patch_count[var] = 2;
      var_sample_counts[var] = 1;
    }
    else if(var % 4 == 2)
    {
      var_patch_count[var] = 4;
      var_sample_counts[var] = 3;
    }
    else if(var % 4 == 3)
    {
      var_patch_count[var] = 8;
      var_sample_counts[var] = 4;
    }
  }

  //The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
    
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

  var_count = malloc(sizeof(int**) * variable_count);
  var_offset = malloc(sizeof(int**) * variable_count);
  
  for(var = 0; var < variable_count; var++)
  {
    var_count[var] = malloc(sizeof(int*) * var_patch_count[var]);
    var_offset[var] = malloc(sizeof(int*) * var_patch_count[var]);
    for(i = 0; i < var_patch_count[var] ; i++)
    {
      var_count[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
      var_offset[var][i] = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);    
    }
  }
  
  for(var = 0; var < variable_count; var++)
  {
    if(var % 4 == 0)                                                       // One patch for this variable
    {
      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
	var_count[var][0][d] = args.count_local[d];
	var_offset[var][0][d] = offset_local[d];
      }
    }
    else if(var % 4 == 1)                                                  // two patches for this variable
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
    else if(var % 4 == 2)
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
    else if(var % 4 == 3)
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
    
  for (ts = 0; ts < args.time_step; ts++) 
  {
    //printf("Variable Count %d\n", variable_count);
    for(var = 0; var < variable_count; var++)
    {
      if(var % 4 == 0)
	var1_double_scalar_data = malloc(sizeof(double*) * var_patch_count[var]);
      if(var % 4 == 1)
	var2_int_scalar_data = malloc(sizeof(int*) * var_patch_count[var]);
      if(var % 4 == 2)
	var3_float_scalar_data = malloc(sizeof(float*) * var_patch_count[var]);
      if(var % 4 == 3)
	var4_char_scalar_data = malloc(sizeof(int*) * var_patch_count[var]);
	
      for(p = 0 ; p < var_patch_count[var] ; p++)
      {
	if(var % 4 == 0)
	{
	  var1_double_scalar_data[p] = malloc(sizeof (double) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * var_sample_counts[var]);
	  for (v = 0; v < var_count[var][p][4]; v++)
	    for (u = 0; u < var_count[var][p][3]; u++)
	      for (k = 0; k < var_count[var][p][2]; k++)
		for (j = 0; j < var_count[var][p][1]; j++)
		  for (i = 0; i < var_count[var][p][0]; i++) 
		  {
		    long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				      (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				      (var_count[var][p][0] * j) + i;
		    for (spv = 0; spv < var_sample_counts[var]; spv++)
		      var1_double_scalar_data[p][index * var_sample_counts[var] + spv] = (100 + 
			((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + 
			(args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + 
			(args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i)));
		  }
	}
	else if(var % 4 == 1)
	{
	  var2_int_scalar_data[p] = malloc(sizeof (int) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * var_sample_counts[var]);
	  for (v = 0; v < var_count[var][p][4]; v++)
	    for (u = 0; u < var_count[var][p][3]; u++)
	      for (k = 0; k < var_count[var][p][2]; k++)
		for (j = 0; j < var_count[var][p][1]; j++)
		  for (i = 0; i < var_count[var][p][0]; i++) 
		  {
		    long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				      (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				      (var_count[var][p][0] * j) + i;
		    for (spv = 0; spv < var_sample_counts[var]; spv++)
		      var2_int_scalar_data[p][index * var_sample_counts[var] + spv] = rank; //100 + ((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + (args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + (args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i));
		  }
	}
	else if(var % 4 == 2)
	{
	  var3_float_scalar_data[p] = malloc(sizeof (float) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * var_sample_counts[var]);
	  for (v = 0; v < var_count[var][p][4]; v++)
	    for (u = 0; u < var_count[var][p][3]; u++)
	      for (k = 0; k < var_count[var][p][2]; k++)
		for (j = 0; j < var_count[var][p][1]; j++)
		  for (i = 0; i < var_count[var][p][0]; i++) 
		  {
		    long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				      (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				      (var_count[var][p][0] * j) + i;
		    for (spv = 0; spv < var_sample_counts[var]; spv++)
		      var3_float_scalar_data[p][index * var_sample_counts[var] + spv] =  ((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + (args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + (args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i));
		  }
	}
	else if(var % 4 == 3)
	{
	  var4_char_scalar_data[p] = malloc(sizeof (int) * var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * var_count[var][p][4] * var_sample_counts[var]);
	  for (v = 0; v < var_count[var][p][4]; v++)
	    for (u = 0; u < var_count[var][p][3]; u++)
	      for (k = 0; k < var_count[var][p][2]; k++)
		for (j = 0; j < var_count[var][p][1]; j++)
		  for (i = 0; i < var_count[var][p][0]; i++) 
		  {
		    long long index = (long long) (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * var_count[var][p][3] * v) + 
				      (var_count[var][p][0] * var_count[var][p][1] * var_count[var][p][2] * u) + (var_count[var][p][0] * var_count[var][p][1] * k) + 
				      (var_count[var][p][0] * j) + i;
		    for (spv = 0; spv < var_sample_counts[var]; spv++)
		      var4_char_scalar_data[p][index * var_sample_counts[var] + spv] = p;//((args.extents[0] * args.extents[1] * args.extents[2] * args.extents[3] * (var_offset[var][p][4] + v)) + (args.extents[0] * args.extents[1] * args.extents[2] * (var_offset[var][p][3] + u)) + (args.extents[0] * args.extents[1]*(var_offset[var][p][2] + k))+(args.extents[0]*(var_offset[var][p][1] + j)) + (var_offset[var][p][0] + i));
		  }
	}
      }
    }
    
    pidx_ptr = PIDX_file_create(output_file, gextent, 0);
    PIDX_set_current_time_step(pidx_ptr, ts);
    PIDX_set_block_size(pidx_ptr, bits_per_block);
    PIDX_set_block_count(pidx_ptr, blocks_per_file);
    
    variable_ptr[0] = PIDX_variable_create(pidx_ptr, "var1_double_scalar_data", var_sample_counts[0], sizeof(double) * 8, "float64");
    for(p = 0 ; p < var_patch_count[0] ; p++)
      PIDX_write_variable(variable_ptr[0], var_offset[0][p], var_count[0][p], var1_double_scalar_data[p], "row");
    
    //PIDX_sync(pidx_ptr);
    
    variable_ptr[1] = PIDX_variable_create(pidx_ptr, "var2_int_scalar_data", var_sample_counts[1], sizeof(int) * 8, "int32");
    for(p = 0 ; p < var_patch_count[1] ; p++)
      PIDX_write_variable(variable_ptr[1], var_offset[1][p], var_count[1][p], var2_int_scalar_data[p], "row");
   
    //PIDX_sync(pidx_ptr);
    
    variable_ptr[2] = PIDX_variable_create(pidx_ptr, "var3_float_scalar_data", var_sample_counts[2], sizeof(float) * 8, "float32");
    for(p = 0 ; p < var_patch_count[2] ; p++)
      PIDX_write_variable(variable_ptr[2], var_offset[2][p], var_count[2][p], var3_float_scalar_data[p], "row");
   
    //PIDX_sync(pidx_ptr);
    
    variable_ptr[3] = PIDX_variable_create(pidx_ptr, "var4_int_scalar_data", var_sample_counts[3], sizeof(int) * 8, "int32");
    for(p = 0 ; p < var_patch_count[3] ; p++)
      PIDX_write_variable(variable_ptr[3], var_offset[3][p], var_count[3][p], var4_char_scalar_data[p], "row");
        
    //PIDX_sync(pidx_ptr);
    
    PIDX_close(pidx_ptr);
    
    for(p = 0; p < var_patch_count[0]; p++)
      free(var1_double_scalar_data[p]);
    for(p = 0; p < var_patch_count[1]; p++)
      free(var2_int_scalar_data[p]); 
    for(p = 0; p < var_patch_count[2]; p++)
      free(var3_float_scalar_data[p]);
    for(p = 0; p < var_patch_count[3]; p++)
      free(var4_char_scalar_data[p]);
    
    free(var1_double_scalar_data);
    free(var2_int_scalar_data); 
    free(var3_float_scalar_data);
    free(var4_char_scalar_data);
  }
  
  free(args.output_file_name);
  free(variable_ptr);
  variable_ptr = 0;
  free(var_sample_counts);
  var_sample_counts = 0;
  
  */
  
  return 0;
}

/*   prints usage instructions   */
void usage_writer(void) 
{
  printf("Usage: test-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}


