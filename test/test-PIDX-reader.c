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

int test_reader(struct Args args, int rank, int nprocs) 
{
#if PIDX_HAVE_MPI
  int ts, var;
  //int spv, i = 0, j = 0, k = 0;
  int slice;
  int sub_div[3], local_offset[3];

  /// IDX file descriptor
  PIDX_file file;
  
  /// IDX File Name
  const char *output_file;
  
  /// IDX File variable counts
  int variable_count;
  
  /// Total number of samples in each block = 2 ^ bits_per_block
  int bits_per_block;   
  
  /// Total number of blocks per file
  int blocks_per_file;                                               
  
  int *values_per_sample;
  
  PIDX_variable* variable;
  double     **double_data;
  
  /// The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.idx_count, 3, MPI_INT, 0, MPI_COMM_WORLD);
  
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

  output_file = args.output_file_name;    
  
  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;

  PIDX_set_point_5D(global_bounding_box, args.extents[0], args.extents[1], args.extents[2], 1, 1);
  PIDX_set_point_5D(local_offset_point, local_offset[0], local_offset[1], local_offset[2], 0, 0);
  PIDX_set_point_5D(local_box_count_point, args.count_local[0], args.count_local[1], args.count_local[2], 1, 1);
  
  for (ts = 0; ts < args.time_step; ts++) 
  {
    PIDX_access access;
    PIDX_create_access(&access);

#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, args.idx_count[0], args.idx_count[1], args.idx_count[2], MPI_COMM_WORLD);
    PIDX_set_process_extent(access, sub_div[0], sub_div[1], sub_div[2]);
#else
    PIDX_set_default_access(access);
#endif
      
    PIDX_file_open(output_file, PIDX_file_rdonly, access, &file);
    
    PIDX_get_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, ts);
    PIDX_get_block_size(file, &bits_per_block);
    PIDX_get_block_count(file, &blocks_per_file);
    PIDX_get_variable_count(file, &variable_count);
    
    //if(rank == 0)
    //  printf("Block: %d %d and Variable Count = %d\n", bits_per_block, blocks_per_file, variable_count);
    
    ///
    variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);
    
    ///
    double_data = (double**)malloc(sizeof(*double_data) * variable_count);
    memset(double_data, 0, sizeof(*double_data) * variable_count);
    
    
    values_per_sample = (int*)malloc(sizeof(*values_per_sample) * variable_count);
    memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count); 
    
    for (var = 0; var < variable_count; var++)
    {
      PIDX_get_next_variable(file, &variable[var]);
      //if(rank == 0)
      //	printf("name = %s and size = %d\n", variable[var]->var_name, variable[var]->values_per_sample);
      
      values_per_sample[var] = variable[var]->values_per_sample;
      
      double_data[var] = (double*)malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * variable[var]->values_per_sample);
      memset(double_data[var], 0, sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * variable[var]->values_per_sample);
      
      PIDX_read_next_variable(variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);
    }
    
    PIDX_close(file);
    PIDX_close_access(access);
    int equal_count = 0;
    for(var = 0; var < variable_count; var++)
    {
      int i, j, k, spv;
      for (k = 0; k < args.count_local[2]; k++)
        for (j = 0; j < args.count_local[1]; j++)
          for (i = 0; i < args.count_local[0]; i++)
          {
            int64_t index = (int64_t) (args.count_local[0] * args.count_local[1] * k) + (args.count_local[0] * j) + i;
            for (spv = 0; spv < values_per_sample[var]; spv++)
            {
              if (double_data[var][index * values_per_sample[var] + spv] != 100 + var + ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i)))
              {
                
              }
              else
                equal_count++;
            }
          }
      free(double_data[var]);
      double_data[var] = 0;
    }
    
    printf("[%d] equal element count = %d\n", rank, equal_count);
    
    free(double_data);
    double_data = 0;
    
    free(variable);
    variable = 0;
    
    free(values_per_sample);
    values_per_sample = 0;
  }
  free(args.output_file_name);
#endif
  
  return 0;
}

/*   prints usage instructions   */
void usage_reader(void) 
{
  printf("Usage: test-multi-var-PIDX-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}
