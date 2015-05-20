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

int serial_reader(struct Args args)
{
#if 0
  int i = 0, j = 0, k = 0;
  int spv = 0, var = 0;
  int ts, variable_count;

  /* IDX file */
  PIDX_file file = 0;                                                // IDX file descriptor
  const char *output_file;                                                      // IDX File Name
  
  /* IDX variables */
  PIDX_variable *variable;                                       // variable descriptor
  int *values_per_sample;
  double** double_data;
  
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
    PIDX_access access;
    PIDX_create_access(&access);
    //PIDX_set_default_access(access);
    
    /// PIDX mandatory calls
    //PIDX_file_open(output_file, PIDX_file_rdonly, access, &file);
    
    /// PIDX calls to set different parameters (optional)
    PIDX_set_current_time_step(file, ts);
    PIDX_get_variable_count(file, &variable_count);
    
    variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);
    values_per_sample = (int*) malloc(sizeof(*values_per_sample) * variable_count);
    memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);      
    double_data = (double**)malloc(sizeof(*double_data) * variable_count);
    memset(double_data, 0, sizeof(*double_data) * variable_count);
    
    for(var = 0; var < variable_count; var++)
    {
      PIDX_get_next_variable(file, &variable[var]);
      values_per_sample[var] = variable[var]->values_per_sample;
      
      double_data[var] = (double*)malloc(sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * variable[var]->values_per_sample);
      memset(double_data[var], 0, sizeof (double) * args.count_local[0] * args.count_local[1] * args.count_local[2]  * variable[var]->values_per_sample);
      
      PIDX_read_next_variable(variable[var], local_offset_point, local_box_count_point, double_data[var], PIDX_row_major);
    }
    
    //PIDX_variable_set_box_metadata_on(variable);
    
    PIDX_close(file);
    
    int elem_count = 0, non_elem_count = 0;
    for (var = 0; var < variable_count; var++)
    {
      printf("extents %lld %lld %lld\n", (long long)args.extents[0], (long long)args.extents[1], (long long)args.extents[2]);
      for (k = 0; k < args.extents[2]; k++)
        for (j = 0; j < args.extents[1]; j++)
          for (i = 0; i < args.extents[0]; i++) 
          {
            int64_t index = (int64_t) (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i;
            for (spv = 0; spv < values_per_sample[var]; spv++)
            {
              if (double_data[var][index * values_per_sample[var] + spv] == (double)100 + (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i)
                elem_count++;
              else
              {
                //printf("elements %f %f\n", double_data[var][index * values_per_sample[var] + spv], (double)100 + (args.extents[0] * args.extents[1] * k) + (args.extents[0] * j) + i);
                non_elem_count++;
              }
            }
          }
    }
    printf("element count %d not equal count %d\n", elem_count, non_elem_count);
    //free(var1_double_scalar_data);
  }
  PIDX_time_step_caching_OFF();
  
  free(args.output_file_name);
#endif
  return 0;
}

/*   prints usage instructions   */
void usage_serial_reader(void) 
{
  printf("Usage: test-PIDX-one-var-writer -g 4x4x4 -l 2x2x2 -f Filename_ -t 4\n");
  printf("  -g: global dimensions\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -f: IDX Filename\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}
