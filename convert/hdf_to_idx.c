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

#include <PIDX.h>
#include <Generic_data_structs.h>

static char output_file_template[512];
static char *output_file_name;
static int time_step;

static double* buffer;

static void dump_data(char* filename);

static void dump_data(char* filename)
{
  FILE  *stream;
  int    numread;
  size_t nelmts = 512*256*256;
  char  *fname=NULL;

  if (filename == NULL)
  {
    usage();
    exit(1);
  }

  fname = strdup(filename);

  int i = 0;
  if( (stream = fopen(fname, "rb" )) != NULL )
  {
    numread = fread( buffer, sizeof( double ), nelmts, stream );
    printf( "Number of items read = %d\n", numread );

    //for (i = 0; i < nelmts; i++)
    //{
      //printf("%f ",buffer[i]);
    //}
    //printf("\n");

    fclose( stream );
  }
  else
    printf( "File %s could not be opened\n",fname );

  free(fname);
}

int main(int argc, char **argv) 
{
  int ret;
  PIDX_file file;                                                
  const char *output_file;                                       
  const int bits_per_block = 15;                                 
  const int blocks_per_file = 256;                               
  
  PIDX_variable variable;                                       
  
  buffer = (double*)malloc(sizeof(double) * 512 * 256 * 256);
  memset(buffer, 0, sizeof(double) * 512 * 256 * 256);
  
  
  char EXE_command[2048];
  char EXE_command_delete[2048];
  int time = 0;
  char binary_file[1024];
  
  output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(output_file_name, "%s%s", "V", ".idx");
    
  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
    
  PIDX_set_point_5D(512, 256, 256, 1, 1, global_bounding_box);
  PIDX_set_point_5D(0, 0, 0, 0, 0, local_offset_point);
  PIDX_set_point_5D(512, 256, 256, 1, 1, local_box_count_point);
  output_file = output_file_name;
    
  for(time = 0; time < 100; time++)
  {
    sprintf(binary_file, "/project/v10/DATA/%02d.bin", time);
    sprintf(EXE_command, "h5dump --binary=LE --dataset=/Species/O2 --output=%s /project/v10/K3DR2/Post2D_h2air_flame_1.%02d00E-03/Solution.h5", binary_file, time);
    printf("%s\n", EXE_command);
   
    system(EXE_command);
    dump_data(binary_file);  
    
    printf("Finished reading Data %s\n", binary_file);
    
    PIDX_access access;
    PIDX_create_access(&access);

#if PIDX_HAVE_MPI
    PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#else
    PIDX_set_default_access(access);
#endif
      
    PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, 0);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    
    PIDX_variable_create(file, "V", sizeof(double) * 8, "1*float64", &variable);
    PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
    
    PIDX_close(&file);
    PIDX_close_access(&access);
    
    sprintf(EXE_command_delete, "rm -rf %s\n", binary_file);
    system(EXE_command_delete);
  }
  
  free(buffer);
  buffer = 0;
  
  return 0;
}