/***************************************************
 ** ViSUS Visualization Project                    **
 ** Copyright (c) 2010 University of Utah          **
 ** Scientific Computing and Imaging Institute     **
 ** 72 S Central Campus Drive, Room 3750           **
 ** Salt Lake City, UT 84112                       **
 **                                                **
 ** For information about this project see:        **
 ** http://www.pascucci.org/visus/                 **
 **                                                **
 **      or contact: pascucci@sci.utah.edu         **
 **                                                **
 ****************************************************/

#define TRUE  1
#define FALSE 0

#define SUCCESS  1
#define FAIL     0

#define MAX_LENGTH 256
#define _DEBUG_MSG_ 1

//===================================================

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "PIDX.h"
#include "mpi.h"

static int parse_args(int argc, char **argv);
static void usage(void);
static int print_error(char *error_message, char* file, int line);

/* global dimensions of 3D volume */
static int extents[5] = {0, 0, 0, 0, 0};

/* per-process dimensions of each sub-block */
static int offset_local[5] = {0, 0, 0, 0, 0};
static int count_local[5] = {0, 0, 0, 0, 0};
static int sub_div[5] = {0, 0, 0, 0, 0};

/* Number of time-steps */
static int nx = 0;
static int ny = 0;
static int nz = 0;

static char *sph_filename;

static char idx_file_template[512] = {0};
static char *idx_file_name;

static double convert_start = 0, convert_end = 0, total_time = 0, max_time = 0, loading_time = 0;

float *alldata;

int main(int argc, char **argv) 
{
    
  int ret;
  int rank, nprocs;

  FILE *sph_file;

  //==========================
  // SPH Data

  int endian_check;       
  int is_big_en;

  int is_scalar;
  int is_vector;
  int is_float;
  int is_double;

  int i_max;
  int j_max;
  int k_max;
  
  float x_orig;
  float y_orig;
  float z_orig;
  
  float time_steps;
  float time_step;
  float time;

  long   data_size;
  float  *float_data, *float_ptr;
  double *double_data, *double_ptr;

  int        int_val;
  int64_t  long_val;
  float      float_val;
  double     double_val;
  //==========================
    
  PIDX_file file;       
  PIDX_variable* variable;                                       // variable descriptor
  char var_name[512];

  int shift;
  int shift_slice;
  int slice = 0;
  #if defined (_SLICE_)
    int slice_dimension;
  #endif
        
  int i = 0, j = 0, k= 0;
  float *alldata_ptr;
        
  const int *gextent; /*Global Extensions of the dataset (64 64 64 0 0)*/
  const int *count; /*Local extents of each process*/
  const int *offset; /*Local counts of each process*/
  const int blocks_per_file = 128; /*Total number of Blocks per file*/
  const int bits_per_block = 15; /*Total number of samples in each block*/

  
  ///SPH Data Reading BEGIN
  
  /// Test Endianess

  endian_check = 1; // 0x00000001

  if (*(char*) & endian_check ) 
  {
    /* little endian. memory image 01 00 00 00 */
    #if defined (_DEBUG_MSG_) 
      printf ("\n <<< LITTLE ENDIAN MACHINE >>> \n\n");
    #endif
  } 
  else
  {
    /* big endian. memory image 00 00 00 01 */
    is_big_en = TRUE;       
    #if defined (_DEBUG_MSG_) 
      printf ("\n <<< BIG ENDIAN MACHINE >>> \n\n");
    #endif
  }
  
  /// Open SPH Binary File
  printf("Starting to load files\n");

  sph_filename = (char *)malloc ( MAX_LENGTH * sizeof(char) ); 
  if ( sph_filename == NULL) 
  {  
    printf("<<< ERROR >> Cannot allocate memory for SPH Filename \n" );
    exit (FAIL);
  }               

  printf("Parsing Arguments\n");
  ret = parse_args(argc, argv);
  if (ret < 0) 
  {
    usage();
    print_error("ret error", __FILE__, __LINE__);
  }

  if ( (sph_file = fopen( sph_filename, "rb")) == NULL ) 
  {
    printf("<<< ERROR >>> Cannot open \"%s\" for reading \n", sph_filename );
    exit(FAIL);
  }

  /// Read Data Attribute
  is_scalar = FALSE;
  is_vector = FALSE;
  is_float  = FALSE;
  is_double = FALSE;
  is_big_en = FALSE;

  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }

  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }
      
  if ( int_val == 1 ) 
    is_scalar = TRUE;
  else if ( int_val == 2 )
    is_vector = TRUE;
  else 
  {
    printf( "<<< ERROR >>> Unrecognized Data Type (1:SCALAR; 2:VECTOR) : %d \n", int_val );
    exit(FAIL);
  }

  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }
      
  if ( int_val == 1 ) 
    is_float = TRUE;
  else if ( int_val == 2 )
    is_double = TRUE;
  else 
    printf( "<<< ERROR >>> Unrecognized Data Type (1:FLOAT; 2:DOUBLE) : %d \n", int_val );

  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }

  #if defined (_DEBUG_MSG_)
    if ( is_scalar )
      printf ("DATA: SCALAR " );
    if ( is_vector )
      printf ("DATA: VECTOR " );
    if ( is_float )
      printf ("(FLOAT)\n" );
    if ( is_double )
      printf ("(DOUBLE)\n\n" );
  #endif

  /// Read Data Size
  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }

  if ( is_float ) 
  {
    if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    i_max = ( int )int_val;
    
    if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    j_max = ( int )int_val;
    
    if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    k_max = ( int )int_val;
  }
  else if ( is_double ) 
  {          
    if ( fread( &long_val, sizeof(int64_t), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    i_max = ( int )long_val;
    
    if ( fread( &long_val, sizeof(int64_t), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    j_max = ( int )long_val;
    
    if ( fread( &long_val, sizeof(int64_t), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    k_max = ( int )long_val;
  }

  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }
      
  #if defined (_DEBUG_MSG_)
    printf ("SIZE: %d x %d x %d \n\n", (int)i_max, (int)j_max, (int)k_max );
  #endif

  /// Read Origin Coordinate of Data
  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }
      
  if ( is_float ) 
  {      
    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    x_orig = ( float )float_val;
    
    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    y_orig = ( float )float_val;

    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    z_orig = ( float )float_val;
  }
  else if ( is_double ) 
  {          
    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    x_orig = ( float )double_val;
    
    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    y_orig = ( float )double_val;

    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    z_orig = ( float )double_val;
  }
              
  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }
      
  #if defined (_DEBUG_MSG_)
    printf ("Origin:  ( %f , %f , %f ) \n\n", (float)x_orig, (float)y_orig, (float)z_orig );
  #endif

  /// Read Data Pitch (Time Step)
  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }

  if ( is_float ) 
  {
    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time_steps = (float)float_val;

    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time_step = (float)float_val;

    if ( fread( &float_val, sizeof(float), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time = (float)float_val;
  }
  else if ( is_double ) 
  {
    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time_steps = (float)double_val;

    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time_step = (float)double_val;
    
    if ( fread( &double_val, sizeof(double), 1, sph_file ) != SUCCESS ) 
    {
      printf("<<< ERROR >>> Cannot read data information \n" );
      exit(FAIL);
    }
    time = (float)double_val;
  }
      
  if ( fread( &int_val, sizeof(int), 1, sph_file ) != SUCCESS ) 
  {
    printf("<<< ERROR >>> Cannot read data information \n" );
    exit(FAIL);
  }

  #if defined (_DEBUG_MSG_)
    printf ("TIME: %f \n\n", (float)time_step );
  #endif

  /// Read Data
  if ( is_scalar ) 
  {
    if ( is_float ) 
    {
      data_size = i_max * j_max * k_max;
      float_data = (float *)malloc ( data_size  * sizeof(float) ); 
      if ( float_data == NULL) 
      { 
        printf("<<< ERROR >> Cannot allocate memory for input data \n" );
        exit (FAIL);
      }
      fread ( float_data, sizeof(float), data_size, sph_file );
    }
    else if ( is_double ) 
    {
      data_size = i_max * j_max * k_max;
      double_data = (double *)malloc ( data_size  * sizeof(double) ); 
      if ( double_data == NULL) 
      { 
        printf("<<< ERROR >> Cannot allocate memory for input data \n" );
        exit (FAIL);
      }
      fread ( double_data, sizeof(double), data_size, sph_file );
    }
  }
  else if ( is_vector )
  {            
    if ( is_float ) 
    {
      data_size = i_max * j_max * k_max * 3;
      float_data = (float *)malloc ( data_size  * sizeof(float) ); 
      if ( float_data == NULL) 
      { 
        printf("<<< ERROR >> Cannot allocate memory for input data \n" );
        exit (FAIL);
      }
      fread ( float_data, sizeof(float), data_size, sph_file );
    }
    else if ( is_double ) 
    {
      data_size = i_max * j_max * k_max * 3 ;
      double_data = (double *)malloc ( data_size * sizeof(double) ); 
      if ( double_data == NULL) 
      { 
        printf("<<< ERROR >> Cannot allocate memory for input data \n" );
        exit (FAIL);
      }
      fread ( double_data, sizeof(double), data_size, sph_file );
    }
  }

  /// SPH Data Reading END
      
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      //================================================
  convert_start = MPI_Wtime();

  /*Rank 0 parses the command Line Arguments*/
  if (rank == 0) 
  {
    printf("Parsing Arguments\n");
    ret = parse_args(argc, argv);
    if (ret < 0) 
    {
      usage();
      print_error("ret error", __FILE__, __LINE__);
    }
    if( (extents[0]/count_local[0]) * (extents[1]/count_local[1]) * (extents[2]/count_local[2]) != nprocs)
    {
      usage();
      print_error("Wrong number of cores", __FILE__, __LINE__);
      
    }
    printf("Starting to load files\n");
  }

  #if defined (_DEBUG_MSG_)
  if ( rank == 0 ) 
  {
    printf("\n");
    printf("Extents      : %d %d %d %d %d \n", extents[0], extents[1], extents[2], extents[3], extents[4]);
    printf("Count_Local  : %d %d %d %d %d \n", count_local[0], count_local[1], count_local[2], count_local[3], count_local[4]);
    printf("\n");
  }
  #endif

  MPI_Bcast(extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&idx_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

  nx = extents[0] / count_local[0];
  ny = extents[1] / count_local[1];
  nz = extents[2] / count_local[2];

  assert(nx * ny * nz == nprocs);
  extents[3] = 1;
  extents[4] = 1;

  sub_div[0] = (extents[0] / count_local[0]);
  sub_div[1] = (extents[1] / count_local[1]);
  sub_div[2] = (extents[2] / count_local[2]);

  slice = rank % (sub_div[0] * sub_div[1]);

  offset_local[0] = (slice % sub_div[0]) * count_local[0];
  offset_local[1] = (slice / sub_div[0]) * count_local[1];
  offset_local[2] = (rank / (sub_div[0] * sub_div[1])) * count_local[2];
  offset_local[3] = 0;
  offset_local[4] = 0;

  count_local[3] = 1;
  count_local[4] = 1;

  #if defined (_DEBUG_MSG_)
    printf("[%d] Offset_Local : %d %d %d %d %d \n", rank, offset_local[0], offset_local[1], offset_local[2], offset_local[3], offset_local[4]);
    printf("[%d] Sub Div      : %d %d %d %d %d \n", rank, sub_div[0], sub_div[1], sub_div[2], sub_div[3], sub_div[4]);
  #endif

  idx_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(idx_file_name, "%s%s", idx_file_template, ".idx");

  alldata = (float*) malloc(sizeof (float) * count_local[0] * count_local[1] * count_local[2] );
  memset( alldata, 0, (sizeof (float) * count_local[0] * count_local[1] * count_local[2]) );
  
  #if defined (_SLICE_)
    printf("SLICE \n");
    slice_dimension = extents[2] / nprocs;
    if (rank == 0)
            printf("Slice Dimension : %d \n", slice_dimension);

    shift = rank * slice_dimension * extents[0] * extents[1];

    if ( (is_scalar) && ( is_float )) 
    {
      memcpy ( alldata, float_data + shift , sizeof(float) * count_local[0] * count_local[1] * count_local[2] ) ;
    }
    else if ( (is_scalar) && ( is_double )) 
    {
      alldata_ptr = alldata;
      double_ptr = double_data + shift;
      for ( k = 0; k < count_local[2]; k++) 
      {
        for ( j = 0; j < count_local[1]; j++) 
        {
          for ( i = 0; i < count_local[0]; j++) 
          {
            *alldata_ptr = (float)*double_ptr;
            alldata_ptr++;;
            double_ptr++;;
          }
        }
      }
    }
  #else 
    printf("BRICKS \n");

    if ( (is_scalar) && ( is_float )) 
    {
      shift_slice = extents[0] * extents[1];
      shift = (offset_local[2] * shift_slice) + (offset_local[1] * extents[0]) + offset_local[0];
      
      alldata_ptr = alldata;
      float_ptr = float_data + shift;

      for ( k = 0; k < count_local[2]; k++) 
      {
        for ( j = 0; j < count_local[1]; j++) 
        {
          memcpy ( alldata_ptr, float_ptr, sizeof(float) * count_local[0]);
          alldata_ptr += count_local[0];
          float_ptr += extents[0];
        }
        float_ptr = float_data + (shift_slice * (k + 1)) + shift;
      }
    }
    else if ( (is_scalar) && ( is_double )) 
    {
      shift_slice = extents[0] * extents[1];
      shift = (offset_local[2] * shift_slice) + (offset_local[1] * extents[0]) + offset_local[0];
      
      alldata_ptr = alldata;
      double_ptr = double_data + shift;

      for ( k = 0; k < count_local[2]; k++) 
      {
        for ( j = 0; j < count_local[1]; j++) 
        {
          for ( i = 0; i < count_local[0]; j++) 
          {
            *alldata_ptr = (float)*double_ptr;
            alldata_ptr++;;
            double_ptr++;;
          }
        }
        double_ptr = double_data + (shift_slice * (k + 1)) + shift;
      }
    }
  #endif

  offset = offset_local;
  count = count_local;
  gextent = extents;

  loading_time = MPI_Wtime();
  //MPI_Barrier(MPI_COMM_WORLD);

  if(rank == 0)
      printf("Done Loading All Files !!!!!!!!!\n");

      //Feeding the loaded chunks to PIDX.
  variable = malloc(sizeof (*variable) * 1);
  if (!variable)
      print_error("Error Allocating Variable pointer", __FILE__, __LINE__);
  memset(variable, 0, sizeof (*variable) * 1);

  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  PIDX_set_point_5D(gextent[0], gextent[1], gextent[2], 1, 1, global_bounding_box);
  PIDX_set_point_5D(offset[0], offset[1], offset[2], 0, 0, local_offset_point);
  PIDX_set_point_5D(count[0], count[1], count[2], 1, 1, local_box_count_point);
  
  PIDX_access access;
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#else
  PIDX_set_default_access(access);
#endif
    
  PIDX_file_create(idx_file_name, PIDX_file_trunc, access, &file);
  PIDX_set_dims(file, global_bounding_box);
  PIDX_set_current_time_step(file, 0);
  PIDX_set_block_size(file, bits_per_block);
  PIDX_set_block_count(file, blocks_per_file);
  PIDX_set_variable_count(file, 1);
  
  char data_type[512] = "1*float32";
  
  sprintf(var_name, "var_%d", 0);
  PIDX_variable_create(file, var_name, sizeof(float) * 8, data_type, &variable[0]);
  PIDX_append_and_write_variable(variable[0], local_offset_point, local_box_count_point, alldata, PIDX_row_major);
  
  PIDX_close(file);
  PIDX_close_access(access);
  
  if(rank == 0)
    printf("Done Writing PIDX !!!!!!!!!\n");

  free(alldata);
  alldata = 0;

  if(rank == 0)
      printf("Done with Cleanups !!!!!!!!!\n");
  
  convert_end = MPI_Wtime();
  total_time = convert_end - convert_start;
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(total_time == max_time)
    printf("Total Conversion Time [Loading Data + IDX Creation] %f = %f + %f\n", total_time, (loading_time - convert_start), (convert_end - loading_time));

  MPI_Finalize();

  return 0;
}

static int parse_args(int argc, char **argv) {
    char flags[] = "g:l:i:o:";
    int one_opt = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        /* postpone error checking for after while loop */
        switch (one_opt) {
            case('g'):
                sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
                break;
            case('l'):
                sscanf(optarg, "%dx%dx%d", &count_local[0], &count_local[1], &count_local[2]);
                break;
            case('i'):
                sprintf(sph_filename, "%s", optarg);
                break;
            case('o'):
                sprintf(idx_file_template, "%s", optarg);
                break;
            case('?'):
                return (-1);
        }
        if(extents[0] < 0 || extents[1] < 0 || extents[2] < 0)
        {
          usage();
          print_error("Wrong Aruments", __FILE__, __LINE__);
        }
        if(count_local[0] < 0 || count_local[1] < 0 || count_local[2] < 0)
        {
          usage();
          print_error("Wrong Arguments", __FILE__, __LINE__);
        }
        if(count_local[0] > extents[0] || count_local[1] > extents[1] || count_local[2] > extents[2])
        {
          usage();
          print_error("Wrong Arguments", __FILE__, __LINE__);
        }
    }
    /* need positive dimensions */

    return (0);
}

/* prints usage instructions */
static void usage(void) 
{
    printf("Usage: hydro_to_idx -g 2048x2048x2048 -l 256x256x256 -i IN_FILENAME -o OUT_FILENAME \n");
    printf("  -g: global dimensions (ex. 2048x2048x2048)\n");
    printf("  -l: local (per-process) dimensions (For 512 cores : 256x256x256)\n");
    printf("  -i: Path to SPH file \n");
    printf("  -o: Path to IDX file (Where you want to create)\n");
    printf("\n");
    return;
}

static int print_error(char *error_message, char* file, int line) {
    fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);

    MPI_Abort(MPI_COMM_WORLD, -1);
        
        return (-1);
}
