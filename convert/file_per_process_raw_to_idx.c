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

#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

#include <PIDX.h>


static int parse_args(int argc, char **argv);
static void usage(void);

/* global dimensions of 3D volume */
static int extents[5] = {0, 0, 0, 0, 0};

/* per-process dimensions of each sub-block of the existing fortran blocks */
static int fortran_count[5] = {0, 0, 0, 0, 0};

/* per-process dimensions of each sub-block of the PIDX blocks */
static int pidx_count[5] = {0, 0, 0, 0, 0};

/* Number of time-steps */
static int time_step = 0;

/*Number of Fortran Cores/data-blocks*/
static int fortran_core_count;


/* Fortran file Name Template*/
static char fortran_file_template[512] = {0};

/* PIDX file Name Template*/
static char pidx_file_template[512] = {0};

static char output_file_template[512] = {0};
static char **output_file_name;

struct NDimension_chunk_bound 
{
  int lower_bound[PIDX_MAX_DIMENSIONS];
  int upper_bound[PIDX_MAX_DIMENSIONS];
};
typedef struct NDimension_chunk_bound ND_chunk;

int intersect_ND_Chunk(ND_chunk* A, ND_chunk* B) 
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
  {
    check_bit = check_bit || A->upper_bound[d] < B->lower_bound[d] || B->upper_bound[d] < A->lower_bound[d];
  }
  return !(check_bit);
}


int main(int argc, char **argv) 
{
  int fp;
  char filename[512];
  int i = 0, j = 0, k = 0, t = 0, var = 0, spv = 0;
  int pidx_core_count, rank; /*process count and rank*/
  int slice, ret;
  int sub_div[5];
  int write_raw = 0;
  int raw_to_idx = 0; 
  
  //Fortran Data related variables
  char **fortran_file_name;
  int **fortran_offset;
  double ***raw_read_buffer;
  double **raw_write_buffer;
  
  //PIDX Data related variables  
  char *pidx_file_name;
  double **pidx_buffer;
  int pidx_offset[5];
  
  //Parameters common to both PIDX and Fortran
  int *fortran_pidx_index;
  int **fortran_pidx_offset;
  
  //Variables/Fields and Number of samples per field
  int number_of_variables;
  int *sample_per_variable;
  char var_name[512];
  char data_type[512];
  
  //PIDX API related variables
  PIDX_file file;                                                // IDX file descriptor
  const int bits_per_block = 15;                                                // Total number of samples in each block = 2 ^ bits_per_block
  const int blocks_per_file = 256;                                               // Total number of blocks per file
  const int *gextent; /*Global Extensions of the dataset (64 64 64 0 0)*/  
  PIDX_variable* variable;                                       // variable descriptor


  //MPI Initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &pidx_core_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  //TODO: Populate this Based on the Fortran input format
  number_of_variables = 1;
  sample_per_variable = (int*) malloc(sizeof (int) * number_of_variables);
  if (!sample_per_variable)
  {}

  sample_per_variable[0] = 1; //Variable 1
  //sample_per_variable[1] = 1; //Variable 2
  //sample_per_variable[2] = 1; //Variable 3
  //sample_per_variable[3] = 1; //Variable 4
    
  //Rank 0 parses the command Line Arguments
  if (rank == 0) 
  {
    ret = parse_args(argc, argv);
    if (ret < 0) 
    {
      usage();
      {}
    }
    if (pidx_count[0] == 0 || pidx_count[1] == 0 || pidx_count[2] == 0 || fortran_count[0] == 0 || fortran_count[1] == 0 || fortran_count[2] == 0) 
    {
      usage();
      {}
    } 
    else 
    {
      if ((extents[0] / pidx_count[0]) * (extents[1] / pidx_count[1]) * (extents[2] / pidx_count[2]) != pidx_core_count) 
      {
        usage();
        {}
      }
    }
  }
  
  // Broadcasting the command line arguments to all processes.
  MPI_Bcast(extents, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(pidx_count, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fortran_count, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&fortran_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);


  // Number of Fortran Cores
  fortran_core_count = (extents[0]/fortran_count[0]) * (extents[1]/fortran_count[1]) * (extents[2]/fortran_count[2]);
  
  // Creating Fortran Folder Names
  //Based on the code fortran-rst.c
  fortran_file_name = (char**) malloc(sizeof (char*) * time_step);
  for (i = 0; i < time_step; i++) 
  {
    fortran_file_name[i] = (char*) malloc(sizeof (char) * 512);
    sprintf(fortran_file_name[i], "%s%04d", fortran_file_template, i);
  }
  
  // Creating Fortran Folder Names
  pidx_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(pidx_file_name, "%s%s", pidx_file_template, ".idx");

    
  /// Figuring out the Fortran Blocks that need to be combined by every PIDX Process.
  /// Calculating offset and count of PIDX Block for process with Rank = rank
  gextent = extents;
  sub_div[0] = (extents[0] / pidx_count[0]);
  sub_div[1] = (extents[1] / pidx_count[1]);
  sub_div[2] = (extents[2] / pidx_count[2]);
  pidx_offset[2] = (rank / (sub_div[0] * sub_div[1])) * pidx_count[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  pidx_offset[1] = (slice / sub_div[0]) * pidx_count[1];
  pidx_offset[0] = (slice % sub_div[0]) * pidx_count[0];

  pidx_offset[3] = 0;
  pidx_offset[4] = 0;
  pidx_count[3] = 1;
  pidx_count[4] = 1;

  assert(pidx_offset[0] < extents[0]);
  assert(pidx_offset[1] < extents[1]);
  assert(pidx_offset[2] < extents[2]);
  assert(pidx_offset[3] < extents[3]);
  assert(pidx_offset[4] < extents[4]);
 
  
  if (write_raw == 1)
  {
    MPI_File fh;
    char filename[512];
    MPI_Status status;
    
    output_file_name = (char**) malloc(sizeof (char*) * time_step);
    for (i = 0; i < time_step; i++) 
    {
      output_file_name[i] = (char*) malloc(sizeof (char) * 512);
      sprintf(output_file_name[i], "%s%04d", output_file_template, i);
    }
    
    double *file_create_start, *file_create_end, *file_open_start, *file_open_end, **var_start, **var_end, *file_close_start, *file_close_end;
    var_start = (double**) malloc(sizeof (double*) * time_step);
    var_end = (double**) malloc(sizeof (double*) * time_step);
    memset(var_start, 0, (sizeof (double*) * time_step));
    memset(var_end, 0, (sizeof (double*) * time_step));
    
    for (i = 0; i < time_step; i++) 
    {
      var_start[i] = (double*) malloc(sizeof (double) * number_of_variables);
      var_end[i] = (double*) malloc(sizeof (double) * number_of_variables);
      memset(var_start[i], 0, (sizeof (double) * number_of_variables));
      memset(var_end[i], 0, (sizeof (double) * number_of_variables));
    }
    
    for (t = 0; t < time_step; t++) 
    {
      //WRITE BUFFER
      raw_write_buffer = (double**) malloc(sizeof (double*) * number_of_variables);
      if (!raw_write_buffer)
      {}
      memset(raw_write_buffer, 0, sizeof (double*) * number_of_variables);

      for (var = 0; var < number_of_variables; var++) 
      {
        raw_write_buffer[var] = (double*) malloc(sizeof (double) * fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[var]);
        if (!raw_write_buffer[var])
        {}
        for (k = 0; k < fortran_count[2]; k++)
          for (j = 0; j < fortran_count[1]; j++)
            for (i = 0; i < fortran_count[0]; i++) 
            {
              int64_t index = (int64_t) (fortran_count[0] * fortran_count[1] * k) + (fortran_count[0] * j) + i;
              for (spv = 0; spv < sample_per_variable[var]; spv++)
                raw_write_buffer[var][index * sample_per_variable[var] + spv] = 100 + var + spv + (extents[0] * extents[1]*(pidx_offset[2] + k))+(extents[0]*(pidx_offset[1] + j)) + (pidx_offset[0] + i);
            }
      }
        
      file_create_start[t] = MPI_Wtime();
      if (rank == 0) 
        mkdir(output_file_name[t], 0770);    
      
      MPI_Barrier(MPI_COMM_WORLD);
      file_create_end[t] = MPI_Wtime();

      file_open_start[t] = MPI_Wtime();
      sprintf(filename, "%s%srank%04d", output_file_name[t], "/", rank);
      MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      file_open_end[t] = MPI_Wtime();

      off_t data_offset = 0;
      for (var = 0; var < number_of_variables; var++) 
      {
        var_start[t][var] = MPI_Wtime();
        MPI_File_write_at(fh, data_offset, raw_write_buffer[var], (fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[var]), MPI_DOUBLE, &status);
        data_offset = data_offset + (fortran_count[0] * fortran_count[1] * fortran_count[2]) * sample_per_variable[var] * sizeof (double);
        var_end[t][var] = MPI_Wtime();
      }
      file_close_start[t] = MPI_Wtime();
      MPI_File_close(&fh);
      file_close_end[t] = MPI_Wtime();

      for (var = 0; var < number_of_variables; var++)
      {
        free(raw_write_buffer[var]);
        raw_write_buffer[var] = 0;
      }
      free(raw_write_buffer);
      raw_write_buffer = 0;
    }
  }
  
  
  
  
  if (raw_to_idx == 1)
  {
    //Creating buffer to hold pidx data
    pidx_buffer = (double**)malloc(number_of_variables * sizeof(double*));
    for(i = 0 ; i < number_of_variables ; i++) 
    {
      pidx_buffer[i] = (double*)malloc(pidx_count[0] * pidx_count[1] * pidx_count[2] * sample_per_variable[i] * sizeof(double));
    }
        
    // Calculating offset and count of Every FORTRAN Block
    fortran_offset = (int**)malloc(fortran_core_count * sizeof(int*));
    memset(fortran_offset, 0, fortran_core_count * sizeof(int*));
    for(i = 0 ; i < fortran_core_count ; i++)
    {
      fortran_offset[i] = (int*)malloc(5 * sizeof(int));
      memset(fortran_offset[i], 0, 5 * sizeof(int));
    }
    
    for(i = 0 ; i < fortran_core_count ; i++)    
    {
      sub_div[0] = (extents[0] / fortran_count[0]);
      sub_div[1] = (extents[1] / fortran_count[1]);
      sub_div[2] = (extents[2] / fortran_count[2]);
      fortran_offset[i][2] = (i / (sub_div[0] * sub_div[1])) * fortran_count[2];
      slice = i % (sub_div[0] * sub_div[1]);
      fortran_offset[i][1] = (slice / sub_div[0]) * fortran_count[1];
      fortran_offset[i][0] = (slice % sub_div[0]) * fortran_count[0];

      fortran_offset[i][3] = 0;
      fortran_offset[i][4] = 0;
      fortran_count[3] = 1;
      fortran_count[4] = 1;
      
      assert(fortran_offset[i][0] < extents[0]);
      assert(fortran_offset[i][1] < extents[1]);
      assert(fortran_offset[i][2] < extents[2]);
      assert(fortran_offset[i][3] < extents[3]);
      assert(fortran_offset[i][4] < extents[4]);
      //printf("[%d] Fortran Offset: %d %d %d %d %d :: %d %d %d %d %d\n", i, fortran_offset[i][0], fortran_offset[i][1], fortran_offset[i][2], fortran_offset[i][3], fortran_offset[i][4], fortran_count[0], fortran_count[1], fortran_count[2], fortran_count[3], fortran_count[4]);
    }
    
    //Creating buffer to hold fortran data
    raw_read_buffer = (double***)malloc((fortran_core_count/pidx_core_count) * sizeof(double**));
    memset(raw_read_buffer, 0, (fortran_core_count/pidx_core_count) * sizeof(double**));
    for(i = 0 ; i < (fortran_core_count/pidx_core_count) ; i++)    
    {
      raw_read_buffer[i] = (double**)malloc(number_of_variables * sizeof(double*));
      memset(raw_read_buffer[i], 0, number_of_variables * sizeof(double*));
      for(j = 0 ; j < (number_of_variables) ; j++)    
      {
        raw_read_buffer[i][j] = (double*)malloc(fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[j] * sizeof(double));
        memset(raw_read_buffer[i][j], 0, fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[j] * sizeof(double));
      }
    }
      
    //Creating memory buffers to hold block information common to PIDX and Fortran
    fortran_pidx_offset = (int**)malloc((fortran_core_count/pidx_core_count) * sizeof(int*));
    memset(fortran_pidx_offset, 0, (fortran_core_count/pidx_core_count) * sizeof(int*));
    for(i = 0 ; i < (fortran_core_count/pidx_core_count) ; i++)    
    {
      fortran_pidx_offset[i] = (int*)malloc(5 * sizeof(int));
      memset(fortran_pidx_offset[i], 0, 5 * sizeof(int));
    }
    fortran_pidx_index = (int*)malloc((fortran_core_count/pidx_core_count) * sizeof(int));
    memset(fortran_pidx_index, 0, (fortran_core_count/pidx_core_count) * sizeof(int));
    
      
    ND_chunk* PIDX_Chunk = (ND_chunk*) malloc(sizeof (ND_chunk));
    memset(PIDX_Chunk, 0, sizeof (*PIDX_Chunk));
    PIDX_Chunk->lower_bound[0] = pidx_offset[0];
    PIDX_Chunk->lower_bound[1] = pidx_offset[1];
    PIDX_Chunk->lower_bound[2] = pidx_offset[2];
    PIDX_Chunk->lower_bound[3] = 0;
    PIDX_Chunk->lower_bound[4] = 0;
    PIDX_Chunk->upper_bound[0] = pidx_offset[0] + pidx_count[0] - 1;
    PIDX_Chunk->upper_bound[1] = pidx_offset[1] + pidx_count[1] - 1;
    PIDX_Chunk->upper_bound[2] = pidx_offset[2] + pidx_count[2] - 1;
    PIDX_Chunk->upper_bound[3] = 0;
    PIDX_Chunk->upper_bound[4] = 0;
      
    int counter = 0;
    for(i = 0 ; i < fortran_core_count ; i++)
    {
      ND_chunk* FORTRAN_Chunk = (ND_chunk*) malloc(sizeof (ND_chunk));
      if (!FORTRAN_Chunk) {}
      memset(FORTRAN_Chunk, 0, sizeof (*FORTRAN_Chunk));
      FORTRAN_Chunk->lower_bound[0] = fortran_offset[i][0];
      FORTRAN_Chunk->lower_bound[1] = fortran_offset[i][1];
      FORTRAN_Chunk->lower_bound[2] = fortran_offset[i][2];
      FORTRAN_Chunk->lower_bound[3] = 0;
      FORTRAN_Chunk->lower_bound[4] = 0;
      FORTRAN_Chunk->upper_bound[0] = fortran_offset[i][0] + fortran_count[0] - 1;
      FORTRAN_Chunk->upper_bound[1] = fortran_offset[i][1] + fortran_count[1] - 1;
      FORTRAN_Chunk->upper_bound[2] = fortran_offset[i][2] + fortran_count[2] - 1;
      FORTRAN_Chunk->upper_bound[3] = 0;
      FORTRAN_Chunk->upper_bound[4] = 0;
      if (intersect_ND_Chunk(FORTRAN_Chunk, PIDX_Chunk)) 
      {
        fortran_pidx_offset[counter][0] = FORTRAN_Chunk->lower_bound[0];
        fortran_pidx_offset[counter][1] = FORTRAN_Chunk->lower_bound[1];
        fortran_pidx_offset[counter][2] = FORTRAN_Chunk->lower_bound[2];
        fortran_pidx_offset[counter][3] = FORTRAN_Chunk->lower_bound[3];
        fortran_pidx_offset[counter][4] = FORTRAN_Chunk->lower_bound[4];
        fortran_pidx_index[counter] = i;
        counter++;
        
        if(counter == fortran_core_count/pidx_core_count)
          break;
      }
      
      free(FORTRAN_Chunk);
      FORTRAN_Chunk = 0;
    }
    
    assert(counter == fortran_core_count/pidx_core_count);
    free(PIDX_Chunk);
    PIDX_Chunk = 0;
      
    
    variable = malloc(sizeof(*variable) * number_of_variables);
    memset(variable, 0, sizeof(*variable) * number_of_variables);
    
    PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
    
    PIDX_set_point_5D((int64_t)gextent[0], (int64_t)gextent[1], (int64_t)gextent[2], 1, 1, global_bounding_box);
    PIDX_set_point_5D((int64_t)pidx_offset[0], (int64_t)pidx_offset[1], (int64_t)pidx_offset[2], 0, 0, local_offset_point);
    PIDX_set_point_5D((int64_t)pidx_count[0], (int64_t)pidx_count[1], (int64_t)pidx_count[2], 1, 1, local_box_count_point);
    
    for(t = 0 ; t < time_step ; t++)
    {
      for(i = 0 ; i < (fortran_core_count/pidx_core_count) ; i++)
      {
        off_t data_offset = 0;
        sprintf(filename, "%s%srank%04d", fortran_file_name[0], "/", fortran_pidx_index[i]);
        fp = open(filename, O_RDONLY);
        
        for(var = 0 ; var < number_of_variables ; var++)
        {
          ret = pread(fp, raw_read_buffer[i][var], (fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[var] * sizeof(double)), data_offset);
          assert(ret == (fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[var] * sizeof(double)));
                      
          data_offset = data_offset + (fortran_count[0] * fortran_count[1] * fortran_count[2] * sample_per_variable[var]);
        }
        data_offset = 0;
        close(fp);
      }
      
      int total_send = 0;
      for(var = 0 ; var < number_of_variables ; var++)
      {
        int k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
        for (r = 0; r < (fortran_core_count/pidx_core_count); r++) 
        {
          for (k1 = fortran_pidx_offset[r][2]; k1 < fortran_pidx_offset[r][2] + fortran_count[2]; k1++)   
          {
            for (j1 = fortran_pidx_offset[r][1]; j1 < fortran_pidx_offset[r][1] + fortran_count[1]; j1++)
            {
              for (i1 = fortran_pidx_offset[r][0]; i1 < fortran_pidx_offset[r][0] + fortran_count[0]; i1 = i1 + fortran_count[0]) 
              {
                index = ((fortran_count[0])* (fortran_count[1]) * (k1 - fortran_pidx_offset[r][2])) + ((fortran_count[0]) * (j1 - fortran_pidx_offset[r][1])) + (i1 - fortran_pidx_offset[r][0]);
                send_o = index * sample_per_variable[var];
                send_c = (fortran_count[0]) * sample_per_variable[var];
                recv_o = ((pidx_count[0]) * (pidx_count[1]) * (k1 - pidx_offset[2])) + ((pidx_count[0])* (j1 - pidx_offset[1])) + (i1 - pidx_offset[0]);
                memcpy(pidx_buffer[var] + (recv_o * sample_per_variable[var]), raw_read_buffer[r][var] + send_o, send_c * sizeof (double));
                total_send = total_send + send_c;
              }
            }
          }
        }
      }
      
      PIDX_access access;
      PIDX_create_access(&access);

  #if PIDX_HAVE_MPI
      PIDX_set_mpi_access(access, MPI_COMM_WORLD);
  #else
      PIDX_set_default_access(access);
  #endif
      
      PIDX_file_create(pidx_file_name, PIDX_file_trunc, access, &file);
      PIDX_set_dims(file, global_bounding_box);
      PIDX_set_current_time_step(file, t);
      PIDX_set_block_size(file, bits_per_block);
      PIDX_set_block_count(file, blocks_per_file);
      PIDX_set_variable_count(file, number_of_variables);
        
      for (var = 0; var < number_of_variables; var++) 
      {
        sprintf(var_name, "var_%d", var);
        sprintf(data_type, "%d*float64", sample_per_variable[var]);
        PIDX_variable_create(file, var_name, sample_per_variable[var] * sizeof(double) * 8, data_type, &variable[var]);
        PIDX_append_and_write_variable(variable[var], local_offset_point, local_box_count_point, pidx_buffer[var], PIDX_row_major);
      }
      PIDX_close(file);
      PIDX_close_access(access);
    }
      
    free(fortran_pidx_index);
    fortran_pidx_index = 0;
      
    for(i = 0 ; i < fortran_core_count ; i++)
    {
      free(fortran_offset[i]);
      fortran_offset[i] = 0;
    }
    free(fortran_offset);
    fortran_offset = 0;
        
    for (var = 0; var < number_of_variables; var++) 
    {
        free(pidx_buffer[var]);
        pidx_buffer[var] = 0;
    }
    free(pidx_buffer);
    pidx_buffer = 0;
  }
  
  for(i = 0 ; i < (fortran_core_count/pidx_core_count) ; i++)
  {
    for(j = 0 ; j < number_of_variables ; j++)
    {
      free(raw_read_buffer[i][j]);
      raw_read_buffer[i][j] = 0;
    }
    free(raw_read_buffer[i]);
    raw_read_buffer[i] = 0;
    free(fortran_pidx_offset[i]);
    fortran_pidx_offset[i] = 0;
  }
  free(raw_read_buffer);
  raw_read_buffer = 0;
  free(fortran_pidx_offset);
  fortran_pidx_offset = 0;
  
  free(variable);
  variable = 0;
  free(sample_per_variable);
  sample_per_variable = 0;
  
  MPI_Finalize();
    
  return 0;
}

/////////////////////////////////////////////////
static int parse_args(int argc, char **argv) 
{
  char flags[] = "g:f:p:i:o:t:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF) 
  {
    /* postpone error checking for after while loop */
    switch (one_opt) 
    {
      case('g'):
          sscanf(optarg, "%dx%dx%d", &extents[0], &extents[1], &extents[2]);
          break;
      case('f'):
          sscanf(optarg, "%dx%dx%d", &fortran_count[0], &fortran_count[1], &fortran_count[2]);
          break;
      case('p'):
          sscanf(optarg, "%dx%dx%d", &pidx_count[0], &pidx_count[1], &pidx_count[2]);
          break;
      case('i'):
          sprintf(fortran_file_template, "%s", optarg);
          break;
      case('o'):
          sprintf(pidx_file_template, "%s", optarg);
          break;
      case('t'):
          sscanf(optarg, "%d", &time_step);
          break;
      case('?'):
          return (-1);
    }
  }
  extents[3] = 1;
  extents[4] = 1;
  /* need positive dimensions */
  if (extents[0] < 1 || extents[1] < 1 || extents[2] < 1 || fortran_count[0] < 1 || fortran_count[1] < 1 || fortran_count[2] < 1) 
  {
    printf("Error: bad dimension specification.\n");
    return (-1);
  }

  /* need global dimension to be larger than the local */
  if (extents[0] < fortran_count[0] || extents[1] < fortran_count[1] || extents[2] < fortran_count[2]) 
  {
    printf("Error: global dimensions and local dimensions aren't evenly divisible\n");
    return (-1);
  }
  return (0);
}

/* prints usage instructions */
static void usage(void) 
{
  printf("Usage: fortran-to-idx -g 4x4x4 -l 2x2x2 -n 8 -f Filename_ -t 4\n");
  printf("  -g: global dimensions of dataset\n");
  printf("  -f: local (per-process) dimensions of dataset\n");
  printf("  -n: number of fortran files\n");
  printf("  -f: fortran file name template\n");
  printf("  -t: number of timesteps\n");
  printf("\n");
  return;
}
