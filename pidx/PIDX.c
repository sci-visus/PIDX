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

/**
 * \file PIDX.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX.h
 *
 */

// TODO: PIDX reader
// (global box)
// 124 128 132
// (local box)
// 31 32 33

#include "PIDX.h"
#define PIDX_DEBUG_OUTPUT 0

static int vp = 0;
static int hp = 0;
static double file_create_time = 0;
static double populate_idx_start_time = 0;
static double populate_idx_end_time = 0;
static double sim_start = 0, sim_end = 0;
static double *write_init_start = 0, *write_init_end = 0;
static double *startup_start = 0, *startup_end = 0;
static double *init_start, *init_end;
static double *rst_start, *rst_end;
static double *hz_start, *hz_end;
static double **agg_start, **agg_end;
static double **io_per_process_start, **io_per_process_end;
static double **io_start, **io_end;
static double *cleanup_start, *cleanup_end;
static double *finalize_start, *finalize_end;
static double *chunk_start, *chunk_end;
static double *buffer_start, *buffer_end;
static double *compression_start, *compression_end;

static int caching_state = 0;
static int time_step_caching = 0;
static uint32_t* cached_header_copy;

static int static_var_counter = 0;

static void PIDX_init_timming_buffers1();
static void PIDX_init_timming_buffers2(PIDX_file file);
static PIDX_return_code PIDX_write(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_read(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_raw_write(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_raw_read(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* file_name, int current_time_step);
static PIDX_return_code PIDX_parameter_validate(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_validate(PIDX_file file);
static PIDX_return_code populate_idx_dataset(PIDX_file file);

/// PIDX File descriptor (equivalent to the descriptor returned by)
/// POSIX or any other IO framework
struct PIDX_file_descriptor
{
  int flags;

#if PIDX_HAVE_MPI
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
  MPI_Comm global_comm;                        ///< MPI super-communicator (includes all processes)
#endif

  PIDX_access access;                          ///< serial or parallel access

  int idx_count[PIDX_MAX_DIMENSIONS];          ///< Number of idx files in each dimensions

  PIDX_header_io_id header_io_id;              ///< IDX metadata id
  PIDX_rst_id rst_id;                          ///< Restructuring phase id
  PIDX_chunk_id chunk_id;              ///< Block restructuring id (prepration for compression)
  PIDX_comp_id comp_id;          ///< Compression (lossy and lossless) id
  PIDX_hz_encode_id hz_id;                     ///< HZ encoding phase id
  PIDX_agg_id agg_id;                          ///< Aggregation phase id
  PIDX_agg_id** tagg_id;                          ///< Aggregation phase id
  PIDX_io_id io_id;                            ///< IO phase id
  PIDX_io_id** tio_id;                            ///< IO phase id

  int local_variable_index;                    ///<
  int local_variable_count;                    ///<
  int var_pipe_length;                         ///<

  int flush_used;
  int write_on_close;                          ///< HPC Writes
  int one_time_initializations;                ///<

  int ROI_writes;

  int debug_rst;                               ///< Debug restructuring phase, works only on the test application
  int debug_chunk;                             ///< Debug chunking phase, works only on the test application
  int debug_compress;                          ///< Debug compression phase, works only on the test application
  int debug_hz;                                ///< Debug HZ encoding phase, works only on the test application
  int debug_agg;                               ///< Debug aggregation phase, works only on the test application

  /// Flags set by user
  int debug_do_rst;                            ///< User controlled flag to activate/deactivate restructuring phase
  int debug_do_chunk;                          ///< User controlled flag to activate/deactivate chunking phase
  int debug_do_compress;                       ///< User controlled flag to activate/deactivate compression
  int debug_do_hz;                             ///< User controlled flag to activate/deactivate hz encoding phase
  int debug_do_agg;                            ///< User controlled flag to activate/deactivate aggregation phase
  int debug_do_io;                             ///< User controlled flag to activate/deactivate I/O phase

  int small_agg_comm;

  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  int enable_raw_dump;

  int debug_output;

  int agg_type;

  int layout_count;
  int layout_start_index;
  int layout_end_index;

  int reduced_res_from;
  int reduced_res_to;
};


/// Returns elapsed time
double PIDX_get_time()
{
#if PIDX_HAVE_MPI
  //double time = MPI_Wtime();
  //printf("Time = %f\n", time);

  return MPI_Wtime();
#else
  struct timeval temp;
  gettimeofday(&temp, NULL);
  return (double)(temp.tv_sec) + (double)(temp.tv_usec)/1000000.0;
#endif
}



PIDX_return_code PIDX_time_step_caching_ON()
{
  caching_state = 1;
  time_step_caching = 1;

  return PIDX_success;
}



PIDX_return_code PIDX_time_step_caching_OFF()
{
  free(cached_header_copy);
  cached_header_copy = 0;

  return PIDX_success;
}



/// Function to create IDX file descriptor (based on flags and access)
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  PIDX_init_timming_buffers1();
  
  sim_start = PIDX_get_time();

  if (flags != PIDX_MODE_CREATE && flags != PIDX_MODE_EXCL)
    return PIDX_err_unsupported_flags;
    
  if (flags == PIDX_MODE_EXCL)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }
  
  uint64_t i = 0, j = 0, k = 0;
  int ret;
  char file_name_skeleton[1024];
  int rank = 0;
  
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
    
  //
  *file = malloc(sizeof (*(*file)));
  memset(*file, 0, sizeof (*(*file)));

  //
  (*file)->idx = malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));

  // IDX derived pointer (everything derived from idx)
  (*file)->idx_d = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_d)));
  memset((*file)->idx_d, 0, sizeof (*((*file)->idx_d)));
  
  (*file)->flags = flags;

  (*file)->access = access_type;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    (*file)->idx_count[i] = 1;

  (*file)->debug_do_rst = 1;
  (*file)->debug_do_chunk = 1;
  (*file)->debug_do_compress = 1;
  (*file)->debug_do_hz = 1;
  (*file)->debug_do_agg = 1;
  (*file)->debug_do_io = 1;

  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->var_pipe_length = -1;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;
  (*file)->agg_type = 0;

  (*file)->debug_rst = 0;
  (*file)->debug_hz = 0;
  (*file)->idx_d->color = 0;

  (*file)->small_agg_comm = 0;

  (*file)->enable_raw_dump = 0;

  (*file)->debug_output = 0;

  (*file)->reduced_res_from = 0;
  (*file)->reduced_res_to = 0;

  (*file)->idx_d->parallel_mode = access_type->parallel;

#if PIDX_HAVE_MPI
  //unsigned int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice = 0;
  int *colors;
  if (access_type->parallel)
  {
    MPI_Comm_rank(access_type->comm, &rank);
    if (access_type->idx_count[0] != 1 || access_type->idx_count[1] != 1 || access_type->idx_count[2] != 1 )
    {
      memcpy ((*file)->idx_count, access_type->idx_count, sizeof(int) * PIDX_MAX_DIMENSIONS);

      colors = malloc(sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      memset(colors, 0, sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);

      for (k = 0; k < (*file)->idx_count[2]; k++)
        for (j = 0; j < (*file)->idx_count[1]; j++)
          for (i = 0; i < (*file)->idx_count[0]; i++)
            colors[((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i] = ((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i;

      int index_x = 0, index_y = 0, index_z = 0;
      for (i = 0; i < access_type->sub_div[0]; i = i + (access_type->sub_div[0] / (*file)->idx_count[0]))
      {
        if (access_type->rank_component[0] >= i && access_type->rank_component[0] < i + (access_type->sub_div[0] / (*file)->idx_count[0]))
        {
          index_x = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[1]; i = i + (access_type->sub_div[1] / (*file)->idx_count[1]))
      {
        if (access_type->rank_component[1] >= i && access_type->rank_component[1] < i + (access_type->sub_div[1] / (*file)->idx_count[1]))
        {
          index_y = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[2]; i = i + (access_type->sub_div[2] / (*file)->idx_count[2]))
      {
        if (access_type->rank_component[2] >= i && access_type->rank_component[2] < i + (access_type->sub_div[2] / (*file)->idx_count[2]))
        {
          index_z = i;
          break;
        }
      }

      (*file)->idx_d->color = colors[((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]))];

      free(colors);
      MPI_Comm_split(access_type->comm, (*file)->idx_d->color, rank, &((*file)->comm));
      MPI_Comm_dup(access_type->comm, &((*file)->global_comm));
    }
    else
      MPI_Comm_dup(access_type->comm, &((*file)->comm));
  }
#endif

  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx->variable_index_tracker = 0;

  (*file)->idx->enable_rst = 1;
  (*file)->idx->enable_agg = 1;
  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';

  if ((*file)->idx_count[0] == 1 && (*file)->idx_count[1] == 1 && (*file)->idx_count[2] == 1)
    sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  else
    sprintf((*file)->idx->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_d->color);

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->bounds[i]=65535;

  //initialize logic_to_physic transform to identity
  (*file)->idx->transform[0]  = 1.0;
  (*file)->idx->transform[5]  = 1.0;
  (*file)->idx->transform[10] = 1.0;
  (*file)->idx->transform[15] = 1.0;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);
  memset((*file)->idx->reg_patch_size, 0, sizeof(int64_t) * PIDX_MAX_DIMENSIONS);

  (*file)->idx->compression_factor = 1;
  (*file)->idx->compression_bit_rate = 64;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx_d->dimension = 0;
  (*file)->idx_d->samples_per_block = pow(2, PIDX_default_bits_per_block);
  (*file)->idx_d->maxh = 0;
  (*file)->idx_d->max_file_count = 0;
  (*file)->idx_d->fs_block_size = 0;
  (*file)->idx_d->start_fs_block = 0;
  //(*file)->idx_d->agg_buffer->aggregation_factor = 1;
  (*file)->idx_d->dump_agg_info = 0;
  (*file)->idx_d->dump_io_info = 0;
  memset((*file)->idx_d->agg_dump_dir_name, 0, 512*sizeof(char));
  memset((*file)->idx_d->io_dump_dir_name, 0, 512*sizeof(char));
  (*file)->idx_d->staged_aggregation = 0;

  if (rank == 0)
  {
    //TODO: close and delete the file (there is a way to do this automatically by fopen...)
    struct stat stat_buf;
    FILE *dummy = fopen("dummy.txt", "w");
    fclose(dummy);
    ret = stat("dummy.txt", &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    (*file)->idx_d->fs_block_size = stat_buf.st_blksize;
  }

#if PIDX_HAVE_MPI
  if ((*file)->idx_d->parallel_mode == 1)
    MPI_Bcast(&((*file)->idx_d->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
    
  file_create_time = PIDX_get_time();

  return PIDX_success;
}


/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  PIDX_init_timming_buffers1();
  
  sim_start = PIDX_get_time();
  
  int i;
  //int ret;
  char file_name_skeleton[1024];
  int rank = 0;

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
  
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  (*file)->flags = flags;
  
  
  (*file)->idx = (idx_dataset)malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));
  
  (*file)->idx_d = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_d)));
  memset((*file)->idx_d, 0, sizeof (*((*file)->idx_d)));


  (*file)->access = access_type;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    (*file)->idx_count[i] = 1;

  (*file)->debug_do_rst = 1;
  (*file)->debug_do_chunk = 1;
  (*file)->debug_do_compress = 1;
  (*file)->debug_do_hz = 1;
  (*file)->debug_do_agg = 1;
  (*file)->debug_do_io = 1;

  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->var_pipe_length = -1;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;

  (*file)->reduced_res_from = 0;
  (*file)->reduced_res_to = 0;

  (*file)->agg_type = 1;

  (*file)->small_agg_comm = 0;
  (*file)->debug_output = 0;

  (*file)->debug_rst = 0;
  (*file)->debug_hz = 0;
  (*file)->idx_d->color = 0;

  (*file)->enable_raw_dump = 0;

  (*file)->idx_d->parallel_mode = access_type->parallel;


#if PIDX_HAVE_MPI
  //unsigned int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice = 0;
  int *colors;
  if (access_type->parallel)
  {
    MPI_Comm_rank(access_type->comm, &rank);
    if (access_type->idx_count[0] != 1 || access_type->idx_count[1] != 1 || access_type->idx_count[2] != 1 )
    {
      uint64_t i = 0, j = 0, k = 0;
      memcpy ((*file)->idx_count, access_type->idx_count, sizeof(int) * PIDX_MAX_DIMENSIONS);

      colors = malloc(sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      memset(colors, 0, sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);

      for (k = 0; k < (*file)->idx_count[2]; k++)
        for (j = 0; j < (*file)->idx_count[1]; j++)
          for (i = 0; i < (*file)->idx_count[0]; i++)
            colors[((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i] = ((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i;

      int index_x = 0, index_y = 0, index_z = 0;
      for (i = 0; i < access_type->sub_div[0]; i = i + (access_type->sub_div[0] / (*file)->idx_count[0]))
      {
        if (access_type->rank_component[0] >= i && access_type->rank_component[0] < i + (access_type->sub_div[0] / (*file)->idx_count[0]))
        {
          index_x = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[1]; i = i + (access_type->sub_div[1] / (*file)->idx_count[1]))
      {
        if (access_type->rank_component[1] >= i && access_type->rank_component[1] < i + (access_type->sub_div[1] / (*file)->idx_count[1]))
        {
          index_y = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[2]; i = i + (access_type->sub_div[2] / (*file)->idx_count[2]))
      {
        if (access_type->rank_component[2] >= i && access_type->rank_component[2] < i + (access_type->sub_div[2] / (*file)->idx_count[2]))
        {
          index_z = i;
          break;
        }
      }

      (*file)->idx_d->color = colors[((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]))];

      free(colors);
      MPI_Comm_split(access_type->comm, (*file)->idx_d->color, rank, &((*file)->comm));
      MPI_Comm_dup(access_type->comm, &((*file)->global_comm));
    }
    else
      MPI_Comm_dup(access_type->comm, &((*file)->comm));
  }
#endif


  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx->variable_index_tracker = 0;

  (*file)->idx->enable_rst = 1;
  (*file)->idx->enable_agg = 1;
  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';

  if ((*file)->idx_count[0] == 1 && (*file)->idx_count[1] == 1 && (*file)->idx_count[2] == 1)
    sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  else
    sprintf((*file)->idx->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_d->color);

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->bounds[i]=65535;

  //initialize logic_to_physic transform to identity
  (*file)->idx->transform[0]  = 1.0;
  (*file)->idx->transform[5]  = 1.0;
  (*file)->idx->transform[10] = 1.0;
  (*file)->idx->transform[15] = 1.0;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);
  memset((*file)->idx->reg_patch_size, 0, sizeof(int64_t) * PIDX_MAX_DIMENSIONS);

  (*file)->idx->compression_bit_rate = 64;
  (*file)->idx->compression_factor = 1;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx_d->dimension = 0;
  (*file)->idx_d->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);
  (*file)->idx_d->maxh = 0;
  (*file)->idx_d->max_file_count = 0;
  (*file)->idx_d->fs_block_size = 0;
  (*file)->idx_d->start_fs_block = 0;
  //(*file)->idx_d->agg_buffer->aggregation_factor = 1;
  (*file)->idx_d->dump_agg_info = 0;
  (*file)->idx_d->dump_io_info = 0;
  memset((*file)->idx_d->agg_dump_dir_name, 0, 512*sizeof(char));
  memset((*file)->idx_d->io_dump_dir_name, 0, 512*sizeof(char));

  int var = 0, variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];

  if (rank == 0)
  {
    FILE *fp = fopen((*file)->idx->filename, "r");
    if (fp == NULL)
    {
      fprintf(stdout, "Error Opening %s\n", (*file)->idx->filename);
      return PIDX_err_file;
    }

    while (fgets(line, sizeof (line), fp) != NULL) 
    {
      //printf("%s", line);
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "(box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL) 
        {
          if (count % 2 == 1)
            (*file)->idx->bounds[count / 2] = atoi(pch) + 1;
          count++;
          pch = strtok(NULL, " ");
        }
      }
      if (strcmp(line, "(fields)") == 0) 
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        count = 0;
        variable_counter = 0;
        
        while (line[0] != '(')
        //while (strcmp(line, "(version)") != 0 && strcmp(line, "(box)") != 0 && strcmp(line, "(bits)") && strcmp(line, "(bitsperblock)") != 0 && strcmp(line, "(blocksperfile)") != 0 && strcmp(line, "(filename_template)") != 0 && strcmp(line, "(time)") != 0)
        {
          (*file)->idx->variable[variable_counter] = malloc(sizeof (*((*file)->idx->variable[variable_counter])));
          memset((*file)->idx->variable[variable_counter], 0, sizeof (*((*file)->idx->variable[variable_counter])));
          
          pch1 = strtok(line, " +");
          while (pch1 != NULL)
          {
            if (count == 0)
              strcpy((*file)->idx->variable[variable_counter]->var_name, strdup(pch1));
              //strcpy((*file)->idx->variable[variable_counter]->var_name, pch1);

            //if (count == 1)
            //  (*file)->idx->variable[variable_counter]->values_per_sample = atoi(pch1);

              
            if (count == 1)
            {
              len = strlen(pch1) - 1;
              if (pch1[len] == '\n')
                pch1[len] = 0;
              //if (strcmp(pch1, "float64") == 0 && (*file)->idx->variable[variable_counter]->values_per_sample == 1)
              //{
                //strcpy((*file)->idx->variable[variable_counter]->type_name, "1*float64");
                strcpy((*file)->idx->variable[variable_counter]->type_name, pch1);
                int ret;
                int bits_per_sample = 0;
                ret = PIDX_default_bits_per_datatype((*file)->idx->variable[variable_counter]->type_name, &bits_per_sample);
                if (ret != PIDX_success)  return PIDX_err_file;

                (*file)->idx->variable[variable_counter]->bits_per_value = bits_per_sample;
                (*file)->idx->variable[variable_counter]->values_per_sample = 1;
              //}
              /*
              else if (strcmp(pch1, "float32") == 0 && (*file)->idx->variable[variable_counter]->values_per_sample == 1)
              {
                strcpy((*file)->idx->variable[variable_counter]->type_name, "1*float32");

                int ret;
                int bits_per_sample = 0;
                ret = PIDX_default_bits_per_datatype((*file)->idx->variable[variable_counter]->type_name, &bits_per_sample);
                if (ret != PIDX_success)  return PIDX_err_file;

                (*file)->idx->variable[variable_counter]->bits_per_value = bits_per_sample;
                (*file)->idx->variable[variable_counter]->values_per_sample = 1;
              }
              else if (strcmp(pch1, "float64") == 0 && (*file)->idx->variable[variable_counter]->values_per_sample == 3)
              {
                strcpy((*file)->idx->variable[variable_counter]->type_name, "3*float64");

                int ret;
                int bits_per_sample = 0;
                ret = PIDX_default_bits_per_datatype((*file)->idx->variable[variable_counter]->type_name, &bits_per_sample);
                if (ret != PIDX_success)  return PIDX_err_file;

                (*file)->idx->variable[variable_counter]->bits_per_value = bits_per_sample;
                (*file)->idx->variable[variable_counter]->values_per_sample = 1;
              }
              */
            }
            count++;
            pch1 = strtok(NULL, " +");
          }
          count = 0;

          if( fgets(line, sizeof line, fp) == NULL)
            return PIDX_err_file;
          line[strcspn(line, "\r\n")] = 0;
          variable_counter++;
        }
        (*file)->idx->variable_count = variable_counter;
      }

      if (strcmp(line, "(bits)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
      }

      if (strcmp(line, "(bitsperblock)") == 0) 
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        (*file)->idx->bits_per_block = atoi(line);
        (*file)->idx_d->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);

        //if (((*file)->idx->chunk_size[0] == 4) && ((*file)->idx->chunk_size[1] == 4) && ((*file)->idx->chunk_size[2] == 4) && ((*file)->idx->bounds[0] == 4) && ((*file)->idx->bounds[1] == 4) && ((*file)->idx->bounds[2] == 4) && (*file)->idx->bits_per_block == 1)
        //{

        //}
      }

      if (strcmp(line, "(compression type)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        (*file)->idx->compression_type = atoi(line);
      }

      if (strcmp(line, "(compressed box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          (*file)->idx->chunk_size[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }

        if((*file)->idx->chunk_size[0] < 0 || (*file)->idx->chunk_size[1] < 0 || (*file)->idx->chunk_size[2] < 0 || (*file)->idx->chunk_size[3] < 0 || (*file)->idx->chunk_size[4] < 0)
          return PIDX_err_box;

        /*
        int reduce_by_sample = 1;
        uint64_t total_chunk_size = (*file)->idx->chunk_size[0] * (*file)->idx->chunk_size[1] * (*file)->idx->chunk_size[2] * (*file)->idx->chunk_size[3] * (*file)->idx->chunk_size[4];
        if (reduce_by_sample == 1)
        {
          (*file)->idx->bits_per_block = (*file)->idx->bits_per_block - (int)log2(total_chunk_size);
          (*file)->idx_d->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);

          if ((*file)->idx->bits_per_block <= 0)
          {
            (*file)->idx->bits_per_block = 1;
            (*file)->idx_d->samples_per_block = 1;
          }
        }
        else
          (*file)->idx->blocks_per_file = (*file)->idx->blocks_per_file / total_chunk_size;
        */
      }

      if (strcmp(line, "(compression bit rate)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        (*file)->idx->compression_bit_rate = atoi(line);
      }

      if (strcmp(line, "(blocksperfile)") == 0) 
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        (*file)->idx->blocks_per_file= atoi(line);
      }

      if (strcmp(line, "(filename_template)") == 0) 
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
      }

      if (strcmp(line, "(time)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
      }
    }
    fclose(fp);
  }

#if PIDX_HAVE_MPI
  if ((*file)->idx_d->parallel_mode == 1)
  {
    MPI_Bcast((*file)->idx->bounds, 5, MPI_LONG_LONG, 0, (*file)->comm);
    MPI_Bcast((*file)->idx->chunk_size, 5, MPI_LONG_LONG, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->blocks_per_file), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->bits_per_block), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->variable_count), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->compression_bit_rate), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->compression_type), 1, MPI_INT, 0, (*file)->comm);

    if ((*file)->idx->compression_type == PIDX_CHUNKING_ZFP)
    {
      if ((*file)->idx->compression_bit_rate == 32)
        (*file)->idx->compression_factor = 2;
      if ((*file)->idx->compression_bit_rate == 16)
        (*file)->idx->compression_factor = 4;
      if ((*file)->idx->compression_bit_rate == 8)
        (*file)->idx->compression_factor = 8;
      if ((*file)->idx->compression_bit_rate == 4)
        (*file)->idx->compression_factor = 16;
      if ((*file)->idx->compression_bit_rate == 2)
        (*file)->idx->compression_factor = 32;
      if ((*file)->idx->compression_bit_rate == 1)
        (*file)->idx->compression_factor = 64;
    }
  }

  (*file)->idx_d->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);
  
  if(rank != 0)
  {
    for (var = 0; var < (*file)->idx->variable_count; var++)
    {
      (*file)->idx->variable[var] = malloc(sizeof (*((*file)->idx->variable[var])));
      memset((*file)->idx->variable[var], 0, sizeof (*((*file)->idx->variable[var])));
    }
  }
#endif
  
  for (var = 0; var < (*file)->idx->variable_count; var++)
  {
#if PIDX_HAVE_MPI
    if ((*file)->idx_d->parallel_mode == 1)
    {
      MPI_Bcast(&((*file)->idx->variable[var]->bits_per_value), 1, MPI_INT, 0, (*file)->comm);
      MPI_Bcast(&((*file)->idx->variable[var]->values_per_sample), 1, MPI_INT, 0, (*file)->comm);
      MPI_Bcast((*file)->idx->variable[var]->var_name, 512, MPI_CHAR, 0, (*file)->comm);
      MPI_Bcast((*file)->idx->variable[var]->type_name, 512, MPI_CHAR, 0, (*file)->comm);
    }
#endif
    (*file)->idx->variable[var]->sim_patch_count = 0;
  }

  //printf("%d %d %d %d %d\n", (*file)->idx->chunk_size[0], (*file)->idx->chunk_size[1], (*file)->idx->chunk_size[2], (*file)->idx->chunk_size[3], (*file)->idx->chunk_size[4]);

  if (rank == 0)
  {
    int ret;
    struct stat stat_buf;
    ret = stat((*file)->idx->filename, &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    (*file)->idx_d->fs_block_size = stat_buf.st_blksize;
  }

#if PIDX_HAVE_MPI
  if ((*file)->idx_d->parallel_mode == 1)
    MPI_Bcast(&((*file)->idx_d->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
  
  file_create_time = PIDX_get_time();
  return PIDX_success;
}


/// validate dims/blocksize
static PIDX_return_code PIDX_validate(PIDX_file file)
{
  int64_t dims;
  int64_t adjusted_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_bounds[0] = file->idx->bounds[0] / file->idx->chunk_size[0];
  adjusted_bounds[1] = file->idx->bounds[1] / file->idx->chunk_size[1];
  adjusted_bounds[2] = file->idx->bounds[2] / file->idx->chunk_size[2];
  adjusted_bounds[3] = file->idx->bounds[3] / file->idx->chunk_size[3];
  adjusted_bounds[4] = file->idx->bounds[4] / file->idx->chunk_size[4];
  
  //if (PIDX_inner_product(&dims, file->idx->bounds))
  if (PIDX_inner_product(&dims, adjusted_bounds))
    return PIDX_err_size;
  if (dims < file->idx_d->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx_d->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx->bits_per_block = getNumBits(file->idx_d->samples_per_block) - 1;
    //file->idx->bits_per_block = getNumBits(file->idx_d->samples_per_block);
  }
 
  if (file->idx->bits_per_block == 0)
  {
    file->idx->bits_per_block = 1;
    file->idx_d->samples_per_block = 2;
  }
  
  // other validations...
  // TODO
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_restructuring_box(PIDX_file file, PIDX_point reg_patch_size)
{
  if(reg_patch_size[0] < 0 || reg_patch_size[1] < 0 || reg_patch_size[2] < 0 || reg_patch_size[3] < 0 || reg_patch_size[4] < 0)
    return PIDX_err_box;
  
  if (file == NULL)
    return PIDX_err_file;
  
  memcpy(file->idx->reg_patch_size, reg_patch_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_dims(PIDX_file file, PIDX_point dims)
{
  if(dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_box;
  
  if(file == NULL)
    return PIDX_err_file;
  
  memcpy(file->idx->bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  file->idx->bounds[0] = file->idx->bounds[0] / file->idx_count[0];
  file->idx->bounds[1] = file->idx->bounds[1] / file->idx_count[1];
  file->idx->bounds[2] = file->idx->bounds[2] / file->idx_count[2];
  
  return PIDX_validate(file);
}



PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims)
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(dims, file->idx->bounds, (sizeof(int64_t) * PIDX_MAX_DIMENSIONS));
  
  return PIDX_success;
}


/*
PIDX_return_code PIDX_set_aggregation_factor(PIDX_file file, int agg_factor)
{
  if(!file)
    return PIDX_err_file;

  file->idx_d->agg_buffer->aggregation_factor = agg_factor;
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_aggregation_factor(PIDX_file file, int *agg_factor)
{
  if(!file)
    return PIDX_err_file;
  
  *agg_factor = file->idx_d->agg_buffer->aggregation_factor;
  
  return PIDX_success;
}
*/



PIDX_return_code PIDX_debug_output(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->debug_output = 1;

  return PIDX_success;
}


PIDX_return_code PIDX_set_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(file->idx->transform, transform, (sizeof(double) * 16));
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(transform, file->idx->transform, (sizeof(double) * 16));
  
  return PIDX_success;
}


/// TODO: get rid of this function
PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* filename, int current_time_step)
{
  int N;
  char dirname[1024], basename[1024];
  int nbits_blocknumber;
  char *directory_path;
  char *data_set_path;
  
  data_set_path = malloc(sizeof(*data_set_path) * 1024);
  memset(data_set_path, 0, sizeof(*data_set_path) * 1024);
  
  directory_path = malloc(sizeof(*directory_path) * 1024);
  memset(directory_path, 0, sizeof(*directory_path) * 1024);
  
  strncpy(directory_path, filename, strlen(filename) - 4);  
  sprintf(data_set_path, "%s/time%09d.idx", directory_path, current_time_step);
  free(directory_path);
  
  nbits_blocknumber = (file->idx_d->maxh - file->idx->bits_per_block - 1);
  VisusSplitFilename(data_set_path, dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--) 
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

#if 0
  //if i put . as the first character, if I move files VisusOpen can do path remapping
  sprintf(pidx->filename_template, "./%s", basename);
#endif
  //pidx does not do path remapping 
  strcpy(file->idx->filename_template, data_set_path);
  for (N = strlen(file->idx->filename_template) - 1; N >= 0; N--)
  {
    int ch = file->idx->filename_template[N];
    file->idx->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
    strcat(file->idx->filename_template, "/%01x.bin");
   
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
      strcat(file->idx->filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12) 
      strcat(file->idx->filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16) 
      strcat(file->idx->filename_template, "/%04x.bin"); //no directories, 65536  files
    else 
    {
      while (nbits_blocknumber > 16) 
      {
    strcat(file->idx->filename_template, "/%02x"); //256 subdirectories
	nbits_blocknumber -= 8;
      }
      strcat(file->idx->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }
  
  free(data_set_path);
  return PIDX_success;
}



PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  if(current_time_step < 0)
    return PIDX_err_time;
   
  file->idx->current_time_step = current_time_step;
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  *current_time_step = file->idx->current_time_step;
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_block_size(PIDX_file file, const int bits_per_block)
{
  if(!file)
    return PIDX_err_file;
  
  if(bits_per_block <= 0)
    return PIDX_err_block;
   
  file->idx->bits_per_block = bits_per_block;
  file->idx_d->samples_per_block = (int)pow(2, bits_per_block);
  
  return PIDX_validate(file);
}



PIDX_return_code PIDX_get_block_size(PIDX_file file, int* bits_per_block)
{ 
  if(!file)
    return PIDX_err_file;
  
  *bits_per_block = file->idx->bits_per_block;
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_block_count(PIDX_file file, const int blocks_per_file)
{
  if(!file)
    return PIDX_err_file;
  
  if(blocks_per_file <= 0)
    return PIDX_err_block;
   
  file->idx->blocks_per_file = blocks_per_file;
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_block_count(PIDX_file file, int* blocks_per_file)
{ 
  if(!file)
    return PIDX_err_file;
  
  *blocks_per_file = file->idx->blocks_per_file;
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_count <= 0)
    return PIDX_err_count;
   
  file->idx->variable_count = variable_count;
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count)
{ 
  if(!file)
    return PIDX_err_file;
  
  *variable_count = file->idx->variable_count;
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type)
{
  if(!file)
    return PIDX_err_file;

  if (compression_type != PIDX_NO_COMPRESSION && compression_type != PIDX_CHUNKING_ONLY && compression_type != PIDX_CHUNKING_ZFP)
    return PIDX_err_unsupported_compression_type;

  file->idx->compression_type = compression_type;

  if (file->idx->compression_type == PIDX_NO_COMPRESSION)
    return PIDX_success;
  else if (file->idx->compression_type == PIDX_CHUNKING_ONLY || file->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    PIDX_point chunk_size;

    chunk_size[0] = 4;
    chunk_size[1] = 4;
    chunk_size[2] = 4;
    chunk_size[3] = 1;
    chunk_size[4] = 1;

    if(chunk_size[0] < 0 || chunk_size[1] < 0 || chunk_size[2] < 0 || chunk_size[3] < 0 || chunk_size[4] < 0)
      return PIDX_err_box;

    memcpy(file->idx->chunk_size, chunk_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

    int reduce_by_sample = 1;
    uint64_t total_chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2] * file->idx->chunk_size[3] * file->idx->chunk_size[4];
    if (reduce_by_sample == 1)
    {
      file->idx->bits_per_block = file->idx->bits_per_block - (int)log2(total_chunk_size);
      file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

      if (file->idx->bits_per_block <= 0)
      {
        file->idx->bits_per_block = 1;
        file->idx_d->samples_per_block = 2;
      }
    }
    else
      file->idx->blocks_per_file = file->idx->blocks_per_file / total_chunk_size;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type)
{
  if(!file)
    return PIDX_err_file;

   *compression_type = file->idx->compression_type;

  return PIDX_success;
}



PIDX_return_code PIDX_set_lossy_compression_bit_rate(PIDX_file file, int compression_bit_rate)
{
  if (!file)
    return PIDX_err_file;

  if (file->idx->compression_type != PIDX_CHUNKING_ZFP)
    return PIDX_err_unsupported_compression_type;

  file->idx->compression_bit_rate = compression_bit_rate;
  
  if (file->idx->compression_bit_rate == 32)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 1;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
    file->idx->compression_factor = 2;
  }
  if (file->idx->compression_bit_rate == 16)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 2;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 4;
  }
  if (file->idx->compression_bit_rate == 8)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 3;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 8;
  }
  if (file->idx->compression_bit_rate == 4)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 4;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 16;
  }
  if (file->idx->compression_bit_rate == 2)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 5;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 32;
  }
  if (file->idx->compression_bit_rate == 1)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 6;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);

    file->idx->compression_factor = 64;
  }
  
  return PIDX_success;
}



PIDX_return_code PIDX_disable_rst(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_agg = 0;

  return PIDX_success;
}




PIDX_return_code PIDX_get_lossy_compression_bit_rate(PIDX_file file, int *compression_bit_rate)
{
  if(!file)
    return PIDX_err_file;

  *compression_bit_rate = file->idx->compression_bit_rate;

  return PIDX_success;
}



PIDX_return_code PIDX_variable_create(char* variable_name, unsigned int bits_per_sample, PIDX_data_type type_name, PIDX_variable* variable)
{
  if (!variable_name)
    return PIDX_err_name;

  if (bits_per_sample <= 0)
    return PIDX_err_size;

  if (!type_name)
    return PIDX_err_type;

  *variable = malloc(sizeof *(*variable));
  memset(*variable, 0, sizeof *(*variable));

  int bits = 0;
  PIDX_default_bits_per_datatype(type_name, &bits);
  if (bits !=0 && bits != bits_per_sample)
    return PIDX_err_type;

  (*variable)->values_per_sample = 1;
  (*variable)->bits_per_value = (bits_per_sample/1);

  /*
  if (strcmp(type_name, FLOAT64)  == 0)
  {
    (*variable)->values_per_sample = 1;
    (*variable)->bits_per_value = (bits_per_sample/1);
  }
  else
  {
    (*variable)->values_per_sample = 3;
    (*variable)->bits_per_value = (bits_per_sample/3);
  }
  */

  strcpy((*variable)->type_name, type_name);
  strcpy((*variable)->var_name, variable_name);

  (*variable)->sim_patch_count = 0;
  (*variable)->dump_meta_data_ON = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_get_next_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;

  *variable = file->idx->variable[file->idx->variable_index_tracker];

  return PIDX_success;
}



PIDX_return_code PIDX_reset_variable_counter(PIDX_file file)
{
  if (!file)
    return PIDX_err_file;

  file->idx->variable_index_tracker = 0;
  file->local_variable_count = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_append_and_write_variable(PIDX_file file, PIDX_variable variable/*, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout*/)
{

  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  if (file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_variable;

  variable->io_state = 1;

  file->idx->variable[file->idx->variable_index_tracker] = variable;

  file->idx->variable_index_tracker++;
  file->local_variable_count++;

  return PIDX_success;
}



PIDX_return_code PIDX_read_next_variable(PIDX_file file, PIDX_variable variable)
{
  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  if (file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_variable;

  variable->io_state = 0;

  variable = file->idx->variable[file->idx->variable_index_tracker];

  file->idx->variable_index_tracker++;
  file->local_variable_count++;

  return PIDX_success;
}


PIDX_return_code PIDX_variable_write_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;

  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;

  if(!dims || dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_count;

#if !SIMULATE_IO
//  if (read_from_this_buffer == NULL)
//    return PIDX_err_block;
#endif

  const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

  temp_buffer = read_from_this_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = (unsigned char*)temp_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  return PIDX_success;
}



PIDX_return_code PIDX_variable_read_data_layout(PIDX_variable variable, PIDX_point offset, PIDX_point dims, void* write_to_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;

  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;

  if(!dims || dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_count;

  //const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

  //write_to_this_buffer = temp_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = write_to_this_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  return PIDX_success;
}


static PIDX_return_code populate_idx_file_structure(PIDX_file file)
{
  PointND bounds_point;
  int d = 0, i = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    if (file->idx->bounds[d] % file->idx->chunk_size[d] == 0)
      file->idx->chunked_bounds[d] = (int) file->idx->bounds[d] / file->idx->chunk_size[d];
    else
      file->idx->chunked_bounds[d] = (int) (file->idx->bounds[d] / file->idx->chunk_size[d]) + 1;
  }

  int64_t* cb = file->idx->chunked_bounds;
  bounds_point.x = (int) cb[0];
  bounds_point.y = (int) cb[1];
  bounds_point.z = (int) cb[2];
  bounds_point.u = (int) cb[3];
  bounds_point.v = (int) cb[4];
  GuessBitmaskPattern(file->idx->bitSequence, bounds_point);
  file->idx_d->maxh = strlen(file->idx->bitSequence);

  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

  int64_t total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]) * getPowerOf2(cb[3]) * getPowerOf2(cb[4]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int64_t max_sample_per_file = (uint64_t) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

  return PIDX_success;
}


static PIDX_return_code populate_idx_layout(PIDX_file file, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level)
{
  int i, j;
  int p = 0, ctr = 1;
  PIDX_return_code ret_code;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  int lvi = file->local_variable_index;

  if (file->idx_d->parallel_mode == 1 && file->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    PIDX_block_layout all_patch_local_block_layout = malloc(sizeof (*all_patch_local_block_layout));
    memset(all_patch_local_block_layout, 0, sizeof (*all_patch_local_block_layout));
    ret_code = PIDX_blocks_initialize_layout(all_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (p = 0 ; p < file->idx->variable[lvi]->sim_patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
        bounding_box[0][i] = file->idx->variable[lvi]->sim_patch[p]->offset[i];
        bounding_box[1][i] = file->idx->variable[lvi]->sim_patch[p]->size[i] + file->idx->variable[lvi]->sim_patch[p]->offset[i];

        bounding_box[0][i] = (bounding_box[0][i] / file->idx->chunk_size[i]);

        if (bounding_box[1][i] % file->idx->chunk_size[i] == 0)
          bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]);
        else
          bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]) + 1;
      }

      PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
      memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
      ret_code = PIDX_blocks_initialize_layout(per_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bitPattern, per_patch_local_block_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      if (all_patch_local_block_layout->resolution_from <= all_patch_local_block_layout->bits_per_block)
      {
        for (i = all_patch_local_block_layout->resolution_from ; i <= all_patch_local_block_layout->bits_per_block ; i++)
        {
          if (per_patch_local_block_layout->hz_block_number_array[i][0] == 0)
          {
            all_patch_local_block_layout->hz_block_number_array[i][0] = per_patch_local_block_layout->hz_block_number_array[i][0];
            break;
          }
        }

        ctr = 1;
        for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
        {
          for (j = 0 ; j < ctr ; j++)
          {
            if(per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
              all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
          }
          ctr = ctr * 2;
        }
      }
      else
      {
        ctr = 1;
        for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
        {
          if (i >= all_patch_local_block_layout->resolution_from)
          {
            for (j = 0 ; j < ctr ; j++)
            {
              if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
                all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
            }
          }
          ctr = ctr * 2;
        }
      }

      PIDX_blocks_free_layout(per_patch_local_block_layout);
      free(per_patch_local_block_layout);
      per_patch_local_block_layout = 0;
    }

    if (block_layout->resolution_from <= block_layout->bits_per_block)
    {
      int level_count = 1;
      for (i = block_layout->resolution_from; i <= block_layout->bits_per_block; i++)
      {
#if PIDX_HAVE_MPI
        if (file->idx_d->parallel_mode == 1)
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
        else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
        memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
      }

      for (i = block_layout->bits_per_block + 1; i < (block_layout->resolution_to); i++)
      {
#if PIDX_HAVE_MPI
        if (file->idx_d->parallel_mode == 1)
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
        else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
        memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
        level_count = level_count * 2;
      }
    }
    else
    {
      int level_count = 1;
      for (i = block_layout->bits_per_block + 1; i < (block_layout->resolution_to); i++)
      {
        if (i >= block_layout->resolution_from)
        {
#if PIDX_HAVE_MPI
          if (file->idx_d->parallel_mode == 1)
            MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
          else
            memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
        }
        level_count = level_count * 2;
      }
    }

    PIDX_blocks_free_layout(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
  }
  else
  {
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      bounding_box[0][i] = 0;
      bounding_box[1][i] = file->idx->chunked_bounds[i];
    }

    ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bitPattern, block_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }


  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  int file_number = 0;
  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
        break;
      }
    }

    ctr = 1;
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->block_count_per_file[file_number]++;
        }
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }


  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;

      count++;
    }
  }

  return PIDX_success;
}


static PIDX_return_code delete_idx_dataset(PIDX_file file)
{
  int lvi = file->local_variable_index;
  PIDX_variable var = file->idx->variable[lvi];

  PIDX_blocks_free_layout(file->idx->variable[lvi]->global_block_layout);
  PIDX_free_layout(file->idx->variable[lvi]->global_block_layout);

  free(file->idx->variable[lvi]->global_block_layout);
  file->idx->variable[lvi]->global_block_layout = 0;

  int i = 0;
  for (i = 0; i < file->layout_count ; i++)
  {
    PIDX_blocks_free_layout(var->block_layout_by_level[i]);
    PIDX_free_layout(var->block_layout_by_level[i]);

    free(var->block_layout_by_level[i]);
    var->block_layout_by_level[i] = 0;
  }
  free(var->block_layout_by_level);
  var->block_layout_by_level = 0;

  return PIDX_success;
}


static PIDX_return_code populate_idx_dataset(PIDX_file file)
{
  int rank = 0;
  int i = 0, j = 0, ctr;
  int file_number = 0;

  PIDX_return_code ret_code;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
    MPI_Comm_rank(file->comm, &rank);
#endif

  ret_code = populate_idx_file_structure(file);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int lvi = file->local_variable_index;
  int lower_hz_level = 0, higher_hz_level = 0;

  PIDX_variable var = file->idx->variable[lvi];
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  file->idx->variable[lvi]->global_block_layout = malloc(sizeof (*file->idx->variable[lvi]->global_block_layout));
  memset(file->idx->variable[lvi]->global_block_layout, 0, sizeof (*file->idx->variable[lvi]->global_block_layout));
  PIDX_block_layout block_layout = file->idx->variable[lvi]->global_block_layout;

  lower_hz_level = file->reduced_res_from;
  higher_hz_level = file->idx_d->maxh - file->reduced_res_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

#if 1
  file->layout_count = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->layout_count <= 0)
    file->layout_count = 1;
#else
  file->layout_count = 1;
#endif

  //file->layout_start_index = 0;
  //file->layout_end_index = file->layout_count;

  var->block_layout_by_level = malloc(sizeof(*(var->block_layout_by_level)) * file->layout_count);
  memset(var->block_layout_by_level, 0, sizeof(*(var->block_layout_by_level)) * file->layout_count);
  for (i = 0; i < file->layout_count ; i++)
  {
    var->block_layout_by_level[i] = malloc(sizeof(*(var->block_layout_by_level[i])));
    memset(var->block_layout_by_level[i], 0, sizeof(*(var->block_layout_by_level[i])));
  }

#if 1
  lower_level_low_layout = lower_hz_level;
  higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

  //if (higher_level_low_layout >= file->idx_d->maxh)
  //  higher_level_low_layout = file->idx_d->maxh;

  if (higher_level_low_layout >= higher_hz_level)
    higher_level_low_layout = higher_hz_level;

  ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret_code = populate_idx_layout(file, file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
    memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int));

  ctr = 1;
  int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
  if (temp_level >= higher_hz_level)
    temp_level = higher_hz_level;
  for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
  {
    memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int) * ctr);
    ctr = ctr * 2;
  }
#else
  lower_level_low_layout = 0;
  higher_level_low_layout = file->idx_d->maxh;

  if (higher_level_low_layout >= file->idx_d->maxh)
    higher_level_low_layout = file->idx_d->maxh;

  ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  ret_code = populate_idx_layout(file, file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (j = 0 ; j < file->idx->bits_per_block + 1 ; j++)
    memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int));

  ctr = 1;
  int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
  if (temp_level >= file->idx_d->maxh)
    temp_level = file->idx_d->maxh;

  for (j = file->idx->bits_per_block + 1 ; j < file->idx_d->maxh ; j++)
  {
    memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int) * ctr);
    ctr = ctr * 2;
  }
#endif


  for (i = 1; i < file->layout_count; i++)
  {
    lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
    higher_level_higher_layout = lower_level_higher_layout + 1;

    //printf("Y [%d] %d %d\n", file->layout_count, lower_level_higher_layout, higher_level_higher_layout);
    ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    ret_code = populate_idx_layout(file, file->idx->variable[lvi]->block_layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], file->idx->variable[lvi]->block_layout_by_level[i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
    ctr = ctr * 2;
  }



  //for (j = temp_level ; j < file->idx_d->maxh ; j++)
  //{

  //}


  /*
  //lower_hz_level = 0;
  //higher_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

  //lower_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
  //higher_hz_level = file->idx_d->maxh;

  //lower_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
  //higher_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 2;

  //lower_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 2;
  //higher_hz_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 3;

  lower_hz_level = 0;
  higher_hz_level = file->idx_d->maxh;

  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  ret_code = populate_idx_layout(file, block_layout, lower_hz_level, higher_hz_level);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  */

  /*
  if (rank == 0)
  {
    printf("[A] Final Block Bitmap [%d %d]\n", file->idx->variable[lvi]->block_layout_by_level[0]->resolution_from, file->idx->variable[lvi]->block_layout_by_level[0]->resolution_to);
    PIDX_blocks_print_layout(file->idx->variable[lvi]->block_layout_by_level[0]);
    printf("[B] Final Block Bitmap\n");
    PIDX_blocks_print_layout(block_layout);
  }
  */

  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
        break;
      }
    }

    ctr = 1;
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->block_count_per_file[file_number]++;
        }
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }


  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      //if (rank == 0)
      //  printf("BPF %d = %d FI = %d\n", i, file->idx->variable[lvi]->block_count_per_file[i], file->idx->variable[lvi]->file_index[i]);
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  //if (rank == 0)
  //{
  //  for (i = 0; i < file->idx_d->max_file_count; i++)
  //  printf("[X] i(%d) = count(%d)\n", i, block_layout->inverse_existing_file_index[i]);
  //}

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->debug_do_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_do_chunk = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_compression(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->debug_compress = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_hz(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_do_hz = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_do_agg = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_io(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_do_io = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_rst = debug_rst;
  
  return PIDX_success;
}



PIDX_return_code PIDX_dump_agg_info(PIDX_file file, int dump_agg_info)
{
  if(!file)
    return PIDX_err_file;
  
  char filename_skeleton[512];
  file->idx_d->dump_agg_info = dump_agg_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_d->agg_dump_dir_name, "%s_agg_dump", filename_skeleton);
  
  return PIDX_success;
}



PIDX_return_code PIDX_dump_io_info(PIDX_file file, int dump_io_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_d->dump_io_info = dump_io_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_d->io_dump_dir_name, "%s_io_dump", filename_skeleton);

  return PIDX_success;
}



PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_hz = debug_hz;
  
  return PIDX_success;
}




PIDX_return_code PIDX_set_resolution(PIDX_file file, int hz_from, int hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  file->reduced_res_from = hz_from;
  file->reduced_res_to = hz_to;
  
  return PIDX_success;
}


PIDX_return_code PIDX_get_resolution(PIDX_file file, int *hz_from, int *hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  *hz_from = file->reduced_res_from;
  *hz_to = file->reduced_res_to;
  
  return PIDX_success;
}




PIDX_return_code PIDX_enable_raw_io(PIDX_file file)
{
  if(file == NULL)
    return PIDX_err_file;

  file->enable_raw_dump = 1;

  return PIDX_success;
}




static PIDX_return_code PIDX_parameter_validate(PIDX_file file, int start_var_index, int end_var_index)
{
  int d, nprocs = 1;

#if PIDX_HAVE_MPI
  int local_patch_count = 0, total_patch_count = 0;
  int ret;
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm, &nprocs);
    PIDX_variable var0 = file->idx->variable[start_var_index];
    if (file->idx->enable_rst == 1)
    {
      /// Calculate the total number of patches across all processes.
      local_patch_count = var0->sim_patch_count;
      ret = MPI_Allreduce(&local_patch_count, &total_patch_count, 1, MPI_INT, MPI_SUM, file->comm);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI error\n", __FILE__, __LINE__);
        return PIDX_err_mpi;
      }

      /// If total patch count != total number of processes
      /// Then NO restructuring
      if (total_patch_count != nprocs)
        file->idx->enable_rst = 0;
    }
  }
  else
  {
    file->idx->enable_rst = 0;
    file->idx->enable_agg = 0;
  }
#else
  file->idx->enable_rst = 0;
#endif

  if (file->idx->compression_type == PIDX_CHUNKING_ONLY || file->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    // No Chunking and compression without restructuring
    if (file->idx->enable_rst != 1)
    {
      file->idx->compression_type = PIDX_NO_COMPRESSION;
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        file->idx->chunk_size[d] = 1;
      file->idx->compression_bit_rate = 64;
    }
    else
    {
      /*
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        if (file->idx->bounds[d] % file->idx->chunk_size[d] != 0)
        {
          file->idx->compression_type = PIDX_NO_COMPRESSION;
          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            file->idx->chunk_size[d] = 1;
          file->idx->compression_bit_rate = 64;

          break;
        }
      }
      */
    }
  }
  //file->var_pipe_length = file->idx->variable_count - 1;

  /*
  if (file->idx->variable_count == file->idx->variable_index_tracker)
  {
    if (file->var_pipe_length >= file->idx->variable_count || file->var_pipe_length == -1)
      file->var_pipe_length = file->idx->variable_count - 1;

    if (file->idx->enable_agg != 0)
    {
      int v, no_of_aggregators = 0;
      for (v = start_var_index; v < start_var_index + file->var_pipe_length ; v++)
        no_of_aggregators = no_of_aggregators + file->idx->variable[v]->values_per_sample * file->idx->variable[file->local_variable_index]->global_block_layout->existing_file_count;

      //for (v = 0; v < file->idx->variable_count ; v++)
      //  no_of_aggregators = no_of_aggregators + file->idx->variable[v]->values_per_sample * file->idx->variable[file->local_variable_index]->existing_file_count;

      if (nprocs < no_of_aggregators)
      {
        if ((file->idx_d->max_file_count == 1 || nprocs == 1))
          file->idx->enable_agg = 0;
        else
          file->idx->enable_agg = 1;
      }
      else
        file->idx->enable_agg = 2;
    }

  }
  else if (file->idx->variable_count > file->idx->variable_index_tracker)
  {
    if (file->var_pipe_length >= file->idx->variable_count || file->var_pipe_length == -1)
      file->var_pipe_length = file->idx->variable_count - 1;

    if (file->idx->enable_agg != 0)
    {
      int v, no_of_aggregators = 0;
      for (v = start_var_index; v < end_var_index ; v++)
        no_of_aggregators = no_of_aggregators + file->idx->variable[v]->values_per_sample * file->idx->variable[file->local_variable_index]->global_block_layout->existing_file_count;

      //for (v = 0; v < file->idx->variable_count ; v++)
      //  no_of_aggregators = no_of_aggregators + file->idx->variable[v]->values_per_sample * file->idx->variable[file->local_variable_index]->existing_file_count;

      if (nprocs < no_of_aggregators)
      {
        if ((file->idx_d->max_file_count == 1 || nprocs == 1))
          file->idx->enable_agg = 0;
        else
          file->idx->enable_agg = 1;
      }
      else
        file->idx->enable_agg = 2;
    }
  }
  */

  return PIDX_success;
}



static PIDX_return_code PIDX_raw_write(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  file->var_pipe_length = file->idx->variable_count - 1;
  if (file->var_pipe_length == 0)
    file->var_pipe_length = 1;


  PIDX_return_code ret;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm,  &nprocs);
  }
#endif

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);

    file->idx->enable_rst = 0;
  }
#else
  file->idx->enable_rst = 0;
#endif

  populate_idx_start_time = MPI_Wtime();

  PIDX_init_timming_buffers2(file);
  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;

    file->one_time_initializations = 1;
  }

  populate_idx_end_time = MPI_Wtime();

  write_init_start[hp] = PIDX_get_time();

#if !SIMULATE_IO
  if (file->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif

    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
#endif

  write_init_end[hp] = PIDX_get_time();
  hp++;

  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);


    /*--------------------------------------Create RST IDs [start]------------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*----------------------------------------Create RST IDs [end]------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    /* Attaching the communicator to the restructuring phase */
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/



    /*--------------------------------------------RST [start]------------------------------------------------*/
    rst_start[vp] = PIDX_get_time();

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Perform data restructuring */
    if (file->debug_do_rst == 1)
    {
      ret = PIDX_rst_write(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

    if (file->debug_do_io == 1)
    {
      PIDX_rst_buf_aggregate_write(file->rst_id);
      //PIDX_rst_buf_aggregate_read(file->rst_id);
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
    rst_end[vp] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/



    /*-------------------------------------------finalize [start]---------------------------------------------*/
    finalize_start[vp] = PIDX_get_time();

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    finalize_end[vp] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

    vp++;
  }

    free(file->idx_d->rank_r_offset);
    file->idx_d->rank_r_offset = 0;

    free(file->idx_d->rank_r_count);
    file->idx_d->rank_r_count = 0;

  return PIDX_success;
}




static PIDX_return_code PIDX_raw_read(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  file->var_pipe_length = file->idx->variable_count - 1;
  if (file->var_pipe_length == 0)
    file->var_pipe_length = 1;


  PIDX_return_code ret;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm,  &nprocs);
  }
#endif

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  }
#endif

  populate_idx_start_time = MPI_Wtime();

  PIDX_init_timming_buffers2(file);
  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;

    file->one_time_initializations = 1;
  }

  populate_idx_end_time = MPI_Wtime();

  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);


    /*--------------------------------------Create RST IDs [start]------------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*----------------------------------------Create RST IDs [end]------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    /* Attaching the communicator to the restructuring phase */
    if (file->idx_d->parallel_mode == 1)
    {
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/



    /*--------------------------------------------RST [start]------------------------------------------------*/
    rst_start[vp] = PIDX_get_time();

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    if (file->debug_do_io == 1)
    {
      ret = PIDX_rst_buf_aggregate_read(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    /* Perform data restructuring */
    if (file->debug_do_rst == 1)
    {
      ret = PIDX_rst_read(file->rst_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }
    }

    rst_end[vp] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/



    /*-------------------------------------------finalize [start]---------------------------------------------*/
    finalize_start[vp] = PIDX_get_time();

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout, "Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_rst;
    }

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    finalize_end[vp] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

    vp++;
  }

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  return PIDX_success;
}




static PIDX_return_code PIDX_write(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  file->var_pipe_length = file->idx->variable_count - 1;
  if (file->var_pipe_length == 0)
    file->var_pipe_length = 1;

#if PIDX_DEBUG_OUTPUT
  unsigned long long l_populate = 0, g_populate = 0;
  unsigned long long l_filec = 0, g_filec = 0;
  unsigned long long l_init = 0, g_init = 0;
  unsigned long long l_rst_buf = 0, g_rst_buf = 0;
  unsigned long long l_chunk_buf = 0, g_chunk_buf = 0;
  unsigned long long l_rst = 0, g_rst = 0;
  unsigned long long l_chunk = 0, g_chunk = 0;
  unsigned long long l_cmp = 0, g_cmp = 0;
  unsigned long long l_hz_buf = 0, g_hz_buf = 0;
  unsigned long long l_hz = 0, g_hz = 0;
  unsigned long long l_agg_buf = 0, g_agg_buf = 0;
  unsigned long long l_agg = 0, g_agg = 0;
  unsigned long long l_io = 0, g_io = 0;
  unsigned long long l_pidx = 0, g_pidx = 0;
#endif

  int j = 0, p, var = 0, d = 0;
  int total_header_size;
  PIDX_return_code ret;
  //static int header_io = 0;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm,  &nprocs);
  }
#endif

  if (file->idx_count[0] != 1 || file->idx_count[1] != 1 || file->idx_count[2] != 1 )
    for (var = start_var_index; var < end_var_index; var++)
      for (p = 0; p < file->idx->variable[var]->sim_patch_count; p++)
        for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
          for (j = 0; j < file->idx->bounds[d] * file->idx_count[d]; j = j + (file->idx->bounds[d]))
            if (file->idx->variable[var]->sim_patch[p]->offset[d] >= j && file->idx->variable[var]->sim_patch[p]->offset[d] < (j + file->idx->bounds[d]))
            {
              file->idx->variable[var]->sim_patch[p]->offset[d] = file->idx->variable[var]->sim_patch[p]->offset[d] - j;
              break;
            }

  populate_idx_start_time = MPI_Wtime();

  ret = populate_idx_dataset(file);
  if (ret != PIDX_success)
    return PIDX_err_file;

  PIDX_init_timming_buffers2(file);

#if PIDX_DEBUG_OUTPUT
  l_populate = 1;
  MPI_Allreduce(&l_populate, &g_populate, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_populate == nprocs)
    printf("Finished Populating IDX dataset (File Count %d maxh = %d)\n", file->idx_d->max_file_count, file->idx_d->maxh);
#endif

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;
    
    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;
    
    ret = PIDX_parameter_validate(file, start_var_index, end_var_index);
    if (ret != PIDX_success)
      return PIDX_err_file;

    file->one_time_initializations = 1;
  }

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  }
#endif

  /*  STEP 1: create files and folders based on the extents of the variable group
   *  STEP 2: if flush used
   *            when variable count is met, then write header information
   *          else
   *            pass the header buffer to the agg phase when no var pipelining is done (else if pipe, then go to if)
   *  STEP 3: at the end of all IO, write the .idx file
   */

  populate_idx_end_time = MPI_Wtime();

  write_init_start[hp] = PIDX_get_time();

#if !SIMULATE_IO

  if (file->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif
    ret = PIDX_header_io_file_create(file->header_io_id, file->idx->variable[file->local_variable_index]->global_block_layout);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (file->idx->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[file->local_variable_index]->global_block_layout,  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      file->flush_used = 1;
    }

    if (file->idx->variable_index_tracker == file->idx->variable_count)
    {
      // Write the header
      if (file->flush_used == 1 /*|| file->idx->enable_agg == 0 || file->idx->enable_agg == 1*/)
      {
        ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[file->local_variable_index]->global_block_layout, 1);
        if (ret != PIDX_success)
          return PIDX_err_header;
      }
      else if (/*(file->var_pipe_length < file->idx->variable_count - 1) && */ caching_state == 0)
      {
        ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[file->local_variable_index]->global_block_layout, 1);
        if (ret != PIDX_success)
          return PIDX_err_header;
      }
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }

#endif

  write_init_end[hp] = PIDX_get_time();
  hp++;

#if PIDX_DEBUG_OUTPUT
  l_filec = 1;
  MPI_Allreduce(&l_filec, &g_filec, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_filec == nprocs)
    printf("Finished Creating File heirarchy\n");
#endif


  int start_index = 0, end_index = 0;
  int i = 0;

  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    startup_start[vp] = PIDX_get_time();
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);

    int agg_io_level = 0, no_of_aggregators = 0;

    if (file->agg_type != 0)
    {
      for (i = 0; i < file->layout_count ; i++)
      {
        no_of_aggregators = file->idx->variable[file->local_variable_index]->block_layout_by_level[i]->existing_file_count;
        if (no_of_aggregators <= nprocs)
          agg_io_level = i;
      }
      agg_io_level = agg_io_level + 1;
    }
    else
    {
      no_of_aggregators = file->idx->variable[file->local_variable_index]->global_block_layout->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level = file->layout_count;
      else
        agg_io_level = 0;
    }

    if (file->idx->enable_agg == 0)
      agg_io_level = 0;

    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the chunking ID */
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the compression ID */
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the HZ encoding ID */
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the aggregation ID */
    /* Create the I/O ID */
    file->tagg_id = malloc(sizeof(*(file->tagg_id)) * file->idx->variable_count);
    memset(file->tagg_id, 0, sizeof(*(file->tagg_id)) * file->idx->variable_count);

    file->tio_id = malloc(sizeof(*(file->tio_id)) * file->idx->variable_count);
    memset(file->tio_id, 0, sizeof(*(file->tio_id)) * file->idx->variable_count);

    int agg_var_pipe = 0;
    int agg_end_index = 0;

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      file->tagg_id[i] = malloc(sizeof(*(file->tagg_id[i])) * file->layout_count);
      memset(file->tagg_id[i], 0, sizeof(*(file->tagg_id[i])) * file->layout_count);

      file->tio_id[i] = malloc(sizeof(*(file->tio_id[i])) * file->layout_count);
      memset(file->tio_id[i], 0, sizeof(*(file->tio_id[i])) * file->layout_count);
    }

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      agg_end_index = ((i + agg_var_pipe) >= (end_index + 1)) ? (end_index) : (i + agg_var_pipe);
      for(j = 0 ; j < file->layout_count; j++)
      {
        file->tagg_id[i][j] = PIDX_agg_init(file->idx, file->idx_d, start_var_index, i, agg_end_index);
        file->tio_id[i][j] = PIDX_io_init(file->idx, file->idx_d, start_var_index, i, agg_end_index);
      }
    }

    file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
    memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      file->idx_d->agg_buffer[i] = malloc(sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      memset(file->idx_d->agg_buffer[i], 0, sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      for(j = 0 ; j < agg_io_level; j++)
      {
        file->idx_d->agg_buffer[i][j] = malloc(sizeof(*(file->idx_d->agg_buffer[i][j])) );
        memset(file->idx_d->agg_buffer[i][j], 0, sizeof(*(file->idx_d->agg_buffer[i][j])) );

        file->idx_d->agg_buffer[i][j]->file_number = -1;
        file->idx_d->agg_buffer[i][j]->var_number = -1;
        file->idx_d->agg_buffer[i][j]->sample_number = -1;

        file->idx_d->agg_buffer[i][j]->no_of_aggregators = 0;
        file->idx_d->agg_buffer[i][j]->aggregator_interval = 0;
        file->idx_d->agg_buffer[i][j]->aggregation_factor = 1;
      }

      /*
      for (j = 1 ; j < agg_io_level; j++)
        file->idx_d->agg_buffer[i][j]->aggregation_factor = 1;//(int)pow(2, (agg_io_level - j));
      */
    }
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      /* Attaching the communicator to the restructuring phase */
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      /* Attaching the communicator to the chunking phase */
      ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_chunk;

      /* Attaching the communicator to the compression phase */
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_compress;

      /* Attaching the communicator to the HZ encodig phase phase */
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_hz;

      /* Attaching the communicator to the aggregation phase */
      /* Attaching the communicator to the I/O phase */
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = 0 ; j < file->layout_count; j++)
        {
          ret = PIDX_agg_set_communicator(file->tagg_id[i][j], file->comm);
          if (ret != PIDX_success)
            return PIDX_err_agg;

          ret = PIDX_io_set_communicator(file->tio_id[i][j], file->comm);
          if (ret != PIDX_success)
            return PIDX_err_io;
        }
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    //if (file->small_agg_comm == 1)
    //  PIDX_create_local_aggregation_comm(file->agg_id);

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
      for (j = 0 ; j < agg_io_level; j++)
      {
        if (file->agg_type == 0)
          ret = PIDX_agg_meta_data_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], file->idx->variable[file->local_variable_index]->global_block_layout);

        else if (file->agg_type == 1)
          ret = PIDX_local_agg_meta_data_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);

        else
          ret = PIDX_local_agg_local_comm_meta_data_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }
    }
    startup_end[vp] = PIDX_get_time();

#if PIDX_DEBUG_OUTPUT
    l_init = 1;
    MPI_Allreduce(&l_init, &g_init, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_init == nprocs)
      printf("All Modules Initialized (Variable index [%d %d] Variable Pipe length %d)\n", start_index, end_index, file->var_pipe_length);
#endif

    /*--------------------------------------------RST [start]------------------------------------------------*/
    rst_start[vp] = PIDX_get_time();

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

#if PIDX_DEBUG_OUTPUT
    l_rst_buf = 1;
    MPI_Allreduce(&l_rst_buf, &g_rst_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_rst_buf == nprocs)
      printf("[R] Restructuring Buffer Created\n");
#endif

    /* Perform data restructuring */
    if (file->debug_do_rst == 1)
    {
      ret = PIDX_rst_write(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

#if PIDX_DEBUG_OUTPUT
    l_rst = 1;
    MPI_Allreduce(&l_rst, &g_rst, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_rst == nprocs)
      printf("[R] Restructuring Completed\n");
#endif

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
    rst_end[vp] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/


    /*----------------------------------------Chunking [start]------------------------------------------------*/
    chunk_start[vp] = PIDX_get_time();

    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

#if PIDX_DEBUG_OUTPUT
    l_chunk_buf = 1;
    MPI_Allreduce(&l_chunk_buf, &g_chunk_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk_buf == nprocs)
      printf("[C] Chunking Buffer Created\n");
#endif

    /* Perform Chunking */
    if (file->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

#if PIDX_DEBUG_OUTPUT
    l_chunk = 1;
    MPI_Allreduce(&l_chunk, &g_chunk, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk == nprocs)
      printf("[C] Chunking Completed\n");
#endif


    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_rst_buf_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    chunk_end[vp] = PIDX_get_time();
    /*-----------------------------------------Chunking [end]------------------------------------------------*/


    /*----------------------------------------Compression [start]--------------------------------------------*/
    compression_start[vp] = PIDX_get_time();

    /* Perform Compression */
    if (file->debug_do_compress == 1)
    {
#if !SIMULATE_IO
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
        return PIDX_err_compress;
#endif
    }

#if PIDX_DEBUG_OUTPUT
    l_cmp = 1;
    MPI_Allreduce(&l_cmp, &g_cmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_cmp == nprocs)
      printf("[CMP] Compression Completed\n");
#endif

    compression_end[vp] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/



    /*---------------------------------------------HZ [start]------------------------------------------------*/
    hz_start[vp] = PIDX_get_time();

    /* Creating the buffers required for HZ encoding */
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

#if PIDX_DEBUG_OUTPUT
    l_hz_buf = 1;
    MPI_Allreduce(&l_hz_buf, &g_hz_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz_buf == nprocs)
      printf("[H] HZ Buffer Created\n");
#endif

    /* Perform HZ encoding */
    if (file->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

#if PIDX_DEBUG_OUTPUT
    l_hz = 1;
    MPI_Allreduce(&l_hz, &g_hz, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz == nprocs)
      printf("[H] HZ Encoding Finished\n");
#endif

    /* Verify the HZ encoding */
    if(file->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

    /* Destroy buffers allocated during chunking phase */
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

    hz_end[vp] = PIDX_get_time();
    /*----------------------------------------------HZ [end]-------------------------------------------------*/


    /*----------------------------------------------Agg [staart]-----------------------------------------------*/


    /* Creating the buffers required for Aggregation */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
      for (j = 0 ; j < agg_io_level; j++)
      {
        /* Creating the buffers required for Aggregation */
        if (file->agg_type == 0)
        {
          if (no_of_aggregators * file->idx->variable_count <= nprocs)
            ret = PIDX_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], file->idx->variable[file->local_variable_index]->global_block_layout, i, j);
          else
            ret = PIDX_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], file->idx->variable[file->local_variable_index]->global_block_layout, 0, j);
        }

        else if (file->agg_type == 1)
        {
          if (no_of_aggregators * file->idx->variable_count <= nprocs)
            ret = PIDX_local_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], j/*0*/);
          else
            ret = PIDX_local_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], 0);
        }

        else
          ret = PIDX_local_agg_local_comm_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], j/*0*/);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_agg_buf = 1;
    MPI_Allreduce(&l_agg_buf, &g_agg_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg_buf == nprocs)
      printf("[A] Aggregation Buffer Created\n");
#endif



    /* Perform Aggregation */
    if (file->debug_do_agg == 1)
    {
      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
        for(j = 0 ; j < agg_io_level; j++)
        {
           agg_start[static_var_counter][j] = PIDX_get_time();

           if (file->agg_type == 0)
             ret = PIDX_agg_write(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], j, file->idx->variable[start_var_index]->block_layout_by_level[j], PIDX_WRITE);

           else if (file->agg_type == 1)
             ret = PIDX_local_agg(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], j, file->idx->variable[start_var_index]->block_layout_by_level[j], PIDX_WRITE);

           else
             ret = PIDX_local_agg_local_comm(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], j, file->idx->variable[start_var_index]->block_layout_by_level[j], PIDX_WRITE);

           if (ret != PIDX_success)
           {
             fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
             return PIDX_err_rst;
           }

           agg_end[static_var_counter][j] = PIDX_get_time();
        }
        static_var_counter++;
      }
    }


    if (file->debug_do_io == 1)
    {
      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = agg_io_level; j < file->layout_count; j++)
        {
          io_per_process_start[static_var_counter][j] = PIDX_get_time();

          ret = PIDX_io_per_process(file->tio_id[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], PIDX_WRITE);
          if (ret != PIDX_success)
            return PIDX_err_io;

          io_per_process_end[static_var_counter][j] = PIDX_get_time();
        }
        static_var_counter++;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_agg = 1;
    MPI_Allreduce(&l_agg, &g_agg, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg == nprocs)
      printf("[A] Aggregation Completed\n");
#endif

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    /*--------------------------------------------Agg [end]--------------------------------------------------*/




    /*--------------------------------------------IO [start]--------------------------------------------------*/


    /* Initialization ONLY ONCE for all TIME STEPS (caching across time) */
    /*
    if (file->idx->enable_agg == 2)
    {
      if (caching_state == 1 && file->idx_d->agg_buffer->var_number == 0 && file->idx_d->agg_buffer->sample_number == 0)
      {
        if (file->flush_used == 0)
        {
          PIDX_cache_headers(file);
          caching_state = 0;
        }
      }
    }
    */

    if (file->debug_do_io == 1)
    {
      if (time_step_caching == 1)
      {
        if (file->flush_used == 0)
          PIDX_io_cached_data(cached_header_copy);
      }

      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
        for(j = 0 ; j < agg_io_level; j++)
        {
          io_start[static_var_counter][j] = PIDX_get_time();

          ret = PIDX_aggregated_io(file->tio_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], PIDX_WRITE);
          if (ret != PIDX_success)
            return PIDX_err_io;

          io_end[static_var_counter][j] = PIDX_get_time();

        }
        static_var_counter++;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_io = 1;
    MPI_Allreduce(&l_io, &g_io, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_io == nprocs)
      printf("[I] I/O completed\n");
#endif

    /* Destroy buffers allocated during aggregation phase */
    static_var_counter = 0;
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
      for(j = 0 ; j < agg_io_level; j++)
      {
        ret = PIDX_agg_buf_destroy(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->agg_type);
        if (ret != PIDX_success)
          return PIDX_err_agg;
      }
      static_var_counter++;
    }


    /*----------------------------------------------IO [end]--------------------------------------------------*/


    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      for(j = 0 ; j < agg_io_level; j++)
      {
        free(file->idx_d->agg_buffer[i][j]);
        file->idx_d->agg_buffer[i][j] = 0;
      }
      free(file->idx_d->agg_buffer[i]);
      file->idx_d->agg_buffer[i] = 0;
    }
    free(file->idx_d->agg_buffer);\
    file->idx_d->agg_buffer = 0;


    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
      for (j = 0 ; j < agg_io_level; j++)
      {
        if (file->agg_type == 0)
          ret = PIDX_agg_meta_data_destroy(file->tagg_id[i][j]);

        else if (file->agg_type == 1)
          ret = PIDX_local_agg_meta_data_destroy(file->tagg_id[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);

        else
          ret = PIDX_local_agg_local_comm_meta_data_destroy(file->tagg_id[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }
    }
#if 1
    ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /*-------------------------------------------finalize [start]---------------------------------------------*/
    finalize_start[vp] = PIDX_get_time();


    //if (file->small_agg_comm == 1)
    //  PIDX_destroy_local_aggregation_comm(file->agg_id);

    /* Deleting the I/O ID */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      for(j = 0 ; j < file->layout_count; j++)
        PIDX_io_finalize(file->tio_id[i][j]);

    /* Deleting the aggregation ID */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      for(j = 0 ; j < file->layout_count; j++)
        PIDX_agg_finalize(file->tagg_id[i][j]);

    for(i = 0 ; i < file->idx->variable_count ; i++)
    {
      free(file->tagg_id[i]);
      file->tagg_id[i] = 0;

      free(file->tio_id[i]);
      file->tio_id[i] = 0;
    }
    free(file->tagg_id);
    file->tagg_id = 0;

    free(file->tio_id);
    file->tio_id = 0;

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);

    /* Deleting the compression ID */
    PIDX_compression_finalize(file->comp_id);

    /* Deleting the chunking ID */
    PIDX_chunk_finalize(file->chunk_id);

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    finalize_end[vp] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/
#endif

    //if (rank == 0)
    //  printf("Finished Writing %d variables\n", end_index - start_index + 1);
    vp++;
  }

  delete_idx_dataset(file);


#if PIDX_DEBUG_OUTPUT
  l_pidx = 1;
  MPI_Allreduce(&l_pidx, &g_pidx, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_pidx == nprocs)
    printf("PIDX closing file\n");

#endif

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  return PIDX_success;
}

static PIDX_return_code PIDX_read(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  file->var_pipe_length = file->idx->variable_count - 1;
  if (file->var_pipe_length == 0)
    file->var_pipe_length = 1;

#if PIDX_DEBUG_OUTPUT
  unsigned long long l_populate = 0, g_populate = 0;
  unsigned long long l_filec = 0, g_filec = 0;
  unsigned long long l_init = 0, g_init = 0;
  unsigned long long l_rst_buf = 0, g_rst_buf = 0;
  unsigned long long l_chunk_buf = 0, g_chunk_buf = 0;
  unsigned long long l_rst = 0, g_rst = 0;
  unsigned long long l_chunk = 0, g_chunk = 0;
  unsigned long long l_cmp = 0, g_cmp = 0;
  unsigned long long l_hz_buf = 0, g_hz_buf = 0;
  unsigned long long l_hz = 0, g_hz = 0;
  unsigned long long l_agg_buf = 0, g_agg_buf = 0;
  unsigned long long l_agg = 0, g_agg = 0;
  unsigned long long l_io = 0, g_io = 0;
  unsigned long long l_pidx = 0, g_pidx = 0;
#endif

  int j = 0, p, var = 0, d = 0;
  int total_header_size;
  PIDX_return_code ret;
  //static int header_io = 0;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm,  &nprocs);
  }
#endif

  populate_idx_start_time = MPI_Wtime();

  ret = populate_idx_dataset(file);
  if (ret != PIDX_success)
    return PIDX_err_file;

  PIDX_init_timming_buffers2(file);

#if PIDX_DEBUG_OUTPUT
  l_populate = 1;
  MPI_Allreduce(&l_populate, &g_populate, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_populate == nprocs)
    printf("Finished Populating IDX dataset (File Count %d maxh = %d)\n", file->idx_d->max_file_count, file->idx_d->maxh);
#endif

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_file;

    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;

    ret = PIDX_parameter_validate(file, start_var_index, end_var_index);
    if (ret != PIDX_success)
      return PIDX_err_file;

    file->one_time_initializations = 1;
  }

  if (file->idx_count[0] != 1 || file->idx_count[1] != 1 || file->idx_count[2] != 1 )
    for (var = start_var_index; var < end_var_index; var++)
      for (p = 0; p < file->idx->variable[var]->sim_patch_count; p++)
        for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
          for (j = 0; j < file->idx->bounds[d] * file->idx_count[d]; j = j + (file->idx->bounds[d]))
            if (file->idx->variable[var]->sim_patch[p]->offset[d] >= j && file->idx->variable[var]->sim_patch[p]->offset[d] < (j + file->idx->bounds[d]))
            {
              file->idx->variable[var]->sim_patch[p]->offset[d] = file->idx->variable[var]->sim_patch[p]->offset[d] - j;
              break;
            }

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI  
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS * sizeof(long long));
  }
#else
  memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
  memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS * sizeof(long long));
#endif

#if PIDX_DEBUG_OUTPUT
  l_filec = 1;
  MPI_Allreduce(&l_filec, &g_filec, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_filec == nprocs)
    printf("Finished Creating File heirarchy\n");
#endif

  populate_idx_end_time = MPI_Wtime();


  int start_index = 0, end_index = 0;
  int i = 0;
  //int agg_by_level_file = 2;//(file->idx_d->maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));

  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);

    int agg_io_level = 0, no_of_aggregators = 0;
    //file->idx->variable[file->local_variable_index]->block_layout_by_level[j];

    for (i = 0; i < file->layout_count ; i++)
    {
      no_of_aggregators = file->idx->variable[file->local_variable_index]->block_layout_by_level[i]->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level = i;
    }
    agg_io_level = agg_io_level + 1;

    if (file->idx->enable_agg == 0)
      agg_io_level = 0;

    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
    /* Create the restructuring ID */
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the chunking ID */
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the compression ID */
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the HZ encoding ID */
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the aggregation ID */
    /* Create the I/O ID */

    file->tagg_id = malloc(sizeof(*(file->tagg_id)) * file->idx->variable_count);
    memset(file->tagg_id, 0, sizeof(*(file->tagg_id)) * file->idx->variable_count);

    file->tio_id = malloc(sizeof(*(file->tio_id)) * file->idx->variable_count);
    memset(file->tio_id, 0, sizeof(*(file->tio_id)) * file->idx->variable_count);

    //file->tagg_id = malloc(sizeof(*(file->tagg_id)) * (end_index - start_index + 1));
    //memset(file->tagg_id, 0, sizeof(*(file->tagg_id)) * (end_index - start_index + 1));

    //file->tio_id = malloc(sizeof(*(file->tio_id)) * (end_index - start_index + 1));
    //memset(file->tio_id, 0, sizeof(*(file->tio_id)) * (end_index - start_index + 1));

    //printf("COUNT [%d %d] : %d\n", start_index, end_index, (end_index - start_index + 1));
    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      file->tagg_id[i] = malloc(sizeof(*(file->tagg_id[i])) * file->layout_count);
      memset(file->tagg_id[i], 0, sizeof(*(file->tagg_id[i])) * file->layout_count);

      file->tio_id[i] = malloc(sizeof(*(file->tio_id[i])) * file->layout_count);
      memset(file->tio_id[i], 0, sizeof(*(file->tio_id[i])) * file->layout_count);
    }

    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for(j = 0 ; j < file->layout_count; j++)
      {
        file->tagg_id[i][j] = PIDX_agg_init(file->idx, file->idx_d, start_var_index, i, i);
        file->tio_id[i][j] = PIDX_io_init(file->idx, file->idx_d, start_var_index, i, i);
      }
    }

    //file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * (end_index - start_index + 1));
    //memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * (end_index - start_index + 1));

    file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
    memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      file->idx_d->agg_buffer[i] = malloc(sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      memset(file->idx_d->agg_buffer[i], 0, sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      for(j = 0 ; j < agg_io_level; j++)
      {
        file->idx_d->agg_buffer[i][j] = malloc(sizeof(*(file->idx_d->agg_buffer[i][j])) );
        memset(file->idx_d->agg_buffer[i][j], 0, sizeof(*(file->idx_d->agg_buffer[i][j])) );

        file->idx_d->agg_buffer[i][j]->file_number = -1;
        file->idx_d->agg_buffer[i][j]->var_number = -1;
        file->idx_d->agg_buffer[i][j]->sample_number = -1;

        file->idx_d->agg_buffer[i][j]->no_of_aggregators = 0;
        file->idx_d->agg_buffer[i][j]->aggregator_interval = 0;
        file->idx_d->agg_buffer[i][j]->aggregation_factor = 1;
      }

      /*
      int start_agg_index;
      if (agg_io_level == 1)
         start_agg_index = 1;
      if (agg_io_level == 2)
        start_agg_index = 2;


      file->idx_d->agg_buffer[i][0]->aggregation_factor = 1;//(int)pow(2, (agg_io_level - 1));
      for (j = 1 ; j < agg_io_level; j++)
      {
        file->idx_d->agg_buffer[i][j]->aggregation_factor = 1;//(int)pow(2, (agg_io_level - j));
        //if (rank == 0)
        //  printf("[%d %d] factor at layout %d = %d\n", agg_io_level, file->layout_count, j, file->idx_d->agg_buffer[i][j]->aggregation_factor);
      }
      */
    }
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    /* Attaching the communicator to the restructuring phase */
    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Attaching the communicator to the chunking phase */
    ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

    /* Attaching the communicator to the compression phase */
    ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_compress;

    /* Attaching the communicator to the HZ encodig phase phase */
    ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    /* Attaching the communicator to the aggregation phase */
    /* Attaching the communicator to the I/O phase */
    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for(j = 0 ; j < file->layout_count; j++)
      {
        ret = PIDX_agg_set_communicator(file->tagg_id[i][j], file->comm);
        if (ret != PIDX_success)
          return PIDX_err_agg;

        ret = PIDX_io_set_communicator(file->tio_id[i][j], file->comm);
        if (ret != PIDX_success)
          return PIDX_err_io;
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/


    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    //if (file->small_agg_comm == 1)
    //  PIDX_create_local_aggregation_comm(file->agg_id);

    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for(j = 0 ; j < agg_io_level; j++)
      {
        //file->local_variable_index
        //ret = PIDX_agg_meta_data_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], file->idx->variable[file->local_variable_index]->global_block_layout);
        ret = PIDX_local_agg_meta_data_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);
        if (ret != PIDX_success)
          return PIDX_err_rst;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_init = 1;
    MPI_Allreduce(&l_init, &g_init, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_init == nprocs)
      printf("All Modules Initialized (Variable index [%d %d] Variable Pipe length %d)\n", start_index, end_index, file->var_pipe_length);
#endif


    /*------------------------------------------Agg and IO [staart]-------------------------------------------*/


    /* Creating the buffers required for Aggregation */
    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for (j = 0 ; j < agg_io_level; j++)
      {
        /* Creating the buffers required for Aggregation */
        //ret = PIDX_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], file->idx->variable[file->local_variable_index]->global_block_layout, i, j);
        ret = PIDX_local_agg_buf_create(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], /*j*/0);
        if (ret != PIDX_success)
          return PIDX_err_agg;
      }
    }

    if (file->debug_do_io == 1)
    {
      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i++)
      {
        for(j = 0 ; j < agg_io_level; j++)
        {
          io_start[static_var_counter][j] = PIDX_get_time();

          ret = PIDX_aggregated_io(file->tio_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], PIDX_READ);
          if (ret != PIDX_success)
            return PIDX_err_io;

          io_end[static_var_counter][j] = PIDX_get_time();
        }
        static_var_counter++;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_agg_buf = 1;
    MPI_Allreduce(&l_agg_buf, &g_agg_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg_buf == nprocs)
      printf("[A] Aggregation Buffer Created\n");
#endif

    /* Creating the buffers required for HZ encoding */
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

#if PIDX_DEBUG_OUTPUT
    l_hz_buf = 1;
    MPI_Allreduce(&l_hz_buf, &g_hz_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz_buf == nprocs)
      printf("[H] HZ Buffer Created\n");
#endif

    /* Perform Aggregation */
    if (file->debug_do_agg == 1)
    {
      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i++)
      {
        for(j = 0 ; j < agg_io_level; j++)
        {
          agg_start[static_var_counter][j] = PIDX_get_time();

          ret = PIDX_local_agg(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], j, file->idx->variable[start_var_index]->block_layout_by_level[j], PIDX_READ);
          //ret = PIDX_agg_write(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], file->idx->variable[start_var_index]->block_layout_by_level[j]);
          if (ret != PIDX_success)
            return PIDX_err_agg;

          agg_end[static_var_counter][j] = PIDX_get_time();
        }
        static_var_counter++;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_io = 1;
    MPI_Allreduce(&l_io, &g_io, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_io == nprocs)
      printf("[I] I/O completed\n");
#endif

    /* Destroy buffers allocated during aggregation phase */
    static_var_counter = 0;
    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for(j = 0 ; j < agg_io_level; j++)
      {
        ret = PIDX_agg_buf_destroy(file->tagg_id[i][j], file->idx_d->agg_buffer[i][j], 0);
        if (ret != PIDX_success)
          return PIDX_err_agg;
      }
      static_var_counter++;
    }

    if (file->debug_do_io == 1)
    {
      static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i++)
      {
        for(j = agg_io_level; j < file->layout_count; j++)
        {
          io_per_process_start[static_var_counter][j] = PIDX_get_time();

          ret = PIDX_io_per_process(file->tio_id[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j], PIDX_READ);
          if (ret != PIDX_success)
            return PIDX_err_io;

          io_per_process_end[static_var_counter][j] = PIDX_get_time();
        }
        static_var_counter++;
      }
    }

    /*-------------------------------------------Agg and IO [end]---------------------------------------------*/


    /*---------------------------------------------HZ [start]------------------------------------------------*/
    hz_start[vp] = PIDX_get_time();
#if PIDX_DEBUG_OUTPUT
    l_agg = 1;
    MPI_Allreduce(&l_agg, &g_agg, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg == nprocs)
      printf("[A] Aggregation Completed\n");
#endif

    /* Verify the HZ encoding */
    if (file->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }


    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

#if PIDX_DEBUG_OUTPUT
    l_chunk_buf = 1;
    MPI_Allreduce(&l_chunk_buf, &g_chunk_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk_buf == nprocs)
      printf("[C] Chunking Buffer Created\n");
#endif




    /* Perform HZ encoding */
    if (file->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_read(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

#if PIDX_DEBUG_OUTPUT
    l_hz = 1;
    MPI_Allreduce(&l_hz, &g_hz, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz == nprocs)
      printf("[H] HZ Encoding Finished\n");
#endif


    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    hz_end[vp] = PIDX_get_time();
    /*----------------------------------------------HZ [end]-------------------------------------------------*/



    /*----------------------------------------Compression [start]--------------------------------------------*/
    compression_start[vp] = PIDX_get_time();

    /* Perform Compression */
    if (file->debug_do_compress == 1)
    {
#if !SIMULATE_IO
      ret = PIDX_decompression(file->comp_id);
      if (ret != PIDX_success)
        return PIDX_err_compress;
#endif
    }

#if PIDX_DEBUG_OUTPUT
    l_cmp = 1;
    MPI_Allreduce(&l_cmp, &g_cmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_cmp == nprocs)
      printf("[CMP] Compression Completed\n");
#endif

    compression_end[vp] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/




    /*----------------------------------------Chunking [start]------------------------------------------------*/
    chunk_start[vp] = PIDX_get_time();

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

#if PIDX_DEBUG_OUTPUT
    l_rst_buf = 1;
    MPI_Allreduce(&l_rst_buf, &g_rst_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_rst_buf == nprocs)
      printf("[R] Restructuring Buffer Created\n");
#endif

    /* Perform Chunking */
    if (file->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_READ);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
//#if 0
#if PIDX_DEBUG_OUTPUT
    l_chunk = 1;
    MPI_Allreduce(&l_chunk, &g_chunk, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk == nprocs)
      printf("[C] Chunking Completed\n");
#endif

    /* Destroy buffers allocated during chunking phase */
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

    chunk_end[vp] = PIDX_get_time();
    /*-----------------------------------------Chunking [end]------------------------------------------------*/




    /*--------------------------------------------RST [start]------------------------------------------------*/
    rst_start[vp] = PIDX_get_time();

    /* Perform data restructuring */
    if (file->debug_do_rst == 1)
    {
      ret = PIDX_rst_read(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

#if PIDX_DEBUG_OUTPUT
    l_rst = 1;
    MPI_Allreduce(&l_rst, &g_rst, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_rst == nprocs)
      printf("[R] Restructuring Completed\n");
#endif

    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_rst_buf_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    rst_end[vp] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/


    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for(j = 0 ; j < agg_io_level; j++)
      {
        free(file->idx_d->agg_buffer[i][j]);
        file->idx_d->agg_buffer[i][j] = 0;
      }
      free(file->idx_d->agg_buffer[i]);
      file->idx_d->agg_buffer[i] = 0;
    }
    free(file->idx_d->agg_buffer);\
    file->idx_d->agg_buffer = 0;


    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      for (j = 0 ; j < agg_io_level; j++)
      {
        //ret = PIDX_agg_meta_data_destroy(file->tagg_id[i][j]);
        ret = PIDX_local_agg_meta_data_destroy(file->tagg_id[i][j], file->idx->variable[file->local_variable_index]->block_layout_by_level[j]);
        if (ret != PIDX_success)
          return PIDX_err_rst;
      }
    }

    ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /*-------------------------------------------finalize [start]---------------------------------------------*/
    finalize_start[vp] = PIDX_get_time();

    //if (file->small_agg_comm == 1)
    //  PIDX_destroy_local_aggregation_comm(file->agg_id);

    /* Deleting the I/O ID */
    for(i = start_index ; i < (end_index + 1) ; i++)
      for(j = 0 ; j < file->layout_count; j++)
        PIDX_io_finalize(file->tio_id[i][j]);

    /* Deleting the aggregation ID */
    for(i = start_index ; i < (end_index + 1) ; i++)
      for(j = 0 ; j < file->layout_count; j++)
        PIDX_agg_finalize(file->tagg_id[i][j]);

    for(i = start_index ; i < (end_index + 1) ; i++)
    {
      free(file->tagg_id[i]);
      file->tagg_id[i] = 0;

      free(file->tio_id[i]);
      file->tio_id[i] = 0;
    }
    free(file->tagg_id);
    file->tagg_id = 0;

    free(file->tio_id);
    file->tio_id = 0;

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);

    /* Deleting the compression ID */
    PIDX_compression_finalize(file->comp_id);

    /* Deleting the chunking ID */
    PIDX_chunk_finalize(file->chunk_id);

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);

    finalize_end[vp] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/


    //if (rank == 0)
    //  printf("Finished Writing %d variables\n", end_index - start_index + 1);
    vp++;
  }

  delete_idx_dataset(file);

#if PIDX_DEBUG_OUTPUT
  l_pidx = 1;
  MPI_Allreduce(&l_pidx, &g_pidx, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_pidx == nprocs)
    printf("PIDX closing file\n");

#endif

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  return PIDX_success;

}


PIDX_return_code PIDX_flush(PIDX_file file)
{
  int i, p;
  int ret;
  if (file->idx->variable_count <= 0)
    return PIDX_err_variable;

  if (file->flags == PIDX_MODE_CREATE)
  {
    if (file->enable_raw_dump == 0)
    {
      ret = PIDX_write(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
    else
    {
      ret = PIDX_raw_write(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
  }

  else if (file->flags == PIDX_MODE_RDONLY)
  {
    if (file->enable_raw_dump == 0)
    {
      ret = PIDX_read(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
    else
    {
      ret = PIDX_raw_read(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
  }

  else if (file->flags == PIDX_MODE_RDWR)
  {
    int state = file->idx->variable[file->local_variable_index]->io_state;
    int state_index = file->local_variable_index;
    int new_state, same_state_count = 0;

    for (i = file->local_variable_index; i < file->local_variable_index + file->local_variable_count; i++)
    {
      new_state = file->idx->variable[i]->io_state;
      if (state == new_state)
      {
        same_state_count++;
        if (i == file->local_variable_index + file->local_variable_count - 1)
        {
          if (state == 1)
            PIDX_write(file, state_index, state_index + same_state_count);
          else if (state == 0)
            PIDX_read(file, state_index, state_index + same_state_count);
        }
      }
      else
      {
        if (state == 1)
          PIDX_write(file, state_index, state_index + same_state_count);
        else if (state == 0)
          PIDX_read(file, state_index, state_index + same_state_count);

        state = new_state;
        state_index = i;
        same_state_count = 1;
      }
    }
  }

  for (i = file->local_variable_index; i < file->local_variable_index + file->local_variable_count; i++)
  {
    for(p = 0; p < file->idx->variable[i]->sim_patch_count; p++)
    {
      free(file->idx->variable[i]->sim_patch[p]);
      file->idx->variable[i]->sim_patch[p] = 0;
    }
  }

  file->local_variable_index = file->idx->variable_index_tracker;
  file->local_variable_count = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_close(PIDX_file file) 
{
  int ret;
  int i = 0;
  file->write_on_close = 1;

  ret = PIDX_flush(file);
  if (ret != PIDX_success)
    return PIDX_err_close;
  
  sim_end = PIDX_get_time();

  if (file->debug_output == 1)
  {

    double total_time = sim_end - sim_start;
    double max_time = total_time;
    int sample_sum = 0, var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode == 1)
    {
      MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
      MPI_Comm_rank(file->comm, &rank);
      MPI_Comm_size(file->comm, &nprocs);
    }
    else
      total_time = max_time;
#endif

    if (file->idx_count[0] *  file->idx_count[1] * file->idx_count[2] == 1)
    {
      if (max_time == total_time)
      {
        for (var = 0; var < file->idx->variable_count; var++)
          sample_sum = sample_sum + file->idx->variable[var]->values_per_sample;
      
        int64_t total_data = file->idx->bounds[0] * file->idx->bounds[1] * file->idx->bounds[2] * file->idx->bounds[3] * file->idx->bounds[4] * sample_sum * 8;
        fprintf(stdout, "\n==========================================================================================================\n");
        fprintf(stdout, "[%d] Time step %d File name %s\n", rank, file->idx->current_time_step, file->idx->filename);
        fprintf(stdout, "Cores %d Global Data %lld %lld %lld Variables %d IDX count %d = %d x %d x %d\n", nprocs, (long long) file->idx->bounds[0], (long long) file->idx->bounds[1], (long long) file->idx->bounds[2], file->idx->variable_count, file->idx_count[0] * file->idx_count[1] * file->idx_count[2], file->idx_count[0], file->idx_count[1], file->idx_count[2]);
        fprintf(stdout, "Rst = %d Comp = %d\n", file->idx->enable_rst, file->idx->compression_type);//, file->idx->enable_agg, file->small_agg_comm);
        fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);//, file->idx->variable[0]->existing_file_count);
        fprintf(stdout, "Chunk Size %d %d %d %d %d\n", (int)file->idx->chunk_size[0], (int)file->idx->chunk_size[1], (int)file->idx->chunk_size[2], (int)file->idx->chunk_size[3], (int)file->idx->chunk_size[4]);
        fprintf(stdout, "Restructuring Box Size %d %d %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2], (int)file->idx->reg_patch_size[3], (int)file->idx->reg_patch_size[4]);
        fprintf(stdout, "Aggregation Type = %d\n", file->agg_type);
        fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", max_time, (float) total_data / (1000 * 1000 * max_time));
        fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
        printf("Block layout creation time %f\n", populate_idx_end_time - populate_idx_start_time);
        fprintf(stdout, "File Create Time: %f Seconds\n", (file_create_time - sim_start));

        double header_io_time = 0;
        for (var = 0; var < hp; var++)
        {
          header_io_time = header_io_time + (write_init_end[var] - write_init_start[var]);
          fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
        }
        double total_time_ai = 0, total_time_a = 0, total_time_i = 0, total_time_pi = 0;
        int p = 0;
        for (var = 0; var < /*file->idx->variable_count*/static_var_counter; var++)
        {
          for (p = 0; p < file->layout_count; p++)
          {
            fprintf(stdout, "[%d %d] Agg time + AGG I/O time + Per-Process I/O time = %f + %f + %f = %f\n", var, p, (agg_end[var][p] - agg_start[var][p]), (io_end[var][p] - io_start[var][p]), (io_per_process_end[var][p] - io_per_process_start[var][p]), (agg_end[var][p] - agg_start[var][p]) + (io_end[var][p] - io_start[var][p]) + (io_per_process_end[var][p] - io_per_process_start[var][p]));

            total_time_a = total_time_a + (agg_end[var][p] - agg_start[var][p]);
            total_time_i = total_time_i + (io_end[var][p] - io_start[var][p]);
            total_time_pi = total_time_pi + (io_per_process_end[var][p] - io_per_process_start[var][p]);
          }
        }
        total_time_ai = total_time_a + total_time_i + total_time_pi;
        fprintf(stdout, "Agg time + AGG I/O time + Per-Process I/O time : %f + %f + %f = %f\n", total_time_a, total_time_i, total_time_pi, total_time_ai);

        int timer_count = 0;
        timer_count = file->idx->variable_count / (file->var_pipe_length + 1);
        if (file->idx->variable_count % (file->var_pipe_length + 1) != 0)
          timer_count = timer_count + 1;


       //printf("static_var_counter = %d\n", static_var_counter);
       double total_time_rch = 0;
        for (var = 0; var < timer_count; var++)
        {
          fprintf(stdout, "[%d] STARTUP + RST + BRST + HZ = %f + %f + %f + %f = %f\n", var, (startup_end[var] - startup_start[var]), (rst_end[var] - rst_start[var]), (chunk_end[var] - chunk_start[var]), (hz_end[var] - hz_start[var]), (startup_end[var] - startup_start[var]) + (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]));
          total_time_rch = total_time_rch + (startup_end[var] - startup_start[var]) + (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]);
      }

        fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f + %f] [%f]\n", total_time_ai + total_time_rch + (file_create_time - sim_start) + (populate_idx_end_time - populate_idx_start_time) + header_io_time, (populate_idx_end_time - populate_idx_start_time), (file_create_time - sim_start), header_io_time, total_time_rch, total_time_ai, max_time);

        fprintf(stdout, "==========================================================================================================\n");
      }
    }
    else
    {
      for (var = 0; var < file->idx->variable_count; var++)
        sample_sum = sample_sum + file->idx->variable[var]->values_per_sample;
    
      int64_t total_data = file->idx->bounds[0] * file->idx->bounds[1] * file->idx->bounds[2] * file->idx->bounds[3] * file->idx->bounds[4] * sample_sum * 8 * file->idx_count[0] * file->idx_count[1] * file->idx_count[2];
      if (max_time == total_time)
        fprintf(stdout, "[EXTRA INFO] %d %s Time %f Seconds Throughput %f MB/sec\n", file->idx->current_time_step, file->idx->filename, max_time, (float) total_data / (1000 * 1000 * max_time));
    
      int global_rank = 0, global_nprocs = 1;
      double global_max_time = 0;
    
#if PIDX_HAVE_MPI
      if (file->idx_d->parallel_mode == 1)
      {
        MPI_Allreduce(&max_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, file->global_comm);
        MPI_Comm_rank(file->global_comm, &global_rank);
        MPI_Comm_size(file->global_comm, &global_nprocs);
      }
      else
        max_time = global_max_time;
#endif
    
      if (global_max_time == total_time)
      {
        fprintf(stdout, "\n==========================================================================================================\n");
        fprintf(stdout, "[%d] Combined Time step %d File name %s\n", global_rank, file->idx->current_time_step, file->idx->filename);
        fprintf(stdout, "Cores %d Global Data %lld %lld %lld Variables %d IDX Count %d = %d x %d x %d\n", global_nprocs, (long long) file->idx->bounds[0] * file->idx_count[0], (long long) file->idx->bounds[1] * file->idx_count[1], (long long) file->idx->bounds[2] * file->idx_count[2], file->idx->variable_count, file->idx_count[0] * file->idx_count[1] * file->idx_count[2], file->idx_count[0], file->idx_count[1], file->idx_count[2]);
        fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);//, file->idx->variable[0]->existing_file_count);
        fprintf(stdout, "Chunk Size %d %d %d %d %d\n", (int)file->idx->chunk_size[0], (int)file->idx->chunk_size[1], (int)file->idx->chunk_size[2], (int)file->idx->chunk_size[3], (int)file->idx->chunk_size[4]);
        fprintf(stdout, "Restructuring Box Size %d %d %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2], (int)file->idx->reg_patch_size[3], (int)file->idx->reg_patch_size[4]);
        fprintf(stdout, "Aggregation Type = %d\n", file->agg_type);
        fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", max_time, (float) total_data / (1000 * 1000 * max_time));
        fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
        printf("Block layout creation time %f\n", populate_idx_end_time - populate_idx_start_time);
        fprintf(stdout, "File Create Time: %f Seconds\n", (file_create_time - sim_start));

        double header_io_time = 0;
        for (var = 0; var < hp; var++)
        {
          header_io_time = header_io_time + (write_init_end[var] - write_init_start[var]);
          fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
        }
        double total_time_ai = 0, total_time_a = 0, total_time_i = 0, total_time_pi = 0;
        int p = 0;
        for (var = 0; var < file->idx->variable_count; var++)
        {
          for (p = 0; p < file->layout_count; p++)
          {
            fprintf(stdout, "[%d %d] Agg time + AGG I/O time + Per-Process I/O time = %f + %f + %f = %f\n", var, p, (agg_end[var][p] - agg_start[var][p]), (io_end[var][p] - io_start[var][p]), (io_per_process_end[var][p] - io_per_process_start[var][p]), (agg_end[var][p] - agg_start[var][p]) + (io_end[var][p] - io_start[var][p]) + (io_per_process_end[var][p] - io_per_process_start[var][p]));

            total_time_a = total_time_a + (agg_end[var][p] - agg_start[var][p]);
            total_time_i = total_time_i + (io_end[var][p] - io_start[var][p]);
            total_time_pi = total_time_pi + (io_per_process_end[var][p] - io_per_process_start[var][p]);
          }
        }
        total_time_ai = total_time_a + total_time_i + total_time_pi;
        fprintf(stdout, "Agg time + AGG I/O time + Per-Process I/O time : %f + %f + %f = %f\n", total_time_a, total_time_i, total_time_pi, total_time_ai);

        int timer_count = 0;
        timer_count = file->idx->variable_count / (file->var_pipe_length + 1);
        if (file->idx->variable_count % (file->var_pipe_length + 1) != 0)
          timer_count = timer_count + 1;

        double total_time_rch = 0;
        for (var = 0; var < timer_count; var++)
        {
          fprintf(stdout, "[%d] STARTUP + RST + BRST + HZ = %f + %f + %f + %f = %f\n", var, (startup_end[var] - startup_start[var]), (rst_end[var] - rst_start[var]), (chunk_end[var] - chunk_start[var]), (hz_end[var] - hz_start[var]), (startup_end[var] - startup_start[var]) + (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]));
          total_time_rch = total_time_rch + (startup_end[var] - startup_start[var]) + (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]);
        }

        fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f + %f] [%f]\n", total_time_ai + total_time_rch + (file_create_time - sim_start) + (populate_idx_end_time - populate_idx_start_time) + header_io_time, (populate_idx_end_time - populate_idx_start_time), (file_create_time - sim_start), header_io_time, total_time_rch, total_time_ai, max_time);

        fprintf(stdout, "==========================================================================================================\n");
      }
    
    }
  }


  vp = 0;
  hp = 0;

  for (i = 0; i < file->idx->variable_count; i++)
  {
    free(agg_start[i]);
    free(agg_end[i]);

    free(io_start[i]);
    free(io_end[i]);

    free(io_per_process_start[i]);
    free(io_per_process_end[i]);
  }


  for (i = 0; i < 1024; i++)
  {
    free(file->idx->variable[i]);
    file->idx->variable[i] = 0;
  }

  file->idx->variable_count = 0;

  //free(file->idx->bounds);        file->idx->bounds = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_free(&(file->comm));
    if (file->idx_count[0] * file->idx_count[1] * file->idx_count[2] != 1)
      MPI_Comm_free(&(file->global_comm));
  }
#endif

  free(file->idx);                  file->idx = 0;
  free(file->idx_d);                file->idx_d = 0;


  free(file);
  
  free(init_start);                     init_start        = 0;
  free(init_end);                       init_end          = 0;
  free(write_init_start);               write_init_start        = 0;
  free(write_init_end);                 write_init_end          = 0;
  free(rst_start);                      rst_start               = 0;
  free(rst_end);                        rst_end                 = 0;
  free(startup_start);                  startup_start           = 0;
  free(startup_end);                    startup_end             = 0;
  free(hz_start);                       hz_start                = 0;
  free(hz_end);                         hz_end                  = 0;
  free(agg_start);                      agg_start               = 0;
  free(agg_end);                        agg_end                 = 0;
  free(io_start);                       io_start                = 0;
  free(io_end);                         io_end                  = 0;
  free(cleanup_start);                  cleanup_start           = 0;
  free(cleanup_end);                    cleanup_end             = 0;
  free(finalize_start);                 finalize_start          = 0;
  free(finalize_end);                   finalize_end            = 0;
  free(buffer_start);                   buffer_start            = 0;
  free(buffer_end);                     buffer_end              = 0;
  free(chunk_start);                chunk_start         = 0;
  free(chunk_end);                  chunk_end           = 0;
  free(compression_start);              compression_start       = 0;
  free(compression_end);                compression_end         = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable)
{
  if(!variable)
    return PIDX_err_variable;
  
  variable->dump_meta_data_ON = 1;
  
  return PIDX_success;
}


PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable)
{
  if(!variable)
    return PIDX_err_variable;
  
  variable->dump_meta_data_ON = 0;
  
  return PIDX_success;
}


PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_get_bits_per_sample(PIDX_data_type type_name, unsigned int bits_per_sample)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_get_patch_count(PIDX_file file, int* patch_count)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point dims)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_get_patch_count_with_rank(PIDX_file file, int MPI_rank, int* patch_count)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, int MPI_rank, PIDX_point offset, PIDX_point dims)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int variable_index)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_index < 0)
    return PIDX_err_size;
  
  if(file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  file->idx->variable_index_tracker = variable_index;
  file->local_variable_count = 1;
  file->local_variable_index = variable_index;

  return PIDX_success;
}



PIDX_return_code PIDX_set_variable_pile_length(PIDX_file file, int var_pipe_length)
{
  if(!file)
    return PIDX_err_file;

  if (var_pipe_length < 0)
    return PIDX_err_size;

  file->var_pipe_length = var_pipe_length;

  return PIDX_success;
}


PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  if(file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  (*variable) = file->idx->variable[file->idx->variable_index_tracker];
  
  return PIDX_success;
}



PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable)
{
  if(!file)
    return PIDX_err_file;

  if(file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;

  file->idx->variable[file->idx->variable_index_tracker] = variable;

  return PIDX_success;
}



PIDX_return_code PIDX_activate_local_aggregation(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->small_agg_comm = 1;

  return PIDX_success;
}



PIDX_return_code PIDX_read_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_write_variable(PIDX_file file, PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if (!file)
    return PIDX_err_file;

  if(!variable)
    return PIDX_err_variable;

  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;

  if(!dims || dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_count;

  if (file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_variable;

  const void *temp_buffer;
  variable->sim_patch[variable->sim_patch_count] = malloc(sizeof(*(variable->sim_patch[variable->sim_patch_count])));
  memset(variable->sim_patch[variable->sim_patch_count], 0, sizeof(*(variable->sim_patch[variable->sim_patch_count])));

  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

  temp_buffer = read_from_this_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = (unsigned char*)temp_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;
  variable->io_state = 1;

  file->idx->variable[file->idx->variable_index_tracker] = variable;

  file->idx->variable_index_tracker++;
  file->local_variable_count++;

  return PIDX_success;
}



int dump_meta_data(PIDX_variable variable 
#if PIDX_HAVE_MPI
                   , MPI_Comm comm
#endif
                  )
{
  int nprocs = 1;
  FILE* meta_data_file;
  int64_t *rank_r_offset, *rank_r_count;
  
  rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  rank_r_count = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  MPI_Comm_size(comm, &nprocs);
  MPI_Allgather(variable->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
  MPI_Allgather(variable->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
#else
  memcpy(rank_r_offset, variable->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(rank_r_count, variable->sim_patch[0]->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
#endif
  
  meta_data_file = fopen("Box_info.txt", "w");
  if (!meta_data_file) 
    return PIDX_err_name;
  
  //for(i = 0; i < nprocs; i++)
  //  fprintf(meta_data_file, "[%d]: %lld %lld %lld %lld %lld - %lld %lld %lld %lld %lld\n", i, rank_r_offset[PIDX_MAX_DIMENSIONS * i + 0], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 1], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 2], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 3], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 4], rank_r_count[PIDX_MAX_DIMENSIONS * i + 0], rank_r_count[PIDX_MAX_DIMENSIONS * i + 1], rank_r_count[PIDX_MAX_DIMENSIONS * i + 2], rank_r_count[PIDX_MAX_DIMENSIONS * i + 3], rank_r_count[PIDX_MAX_DIMENSIONS * i + 4]);
  
  fclose(meta_data_file);
  
  return PIDX_success;
}

static void PIDX_init_timming_buffers2(PIDX_file file)
{
  int timer_count = file->idx->variable_count / (file->var_pipe_length + 1);
  if (file->idx->variable_count % (file->var_pipe_length + 1) != 0)
    timer_count = timer_count + 1;

  startup_start = malloc (sizeof(double) * timer_count);                 memset(startup_start, 0, sizeof(double) * timer_count);
  startup_end = malloc (sizeof(double) * timer_count);                   memset(startup_end, 0, sizeof(double) * timer_count);

  rst_start = malloc (sizeof(double) * timer_count);                     memset(rst_start, 0, sizeof(double) * timer_count);
  rst_end = malloc (sizeof(double) * timer_count);                       memset(rst_end, 0, sizeof(double) * timer_count);
  hz_start = malloc (sizeof(double) * timer_count);                      memset(hz_start, 0, sizeof(double) * timer_count);
  hz_end = malloc (sizeof(double) * timer_count);                        memset(hz_end, 0, sizeof(double) * timer_count);

  cleanup_start = malloc (sizeof(double) * timer_count);                 memset(cleanup_start, 0, sizeof(double) * timer_count);
  cleanup_end = malloc (sizeof(double) * timer_count);                   memset(cleanup_end, 0, sizeof(double) * timer_count);
  finalize_start = malloc (sizeof(double) * timer_count);                memset(finalize_start, 0, sizeof(double) * timer_count);
  finalize_end = malloc (sizeof(double) * timer_count);                  memset(finalize_end, 0, sizeof(double) * timer_count);

  buffer_start = malloc (sizeof(double) * timer_count);                  memset(buffer_start, 0, sizeof(double) * timer_count);
  buffer_end = malloc (sizeof(double) * timer_count);                    memset(buffer_end, 0, sizeof(double) * timer_count);

  chunk_start =  malloc (sizeof(double) * timer_count);                  memset(chunk_start, 0, sizeof(double) * timer_count);
  chunk_end =  malloc (sizeof(double) * timer_count);                    memset(chunk_end, 0, sizeof(double) * timer_count);
  compression_start =  malloc (sizeof(double) * timer_count);            memset(compression_start, 0, sizeof(double) * timer_count);
  compression_end =  malloc (sizeof(double) * timer_count);              memset(compression_end, 0, sizeof(double) * timer_count);

  agg_start = malloc (sizeof(double*) * file->idx->variable_count);       memset(agg_start, 0, sizeof(double*) * file->idx->variable_count);
  agg_end = malloc (sizeof(double*) * file->idx->variable_count);         memset(agg_end, 0, sizeof(double*) * file->idx->variable_count);
  io_start = malloc (sizeof(double*) * file->idx->variable_count);        memset(io_start, 0, sizeof(double*) * file->idx->variable_count);
  io_end = malloc (sizeof(double*) * file->idx->variable_count);          memset(io_end, 0, sizeof(double*) * file->idx->variable_count);
  io_per_process_start = malloc (sizeof(double*) * file->idx->variable_count);          memset(io_per_process_start, 0, sizeof(double*) * file->idx->variable_count);
  io_per_process_end = malloc (sizeof(double*) * file->idx->variable_count);          memset(io_per_process_end, 0, sizeof(double*) * file->idx->variable_count);

  int i = 0;
  for (i = 0; i < file->idx->variable_count; i++)
  {
    agg_start[i] = malloc (sizeof(double) * file->layout_count);       memset(agg_start[i], 0, sizeof(double) * file->layout_count);
    agg_end[i] = malloc (sizeof(double) * file->layout_count);         memset(agg_end[i], 0, sizeof(double) * file->layout_count);
    io_start[i] = malloc (sizeof(double) * file->layout_count);        memset(io_start[i], 0, sizeof(double) * file->layout_count);
    io_end[i] = malloc (sizeof(double) * file->layout_count);          memset(io_end[i], 0, sizeof(double) * file->layout_count);

    io_per_process_start[i] = malloc (sizeof(double) * file->layout_count);        memset(io_per_process_start[i], 0, sizeof(double) * file->layout_count);
    io_per_process_end[i] = malloc (sizeof(double) * file->layout_count);          memset(io_per_process_end[i], 0, sizeof(double) * file->layout_count);
  }
}


static void PIDX_init_timming_buffers1()
{
  init_start = malloc (sizeof(double) * 160);                    memset(init_start, 0, sizeof(double) * 160);
  init_end = malloc (sizeof(double) * 160);                      memset(init_end, 0, sizeof(double) * 160);
  write_init_start = malloc (sizeof(double) * 160);              memset(write_init_start, 0, sizeof(double) * 160);
  write_init_end = malloc (sizeof(double) * 160);                memset(write_init_end, 0, sizeof(double) * 160);
}


PIDX_return_code PIDX_set_ROI_writes(PIDX_file file)
{
  if (file == NULL)
    return PIDX_err_file;

  file->ROI_writes = 1;
  file->idx->enable_rst = 0;
  //file->idx->enable_agg = 0;

  return PIDX_success;
}
