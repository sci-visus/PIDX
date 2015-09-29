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
#define PIDX_DEBUG_OUTPUT 1

static int vp = 0;
static int hp = 0;
static double sim_start = 0, sim_end = 0;
static double *write_init_start = 0, *write_init_end = 0;
static double *init_start, *init_end;
static double *rst_start, *rst_end;
static double *hz_start, *hz_end;
static double *agg_start, *agg_end;
static double *io_start, *io_end;
static double *cleanup_start, *cleanup_end;
static double *finalize_start, *finalize_end;
static double *chunk_start, *chunk_end;
static double *buffer_start, *buffer_end;
static double *compression_start, *compression_end;
static double *agg_1, *agg_2, *agg_3, *agg_4, *agg_5, *agg_6;
static double *block_1, *block_2, *block_3;

static int caching_state = 0;
static int time_step_caching = 0;
static uint32_t* cached_header_copy;

static void PIDX_init_timming_buffers();
static PIDX_return_code PIDX_write(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_read(PIDX_file file, int start_var_index, int end_var_index);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* file_name, int current_time_step);
static PIDX_return_code PIDX_cache_headers(PIDX_file file);
static PIDX_return_code PIDX_parameter_validate(PIDX_file file, int var_index);
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
  PIDX_io_id io_id;                            ///< IO phase id
    
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

  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;    ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated
};


/// Returns elapsed time
double PIDX_get_time()
{
#if PIDX_HAVE_MPI
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
  PIDX_init_timming_buffers();
  
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
  (*file)->var_pipe_length = 0;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;

  (*file)->debug_rst = 0;
  (*file)->debug_hz = 0;
  (*file)->idx_d->color = 0;

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
  (*file)->idx->enable_agg = 2;
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
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx_d->dimension = 0;
  (*file)->idx_d->samples_per_block = pow(2, PIDX_default_bits_per_block);
  (*file)->idx_d->maxh = 0;
  (*file)->idx_d->max_file_count = 0;
  (*file)->idx_d->fs_block_size = 0;
  (*file)->idx_d->start_fs_block = 0;
  (*file)->idx_d->aggregation_factor = 1;
  (*file)->idx_d->dump_agg_info = 0;
  (*file)->idx_d->dump_io_info = 0;
  memset((*file)->idx_d->agg_dump_dir_name, 0, 512*sizeof(char));
  memset((*file)->idx_d->io_dump_dir_name, 0, 512*sizeof(char));
  (*file)->idx_d->res_from = 0;
  (*file)->idx_d->res_to = 0;
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
  MPI_Bcast(&((*file)->idx_d->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
    
  return PIDX_success;
}


/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  PIDX_init_timming_buffers();
  
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
  (*file)->var_pipe_length = 0;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;

  (*file)->debug_rst = 0;
  (*file)->debug_hz = 0;
  (*file)->idx_d->color = 0;

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
  (*file)->idx->enable_agg = 2;
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
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx_d->dimension = 0;
  (*file)->idx_d->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);
  (*file)->idx_d->maxh = 0;
  (*file)->idx_d->max_file_count = 0;
  (*file)->idx_d->fs_block_size = 0;
  (*file)->idx_d->start_fs_block = 0;
  (*file)->idx_d->aggregation_factor = 1;
  (*file)->idx_d->dump_agg_info = 0;
  (*file)->idx_d->dump_io_info = 0;
  memset((*file)->idx_d->agg_dump_dir_name, 0, 512*sizeof(char));
  memset((*file)->idx_d->io_dump_dir_name, 0, 512*sizeof(char));
  (*file)->idx_d->res_from = 0;
  (*file)->idx_d->res_to = 0;

  int var = 0, variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];
  
  if (rank == 0)
  {
    FILE *fp = fopen((*file)->idx->filename, "r");
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
          
          pch1 = strtok(line, " *+");
          while (pch1 != NULL)
          {
            if (count == 0)
              strcpy((*file)->idx->variable[variable_counter]->var_name, strdup(pch1));
              //strcpy((*file)->idx->variable[variable_counter]->var_name, pch1);

            if (count == 1)
              (*file)->idx->variable[variable_counter]->values_per_sample = atoi(pch1);

              
            if (count == 2) 
            {
              len = strlen(pch1) - 1;
              if (pch1[len] == '\n')
                pch1[len] = 0;
              if (strcmp(pch1, "float64") == 0 && (*file)->idx->variable[variable_counter]->values_per_sample == 1)
                strcpy((*file)->idx->variable[variable_counter]->type_name, "1*float64");
              else if (strcmp(pch1, "float64") == 0 && (*file)->idx->variable[variable_counter]->values_per_sample == 3)
                  strcpy((*file)->idx->variable[variable_counter]->type_name, "3*float64");
            }
            count++;
            pch1 = strtok(NULL, " *+");
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
  MPI_Bcast((*file)->idx->bounds, 5, MPI_LONG_LONG, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx->blocks_per_file), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx->bits_per_block), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx->variable_count), 1, MPI_INT, 0, (*file)->comm);

  (*file)->idx_d->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);
  
  if(rank != 0)
  {
    for (var = 0; var < (*file)->idx->variable_count; var++)
      (*file)->idx->variable[var] = malloc(sizeof (*((*file)->idx->variable[var])));
  }
#endif
  
  for (var = 0; var < (*file)->idx->variable_count; var++)
  {
#if PIDX_HAVE_MPI
    MPI_Bcast(&((*file)->idx->variable[var]->values_per_sample), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast((*file)->idx->variable[var]->var_name, 512, MPI_CHAR, 0, (*file)->comm);
    MPI_Bcast((*file)->idx->variable[var]->type_name, 512, MPI_CHAR, 0, (*file)->comm);

    int ret;
    int bits_per_sample = 0;
    ret = PIDX_default_bits_per_datatype((*file)->idx->variable[var]->type_name, &bits_per_sample);
    if (ret != PIDX_success)  return PIDX_err_file;
    (*file)->idx->variable[var]->bits_per_value = bits_per_sample/(*file)->idx->variable[var]->values_per_sample;
    //(*file)->idx->variable[var]->values_per_sample = 1;

#else
    (*file)->idx->variable[var]->values_per_sample = 1;
    (*file)->idx->variable[var]->bits_per_value = (int)sizeof(double) * 8;
    strcpy((*file)->idx->variable[var]->type_name, "1*float64");
#endif
    (*file)->idx->variable[var]->sim_patch_count = 0;
  }

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
  MPI_Bcast(&((*file)->idx_d->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
  
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
    file->idx_d->samples_per_block = 1;
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



PIDX_return_code PIDX_set_aggregation_factor(PIDX_file file, int agg_factor)
{
  if(!file)
    return PIDX_err_file;

  file->idx_d->aggregation_factor = agg_factor;
  
  return PIDX_success;
}



PIDX_return_code PIDX_get_aggregation_factor(PIDX_file file, int *agg_factor)
{
  if(!file)
    return PIDX_err_file;
  
  *agg_factor = file->idx_d->aggregation_factor;
  
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
  }
  if (file->idx->compression_bit_rate == 16)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 2;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
  }
  if (file->idx->compression_bit_rate == 8)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 3;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
  }
  if (file->idx->compression_bit_rate == 4)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 4;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
  }
  if (file->idx->compression_bit_rate == 2)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 5;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
  }
  if (file->idx->compression_bit_rate == 1)
  {
    file->idx->bits_per_block = file->idx->bits_per_block + 6;
    file->idx_d->samples_per_block = (int)pow(2, file->idx->bits_per_block);
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
  memcpy(variable->sim_patch[variable->sim_patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->sim_patch[variable->sim_patch_count]->size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

  //write_to_this_buffer = temp_buffer;
  variable->sim_patch[variable->sim_patch_count]->buffer = write_to_this_buffer;

  variable->data_layout = data_layout;
  variable->sim_patch_count = variable->sim_patch_count + 1;

  return PIDX_success;
}



static PIDX_return_code populate_idx_dataset(PIDX_file file)
{
  int d, rank;
  MPI_Comm_rank(file->comm, &rank);

  int i, j, counter = 0, file_number = 0;
    int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };
  
  PointND bounds_point;
  
  if(file->debug_do_chunk == 1)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      if (file->idx->bounds[d] % file->idx->chunk_size[d] == 0)
        file->idx->chunked_bounds[d] = (int) file->idx->bounds[d] / file->idx->chunk_size[d];
      else
        file->idx->chunked_bounds[d] = (int) (file->idx->bounds[d] / file->idx->chunk_size[d]) + 1;
    }
  }
  else
    memcpy(file->idx->chunked_bounds, file->idx->bounds, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);

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
  int64_t max_sample_per_file = (uint64_t) file->idx_d->samples_per_block * file->idx->blocks_per_file;

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;
  
  file->idx_d->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(file->idx_d->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));
  
  file->idx->variable[file->local_variable_index]->global_block_layout = malloc(sizeof (*file->idx->variable[file->local_variable_index]->global_block_layout));

  PIDX_block_layout block_layout = file->idx->variable[file->local_variable_index]->global_block_layout;
#if 0
  //if (file->ROI_writes == 0)
  //{
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++) 
    {
      bounding_box[0][i] = 0;
      bounding_box[1][i] = cb[i];
    }

    PIDX_blocks_create_layout(bounding_box, file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->maxh, file->idx_d->res_from, file->idx_d->res_to, file->idx->bitPattern, block_layout);
    
    int k = 1;
    for (i = 1; i < (block_layout->levels); i++)
    {
      counter = 0;
      for (j = 0 ; j < k ; j++)
      {
        if(block_layout->hz_block_number_array[i][j] != 0)
        {
          block_layout->hz_block_number_array[i][counter] = block_layout->hz_block_number_array[i][j];
          counter++;
        }
      }
      k = k * 2;
    }
  //}
#endif
#if 1
  //else
  //{
    int p = 0, ctr = 1;
    PIDX_blocks_initialize_layout(block_layout, file->idx_d->maxh, file->idx->bits_per_block);

    PIDX_block_layout all_patch_local_block_layout = malloc(sizeof (*all_patch_local_block_layout));
    PIDX_blocks_initialize_layout(all_patch_local_block_layout, file->idx_d->maxh, file->idx->bits_per_block);

    for (p = 0 ; p < file->idx->variable[file->local_variable_index]->sim_patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
        bounding_box[0][i] = file->idx->variable[file->local_variable_index]->sim_patch[p]->offset[i];
        bounding_box[1][i] = file->idx->variable[file->local_variable_index]->sim_patch[p]->size[i] + file->idx->variable[file->local_variable_index]->sim_patch[p]->offset[i];
      }
      
      PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
      PIDX_blocks_create_layout (bounding_box, file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->maxh, file->idx_d->res_from, file->idx_d->res_to, file->idx->bitPattern, per_patch_local_block_layout);
      
      ctr = 1;
      for (i = 1; i < (all_patch_local_block_layout->levels); i++)
      {
        for (j = 0 ; j < ctr ; j++)
        {
          if(per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
            all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
        }
        ctr = ctr * 2;
      }
      PIDX_blocks_free_layout(per_patch_local_block_layout);
      free(per_patch_local_block_layout);
      per_patch_local_block_layout = 0;
    }
  
    int level_count = 1;
    for (i = 1; i < (block_layout->levels); i++)
    {
#if PIDX_HAVE_MPI
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
#else
      memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
      level_count = level_count * 2;
    }
    
    PIDX_blocks_free_layout(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
    
    ctr = 1;
    for (i = 1; i < (block_layout->levels); i++)
    {
      for (j = 0 ; j < ctr ; j++)
      {
        if(block_layout->hz_block_number_array[i][j] != 0)
          block_layout->hz_block_count_array[i]++;
      }    
      
      counter = 0;
      for (j = 0 ; j < ctr ; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          block_layout->hz_block_number_array[i][counter] = block_layout->hz_block_number_array[i][j];
          counter++;
        }
      }
      memset(block_layout->hz_block_number_array[i]+counter, 0, (ctr-counter)*sizeof(int));
      ctr = ctr * 2;
    }
  //}
#endif
  //if (rank == 0)
  //  PIDX_blocks_print_layout(block_layout);
  file->idx->variable[file->local_variable_index]->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(file->idx->variable[file->local_variable_index]->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));
  
  file->idx->variable[file->local_variable_index]->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(file->idx->variable[file->local_variable_index]->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  file->idx->variable[file->local_variable_index]->file_index[0] = 1;
  file->idx_d->file_bitmap[0] = 1;
  file->idx->variable[file->local_variable_index]->block_count_per_file[0] = 1;
  
  for (i = 1; i < block_layout->levels; i++)
  {
    for (j = 0; j < block_layout->hz_block_count_array[i]; j++)
    {
      file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
      file->idx_d->file_bitmap[file_number] = 1;
      file->idx->variable[file->local_variable_index]->file_index[file_number] = 1;
      file->idx->variable[file->local_variable_index]->block_count_per_file[file_number]++;
    }
  }
  
  //PIDX_blocks_print_layout(block_layout);
  
  file->idx->variable[file->local_variable_index]->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (file->idx->variable[file->local_variable_index]->file_index[i] == 1)
      file->idx->variable[file->local_variable_index]->existing_file_count++;

  file->idx->variable[file->local_variable_index]->existing_file_index = (int*) malloc(file->idx->variable[file->local_variable_index]->existing_file_count * sizeof (int));
  memset(file->idx->variable[file->local_variable_index]->existing_file_index, 0, file->idx->variable[file->local_variable_index]->existing_file_count * sizeof (int));
  
  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (file->idx->variable[file->local_variable_index]->file_index[i] == 1)
    {
      //if (rank == 0)
      //  printf("BPF %d = %d\n", i, file->idx->variable[file->local_variable_index]->block_count_per_file[i]);
      file->idx->variable[file->local_variable_index]->existing_file_index[count] = i;
      count++;
    }
  }

    
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



static PIDX_return_code PIDX_cache_headers(PIDX_file file)
{
#if 1
  if(!file)
    return PIDX_err_file;

  int i = 0, j = 0, k = 0;
  off_t data_offset = 0, base_offset = 0;
  int block_negative_offset = 0;
  int all_scalars = 0;
  int total_header_size;
  int64_t total_chunk_size = (file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2] * file->idx->chunk_size[3] * file->idx->chunk_size[4]) / (64 / file->idx->compression_bit_rate);
  
  total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
  if (total_header_size % file->idx_d->fs_block_size)
    file->idx_d->start_fs_block++;
  
  cached_header_copy = (uint32_t*)malloc(total_header_size);
  memset(cached_header_copy, 0, total_header_size);

  /*
  for (i = 0; i < file->idx->variable_count; i++)
  {
    if (file->idx->variable[i]->values_per_sample != 1)
    {
      all_scalars = 0;
      break;
    }
  }
  */

  for (i = 0; i < file->idx->blocks_per_file; i++)
  {
    if (PIDX_blocks_is_block_present((i + (file->idx->blocks_per_file * file->idx_d->agg_buffer->file_number)), file->idx->variable[file->local_variable_index]->global_block_layout))
    {
      block_negative_offset = PIDX_blocks_find_negative_offset(file->idx->blocks_per_file, (i + (file->idx->blocks_per_file * file->idx_d->agg_buffer->file_number)), file->idx->variable[file->local_variable_index]->global_block_layout);
      
      for (j = 0; j < file->idx->variable_count; j++)
      {
        base_offset = 0;
        if (all_scalars == 0)
        {
          for (k = 0; k < j; k++)
            base_offset = base_offset + (file->idx->variable[file->local_variable_index]->block_count_per_file[file->idx_d->agg_buffer->file_number]) * (file->idx->variable[k]->bits_per_value / 8) * total_chunk_size * file->idx_d->samples_per_block * file->idx->variable[k]->values_per_sample;
        }
        else
          base_offset =  j * (file->idx->variable[file->local_variable_index]->block_count_per_file[file->idx_d->agg_buffer->file_number]) * (file->idx->variable[file->local_variable_index]->bits_per_value / 8) * total_chunk_size * file->idx_d->samples_per_block * file->idx->variable[file->local_variable_index]->values_per_sample;
        
        data_offset = (((i) - block_negative_offset) * file->idx_d->samples_per_block) * (file->idx->variable[j]->bits_per_value / 8) * total_chunk_size * file->idx->variable[j]->values_per_sample;
        
        data_offset = base_offset + data_offset + file->idx_d->start_fs_block * file->idx_d->fs_block_size;

        //TODO
        cached_header_copy[12 + ((i + (file->idx->blocks_per_file * j))*10 )] = htonl(data_offset);
        cached_header_copy[14 + ((i + (file->idx->blocks_per_file * j))*10)] = htonl(file->idx_d->samples_per_block * (file->idx->variable[j]->bits_per_value / 8) * total_chunk_size * file->idx->variable[j]->values_per_sample);
      }
    }
  }
#endif
  return PIDX_success;
}



PIDX_return_code PIDX_set_resolution(PIDX_file file, int hz_from, int hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  file->idx_d->res_from = hz_from;
  file->idx_d->res_to = hz_to;
  
  return PIDX_success;
}


PIDX_return_code PIDX_get_resolution(PIDX_file file, int *hz_from, int *hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  *hz_from = file->idx_d->res_from;
  *hz_to = file->idx_d->res_to;
  
  return PIDX_success;
}


static PIDX_return_code PIDX_parameter_validate(PIDX_file file, int var_index)
{
  int d, nprocs = 1;

#if PIDX_HAVE_MPI
  int local_patch_count = 0, total_patch_count = 0;
  int ret;
  MPI_Comm_size(file->comm, &nprocs);
  PIDX_variable var0 = file->idx->variable[var_index];
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

  if (file->idx->variable_count == file->idx->variable_index_tracker)
  {
    if (file->idx->enable_agg != 0)
    {
      int v, no_of_aggregators = 0;
      for (v = 0; v < file->idx->variable_count ; v++)
        no_of_aggregators = no_of_aggregators + file->idx->variable[v]->values_per_sample * file->idx->variable[file->local_variable_index]->existing_file_count;

      if (nprocs < no_of_aggregators * file->idx_d->aggregation_factor)
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

  file->var_pipe_length = file->idx->variable_count - 1;
  //file->var_pipe_length = 1;

  return PIDX_success;
}


static PIDX_return_code PIDX_write(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  int j = 0, p, var = 0, d = 0;
  int total_header_size;
  PIDX_return_code ret;
  //static int header_io = 0;
  int rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm,  &nprocs);
#endif

  ret = populate_idx_dataset(file);
  if (ret != PIDX_success)
    return PIDX_err_file;

#if PIDX_DEBUG_OUTPUT
  if (rank == 0)
    printf("[1] Finish Populating IDX dataset: File Count %d maxh = %d\n", file->idx_d->max_file_count, file->idx_d->maxh);
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
    
    ret = PIDX_parameter_validate(file, start_var_index);
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

#if PIDX_HAVE_MPI
    file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
#endif

  /*  STEP 1: create files and folders based on the extents of the variable group
   *  STEP 2: if flush used
   *            when variable count is met, then write header information
   *          else
   *            pass the header buffer to the agg phase when no var pipelining is done (else if pipe, then go to if)
   *  STEP 3: at the end of all IO, write the .idx file
   */

  write_init_start[hp] = PIDX_get_time();

#if !SIMULATE_IO

  /* STEP 1 */
  file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
  ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
  if (ret != PIDX_success)
    return PIDX_err_header;
#endif
  ret = PIDX_header_io_file_create(file->header_io_id);
  if (ret != PIDX_success)
    return PIDX_err_header;

  /* STEP 2 */
  if (file->idx->variable_index_tracker < file->idx->variable_count )
  {
    // Create the header
    ret = PIDX_header_io_file_write(file->header_io_id, 0);
    if (ret != PIDX_success)
      return PIDX_err_header;
    file->flush_used = 1;
  }

  if (file->idx->variable_index_tracker == file->idx->variable_count)
  {
    // Write the header
    if (file->flush_used == 1 || file->idx->enable_agg == 0 || file->idx->enable_agg == 1)
    {
      ret = PIDX_header_io_file_write(file->header_io_id, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
    else if (/*(file->var_pipe_length < file->idx->variable_count - 1) && */ caching_state == 0)
    {
      ret = PIDX_header_io_file_write(file->header_io_id, 1);
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

#endif

  write_init_end[hp] = PIDX_get_time();
  hp++;

#if PIDX_DEBUG_OUTPUT
  if (rank == 0)
    printf("[2] Finish Creating File heirarchy\n");
#endif

#if 1
  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);

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
    file->agg_id = PIDX_agg_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the I/O ID */
    file->io_id = PIDX_io_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
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
    ret = PIDX_agg_set_communicator(file->agg_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    /* Attaching the communicator to the I/O phase */
    ret = PIDX_io_set_communicator(file->io_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_io;
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

    ret = PIDX_agg_meta_data_create(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[3] All Modules [R B C H A I] Initialized: Variable index [%d %d] Variable Pipe length %d\n", start_index, end_index, file->var_pipe_length);
#endif

    /*--------------------------------------------RST [start]------------------------------------------------*/
    rst_start[vp] = PIDX_get_time();

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
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
    if (rank == 0)
      printf("[R] Restructuring Phase Completed\n");
#endif

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      ret = HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[R] Restructuring Verified\n");
#endif

    rst_end[vp] = PIDX_get_time();
    /*--------------------------------------------RST [end]---------------------------------------------------*/

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[4] Restructuring Phase Complete\n");
#endif

    /*----------------------------------------Chunking [start]------------------------------------------------*/
    chunk_start[vp] = PIDX_get_time();

    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[4] Restructuring Phase Complete\n");
#endif

    /* Perform Chunking */
    if (file->debug_do_chunk == 1)
    {
      ret = PIDX_chunk_write(file->chunk_id);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

    /* Verifying the correctness of the chunking phase */
    if (file->debug_chunk == 1)
    {
      ret = HELPER_chunking(file->chunk_id);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_rst_buf_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    chunk_end[vp] = PIDX_get_time();
    /*-----------------------------------------Chunking [end]------------------------------------------------*/

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[5] Chunking Phase Complete\n");
#endif
    
    /*----------------------------------------Compression [start]--------------------------------------------*/
    compression_start[vp] = PIDX_get_time();

    /* Perform Compression */
    if (file->debug_do_compress == 1)
    {
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
        return PIDX_err_compress;
    }

    compression_end[vp] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[6] Compression Phase Complete\n");
#endif


    /*---------------------------------------------HZ [start]------------------------------------------------*/
    hz_start[vp] = PIDX_get_time();

    /* Creating the buffers required for HZ encoding */
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    /* Perform HZ encoding */
    if (file->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

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


#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[6] HZ Encoding Phase Complete\n");
#endif


    /*----------------------------------------------Agg [start]-----------------------------------------------*/
    agg_start[vp] = PIDX_get_time();

    /* Creating the buffers required for Aggregation */
    ret = PIDX_agg_buf_create(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    /* Perform Aggregation */
    if (file->debug_do_agg == 1)
    {
      ret = PIDX_agg_write(file->agg_id);
      if (ret != PIDX_success)
        return PIDX_err_agg;
    }

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    agg_end[vp] = PIDX_get_time();
    /*--------------------------------------------Agg [end]--------------------------------------------------*/
    
#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[7] Aggregation Phase Complete\n");
#endif

#if 1

    /*--------------------------------------------IO [start]--------------------------------------------------*/
    io_start[vp] = PIDX_get_time();

    /* Initialization ONLY ONCE for all TIME STEPS (caching across time) */
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

    if (file->debug_do_io == 1)
    {
      if (time_step_caching == 1)
      {
        if (file->flush_used == 0)
          PIDX_io_cached_data(cached_header_copy);
      }

      ret = PIDX_io_aggregated_write(file->io_id);
      if (ret != PIDX_success)
        return PIDX_err_io;
    }

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[8] I/O Phase Complete\n");
#endif

    /* Destroy buffers allocated during aggregation phase */
    ret = PIDX_agg_buf_destroy(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    io_end[vp] = PIDX_get_time();
    /*----------------------------------------------IO [end]--------------------------------------------------*/

    
    ret = PIDX_agg_meta_data_destroy(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

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

    /* Deleting the I/O ID */
    PIDX_io_finalize(file->io_id);

    /* Deleting the aggregation ID */
    PIDX_agg_finalize(file->agg_id);

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

  PIDX_variable var00 = file->idx->variable[file->local_variable_index];
  free(var00->block_count_per_file);
  var00->block_count_per_file = 0;

  PIDX_blocks_free_layout(var00->global_block_layout);
  free(var00->global_block_layout);
  var00->global_block_layout = 0;
  free(var00->existing_file_index);
  var00->existing_file_index = 0;

  free(var00->file_index);

#if PIDX_DEBUG_OUTPUT
    if (rank == 0)
      printf("[9] Complete !!!!\n");
#endif

#endif
  return PIDX_success;
}



static PIDX_return_code PIDX_read(PIDX_file file, int start_var_index, int end_var_index)
{
  if (file->local_variable_index == file->idx->variable_count)
    return PIDX_success;

  int j = 0, p, var = 0, d = 0;
  int total_header_size;
  PIDX_return_code ret;
  //static int header_io = 0;  

#if PIDX_HAVE_MPI
  int rank = 0, nprocs = 1;
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm,  &nprocs);
#endif

  populate_idx_dataset(file);

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);

    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;

    PIDX_parameter_validate(file, start_var_index);

    file->one_time_initializations = 1;
  }

#if PIDX_HAVE_MPI
    file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
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

#if 1
  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->var_pipe_length + 1))
  {
    end_index = ((start_index + file->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->var_pipe_length);

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
    file->agg_id = PIDX_agg_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the I/O ID */
    file->io_id = PIDX_io_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
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
    ret = PIDX_agg_set_communicator(file->agg_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    /* Attaching the communicator to the I/O phase */
    ret = PIDX_io_set_communicator(file->io_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_io;
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

    ret = PIDX_agg_meta_data_create(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;


    /*----------------------------------------------IO [start]-----------------------------------------------*/
    io_start[vp] = PIDX_get_time();

    /* Creating the buffers required for Aggregation */
    ret = PIDX_agg_buf_create(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    // TODO
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;

    if (file->debug_do_io == 1)
    {
      ret = PIDX_io_aggregated_read(file->io_id);
      if (ret != PIDX_success)
        return PIDX_err_io;
    }

    io_end[vp] = PIDX_get_time();
    /*----------------------------------------------IO [End]-----------------------------------------------*/



    /*--------------------------------------------Agg [start]-----------------------------------------------*/
    agg_start[vp] = PIDX_get_time();
    /* Creating the buffers required for HZ encoding */
    // TODO
    //ret = PIDX_hz_encode_buf_create(file->hz_id);
    //if (ret != PIDX_success)
    //  return PIDX_err_hz;

    /* Perform Aggregation */
    if (file->debug_do_agg == 1)
    {
      ret = PIDX_agg_read(file->agg_id);
      if (ret != PIDX_success)
        return PIDX_err_agg;
    }

    /* Verify the HZ encoding */
    if(file->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

    // TODO
    /* Destroy buffers allocated during aggregation phase */
    //ret = PIDX_agg_buf_destroy(file->agg_id);
    //if (ret != PIDX_success)
    //  return PIDX_err_agg;
    agg_end[vp] = PIDX_get_time();
    /*----------------------------------------------Agg [end]-----------------------------------------------*/



    /*----------------------------------------Compression [start]--------------------------------------------*/
    compression_start[vp] = PIDX_get_time();
    /* Perform Compression */

    /*
    if (file->debug_do_compress == 1)
    {
      ret = PIDX_decompression(file->comp_id);
      if (ret != PIDX_success)
        return PIDX_err_compress;
    }
    */

    compression_end[vp] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/



    /*-----------------------------------------HZ encoding [start]--------------------------------------------*/
    hz_start[vp] = PIDX_get_time();

    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;

    /* Perform HZ encoding */
    if (file->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_read(file->hz_id);
      if (ret != PIDX_success)
        return PIDX_err_hz;
    }

    /* Verifying the correctness of the chunking phase */
    if (file->debug_chunk == 1)
    {
      ret = HELPER_chunking(file->chunk_id);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

    // TODO
    ret = PIDX_agg_buf_destroy(file->agg_id);
    if (ret != PIDX_success)
      return PIDX_err_agg;

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
      return PIDX_err_hz;
    hz_end[vp] = PIDX_get_time();
    /*-----------------------------------------HZ encoding [end]-------------------------------------------*/



    /*------------------------------------------Chunking [start]--------------------------------------------*/
    chunk_start[vp] = PIDX_get_time();
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Perform Chunking */
    if (file->debug_do_chunk == 1)
    {
      ret = PIDX_chunk_read(file->chunk_id);
      if (ret != PIDX_success)
        return PIDX_err_chunk;
    }

    /* Verifying the correctness of the restructuring phase */
    if (file->debug_rst == 1)
    {
      HELPER_rst(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

    /* Destroy buffers allocated during chunking phase */
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
      return PIDX_err_chunk;
    chunk_end[vp] = PIDX_get_time();
    /*--------------------------------------------Chunking [end]--------------------------------------------*/



    /*-----------------------------------------------Rst [start]--------------------------------------------*/
    rst_start[vp] = PIDX_get_time();
    /* Perform data restructuring */
    if (file->debug_do_rst == 1)
    {
      ret = PIDX_rst_read(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }

    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_rst_buf_destroy(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;
    rst_end[vp] = PIDX_get_time();
    /*-------------------------------------------------Rst [end]--------------------------------------------*/



    /*-------------------------------------------finalize [start]---------------------------------------------*/
    finalize_start[vp] = PIDX_get_time();

    /* Deleting the I/O ID */
    PIDX_io_finalize(file->io_id);

    /* Deleting the aggregation ID */
    PIDX_agg_finalize(file->agg_id);

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

    vp++;
  }

  PIDX_variable var00 = file->idx->variable[file->local_variable_index];
  free(var00->block_count_per_file);
  var00->block_count_per_file = 0;

  PIDX_blocks_free_layout(var00->global_block_layout);
  free(var00->global_block_layout);
  var00->global_block_layout = 0;
  free(var00->existing_file_index);
  var00->existing_file_index = 0;

  free(var00->file_index);

//#endif
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
    ret = PIDX_write(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
    if (ret != PIDX_success)
      return PIDX_err_flush;
  }

  else if (file->flags == PIDX_MODE_RDONLY)
  {
    ret = PIDX_read(file, file->local_variable_index, file->local_variable_index + file->local_variable_count);
    if (ret != PIDX_success)
      return PIDX_err_flush;
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
  file->write_on_close = 1;

  ret = PIDX_flush(file);
  if (ret != PIDX_success)
    return PIDX_err_close;
  
  sim_end = PIDX_get_time();
  
#if PIDX_DEBUG_OUTPUT
  int i = 0;
  double total_time = sim_end - sim_start;
  double max_time = total_time;
  int sample_sum = 0, var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
#endif

#if PIDX_HAVE_MPI
    free(file->idx_d->rank_r_offset);
    free(file->idx_d->rank_r_count);
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
      fprintf(stdout, "Rst = %d Comp = %d Agg = %d\n", file->idx->enable_rst, file->idx->compression_type, file->idx->enable_agg);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d Aggregation Factor %d Aggregator Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count, file->idx_d->aggregation_factor, file->idx->variable_count * file->idx_d->max_file_count * file->idx_d->aggregation_factor);
      fprintf(stdout, "Chunk Size %d %d %d %d %d\n", (int)file->idx->chunk_size[0], (int)file->idx->chunk_size[1], (int)file->idx->chunk_size[2], (int)file->idx->chunk_size[3], (int)file->idx->chunk_size[4]);
      fprintf(stdout, "Restructuring Box Size %d %d %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2], (int)file->idx->reg_patch_size[3], (int)file->idx->reg_patch_size[4]);
      fprintf(stdout, "Staged Aggregation = %d\n", file->idx_d->staged_aggregation);
      fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", max_time, (float) total_data / (1000 * 1000 * max_time));
      fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
      //printf("File creation time %f\n", write_init_end - write_init_start);
      
      for (var = 0; var < hp; var++)
        fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
      
      for (var = 0; var < vp; var++)
      {
        fprintf(stdout, "------------------------------------------------VG %d (START)----------------------------------------------\n", var);
        fprintf(stdout, "Init time  [RST + BRST + HZ + AGG + IO] %f \n", (init_end[var] - init_start[var]));
        
        fprintf(stdout, "Write time [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_end[var] - rst_start[var]), (chunk_end[var] - chunk_start[var]), (hz_end[var] - hz_start[var]), (agg_end[var] - agg_start[var]), (io_end[var] - io_start[var]), (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]) + (agg_end[var] - agg_start[var]) + (io_end[var] - io_start[var]));
        
        fprintf(stdout, "Block Restructuring time %f = %f + %f\n", (chunk_end[var] - chunk_start[var]), (block_2[var] - block_1[var]), (block_3[var] - block_2[var]));
        
        fprintf(stdout, "Agg time %f = %f + %f + %f + %f + %f\n", (agg_end[var] - agg_start[var]), (agg_2[var] - agg_1[var]), (agg_3[var] - agg_2[var]), (agg_4[var] - agg_3[var]), (agg_5[var] - agg_4[var]), (agg_6[var] - agg_5[var]));
        
        fprintf(stdout, "Cleanup time %f\n", cleanup_end[var] - cleanup_start[var]);
        fprintf(stdout, "-------------------------------------------------VG %d (END)-----------------------------------------------\n", var);
      }
      
      //
      //double total_agg_time = 0, all_time = 0;
      //for (p = 0; p < file->idx->variable[0]->patch_group_count; p++)
      //  for (var = 0; var < file->idx->variable_count; var++)
      //    for (i = file->idx->variable[var]->HZ_patch[p]->HZ_level_from; i < file->idx->variable[var]->HZ_patch[p]->HZ_level_to; i++)
      //    {
      //      printf("Agg Time [Patch %d Var %d Level %d] = %f\n", p, var, i, (file->idx_d->agg_level_end[p][var][i] - file->idx_d->agg_level_start[p][var][i]));
      //      total_agg_time = total_agg_time + (file->idx_d->agg_level_end[p][var][i] - file->idx_d->agg_level_start[p][var][i]);
      //    }
      
      //all_time = total_agg_time + (file->idx_d->win_time_end - file->idx_d->win_time_start) + (file->idx_d->win_free_time_end - file->idx_d->win_free_time_start);
      //printf("Total Agg Time %f = [Network + Win_Create + Win_free] %f + %f + %f\n", all_time, total_agg_time, (file->idx_d->win_time_end - file->idx_d->win_time_start), (file->idx_d->win_free_time_end - file->idx_d->win_free_time_start));
      
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
    MPI_Allreduce(&max_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, file->global_comm);
    MPI_Comm_rank(file->global_comm, &global_rank);
    MPI_Comm_size(file->global_comm, &global_nprocs);
#endif
    
    if (global_max_time == total_time)
    {
      fprintf(stdout, "\n------------- PIDX -------------\n");
      
      fprintf(stdout, "\n==========================================================================================================\n");
      fprintf(stdout, "[%d] Combined Time step %d File name %s\n", global_rank, file->idx->current_time_step, file->idx->filename);
      fprintf(stdout, "Cores %d Global Data %lld %lld %lld Variables %d IDX Count %d = %d x %d x %d\n", global_nprocs, (long long) file->idx->bounds[0] * file->idx_count[0], (long long) file->idx->bounds[1] * file->idx_count[1], (long long) file->idx->bounds[2] * file->idx_count[2], file->idx->variable_count, file->idx_count[0] * file->idx_count[1] * file->idx_count[2], file->idx_count[0], file->idx_count[1], file->idx_count[2]);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d Aggregation Factor %d Aggregator Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count, file->idx_d->aggregation_factor, file->idx->variable_count * file->idx_d->max_file_count * file->idx_d->aggregation_factor );
      fprintf(stdout, "Blocks Restructuring Size %d %d %d %d %d\n", (int)file->idx->chunk_size[0], (int)file->idx->chunk_size[1], (int)file->idx->chunk_size[2], (int)file->idx->chunk_size[3], (int)file->idx->chunk_size[4]);
      fprintf(stdout, "Staged Aggregation = %d\n", file->idx_d->staged_aggregation);
      fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", global_max_time, (float) total_data / (1000 * 1000 * global_max_time));
      fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
      //printf("File creation time %f\n", write_init_end - write_init_start);
      
      for (var = 0; var < hp; var++)
        fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
      
      for (var = 0; var < vp; var++)
      {
        fprintf(stdout, "------------------------------------------------VG %d (START)----------------------------------------------\n", var);

        fprintf(stdout, "Init time  [RST + BRST + HZ + AGG + IO] %f\n", (init_end[var] - init_start[var]));
        
        fprintf(stdout, "Write time [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_end[var] - rst_start[var]), (chunk_end[var] - chunk_start[var]), (hz_end[var] - hz_start[var]), (agg_end[var] - agg_start[var]), (io_end[var] - io_start[var]), (rst_end[var] - rst_start[var]) + (chunk_end[var] - chunk_start[var]) + (hz_end[var] - hz_start[var]) + (agg_end[var] - agg_start[var]) + (io_end[var] - io_start[var]));
        
        fprintf(stdout, "Block Restructuring time %f = %f + %f\n", (chunk_end[var] - chunk_start[var]), (block_2[var] - block_1[var]), (block_3[var] - block_2[var]));
        
        fprintf(stdout, "Agg time %f = %f + %f + %f + %f = %f\n", (agg_end[var] - agg_start[var]), (agg_2[var] - agg_1[var]), (agg_3[var] - agg_2[var]), (agg_4[var] - agg_3[var]), (agg_5[var] - agg_4[var]), (agg_6[var] - agg_5[var]));
        
        fprintf(stdout, "Cleanup time %f\n", cleanup_end[var] - cleanup_start[var]);
        fprintf(stdout, "--------------------------------------------------VG %d (END)-----------------------------------------------\n", var);
      }
      fprintf(stdout, "==========================================================================================================\n");
    }
    
  }
#endif

  /*
  for (p = 0; p < file->idx->variable[0]->patch_group_count; p++)
  {
    for(var = 0; var < file->idx->variable_count; var++)
    {
      free(file->idx_d->agg_level_start[p][var]);
      free(file->idx_d->agg_level_end[p][var]);
    }
    free(file->idx_d->agg_level_start[p]);
    free(file->idx_d->agg_level_end[p]);

  }
  free(file->idx_d->agg_level_start);
  free(file->idx_d->agg_level_end);
  */
#if 1
  for (i = 0; i < 1024; i++)
  {
    free(file->idx->variable[i]);
    file->idx->variable[i] = 0;
  }
  
  file->idx->variable_count = 0;
  
  //free(file->idx->bounds);         file->idx->bounds = 0;
  free(file->idx);                        file->idx = 0;
  free(file->idx_d->file_bitmap);   file->idx_d->file_bitmap = 0;
  free(file->idx_d);                file->idx_d = 0;
#endif
  
#if PIDX_HAVE_MPI
  MPI_Comm_free(&(file->comm));
  if (file->idx_count[0] * file->idx_count[1] * file->idx_count[2] != 1)
  //if (file->access->topology_aware_io != 0)
    MPI_Comm_free(&(file->global_comm));
#endif
  
  free(file);
  
  vp = 0;
  hp = 0;
  
  free(init_start);                     init_start        = 0;
  free(init_end);                       init_end          = 0;
  free(write_init_start);               write_init_start        = 0;
  free(write_init_end);                 write_init_end          = 0;
  free(rst_start);                      rst_start               = 0;
  free(rst_end);                        rst_end                 = 0;
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
  free(block_1);                        block_1                 = 0;
  free(block_2);                        block_2                 = 0;
  free(block_3);                        block_3                 = 0;
  free(agg_1);                          agg_1                   = 0;
  free(agg_2);                          agg_2                   = 0;
  free(agg_3);                          agg_3                   = 0;
  free(agg_4);                          agg_4                   = 0;
  free(agg_5);                          agg_5                   = 0;
  free(agg_6);                          agg_6                   = 0;

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
  
  if(variable_index <= 0)
    return PIDX_err_size;
  
  if(file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  file->idx->variable_index_tracker = variable_index;
  return PIDX_success;
}



PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  if(file->idx->variable_index_tracker >= file->idx->variable_count)
    return PIDX_err_count;
  
  (*variable) = file->idx->variable[file->idx->variable_index_tracker];
  
  return PIDX_err_not_implemented;
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

static void PIDX_init_timming_buffers()
{
  init_start = malloc (sizeof(double) * 64);                    memset(init_start, 0, sizeof(double) * 64);
  init_end = malloc (sizeof(double) * 64);                      memset(init_end, 0, sizeof(double) * 64);
  write_init_start = malloc (sizeof(double) * 64);              memset(write_init_start, 0, sizeof(double) * 64);
  write_init_end = malloc (sizeof(double) * 64);                memset(write_init_end, 0, sizeof(double) * 64);
  rst_start = malloc (sizeof(double) * 64);                     memset(rst_start, 0, sizeof(double) * 64);
  rst_end = malloc (sizeof(double) * 64);                       memset(rst_end, 0, sizeof(double) * 64);
  hz_start = malloc (sizeof(double) * 64);                      memset(hz_start, 0, sizeof(double) * 64);
  hz_end = malloc (sizeof(double) * 64);                        memset(hz_end, 0, sizeof(double) * 64);
  agg_start = malloc (sizeof(double) * 64);                     memset(agg_start, 0, sizeof(double) * 64);
  agg_end = malloc (sizeof(double) * 64);                       memset(agg_end, 0, sizeof(double) * 64);
  io_start = malloc (sizeof(double) * 64);                      memset(io_start, 0, sizeof(double) * 64);
  io_end = malloc (sizeof(double) * 64);                        memset(io_end, 0, sizeof(double) * 64);
  cleanup_start = malloc (sizeof(double) * 64);                 memset(cleanup_start, 0, sizeof(double) * 64);
  cleanup_end = malloc (sizeof(double) * 64);                   memset(cleanup_end, 0, sizeof(double) * 64);
  finalize_start = malloc (sizeof(double) * 64);                memset(finalize_start, 0, sizeof(double) * 64);
  finalize_end = malloc (sizeof(double) * 64);                  memset(finalize_end, 0, sizeof(double) * 64);

  buffer_start = malloc (sizeof(double) * 64);                  memset(buffer_start, 0, sizeof(double) * 64);
  buffer_end = malloc (sizeof(double) * 64);                    memset(buffer_end, 0, sizeof(double) * 64);

  chunk_start =  malloc (sizeof(double) * 64);              memset(chunk_start, 0, sizeof(double) * 64);
  chunk_end =  malloc (sizeof(double) * 64);                memset(chunk_end, 0, sizeof(double) * 64);
  compression_start =  malloc (sizeof(double) * 64);            memset(compression_start, 0, sizeof(double) * 64);
  compression_end =  malloc (sizeof(double) * 64);              memset(compression_end, 0, sizeof(double) * 64);

  block_1 = malloc (sizeof(double) * 64);                       memset(block_1, 0, sizeof(double) * 64);
  block_2 = malloc (sizeof(double) * 64);                       memset(block_2, 0, sizeof(double) * 64);
  block_3 = malloc (sizeof(double) * 64);                       memset(block_3, 0, sizeof(double) * 64);

  agg_1 = malloc (sizeof(double) * 64);                         memset(agg_1, 0, sizeof(double) * 64);
  agg_2 = malloc (sizeof(double) * 64);                         memset(agg_2, 0, sizeof(double) * 64);
  agg_3 = malloc (sizeof(double) * 64);                         memset(agg_3, 0, sizeof(double) * 64);
  agg_4 = malloc (sizeof(double) * 64);                         memset(agg_4, 0, sizeof(double) * 64);
  agg_5 = malloc (sizeof(double) * 64);                         memset(agg_5, 0, sizeof(double) * 64);
  agg_6 = malloc (sizeof(double) * 64);                         memset(agg_6, 0, sizeof(double) * 64);
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
