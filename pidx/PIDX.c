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

static int vp = 0;
static int hp = 0;
static double sim_start = 0, sim_end = 0;
static double *write_init_start = 0, *write_init_end = 0;
static double *rst_init_start, *rst_init_end, *rst_start, *rst_end;
static double *hz_init_start, *hz_init_end, *hz_start, *hz_end;
static double *agg_init_start, *agg_init_end, *agg_start, *agg_end;
static double *io_init_start, *io_init_end, *io_start, *io_end;
static double *var_init_start, *var_init_end;
static double *cleanup_start, *cleanup_end;
static double *finalize_start, *finalize_end;
static double *block_rst_init_start, *block_rst_init_end;
static double *block_rst_start, *block_rst_end;
static double *buffer_start, *buffer_end;
static double *compression_init_start, *compression_init_end;
static double *compression_start, *compression_end;
static double *agg_1, *agg_2, *agg_3, *agg_4, *agg_5, *agg_6;
static double *block_1, *block_2, *block_3;

static int caching_state = 0;
static int time_step_caching = 0;
static int hz_caching = 0;
static uint32_t* cached_header_copy;

static PIDX_return_code PIDX_cleanup(PIDX_file file);
static PIDX_return_code PIDX_write(PIDX_file file);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* file_name, int current_time_step);
static PIDX_return_code PIDX_cache_headers(PIDX_file file);

/// PIDX File descriptor (equivalent to the descriptor returned by)
/// POSIX or any other IO framework
struct PIDX_file_descriptor 
{
  idx_dataset idx_ptr;                                  ///< Contains all relevant IDX file info (Blocks per file, samples per block, bitmask, box, file name template... )
  idx_dataset_derived_metadata idx_derived_ptr;         ///< Contains all derieved IDX file info (number of files, files that are ging to be populated)
  
#if PIDX_HAVE_MPI
  MPI_Comm comm;                                        ///< MPI sub-communicator (including all processes per IDX file)
  MPI_Comm global_comm;                                 ///< MPI super-communicator (includes all processes)
#endif
  
  int IDX_WRITE;                                        ///< 1 for write 0 for read
  int IDX_READ;                                         ///< 0 for write 1 for read
  int flags;
  
  PIDX_access access;                                   ///< serial or parallel access
  
  int idx_count[PIDX_MAX_DIMENSIONS];                   ///< Number of idx files in each dimensions
  
  PIDX_header_io_id header_io_id;                       ///< IDX metadata id
#if PIDX_HAVE_MPI
  PIDX_rst_id rst_id;                                   ///< Restructuring phase id
#endif
  PIDX_block_rst_id block_rst_id;                       ///< Block restructuring id (prepration for compression)
  PIDX_hz_encode_id hz_id;                              ///< HZ encoding phase id
  PIDX_compression_id compression_id;                   ///< Compression (lossy and lossless) id
  PIDX_agg_id agg_id;                                   ///< Aggregation phase id
  PIDX_io_id io_id;                                     ///< IO phase id
  
  int local_variable_index;                             ///<
  int local_variable_count;                             ///<
  int variable_pipelining_factor;                       ///<
  
  int write_on_close;                                   ///< HPC Writes
  int one_time_initializations;                         ///<
  
  int debug_rst;                                        ///< Debug restructuring phase (1 or 0), works only on the test application
  int debug_hz;                                         ///< Debug HZ encoding phase (1 or 0), works only on the test application
  
  int perform_block_rst;                                ///< Counter to activate/deactivate (1/0) block restructuring phase
  int perform_hz;                                       ///< Counter to activate/deactivate (1/0) hz encoding phase
  int perform_agg;                                      ///< Counter to activate/deactivate (1/0) aggregation phase
  int perform_io;                                       ///< Counter to activate/deactivate (1/0) I/O phase
  //int perform_compression;                              ///< Counter to activate/deactivate (1/0) compression
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

void PIDX_init_timming_buffers()
{
  write_init_start = malloc (sizeof(double) * 64);              memset(write_init_start, 0, sizeof(double) * 64);
  write_init_end = malloc (sizeof(double) * 64);                memset(write_init_end, 0, sizeof(double) * 64);
  rst_init_start = malloc (sizeof(double) * 64);                memset(rst_init_start, 0, sizeof(double) * 64);
  rst_init_end = malloc (sizeof(double) * 64);                  memset(rst_init_end, 0, sizeof(double) * 64);
  hz_init_start = malloc (sizeof(double) * 64);                 memset(hz_init_start, 0, sizeof(double) * 64);
  hz_init_end = malloc (sizeof(double) * 64);                   memset(hz_init_end, 0, sizeof(double) * 64);
  agg_init_start = malloc (sizeof(double) * 64);                memset(agg_init_start, 0, sizeof(double) * 64);
  agg_init_end = malloc (sizeof(double) * 64);                  memset(agg_init_end, 0, sizeof(double) * 64);
  io_init_start = malloc (sizeof(double) * 64);                 memset(io_init_start, 0, sizeof(double) * 64);
  io_init_end = malloc (sizeof(double) * 64);                   memset(io_init_end, 0, sizeof(double) * 64);
  var_init_start = malloc (sizeof(double) * 64);                memset(var_init_start, 0, sizeof(double) * 64);
  var_init_end = malloc (sizeof(double) * 64);                  memset(var_init_end, 0, sizeof(double) * 64);
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
  
  block_rst_init_start =  malloc (sizeof(double) * 64);         memset(block_rst_init_start, 0, sizeof(double) * 64);
  block_rst_init_end =  malloc (sizeof(double) * 64);           memset(block_rst_init_end, 0, sizeof(double) * 64);
  compression_init_start =  malloc (sizeof(double) * 64);       memset(compression_init_start, 0, sizeof(double) * 64);
  compression_init_end =  malloc (sizeof(double) * 64);         memset(compression_init_end, 0, sizeof(double) * 64);
  block_rst_start =  malloc (sizeof(double) * 64);              memset(block_rst_start, 0, sizeof(double) * 64);
  block_rst_end =  malloc (sizeof(double) * 64);                memset(block_rst_end, 0, sizeof(double) * 64);
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

/// Function to create IDX file descriptor (based on flags and access)
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  PIDX_init_timming_buffers();
  
  sim_start = PIDX_get_time();

  if(flags != PIDX_file_excl && flags != PIDX_file_trunc)
    return PIDX_err_unsupported_flags;
    
  if(flags == PIDX_file_excl)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }
  
  int i; 
  int ret;
  char file_name_skeleton[1024];
  int rank = 0;
  
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
    
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  (*file)->flags = flags;
  (*file)->idx_ptr = (idx_dataset)malloc(sizeof (*((*file)->idx_ptr)));
  memset((*file)->idx_ptr, 0, sizeof (*((*file)->idx_ptr)));
  
  (*file)->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_derived_ptr)));
  memset((*file)->idx_derived_ptr, 0, sizeof (*((*file)->idx_derived_ptr)));
  
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++) 
    (*file)->idx_ptr->global_bounds[i]=65535;

  //initialize logic_to_physic transform to identity
  (*file)->idx_ptr->transform[0]  = 1.0;
  (*file)->idx_ptr->transform[5]  = 1.0;
  (*file)->idx_ptr->transform[10] = 1.0;
  (*file)->idx_ptr->transform[15] = 1.0;
  
  (*file)->idx_ptr->variable_count = -1;
  (*file)->idx_ptr->variable_index_tracker = 0;

  (*file)->access = access_type;
  (*file)->idx_ptr->current_time_step = 0;
  (*file)->idx_derived_ptr->aggregation_factor = 1;
  (*file)->idx_derived_ptr->color = 0;
  (*file)->idx_count[0] = 1;
  (*file)->idx_count[1] = 1;
  (*file)->idx_count[2] = 1;
  (*file)->idx_count[3] = 1;
  (*file)->idx_count[4] = 1;
  
  (*file)->perform_hz = 1;
  (*file)->perform_agg = 1;
  (*file)->perform_io = 1;
  
#if PIDX_HAVE_MPI
  if (access_type->parallel)
    MPI_Comm_rank(access_type->comm, &rank);
#endif
  
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
    (*file)->idx_derived_ptr->fs_block_size = stat_buf.st_blksize;
  }
  
#if PIDX_HAVE_MPI
  MPI_Bcast(&((*file)->idx_derived_ptr->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
  
#if PIDX_HAVE_MPI
  
  int rank_x = 0, rank_y = 0, rank_z = 0, rank_slice = 0;
  int *colors;
  if (access_type->parallel)
  {
    memcpy ((*file)->idx_count, access_type->idx_count, sizeof(int) * /*PIDX_MAX_DIMENSIONS*/3);
    
    if ((*file)->idx_count[0] != 1 || (*file)->idx_count[1] != 1 || (*file)->idx_count[2] != 1 )
    {
      int i = 0, j = 0, k = 0;
      rank_z = rank / (access_type->sub_div[0] * access_type->sub_div[1]);
      rank_slice = rank % (access_type->sub_div[0] * access_type->sub_div[1]);
      rank_y = (rank_slice / access_type->sub_div[0]);
      rank_x = (rank_slice % access_type->sub_div[0]);
      
      colors = malloc(sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      memset(colors, 0, sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      
      for (k = 0; k < (*file)->idx_count[2]; k++)
        for (j = 0; j < (*file)->idx_count[1]; j++)
          for (i = 0; i < (*file)->idx_count[0]; i++)
            colors[((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i] = ((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i;
            
      int index_x = 0, index_y = 0, index_z = 0;
      for (i = 0; i < access_type->sub_div[0]; i = i + (access_type->sub_div[0] / (*file)->idx_count[0]))
      {
        if (rank_x >= i && rank_x < i + (access_type->sub_div[0] / (*file)->idx_count[0]))
        {
          index_x = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[1]; i = i + (access_type->sub_div[1] / (*file)->idx_count[1]))
      {
        if (rank_y >= i && rank_y < i + (access_type->sub_div[1] / (*file)->idx_count[1]))
        {
          index_y = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[2]; i = i + (access_type->sub_div[2] / (*file)->idx_count[2]))
      {
        if (rank_z >= i && rank_z < i + (access_type->sub_div[2] / (*file)->idx_count[2]))
        {
          index_z = i;
          break;
        }
      }
      
      //((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]));
      (*file)->idx_derived_ptr->color = colors[((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]))];
      
      free(colors);
      
      //printf("%d = %d %d %d :: %d\n", rank, rank_x, rank_y, rank_z, (*file)->idx_derived_ptr->color);
      MPI_Comm_split(access_type->comm, (*file)->idx_derived_ptr->color, rank, &((*file)->comm));
      MPI_Comm_dup(access_type->comm, &((*file)->global_comm));
      
    }
    else
      MPI_Comm_dup( access_type->comm , &((*file)->comm));
  }
#endif
  
  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';
  sprintf((*file)->idx_ptr->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_derived_ptr->color);
  
  //printf("[%d = %d %d %d] Color == %d Filename = %s\n", rank, rank_x, rank_y, rank_z, (*file)->idx_derived_ptr->color, (*file)->idx_ptr->filename);
  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->write_on_close = 0; 
  (*file)->one_time_initializations = 0;
  (*file)->idx_ptr->compression_block_size[0] = 1;
  (*file)->idx_ptr->compression_block_size[1] = 1;
  (*file)->idx_ptr->compression_block_size[2] = 1;
  (*file)->idx_ptr->compression_block_size[3] = 1;
  (*file)->idx_ptr->compression_block_size[4] = 1;
  
  (*file)->idx_ptr->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx_derived_ptr->samples_per_block = pow(2, PIDX_default_bits_per_block);
  (*file)->idx_ptr->blocks_per_file = PIDX_default_blocks_per_file;
    
  return PIDX_success;
}

/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  PIDX_init_timming_buffers();
  
  sim_start = PIDX_get_time();
  
  int rank = 0, ret;
  if(flags != PIDX_file_rdonly)
    return PIDX_err_unsupported_flags;
  
  if(flags != PIDX_file_rdonly)
  {
    //struct stat buffer;
    //if (stat(filename, &buffer) != 0)
    return PIDX_err_file;
  }
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
  
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  (*file)->flags = flags;
  
  
  (*file)->idx_ptr = (idx_dataset)malloc(sizeof (*((*file)->idx_ptr)));
  memset((*file)->idx_ptr, 0, sizeof (*((*file)->idx_ptr)));
  
  (*file)->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_derived_ptr)));
  memset((*file)->idx_derived_ptr, 0, sizeof (*((*file)->idx_derived_ptr)));

  //(*file)->idx_ptr->filename = strdup(filename);
  
  (*file)->idx_ptr->variable_count = 0;
  (*file)->idx_ptr->variable_index_tracker = 0;
  
  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
 
  (*file)->idx_ptr->compression_block_size[0] = 1;
  (*file)->idx_ptr->compression_block_size[1] = 1;
  (*file)->idx_ptr->compression_block_size[2] = 1;
  (*file)->idx_ptr->compression_block_size[3] = 1;
  (*file)->idx_ptr->compression_block_size[4] = 1;
  
  (*file)->idx_derived_ptr->aggregation_factor = 1;
  
  (*file)->idx_ptr->current_time_step = 0;
  (*file)->idx_derived_ptr->aggregation_factor = 1;
  (*file)->idx_derived_ptr->color = 0;
  (*file)->idx_count[0] = 1;
  (*file)->idx_count[1] = 1;
  (*file)->idx_count[2] = 1;
  (*file)->idx_count[3] = 1;
  (*file)->idx_count[4] = 1;
  
  (*file)->perform_hz = 1;
  (*file)->perform_agg = 1;
  (*file)->perform_io = 1;
  
#if PIDX_HAVE_MPI
  if (access_type->parallel)
    MPI_Comm_rank(access_type->comm, &rank);
#endif
  
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
    (*file)->idx_derived_ptr->fs_block_size = stat_buf.st_blksize;
  }
  
#if PIDX_HAVE_MPI
  MPI_Bcast(&((*file)->idx_derived_ptr->fs_block_size), 1, MPI_INT, 0, access_type->comm);
#endif
  
#if PIDX_HAVE_MPI
  if (access_type->parallel)
  {
    int rank_x, rank_y, rank_z, rank_slice;
    int *colors;
    memcpy ((*file)->idx_count, access_type->idx_count, sizeof(int) * /*PIDX_MAX_DIMENSIONS*/3);
    
    //printf("[%d] :: %d %d %d\n", rank, (*file)->idx_count[0], (*file)->idx_count[1], (*file)->idx_count[2]);
    
    if ((*file)->idx_count[0] != 1 || (*file)->idx_count[1] != 1 || (*file)->idx_count[2] != 1 )
    {
      int i = 0, j = 0, k = 0;
      rank_z = rank / (access_type->sub_div[0] * access_type->sub_div[1]);
      rank_slice = rank % (access_type->sub_div[0] * access_type->sub_div[1]);
      rank_y = (rank_slice / access_type->sub_div[0]);
      rank_x = (rank_slice % access_type->sub_div[0]);
      
      colors = malloc(sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      memset(colors, 0, sizeof(*colors) * (*file)->idx_count[0] * (*file)->idx_count[1] * (*file)->idx_count[2]);
      
      for (k = 0; k < (*file)->idx_count[2]; k++)
        for (j = 0; j < (*file)->idx_count[1]; j++)
          for (i = 0; i < (*file)->idx_count[0]; i++)
            colors[((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i] = ((*file)->idx_count[0] * (*file)->idx_count[1] * k) + ((*file)->idx_count[0] * j) + i;
            
      int index_x = 0, index_y = 0, index_z = 0;
      for (i = 0; i < access_type->sub_div[0]; i = i + (access_type->sub_div[0] / (*file)->idx_count[0]))
      {
        if (rank_x >= i && rank_x < i + (access_type->sub_div[0] / (*file)->idx_count[0]))
        {
          index_x = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[1]; i = i + (access_type->sub_div[1] / (*file)->idx_count[1]))
      {
        if (rank_y >= i && rank_y < i + (access_type->sub_div[1] / (*file)->idx_count[1]))
        {
          index_y = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[2]; i = i + (access_type->sub_div[2] / (*file)->idx_count[2]))
      {
        if (rank_z >= i && rank_z < i + (access_type->sub_div[2] / (*file)->idx_count[2]))
        {
          index_z = i;
          break;
        }
      }
      
      //((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]));
      (*file)->idx_derived_ptr->color = colors[((*file)->idx_count[0] * (*file)->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_count[2]))) + ((*file)->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_count[0]))];
      
      free(colors);
      
      MPI_Comm_split(access_type->comm, (*file)->idx_derived_ptr->color, rank, &((*file)->comm));
      MPI_Comm_dup(access_type->comm, &((*file)->global_comm));

    }
    else
      MPI_Comm_dup( access_type->comm , &((*file)->comm));
  }
#endif

  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';
  sprintf((*file)->idx_ptr->filename, "%s_%d.idx", file_name_skeleton, (*file)->idx_derived_ptr->color);
  
  (*file)->write_on_close = 0; 
  (*file)->one_time_initializations = 0;
  
  int var = 0, variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];
  
  if (rank == 0)
  {
    FILE *fp = fopen((*file)->idx_ptr->filename, "r");
    while (fgets(line, sizeof (line), fp) != NULL) 
    {
      //printf("%s", line);
      len = strlen(line) - 1;
      if (line[len] == '\n')
	line[len] = 0;
	
      if (strcmp(line, "(version)") == 0)
	fgets(line, sizeof line, fp);
	
      if (strcmp(line, "(box)") == 0) 
      {
	fgets(line, sizeof line, fp);
	len = strlen(line) - 1;
	if (line[len] == '\n')
	  line[len] = 0;

	pch = strtok(line, " ");
	count = 0;
	while (pch != NULL) 
	{
	  if (count % 2 == 1)
	    (*file)->idx_ptr->global_bounds[count / 2] = atoi(pch) + 1;
	  count++;
	  pch = strtok(NULL, " ");
	}
      }
      if (strcmp(line, "(fields)") == 0) 
      {
	fgets(line, sizeof (line), fp);
	len = strlen(line) - 1;
	if (line[len] == '\n')
	  line[len] = 0;
	count = 0;
	variable_counter = 0;
	
	while (strcmp(line, "(version)") != 0 && strcmp(line, "(box)") != 0 && strcmp(line, "(bits)") && strcmp(line, "(bitsperblock)") != 0 && strcmp(line, "(blocksperfile)") != 0 && strcmp(line, "(filename_template)") != 0 && strcmp(line, "(time)") != 0)
	{
	  (*file)->idx_ptr->variable[variable_counter] = malloc(sizeof (*((*file)->idx_ptr->variable[variable_counter])));
	  
	  pch1 = strtok(line, " *+");
	  while (pch1 != NULL)
	  {
	    //printf("");
	    if (count == 0)
	      strcpy((*file)->idx_ptr->variable[variable_counter]->var_name, strdup(pch1));
	      //strcpy((*file)->idx_ptr->variable[variable_counter]->var_name, pch1);

	    if (count == 1)
	      (*file)->idx_ptr->variable[variable_counter]->values_per_sample = atoi(pch1);
	      
	    if (count == 2) 
	    {
	      len = strlen(pch1) - 1;
	      if (pch1[len] == '\n')
		pch1[len] = 0;
	      if (strcmp(pch1, "float64") == 0)
		strcpy((*file)->idx_ptr->variable[variable_counter]->type_name, "float64");
	    }
	    count++;
	    pch1 = strtok(NULL, " *+");
	  }
	  count = 0;

	  fgets(line, sizeof (line), fp);
	  len = strlen(line) - 1;
	  if (line[len] == '\n')
	    line[len] = 0;
	  variable_counter++;
	}
	(*file)->idx_ptr->variable_count = variable_counter;
      }
      if (strcmp(line, "(bits)") == 0)
	fgets(line, sizeof line, fp);

      if (strcmp(line, "(bitsperblock)") == 0) 
      {
	fgets(line, sizeof line, fp);
	len = strlen(line) - 1;
	if (line[len] == '\n')
	  line[len] = 0;
	(*file)->idx_ptr->bits_per_block = atoi(line);
	(*file)->idx_derived_ptr->samples_per_block = pow(2, (*file)->idx_ptr->bits_per_block);
      }
      if (strcmp(line, "(blocksperfile)") == 0) 
      {
	fgets(line, sizeof line, fp);
	len = strlen(line) - 1;
	if (line[len] == '\n')
	  line[len] = 0;
	(*file)->idx_ptr->blocks_per_file= atoi(line);
      }
      if (strcmp(line, "(filename_template)") == 0) 
      {
	fgets(line, sizeof line, fp);
	len = strlen(line) - 1;
	if (line[len] == '\n')
	  line[len] = 0;
      }
      if (strcmp(line, "(time)") == 0)
	fgets(line, sizeof line, fp);
    }
    fclose(fp);
  }

#if PIDX_HAVE_MPI
  MPI_Bcast((*file)->idx_ptr->global_bounds, 5, MPI_LONG_LONG, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->blocks_per_file), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->bits_per_block), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->variable_count), 1, MPI_INT, 0, (*file)->comm);

  (*file)->idx_derived_ptr->samples_per_block = pow(2, (*file)->idx_ptr->bits_per_block);
  
  if(rank != 0)
  {
    for (var = 0; var < (*file)->idx_ptr->variable_count; var++) 
      (*file)->idx_ptr->variable[var] = malloc(sizeof (*((*file)->idx_ptr->variable[var])));
  }
#endif
  
  for (var = 0; var < (*file)->idx_ptr->variable_count; var++) 
  {
#if PIDX_HAVE_MPI
    MPI_Bcast(&((*file)->idx_ptr->variable[var]->values_per_sample), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast((*file)->idx_ptr->variable[var]->var_name, 512, MPI_CHAR, 0, (*file)->comm);
    MPI_Bcast((*file)->idx_ptr->variable[var]->type_name, 512, MPI_CHAR, 0, (*file)->comm);
#endif
    
    (*file)->idx_ptr->variable[var]->patch_count = 0;
    (*file)->idx_ptr->variable[var]->values_per_sample = 1;
    (*file)->idx_ptr->variable[var]->bits_per_value = (int)sizeof(double) * 8;
    strcpy((*file)->idx_ptr->variable[var]->type_name, "1*float64");
  }
  
  return PIDX_success;
}

/// validate dims/blocksize
PIDX_return_code PIDX_validate(PIDX_file file)
{
  int64_t dims;
  int64_t adjusted_global_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_global_bounds[0] = file->idx_ptr->global_bounds[0] / file->idx_ptr->compression_block_size[0];
  adjusted_global_bounds[1] = file->idx_ptr->global_bounds[1] / file->idx_ptr->compression_block_size[1];
  adjusted_global_bounds[2] = file->idx_ptr->global_bounds[2] / file->idx_ptr->compression_block_size[2];
  adjusted_global_bounds[3] = file->idx_ptr->global_bounds[3] / file->idx_ptr->compression_block_size[3];
  adjusted_global_bounds[4] = file->idx_ptr->global_bounds[4] / file->idx_ptr->compression_block_size[4];
  
  //if (PIDX_inner_product(&dims, file->idx_ptr->global_bounds))
  if (PIDX_inner_product(&dims, adjusted_global_bounds))
    return PIDX_err_size;
  if (dims < file->idx_derived_ptr->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx_derived_ptr->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx_ptr->bits_per_block = getNumBits(file->idx_derived_ptr->samples_per_block) - 1;
    //file->idx_ptr->bits_per_block = getNumBits(file->idx_derived_ptr->samples_per_block);
  }
  
  /*
  int reduce_by_sample = 1;
  if (reduce_by_sample == 1)
  {
    file->idx_ptr->bits_per_block = file->idx_ptr->bits_per_block - log2(file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4]);
    
    file->idx_derived_ptr->samples_per_block = pow(2, file->idx_ptr->bits_per_block);
  }
  else
  {
    file->idx_ptr->blocks_per_file = file->idx_ptr->blocks_per_file /  (file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4]);
  }
  */
  //file->idx_derived_ptr->samples_per_block = file->idx_derived_ptr->samples_per_block * 2;
  //file->idx_ptr->bits_per_block = file->idx_ptr->bits_per_block + 1;
  
  if (file->idx_ptr->bits_per_block == 0)
  {
    file->idx_ptr->bits_per_block = 1;
    file->idx_derived_ptr->samples_per_block = 1;
  }
  

  // other validations...
  // TODO
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_access(PIDX_file file, PIDX_access *access)
{
  if (file == NULL)
    return PIDX_err_file;
  
  if (!access)
    return PIDX_err_access;

  (*access) = file->access;
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_dims(PIDX_file file, PIDX_point dims)
{
  if(dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_box;
  
  if(file == NULL)
    return PIDX_err_file;
  
  memcpy(file->idx_ptr->global_bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  /*
  if (file->idx_count[0] == 1 && file->idx_count[1] == 1 && file->idx_count[2] != 1 )
    file->idx_ptr->global_bounds[2] = file->idx_ptr->global_bounds[2] / file->idx_count[2];
  
  else if (file->idx_count[0] != 1 && file->idx_count[1] == 1 && file->idx_count[2] == 1 )
    file->idx_ptr->global_bounds[0] = file->idx_ptr->global_bounds[0] / file->idx_count[0];
  
  else if (file->idx_count[0] == 1 && file->idx_count[1] != 1 && file->idx_count[2] == 1 )
    file->idx_ptr->global_bounds[1] = file->idx_ptr->global_bounds[1] / file->idx_count[1];
  */
  //else if (file->idx_count[0] != 1 && file->idx_count[1] != 1 && file->idx_count[2] != 1 )
  //{
  //TODO: check this what if idx_count is not set here
  file->idx_ptr->global_bounds[0] = file->idx_ptr->global_bounds[0] / file->idx_count[0];
  file->idx_ptr->global_bounds[1] = file->idx_ptr->global_bounds[1] / file->idx_count[1];
  file->idx_ptr->global_bounds[2] = file->idx_ptr->global_bounds[2] / file->idx_count[2];
  //}
  
  return PIDX_validate(file);
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims)
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(dims, file->idx_ptr->global_bounds, (sizeof(int64_t) * PIDX_MAX_DIMENSIONS));
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_aggregation_factor(PIDX_file file, int agg_factor)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_derived_ptr->aggregation_factor = agg_factor;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_aggregation_factor(PIDX_file file, int *agg_factor)
{
  if(!file)
    return PIDX_err_file;
  
  *agg_factor = file->idx_derived_ptr->aggregation_factor;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(file->idx_ptr->transform, transform, (sizeof(double) * 16));
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_transform(PIDX_file file, double transform[16])
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(transform, file->idx_ptr->transform, (sizeof(double) * 16));
  
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
  sprintf(data_set_path, "%s/time%06d.idx", directory_path, current_time_step);
  free(directory_path);
  
  nbits_blocknumber = (file->idx_derived_ptr->maxh - file->idx_ptr->bits_per_block - 1);
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
  strcpy(file->idx_ptr->filename_template, data_set_path);
  for (N = strlen(file->idx_ptr->filename_template) - 1; N >= 0; N--) 
  {
    int ch = file->idx_ptr->filename_template[N];
    file->idx_ptr->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
    strcat(file->idx_ptr->filename_template, "/%01x.bin");
   
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
      strcat(file->idx_ptr->filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12) 
      strcat(file->idx_ptr->filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16) 
      strcat(file->idx_ptr->filename_template, "/%04x.bin"); //no directories, 65536  files
    else 
    {
      while (nbits_blocknumber > 16) 
      {
	strcat(file->idx_ptr->filename_template, "/%02x"); //256 subdirectories
	nbits_blocknumber -= 8;
      }
      strcat(file->idx_ptr->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      assert(nbits_blocknumber <= 0);
    }
  }
  
  free(data_set_path);
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  if(current_time_step < 0)
    return PIDX_err_time;
   
  file->idx_ptr->current_time_step = current_time_step;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  *current_time_step = file->idx_ptr->current_time_step;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_block_size(PIDX_file file, const int bits_per_block)
{
  if(!file)
    return PIDX_err_file;
  
  if(bits_per_block <= 0)
    return PIDX_err_block;
   
  file->idx_ptr->bits_per_block = bits_per_block;
  file->idx_derived_ptr->samples_per_block = pow(2, bits_per_block);
  
  return PIDX_validate(file);
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_block_size(PIDX_file file, int* bits_per_block)
{ 
  if(!file)
    return PIDX_err_file;
  
  *bits_per_block = file->idx_ptr->bits_per_block;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_block_count(PIDX_file file, const int blocks_per_file)
{
  if(!file)
    return PIDX_err_file;
  
  if(blocks_per_file <= 0)
    return PIDX_err_block;
   
  file->idx_ptr->blocks_per_file = blocks_per_file;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_block_count(PIDX_file file, int* blocks_per_file)
{ 
  if(!file)
    return PIDX_err_file;
  
  *blocks_per_file = file->idx_ptr->blocks_per_file;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_count <= 0)
    return PIDX_err_count;
   
  file->idx_ptr->variable_count = variable_count;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count)
{ 
  if(!file)
    return PIDX_err_file;
  
  *variable_count = file->idx_ptr->variable_count;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_next_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  *variable = file->idx_ptr->variable[file->idx_ptr->variable_index_tracker];
  
  file->idx_ptr->variable_index_tracker++;
  file->local_variable_count++;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_variable_create(PIDX_file file, char* variable_name, unsigned int bits_per_sample, PIDX_type type_name, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  if(!variable_name)
    return PIDX_err_name;

  if(bits_per_sample <= 0)
    return PIDX_err_size;
  
  if(!type_name)
    return PIDX_err_type;
  
  //if(file->idx_ptr->variable_count == -1)
    //return PIDX_err_variable;
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker] = malloc(sizeof (*file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]));
  memset(file->idx_ptr->variable[file->idx_ptr->variable_index_tracker], 0, sizeof (*file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]));
  
  //file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->bits_per_sample = bits_per_sample;
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->values_per_sample = 1;
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->bits_per_value = (bits_per_sample/1);
  
  //file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->type_name = strdup(type_name);    
  strcpy(file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->type_name, type_name);    
  //file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->var_name = strdup(variable_name);
  strcpy(file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->var_name, variable_name);
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->patch_count = 0;
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->dump_meta_data_ON = 0;
  
  *variable = file->idx_ptr->variable[file->idx_ptr->variable_index_tracker];
  
  file->idx_ptr->variable_index_tracker++;
  file->local_variable_count++;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_read_next_variable(PIDX_variable variable, PIDX_point offset, PIDX_point dims, void* write_to_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;
  
  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;
  
  if(!dims || dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_count;
  
  variable->patch[variable->patch_count] = malloc(sizeof(*(variable->patch[variable->patch_count])));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  variable->patch[variable->patch_count]->Ndim_box_buffer = write_to_this_buffer;
  
  variable->data_layout = data_layout;
  variable->patch_count = variable->patch_count + 1;
  
  return PIDX_success;
}


/////////////////////////////////////////////////
PIDX_return_code PIDX_append_and_write_variable(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;
  
  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;
  
  if(!dims || dims[0] < 0 || dims[1] < 0 || dims[2] < 0 || dims[3] < 0 || dims[4] < 0)
    return PIDX_err_count;
  
  const void *temp_buffer;
  variable->patch[variable->patch_count] = malloc(sizeof(*(variable->patch[variable->patch_count])));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_size, dims, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  temp_buffer = read_from_this_buffer;
  variable->patch[variable->patch_count]->Ndim_box_buffer = (unsigned char*)temp_buffer;
  
  variable->data_layout = data_layout;
  variable->patch_count = variable->patch_count + 1;
  return PIDX_success; 
}

/////////////////////////////////////////////////
PIDX_return_code populate_idx_dataset(PIDX_file file)
{
  int i, j, counter = 0, file_number = 0;
    int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };
  
  PointND global_bounds_point;
  
  if(file->perform_block_rst == 1)
  {
    if (file->idx_ptr->global_bounds[0] % file->idx_ptr->compression_block_size[0] == 0)
      file->idx_ptr->compressed_global_bounds[0] = (int) file->idx_ptr->global_bounds[0] / file->idx_ptr->compression_block_size[0];
    else
      file->idx_ptr->compressed_global_bounds[0] = (int) (file->idx_ptr->global_bounds[0] / file->idx_ptr->compression_block_size[0]) + 1;
    
    if (file->idx_ptr->global_bounds[1] % file->idx_ptr->compression_block_size[1] == 0)
      file->idx_ptr->compressed_global_bounds[1] = (int) file->idx_ptr->global_bounds[1] / file->idx_ptr->compression_block_size[1];
    else
      file->idx_ptr->compressed_global_bounds[1] = (int) (file->idx_ptr->global_bounds[1] / file->idx_ptr->compression_block_size[1]) + 1;
    
    if (file->idx_ptr->global_bounds[2] % file->idx_ptr->compression_block_size[2] == 0)
      file->idx_ptr->compressed_global_bounds[2] = (int) file->idx_ptr->global_bounds[2] / file->idx_ptr->compression_block_size[2];
    else
      file->idx_ptr->compressed_global_bounds[2] = (int) (file->idx_ptr->global_bounds[2] / file->idx_ptr->compression_block_size[2]) + 1;
    
    if (file->idx_ptr->global_bounds[3] % file->idx_ptr->compression_block_size[3] == 0)
      file->idx_ptr->compressed_global_bounds[3] = (int) file->idx_ptr->global_bounds[3] / file->idx_ptr->compression_block_size[3];
    else
      file->idx_ptr->compressed_global_bounds[3] = (int) (file->idx_ptr->global_bounds[3] / file->idx_ptr->compression_block_size[3]) + 1;
    
    if (file->idx_ptr->global_bounds[4] % file->idx_ptr->compression_block_size[4] == 0)
      file->idx_ptr->compressed_global_bounds[4] = (int) file->idx_ptr->global_bounds[4] / file->idx_ptr->compression_block_size[4];
    else
      file->idx_ptr->compressed_global_bounds[4] = (int) (file->idx_ptr->global_bounds[4] / file->idx_ptr->compression_block_size[4]) + 1;
  }
  else
    memcpy(file->idx_ptr->compressed_global_bounds, file->idx_ptr->global_bounds, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  
  global_bounds_point.x = (int) file->idx_ptr->compressed_global_bounds[0];
  global_bounds_point.y = (int) file->idx_ptr->compressed_global_bounds[1];
  global_bounds_point.z = (int) file->idx_ptr->compressed_global_bounds[2];
  global_bounds_point.u = (int) file->idx_ptr->compressed_global_bounds[3];
  global_bounds_point.v = (int) file->idx_ptr->compressed_global_bounds[4];
  GuessBitmaskPattern(file->idx_ptr->bitSequence, global_bounds_point);
  file->idx_derived_ptr->maxh = strlen(file->idx_ptr->bitSequence);
  
  for (i = 0; i <= file->idx_derived_ptr->maxh; i++)
    file->idx_ptr->bitPattern[i] = RegExBitmaskBit(file->idx_ptr->bitSequence, i);
  
  file->idx_derived_ptr->max_file_count = (getPowerOf2(file->idx_ptr->compressed_global_bounds[0]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[1]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[2]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[3]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[4])) / ((uint64_t) file->idx_derived_ptr->samples_per_block * (uint64_t) file->idx_ptr->blocks_per_file);
  if ((getPowerOf2(file->idx_ptr->compressed_global_bounds[0]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[1]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[2]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[3]) * getPowerOf2(file->idx_ptr->compressed_global_bounds[4])) % ((uint64_t) file->idx_derived_ptr->samples_per_block * (uint64_t) file->idx_ptr->blocks_per_file))
    file->idx_derived_ptr->max_file_count++;
  
  file->idx_derived_ptr->file_bitmap = (int*) malloc(file->idx_derived_ptr->max_file_count * sizeof (int));
  memset(file->idx_derived_ptr->file_bitmap, 0, file->idx_derived_ptr->max_file_count * sizeof (int));
  
#ifdef PIDX_VAR_SLOW_LOOP
  int var, level_count = 1, p, ctr;
  for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
  {
    block_layout* all_patch_local_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(all_patch_local_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    file->idx_ptr->variable[var]->VAR_global_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(file->idx_ptr->variable[var]->VAR_global_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    for (p = 0 ; p < file->idx_ptr->variable[var]->patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
        bounding_box[0][i] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[i];
        bounding_box[1][i] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_size[i] + file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[i];
      }
      
      block_layout* per_patch_local_block_layout = (block_layout*) malloc(sizeof (block_layout));
      createBlockBitmap (bounding_box, file->idx_ptr->blocks_per_file, file->idx_ptr->bits_per_block, file->idx_derived_ptr->maxh, file->idx_ptr->bitPattern, per_patch_local_block_layout);
      
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
      destroyBlockBitmap(per_patch_local_block_layout);
      free(per_patch_local_block_layout);
      per_patch_local_block_layout = 0;
    }
  
    level_count = 1;
    for (i = 1; i < (file->idx_ptr->variable[var]->VAR_global_block_layout->levels); i++) 
    {
#if PIDX_HAVE_MPI
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i], level_count,
                MPI_INT, MPI_BOR, file->comm);
#else
      memcpy(file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
      level_count = level_count * 2;
    }
    
    destroyBlockBitmap(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
    
    ctr = 1;
    for (i = 1; i < (file->idx_ptr->variable[var]->VAR_global_block_layout->levels); i++)
    {
      for (j = 0 ; j < ctr ; j++)
      {
        if(file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i][j] != 0)
          file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_count_array[i]++;
      }    
      
      counter = 0;
      for (j = 0 ; j < ctr ; j++)
      {
        if (file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i][j] != 0)
        {
          file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i][counter] = file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i][j];
          counter++;
        }
      }
      ctr = ctr * 2;
    }
    
    int *temp_file_index = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
    memset(temp_file_index, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));
    
    file->idx_ptr->variable[var]->VAR_blocks_per_file = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
    memset(file->idx_ptr->variable[var]->VAR_blocks_per_file, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));
  
    temp_file_index[0] = 1;
    file->idx_derived_ptr->file_bitmap[0] = 1;
    file->idx_ptr->variable[var]->VAR_blocks_per_file[0] = 1;
    
    for (i = 1; i < file->idx_ptr->variable[var]->VAR_global_block_layout->levels; i++) 
    {
      for (j = 0; j < file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_count_array[i]; j++) 
      {
        file_number = file->idx_ptr->variable[var]->VAR_global_block_layout->hz_block_number_array[i][j] / file->idx_ptr->blocks_per_file;
        file->idx_derived_ptr->file_bitmap[file_number] = 1;
        temp_file_index[file_number] = 1;
        file->idx_ptr->variable[var]->VAR_blocks_per_file[file_number]++;
      }
    }
    
    file->idx_ptr->variable[var]->VAR_existing_file_count = 0;
    for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
      if (temp_file_index[i] == 1)
        file->idx_ptr->variable[var]->VAR_existing_file_count++;
    
    file->idx_ptr->variable[var]->VAR_existing_file_index = (int*) malloc(file->idx_ptr->variable[var]->VAR_existing_file_count * sizeof (int));
    memset(file->idx_ptr->variable[var]->VAR_existing_file_index, 0, file->idx_ptr->variable[var]->VAR_existing_file_count * sizeof (int));
    
    int count = 0;
    for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
      if (temp_file_index[i] == 1)
      {
        file->idx_ptr->variable[var]->VAR_existing_file_index[count] = i;
        count++;
      }
    free(temp_file_index);
  }
#else

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++) 
  {
    bounding_box[0][i] = 0;
    //bounding_box[1][i] = file->idx_ptr->global_bounds[i];
    bounding_box[1][i] = file->idx_ptr->compressed_global_bounds[i];
  }
  
  file->idx_derived_ptr->global_block_layout =  malloc(sizeof (*file->idx_derived_ptr->global_block_layout));
  PIDX_blocks_create_layout(bounding_box, file->idx_ptr->blocks_per_file, file->idx_ptr->bits_per_block, file->idx_derived_ptr->maxh, file->idx_ptr->bitPattern, file->idx_derived_ptr->global_block_layout);
  
  int k = 1;
  for (i = 1; i < (file->idx_derived_ptr->global_block_layout->levels); i++)
  {
    counter = 0;
    for (j = 0 ; j < k ; j++)
    {
      if(file->idx_derived_ptr->global_block_layout->hz_block_number_array[i][j] != 0)
      {
        file->idx_derived_ptr->global_block_layout->hz_block_number_array[i][counter] = file->idx_derived_ptr->global_block_layout->hz_block_number_array[i][j];
        counter++;
      }
    }
    k = k * 2;
  }
    
  int *temp_file_index = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
  memset(temp_file_index, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));
  
  file->idx_derived_ptr->existing_blocks_index_per_file = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
  memset(file->idx_derived_ptr->existing_blocks_index_per_file, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));

  temp_file_index[0] = 1;
  file->idx_derived_ptr->file_bitmap[0] = 1;
  file->idx_derived_ptr->existing_blocks_index_per_file[0] = 1;
  
  for (i = 1; i < file->idx_derived_ptr->global_block_layout->levels; i++) 
  {
    for (j = 0; j < file->idx_derived_ptr->global_block_layout->hz_block_count_array[i]; j++) 
    {
      file_number = file->idx_derived_ptr->global_block_layout->hz_block_number_array[i][j] / file->idx_ptr->blocks_per_file;
      file->idx_derived_ptr->file_bitmap[file_number] = 1;
      temp_file_index[file_number] = 1;
      file->idx_derived_ptr->existing_blocks_index_per_file[file_number]++;
    }
  }
  
  //PIDX_blocks_print_layout(file->idx_derived_ptr->global_block_layout);
  
  file->idx_derived_ptr->existing_file_count = 0;
  for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
    if (temp_file_index[i] == 1)
      file->idx_derived_ptr->existing_file_count++;
  
  file->idx_derived_ptr->existing_file_index = (int*) malloc(file->idx_derived_ptr->existing_file_count * sizeof (int));
  memset(file->idx_derived_ptr->existing_file_index, 0, file->idx_derived_ptr->existing_file_count * sizeof (int));
  
  int count = 0;
  for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
    if (temp_file_index[i] == 1)
    {
      file->idx_derived_ptr->existing_file_index[count] = i;
      count++;
    }
  free(temp_file_index); 
#endif
    
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_read(PIDX_file file)
{
  if (file->local_variable_index == file->idx_ptr->variable_count)
    return PIDX_success;
    
  int j = 0, p, var = 0;
  int rank = 0;
  int var_used_in_binary_file, total_header_size;
  //static int header_io = 0;
#if PIDX_HAVE_MPI
  MPI_Comm_rank(file->comm, &rank);
#endif
  
  populate_idx_dataset(file);
  
  /// Initialization ONLY ONCE per IDX file
  if(file->one_time_initializations == 0)
  {
    PIDX_file_initialize_time_step(file, file->idx_ptr->filename, file->idx_ptr->current_time_step);
    var_used_in_binary_file = (file->idx_ptr->variable_count < 0) ? 64 : file->idx_ptr->variable_count;
    
    total_header_size = (10 + (10 * file->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * var_used_in_binary_file;
    file->idx_derived_ptr->start_fs_block = total_header_size / file->idx_derived_ptr->fs_block_size;
    if (total_header_size % file->idx_derived_ptr->fs_block_size)
      file->idx_derived_ptr->start_fs_block++;
    
    file->one_time_initializations = 1;
  }
  
  /// HPC version of writes
  /// All variables are written in one go,
  /// basically NO intermittent PIDX_flush involved.
  /// Variable Count maybe Set or not.
  if (file->local_variable_count == file->idx_ptr->variable_index_tracker && file->write_on_close == 1)
  {
    write_init_start[hp] = PIDX_get_time();
    if (file->idx_ptr->variable_count == -1)
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    else if (file->idx_ptr->variable_count < file->idx_ptr->variable_index_tracker)
    {
      if(rank == 0)
        fprintf(stderr, "(Warning !!!!) Variable count set at %d while attempting to write %d variables\n", file->idx_ptr->variable_count, file->idx_ptr->variable_index_tracker);
      
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    }
    write_init_end[hp] = PIDX_get_time();
    hp++;
  }
  
  /// Variable count is set (and all variables have been initialized) 
  /// Used with flush (but all header is written once)
  /// Takes less memory and is IO efficient way of writing data
  //else if (file->idx_ptr->variable_count == file->idx_ptr->variable_index_tracker && header_io == 0)
  //{
  //  file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_count);
  //  PIDX_header_io_set_communicator(file->header_io_id, file->comm);
  //  PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
  //  PIDX_header_io_file_create(file->header_io_id);
  //  PIDX_header_io_finalize(file->header_io_id);
  //  header_io = 1;
  //}
  
  /// All variables are not initialized (variable count might or might not be set)
  /// Used with flush (header is also written in groups)
  /// Takes less memory but is not IO efficient way of writing data
  else
  {
    write_init_start[hp] = PIDX_get_time();
    //file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    write_init_end[hp] = PIDX_get_time();
    hp++;
  }
  
  if (file->idx_count[0] != 1 || file->idx_count[1] != 1 || file->idx_count[2] != 1 )
  {
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j < file->idx_ptr->compressed_global_bounds[0] * file->idx_count[0]; j = j + (file->idx_ptr->compressed_global_bounds[0]))
        {          
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] < (j + file->idx_ptr->compressed_global_bounds[0]))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] - (j);
            break;
          }
        }
        
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j <= file->idx_ptr->compressed_global_bounds[1] * file->idx_count[1]; j = j + (file->idx_ptr->compressed_global_bounds[1]))
        {
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] < j + (file->idx_ptr->compressed_global_bounds[1] ))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] - (j);
            break;
          }
        }
        
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j < file->idx_ptr->compressed_global_bounds[2] * file->idx_count[2]; j = j + (file->idx_ptr->compressed_global_bounds[2]))
        {
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] < j + (file->idx_ptr->compressed_global_bounds[2] /* file->idx_count */))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] - (j);
            break;
          }
        }
  }
    

  int do_agg = 1;
  int local_do_rst = 0, global_do_rst = 0;
  int start_index = 0, end_index = 0;
  
  file->variable_pipelining_factor = 31;
  for (start_index = file->local_variable_index; start_index < file->local_variable_index + file->local_variable_count; start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    ///----------------------------------IO init start-------------------------------------------------///
    io_init_start[vp] = PIDX_get_time();
    
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
    PIDX_io_set_communicator(file->io_id, file->comm);
#endif
    
    io_init_end[vp] = PIDX_get_time();
    ///----------------------------------IO init end---------------------------------------------------///
    
    
    ///-----------------------------------AGG init start-----------------------------------------------///
    agg_init_start[vp] = PIDX_get_time();                                            
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
      PIDX_agg_set_communicator(file->agg_id, file->comm);
#endif
    }
    agg_init_end[vp] = PIDX_get_time();
    ///-----------------------------------AGG init end-------------------------------------------------///
    
    
    ///------------------------------------HZ init start-----------------------------------------------///
    hz_init_start[vp] = PIDX_get_time();                                             
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
    PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
#endif
    hz_init_end[vp] = PIDX_get_time();
    ///------------------------------------HZ init end-------------------------------------------------///
    
    
#if PIDX_HAVE_MPI
    ///----------------------------------- RST init start----------------------------------------------///
    rst_init_start[vp] = PIDX_get_time();
    if (file->idx_ptr->variable[start_index]->patch_count == 1)
      local_do_rst = 1;
    
    MPI_Allreduce(&local_do_rst, &global_do_rst, 1, MPI_INT, MPI_LOR, file->comm);
    global_do_rst = 1;
    if(global_do_rst == 1)
    {
      file->rst_id = PIDX_rst_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_rst_set_communicator(file->rst_id, file->comm);
    }
    
    rst_init_end[vp] = PIDX_get_time();
    ///----------------------------------- RST init end------------------------------------------------///
#endif
    
    ///------------------------------Var buffer init start---------------------------------------------///
    var_init_start[vp] = PIDX_get_time();
#ifdef PIDX_VAR_SLOW_LOOP
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = 0;
#if PIDX_HAVE_MPI
      if(global_do_rst == 1)
        file->idx_ptr->variable[var]->patch_group_count = PIDX_rst_attach_restructuring_box(file->rst_id, 0, NULL);
      else
        file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
#else
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
#endif
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
#else
    file->idx_ptr->variable[start_index]->patch_group_count = 0;
#if PIDX_HAVE_MPI
    if (global_do_rst == 1)
      file->idx_ptr->variable[start_index]->patch_group_count = PIDX_rst_attach_restructuring_box(file->rst_id, 0, NULL);
    else
      file->idx_ptr->variable[start_index]->patch_group_count = file->idx_ptr->variable[start_index]->patch_count;
#else
      file->idx_ptr->variable[start_index]->patch_group_count = file->idx_ptr->variable[start_index]->patch_count;
#endif
    
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[start_index]->patch_group_count;
    
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
#endif
    var_init_end[vp] = PIDX_get_time();
    ///------------------------------Var buffer init end--------------------------------------------------///
    
    
    ///------------------------------------ALL buffer start time-------------------------------------------------///
    buffer_start[vp] = PIDX_get_time();

#if PIDX_HAVE_MPI
    if(global_do_rst == 0)
    {
      for (var = start_index; var <= end_index; var++)
      {
        for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        {
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count = 1;
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type = 0;
          
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count);
          for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
          {
            file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j])));
            memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
          }
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
      }
    }
    else
      PIDX_rst_buf_create(file->rst_id);
#else
    
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count = 1;
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type = 1;
        
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count);
        for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
        {
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j])));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
        }
        memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
      }
    }
    
#endif

    PIDX_hz_encode_buf_create(file->hz_id);
    
    file->idx_derived_ptr->agg_buffer = malloc(sizeof(*file->idx_derived_ptr->agg_buffer));
    PIDX_agg_buf_create(file->agg_id);  
    
    buffer_end[vp] = PIDX_get_time();
    ///--------------------------------------ALL buffer end time---------------------------------------------------///
 
    
    ///------------------------------------IO start time---------------------------------------------------///
    io_start[vp] = PIDX_get_time();
    if (do_agg == 1)
      PIDX_io_aggregated_read(file->io_id);
    else
      PIDX_io_per_process_read(file->io_id);
    io_end[vp] = PIDX_get_time();
    ///---------------------------------------IO end time---------------------------------------------------///
    

    ///---------------------------------------Agg start time--------------------------------------------------///
    agg_start[vp] = PIDX_get_time();  
    if (do_agg == 1)
    {
      PIDX_agg_read(file->agg_id);
      
      PIDX_agg_buf_destroy(file->agg_id);
      free(file->idx_derived_ptr->agg_buffer);
    }
    agg_end[vp] = PIDX_get_time();
    ///---------------------------------------IO end time---------------------------------------------------///
    
    
    ///-------------------------------------HZ start time---------------------------------------------------///
    hz_start[vp] = PIDX_get_time();
    PIDX_hz_encode_read(file->hz_id);
    //HELPER_Hz_encode(file->hz_id);
    PIDX_hz_encode_buf_destroy(file->hz_id);
    hz_end[vp] = PIDX_get_time();
    ///------------------------------------HZ end time------------------------------------------------------///
    
    
    ///---------------------------------------Agg start time--------------------------------------------------///
    rst_start[vp] = PIDX_get_time();  
#if PIDX_HAVE_MPI
    PIDX_rst_read(file->rst_id);
    //HELPER_rst(file->rst_id);
#endif
    rst_end[vp] = PIDX_get_time();  
    ///-----------------------------------------Agg end time--------------------------------------------------///
    
    
    ///--------------------------------------cleanup start time---------------------------------------------///
    cleanup_start[vp] = PIDX_get_time();
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        if(global_do_rst == 0)
        {
          for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
            free(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]);
          free(file->idx_ptr->variable[var]->patch_group_ptr[p]->box);
        }
        free(file->idx_ptr->variable[var]->patch_group_ptr[p]);
      }
      free(file->idx_ptr->variable[var]->patch_group_ptr);
    }
    cleanup_end[vp] = PIDX_get_time();
    ///-------------------------------------cleanup end time------------------------------------------------///
    
    
    ///-------------------------------------finalize start time---------------------------------------------///
    finalize_start[vp] = PIDX_get_time();
    PIDX_io_finalize(file->io_id);
    if(do_agg == 1)
      PIDX_agg_finalize(file->agg_id);
    PIDX_hz_encode_finalize(file->hz_id);
    
#if PIDX_HAVE_MPI
    if(global_do_rst == 1)
      PIDX_rst_finalize(file->rst_id);
#endif
      
    finalize_end[vp] = PIDX_get_time();
    ///------------------------------------finalize end time------------------------------------------------///
    
    vp++;
  }

  return PIDX_success;
}

///
PIDX_return_code PIDX_hz_encoding_caching_ON()
{
  hz_caching = 1;
  return PIDX_success;
}

///
PIDX_return_code PIDX_hz_encoding_caching_OFF()
{
  PIDX_hz_encode_delete_cache_buffers();
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_time_step_caching_ON()
{
  caching_state = 1;
  time_step_caching = 1;
  
  return PIDX_success;
}


/////////////////////////////////////////////////
PIDX_return_code PIDX_time_step_caching_OFF()
{
  //free(cached_header_copy);
  //cached_header_copy = 0;
  
  return PIDX_success;
}


/////////////////////////////////////////////////
PIDX_return_code PIDX_enable_compression(PIDX_file file, int compression)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_ptr->enable_compression = compression;
  
  return PIDX_success;
}


PIDX_return_code PIDX_enable_block_restructuring(PIDX_file file, int brst)
{
  if(!file)
    return PIDX_err_file;
  
  file->perform_block_rst = brst;
  
  return PIDX_success;
}


PIDX_return_code PIDX_enable_hz(PIDX_file file, int hz)
{
  if(!file)
    return PIDX_err_file;
  
  file->perform_hz = hz;
  
  return PIDX_success;
}


PIDX_return_code PIDX_enable_agg(PIDX_file file, int agg)
{
  if(!file)
    return PIDX_err_file;
  
  file->perform_agg = agg;
  
  return PIDX_success;
}


PIDX_return_code PIDX_enable_io(PIDX_file file, int io)
{
  if(!file)
    return PIDX_err_file;
  
  file->perform_io = io;
  
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
  file->idx_derived_ptr->dump_agg_info = dump_agg_info;
  strncpy(filename_skeleton, file->idx_ptr->filename, strlen(file->idx_ptr->filename) - 4);
  filename_skeleton[strlen(file->idx_ptr->filename) - 4] = '\0';
  sprintf(file->idx_derived_ptr->agg_dump_dir_name, "%s_agg_dump", filename_skeleton);
  
  return PIDX_success;
}


PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz)
{
  if(!file)
    return PIDX_err_file;
  
  file->debug_hz = debug_hz;
  
  return PIDX_success;
}


/////////////////////////////////////////////////
static PIDX_return_code PIDX_cache_headers(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  int i = 0, j = 0, k = 0;
  off_t data_offset = 0, base_offset = 0;
  int /* empty_blocks = 0, */ block_negative_offset = 0;
  int all_scalars = 1;
  int var_used_in_binary_file, total_header_size;
  int64_t total_compression_block_size = file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4];
  
  var_used_in_binary_file = (file->idx_ptr->variable_count < 0) ? 64 : file->idx_ptr->variable_count;
    
  total_header_size = (10 + (10 * file->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * var_used_in_binary_file;
  file->idx_derived_ptr->start_fs_block = total_header_size / file->idx_derived_ptr->fs_block_size;
  if (total_header_size % file->idx_derived_ptr->fs_block_size)
    file->idx_derived_ptr->start_fs_block++;
  
  cached_header_copy = (uint32_t*)malloc(total_header_size);
  memset(cached_header_copy, 0, total_header_size);
  /*
  for (i = 0; i < file->idx_ptr->variable_count; i++)
  {
    if (file->idx_ptr->variable[i]->values_per_sample != 1)
    {
      all_scalars = 0;
      break;
    }
  }
  
  */
  //all_scalars = 0;
#ifdef PIDX_VAR_SLOW_LOOP
  for (i = 0; i < file->idx_ptr->blocks_per_file; i++)
  {
    //empty_blocks = find_block_negative_offset(file->idx_ptr->blocks_per_file, ((file->idx_ptr->variable[file->local_variable_index]->VAR_blocks_per_file[file->idx_derived_ptr->agg_buffer->file_number] - 1) + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_ptr->variable[file->local_variable_index]->VAR_global_block_layout);  
    if (is_block_present((i + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_ptr->variable[file->local_variable_index]->VAR_global_block_layout))
    {
      block_negative_offset = find_block_negative_offset(file->idx_ptr->blocks_per_file, (i + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_ptr->variable[file->local_variable_index]->VAR_global_block_layout);
      
      for (j = 0; j < file->idx_ptr->variable_count; j++)
      {
        base_offset = 0;
        if (all_scalars == 0)
        {
          for (k = 0; k < j; k++)
            base_offset = base_offset + (file->idx_ptr->variable[file->local_variable_index]->VAR_blocks_per_file[file->idx_derived_ptr->agg_buffer->file_number] /*- empty_blocks*/) * (file->idx_ptr->variable[k]->bits_per_value / 8) * total_compression_block_size * file->idx_derived_ptr->samples_per_block * file->idx_ptr->variable[k]->values_per_sample;
        }
        else
          base_offset =  j * (file->idx_ptr->variable[file->local_variable_index]->VAR_blocks_per_file[file->idx_derived_ptr->agg_buffer->file_number] /*- empty_blocks*/) * (file->idx_ptr->variable[file->local_variable_index]->bits_per_value / 8) * total_compression_block_size * file->idx_derived_ptr->samples_per_block * file->idx_ptr->variable[file->local_variable_index]->values_per_sample;
        
        data_offset = (((i) - block_negative_offset) * file->idx_derived_ptr->samples_per_block) * (file->idx_ptr->variable[j]->bits_per_value / 8) * total_compression_block_size * file->idx_ptr->variable[j]->values_per_sample;
        data_offset = base_offset + data_offset + file->idx_derived_ptr->start_fs_block * file->idx_derived_ptr->fs_block_size;
        
        cached_header_copy[12 + ((i + (file->idx_ptr->blocks_per_file * j))*10)] = htonl(data_offset);
        cached_header_copy[14 + ((i + (file->idx_ptr->blocks_per_file * j))*10)] = htonl(file->idx_derived_ptr->samples_per_block * (file->idx_ptr->variable[j]->bits_per_value / 8) * total_compression_block_size * file->idx_ptr->variable[j]->values_per_sample);  
      }
    }
  }
#else
  for (i = 0; i < file->idx_ptr->blocks_per_file; i++)
  {
    //empty_blocks = find_block_negative_offset(file->idx_ptr->blocks_per_file, ((file->idx_derived_ptr->existing_blocks_index_per_file[file->idx_derived_ptr->agg_buffer->file_number] - 1) + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_derived_ptr->global_block_layout);
    if (PIDX_blocks_is_block_present((i + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_derived_ptr->global_block_layout))
    {
      block_negative_offset = PIDX_blocks_find_negative_offset(file->idx_ptr->blocks_per_file, (i + (file->idx_ptr->blocks_per_file * file->idx_derived_ptr->agg_buffer->file_number)), file->idx_derived_ptr->global_block_layout);
      
      for (j = 0; j < file->idx_ptr->variable_count; j++)
      {
        base_offset = 0;
        if (all_scalars == 0)
        {
          for (k = 0; k < j; k++)
            base_offset = base_offset + (file->idx_derived_ptr->existing_blocks_index_per_file[file->idx_derived_ptr->agg_buffer->file_number] /*- empty_blocks*/) * (file->idx_ptr->variable[k]->bits_per_value / 8) * total_compression_block_size * file->idx_derived_ptr->samples_per_block * file->idx_ptr->variable[k]->values_per_sample;
        }
        else
          base_offset =  j * (file->idx_derived_ptr->existing_blocks_index_per_file[file->idx_derived_ptr->agg_buffer->file_number] /*- empty_blocks*/) * (file->idx_ptr->variable[file->local_variable_index]->bits_per_value / 8) * total_compression_block_size * file->idx_derived_ptr->samples_per_block * file->idx_ptr->variable[file->local_variable_index]->values_per_sample;
        
        data_offset = (((i) - block_negative_offset) * file->idx_derived_ptr->samples_per_block) * (file->idx_ptr->variable[j]->bits_per_value / 8) * total_compression_block_size * file->idx_ptr->variable[j]->values_per_sample;
        
        //if (j == 1)
        //printf("%d Block %d: [BO %ld DO %ld = %ld] Offset %ld Count %d\n", j, i, 
          //     base_offset, data_offset, (base_offset + data_offset), data_offset, (file->idx_derived_ptr->samples_per_block * (file->idx_ptr->variable[j]->bits_per_value / 8) * file->idx_ptr->variable[j]->values_per_sample));
        
        data_offset = base_offset + data_offset + file->idx_derived_ptr->start_fs_block * file->idx_derived_ptr->fs_block_size;
        //printf("[%d %d]: (%d) O %ld (%ld) C %d\n", j, i, file->idx_derived_ptr->existing_blocks_index_per_file[file->idx_derived_ptr->agg_buffer->file_number], data_offset, base_offset, (file->idx_derived_ptr->samples_per_block * (file->idx_ptr->variable[j]->bits_per_value / 8) * file->idx_ptr->variable[j]->values_per_sample));
        cached_header_copy[12 + ((i + (file->idx_ptr->blocks_per_file * j))*10)] = htonl(data_offset);
        cached_header_copy[14 + ((i + (file->idx_ptr->blocks_per_file * j))*10)] = htonl(file->idx_derived_ptr->samples_per_block * (file->idx_ptr->variable[j]->bits_per_value / 8) * total_compression_block_size * file->idx_ptr->variable[j]->values_per_sample);  
      }
      //printf("\n");
    }
  }
#endif
  return PIDX_success;
}


PIDX_return_code PIDX_set_compression_type(PIDX_file file, int compression_type)
{
  if(compression_type != 1)
    return PIDX_err_unsupported_compression_type;
  
  if(file == NULL)
    return PIDX_err_file;
  
  file->idx_ptr->compression_type = compression_type;
  
  return PIDX_success;
}


PIDX_return_code PIDX_get_compression_type(PIDX_file file, int *compression_type)
{
  if(file == NULL)
    return PIDX_err_file;
  
  *compression_type = file->idx_ptr->compression_type;
  
  return PIDX_success;
}


PIDX_return_code PIDX_set_compression_block_size(PIDX_file file, PIDX_point compression_block_size)
{
  if(compression_block_size[0] < 0 || compression_block_size[1] < 0 || compression_block_size[2] < 0 || compression_block_size[3] < 0 || compression_block_size[4] < 0)
    return PIDX_err_box;
  
  if(file == NULL)
    return PIDX_err_file;
  
  memcpy(file->idx_ptr->compression_block_size, compression_block_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  return PIDX_validate(file);
}


PIDX_return_code PIDX_get_compression_block_size(PIDX_file file, PIDX_point compression_block_size)
{
  if(file == NULL)
    return PIDX_err_file;
  
  memcpy(compression_block_size, file->idx_ptr->compression_block_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  
  return PIDX_success;
}

/////////////////////////////////////////////////
static PIDX_return_code PIDX_write(PIDX_file file)
{
  if (file->local_variable_index == file->idx_ptr->variable_count)
    return PIDX_success;
  
#if defined(BGL) || defined(BGP) || defined(BGQ)
  identity(file->comm);
#endif
  
#if 1
  int j = 0, p, var = 0;
  int rank = 0, nprocs = 1;
  int var_used_in_binary_file, total_header_size;
  //static int header_io = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
#endif
  
  populate_idx_dataset(file);
  /// Initialization ONLY ONCE per IDX file
  if(file->one_time_initializations == 0)
  {
    PIDX_file_initialize_time_step(file, file->idx_ptr->filename, file->idx_ptr->current_time_step);
    var_used_in_binary_file = (file->idx_ptr->variable_count < 0) ? 64 : file->idx_ptr->variable_count;
    
    total_header_size = (10 + (10 * file->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * var_used_in_binary_file;
    file->idx_derived_ptr->start_fs_block = total_header_size / file->idx_derived_ptr->fs_block_size;
    if (total_header_size % file->idx_derived_ptr->fs_block_size)
      file->idx_derived_ptr->start_fs_block++;
    
    file->one_time_initializations = 1;
  }
  
  /// HPC version of writes
  /// All variables are written in one go,
  /// basically NO intermittent PIDX_flush involved.
  /// Variable Count maybe Set or not.
  if (file->local_variable_count == file->idx_ptr->variable_index_tracker && file->write_on_close == 1)
  {
    write_init_start[hp] = PIDX_get_time();
    if (file->idx_ptr->variable_count == -1)
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    else if (file->idx_ptr->variable_count < file->idx_ptr->variable_index_tracker)
    {
      if(rank == 0)
	fprintf(stderr, "(Warning !!!!) Variable count set at %d while attempting to write %d variables\n", file->idx_ptr->variable_count, file->idx_ptr->variable_index_tracker);
      
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    }
    if (file->perform_io == 1)
    {
      file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_count);
#if PIDX_HAVE_MPI
      PIDX_header_io_set_communicator(file->header_io_id, file->comm);
#endif
      PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
      PIDX_header_io_file_create(file->header_io_id);
      //PIDX_header_io_file_write(file->header_io_id);
      PIDX_header_io_finalize(file->header_io_id);
    }
    write_init_end[hp] = PIDX_get_time();
    hp++;
  }
  
  /// Variable count is set (and all variables have been initialized) 
  /// Used with flush (but all header is written once)
  /// Takes less memory and is IO efficient way of writing data
  //else if (file->idx_ptr->variable_count == file->idx_ptr->variable_index_tracker && header_io == 0)
  //{
  //  file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_count);
  //  PIDX_header_io_set_communicator(file->header_io_id, file->comm);
  //  PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
  //  PIDX_header_io_file_create(file->header_io_id);
  //  PIDX_header_io_finalize(file->header_io_id);
  //  header_io = 1;
  //}
  
  /// All variables are not initialized (variable count might or might not be set)
  /// Used with flush (header is also written in groups)
  /// Takes less memory but is not IO efficient way of writing data
  else
  {
    write_init_start[hp] = PIDX_get_time();
    //file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    
    file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_index_tracker);
#if PIDX_HAVE_MPI
    PIDX_header_io_set_communicator(file->header_io_id, file->comm);
#endif
    
    if (hp == 0)
      PIDX_header_io_file_create(file->header_io_id);
    
    if (file->idx_ptr->variable_count == -1 || (file->idx_ptr->variable_count == file->idx_ptr->variable_index_tracker))
    {
      PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
      //PIDX_header_io_file_write(file->header_io_id);
    }
        
    PIDX_header_io_finalize(file->header_io_id);
    
    write_init_end[hp] = PIDX_get_time();
    hp++;
  }
  if (file->idx_count[0] != 1 || file->idx_count[1] != 1 || file->idx_count[2] != 1 )
  {
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j < file->idx_ptr->compressed_global_bounds[0] * file->idx_count[0]; j = j + (file->idx_ptr->compressed_global_bounds[0]))
        {          
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] < (j + file->idx_ptr->compressed_global_bounds[0]))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[0] - 
            (j);
            break;
          }
        }
        
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j <= file->idx_ptr->compressed_global_bounds[1] * file->idx_count[1]; j = j + (file->idx_ptr->compressed_global_bounds[1]))
        {
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] < j + (file->idx_ptr->compressed_global_bounds[1] ))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[1] - 
            (j);
            break;
          }
        }
        
    for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
        for (j = 0; j < file->idx_ptr->compressed_global_bounds[2] * file->idx_count[2]; j = j + (file->idx_ptr->compressed_global_bounds[2]))
        {
          if (file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] >= j && file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] < j + (file->idx_ptr->compressed_global_bounds[2] /* file->idx_count */))
          {
            file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[2] - 
            (j);
            break;
          }
        }
  }
  
  int i;
  int do_agg = 1;
  int local_do_rst = 0, global_do_rst = 0;
  int start_index = 0, end_index = 0;
  file->variable_pipelining_factor = 31;
  
  for (start_index = file->local_variable_index; start_index < file->local_variable_index + file->local_variable_count; start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    ///----------------------------------- RST init start----------------------------------------------///
    rst_init_start[vp] = PIDX_get_time();
#if PIDX_HAVE_MPI
    if (file->idx_ptr->variable[start_index]->patch_count == 1)
      local_do_rst = 1;
    
    MPI_Allreduce(&local_do_rst, &global_do_rst, 1, MPI_INT, MPI_LOR, file->comm);
    if(global_do_rst == 1)
    {
      file->rst_id = PIDX_rst_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_rst_set_communicator(file->rst_id, file->comm);
    }
#endif
    rst_init_end[vp] = PIDX_get_time();
    ///----------------------------------- RST init end------------------------------------------------///
    
    
    ///----------------------------BLOCK restructure init start ---------------------------------------///
    block_rst_init_start[vp] = PIDX_get_time();                                    
    if(file->perform_block_rst == 1)
    {
      file->block_rst_id = PIDX_block_rst_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
      PIDX_block_rst_set_communicator(file->block_rst_id, file->comm);
#endif
    }
    block_rst_init_end[vp] = PIDX_get_time();
    ///----------------------------BLOCK restructure init end -----------------------------------------///
    
    
    ///------------------------------------HZ init start-----------------------------------------------///
    hz_init_start[vp] = PIDX_get_time();                                             
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
    PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
#endif
    hz_init_end[vp] = PIDX_get_time();
    ///------------------------------------HZ init end-------------------------------------------------///
    
#if 0
    ///------------------------------Compression init start -------------------------------------------///
    compression_init_start[vp] = PIDX_get_time();                                    
    if(file->perform_compression == 1)
    {
      file->compression_id = PIDX_compression_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
      PIDX_compression_set_communicator(file->compression_id, file->comm);
#endif
    }
    compression_init_end[vp] = PIDX_get_time();
    ///------------------------------Compression init end---------------------------------------------///
#endif
    
    ///-----------------------------------AGG init start-----------------------------------------------///
    agg_init_start[vp] = PIDX_get_time();                                            
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
      PIDX_agg_set_communicator(file->agg_id, file->comm);
#endif
    }
    agg_init_end[vp] = PIDX_get_time();
    ///-----------------------------------AGG init end-------------------------------------------------///
    
    ///----------------------------------IO init start-------------------------------------------------///
    io_init_start[vp] = PIDX_get_time();
    if (file->perform_io == 1)
    {
      file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
#if PIDX_HAVE_MPI
      PIDX_io_set_communicator(file->io_id, file->comm);
#endif
    }
    io_init_end[vp] = PIDX_get_time();
    ///----------------------------------IO init end---------------------------------------------------///
    
    
    ///------------------------------Var buffer init start---------------------------------------------///
    var_init_start[vp] = PIDX_get_time();
    
#if PIDX_HAVE_MPI
    file->idx_ptr->variable[start_index]->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_ptr->variable[start_index]->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    file->idx_ptr->variable[start_index]->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
    memset(file->idx_ptr->variable[start_index]->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

    MPI_Allgather(file->idx_ptr->variable[start_index]->patch[0]->Ndim_box_offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_ptr->variable[start_index]->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
    
    MPI_Allgather(file->idx_ptr->variable[start_index]->patch[0]->Ndim_box_size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_ptr->variable[start_index]->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
#endif
    
#ifdef PIDX_VAR_SLOW_LOOP
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = 0;

#if PIDX_HAVE_MPI
      if(global_do_rst == 1)
        file->idx_ptr->variable[var]->patch_group_count = PIDX_rst_attach_restructuring_box(file->rst_id, 0, NULL);
      else
        file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
#else
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
#endif
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
#else
#if PIDX_HAVE_MPI
    file->idx_ptr->variable[start_index]->patch_group_count = 0;
    if (global_do_rst == 1)
      file->idx_ptr->variable[start_index]->patch_group_count = PIDX_rst_attach_restructuring_box(file->rst_id, 0, NULL);
    else
      file->idx_ptr->variable[start_index]->patch_group_count = file->idx_ptr->variable[start_index]->patch_count;
#else
    file->idx_ptr->variable[start_index]->patch_group_count = file->idx_ptr->variable[start_index]->patch_count;
#endif
    
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[start_index]->patch_group_count;
    
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
#endif
    var_init_end[vp] = PIDX_get_time();
    ///------------------------------Var buffer init end--------------------------------------------------///
    
    ///------------------------------------RST start time-------------------------------------------------///
    rst_start[vp] = PIDX_get_time();
#if PIDX_HAVE_MPI
    if(global_do_rst == 0)
    {
      for (var = start_index; var <= end_index; var++)
      {
        for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        {
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count = 1;
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type = 0;
          
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count);
          for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
          {
            file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j])));
            memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
          }
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
      }
    }
    else
    {
      PIDX_rst_buf_create(file->rst_id);
      PIDX_rst_write(file->rst_id);
      if(global_do_rst == 1 && file->debug_rst == 1)
        HELPER_rst(file->rst_id);
    }
#else
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count = 1;
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type = 0;
        
        file->idx_ptr->variable[var]->patch_group_ptr[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count);
        for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
        {
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j])));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
        }
        memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
      }
    }
#endif
    rst_end[vp] = PIDX_get_time();
    ///--------------------------------------RST end time---------------------------------------------------///

    
    ///----------------------------BLOCK restructure start time---------------------------------------------///
    block_rst_start[vp] = PIDX_get_time();
    
    for (var = start_index; var <= end_index; var++)
    {
      //TODO
      if ( file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4] != 1)
      {
        file->idx_ptr->variable[var]->post_rst_block = malloc(sizeof(*file->idx_ptr->variable[var]->post_rst_block) * file->idx_ptr->variable[var]->patch_group_count);
        memset(file->idx_ptr->variable[var]->post_rst_block, 0, sizeof(*file->idx_ptr->variable[var]->post_rst_block) * file->idx_ptr->variable[var]->patch_group_count);
        for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        {
          file->idx_ptr->variable[var]->post_rst_block[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p])));
          memset(file->idx_ptr->variable[var]->post_rst_block[p], 0, sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p])));
          
          file->idx_ptr->variable[var]->post_rst_block[p]->box_count = 1;
          file->idx_ptr->variable[var]->post_rst_block[p]->box_group_type = file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type;
          
          file->idx_ptr->variable[var]->post_rst_block[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p]->box)) * file->idx_ptr->variable[var]->post_rst_block[p]->box_count);
          
          file->idx_ptr->variable[var]->post_rst_block[p]->box[0] = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p]->box[0])));
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->box[0]->Ndim_box_offset, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->box[0]->Ndim_box_size, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, PIDX_MAX_DIMENSIONS *       sizeof(int64_t));
          //file->idx_ptr->variable[var]->post_rst_block[p]->box[0]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
          
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
      }
      else
      {
        file->idx_ptr->variable[var]->post_rst_block = malloc(sizeof(*file->idx_ptr->variable[var]->post_rst_block) * file->idx_ptr->variable[var]->patch_group_count);
        memset(file->idx_ptr->variable[var]->post_rst_block, 0, sizeof(*file->idx_ptr->variable[var]->post_rst_block) * file->idx_ptr->variable[var]->patch_group_count);
        for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        {
          file->idx_ptr->variable[var]->post_rst_block[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p])));
          memset(file->idx_ptr->variable[var]->post_rst_block[p], 0, sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p])));
          
          file->idx_ptr->variable[var]->post_rst_block[p]->box_count = file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count;
          file->idx_ptr->variable[var]->post_rst_block[p]->box_group_type = file->idx_ptr->variable[var]->patch_group_ptr[p]->box_group_type;
          
          file->idx_ptr->variable[var]->post_rst_block[p]->box = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p]->box)) * file->idx_ptr->variable[var]->post_rst_block[p]->box_count);
          
          for(j = 0; j < file->idx_ptr->variable[var]->post_rst_block[p]->box_count; j++)
          {
            file->idx_ptr->variable[var]->post_rst_block[p]->box[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->post_rst_block[p]->box[j])));
            memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->box[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->box[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            file->idx_ptr->variable[var]->post_rst_block[p]->box[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]->Ndim_box_buffer;
          }
          
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->enclosing_box_offset, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(file->idx_ptr->variable[var]->post_rst_block[p]->enclosing_box_size, file->idx_ptr->variable[var]->patch_group_ptr[p]->enclosing_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
        
      }
    }
    
    if (file->perform_block_rst == 1)
    {
      block_1[vp] =  PIDX_get_time();
      if ( file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4] != 1)
        PIDX_block_rst_prepare(file->block_rst_id);
      //PIDX_block_rst_compress(file->block_rst_id);
      block_2[vp] =  PIDX_get_time();
      
#if PIDX_HAVE_MPI
      if ( file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4] != 1)
      {
        if(global_do_rst == 1)
          PIDX_rst_buf_destroy(file->rst_id);
        else
        {
          for(var = start_index; var <= end_index; var++)
          {
            for(i = 0; i < file->idx_ptr->variable[start_index]->patch_group_count; i++)
            {
              for(j = 0; j < file->idx_ptr->variable[start_index]->patch_group_ptr[i]->box_count; j++)
              {                
                free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer);
                file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer = 0;
                                
                free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]);
                file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j] = 0;
              }
              free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box);
              file->idx_ptr->variable[var]->patch_group_ptr[i]->box = 0;
            }
          }
        }
      }
#else
      for(var = start_index; var <= end_index; var++)
      {
        for(i = 0; i < file->idx_ptr->variable[start_index]->patch_group_count; i++)
        {
          for(j = 0; j < file->idx_ptr->variable[start_index]->patch_group_ptr[i]->box_count; j++)
          {
            free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer);
            file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer = 0;
            
            free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]);
            file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j] = 0;
          }
          free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box);
          file->idx_ptr->variable[var]->patch_group_ptr[i]->box = 0;
        }
      }
#endif
      block_3[vp] =  PIDX_get_time();
    }
    block_rst_end[vp] = PIDX_get_time();
    ///----------------------------BLOCK restructure end time-----------------------------------------------///
    

    ///-------------------------------------HZ start time---------------------------------------------------///
    hz_start[vp] = PIDX_get_time();
    PIDX_hz_encode_buf_create(file->hz_id);
    
    if (file->perform_hz == 1)
      PIDX_hz_encode_write(file->hz_id);
    
    if(global_do_rst == 1 && file->debug_hz == 1)
      HELPER_Hz_encode(file->hz_id);
    
    
    if ( file->idx_ptr->compression_block_size[0] * file->idx_ptr->compression_block_size[1] * file->idx_ptr->compression_block_size[2] * file->idx_ptr->compression_block_size[3] * file->idx_ptr->compression_block_size[4] != 1)
    {
      for (var = start_index; var <= end_index; var++)
      {
        for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
        {
          
          free(file->idx_ptr->variable[var]->post_rst_block[p]->box[0]->Ndim_box_buffer);
          file->idx_ptr->variable[var]->post_rst_block[p]->box[0]->Ndim_box_buffer = 0;
          
          free(file->idx_ptr->variable[var]->post_rst_block[p]->box[0]);
          file->idx_ptr->variable[var]->post_rst_block[p]->box[0] = 0;
          
          free(file->idx_ptr->variable[var]->post_rst_block[p]->box);
          file->idx_ptr->variable[var]->post_rst_block[p]->box = 0;
                    
          free(file->idx_ptr->variable[var]->post_rst_block[p]);
          file->idx_ptr->variable[var]->post_rst_block[p] = 0;
        }
        free(file->idx_ptr->variable[var]->post_rst_block);
        file->idx_ptr->variable[var]->post_rst_block = 0;
      }
    }
    else
    {
      if(global_do_rst == 1)
        PIDX_rst_buf_destroy(file->rst_id);
      else
      {
        for(var = start_index; var <= end_index; var++)
        {
          for(i = 0; i < file->idx_ptr->variable[start_index]->patch_group_count; i++)
          {
            for(j = 0; j < file->idx_ptr->variable[start_index]->patch_group_ptr[i]->box_count; j++)
            {                
              free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer);
              file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]->Ndim_box_buffer = 0;
                              
              free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j]);
              file->idx_ptr->variable[var]->patch_group_ptr[i]->box[j] = 0;
            }
            free(file->idx_ptr->variable[var]->patch_group_ptr[i]->box);
            file->idx_ptr->variable[var]->patch_group_ptr[i]->box = 0;
          }
        }
      }
    }
    
    hz_end[vp] = PIDX_get_time();
    ///------------------------------------HZ end time------------------------------------------------------///
    
#if 0
    ///---------------------------------Compression start time----------------------------------------------///
    compression_start[vp] = PIDX_get_time();
    if(file->perform_compression == 1)
    {
      PIDX_compression_prepare(file->compression_id);
      PIDX_compression_compress(file->compression_id);
      PIDX_compression_finalize(file->compression_id);
    }
    compression_end[vp] = PIDX_get_time();
    ///----------------------------------Compression end time-----------------------------------------------///
#endif
    
    ///------------------------------------Agg start time---------------------------------------------------///
    agg_start[vp] = PIDX_get_time();
    if(do_agg == 1)
    {
      agg_1[vp] = PIDX_get_time();
      file->idx_derived_ptr->agg_buffer = malloc(sizeof(*file->idx_derived_ptr->agg_buffer));
      
      int p;
      file->idx_derived_ptr->agg_level_start = malloc(sizeof(*file->idx_derived_ptr->agg_level_start) * file->idx_ptr->variable[start_index]->patch_group_count);
      memset(file->idx_derived_ptr->agg_level_start, 0, sizeof(*file->idx_derived_ptr->agg_level_start) * file->idx_ptr->variable[start_index]->patch_group_count);
      file->idx_derived_ptr->agg_level_end = malloc(sizeof(*file->idx_derived_ptr->agg_level_end) * file->idx_ptr->variable[start_index]->patch_group_count);
      memset(file->idx_derived_ptr->agg_level_end, 0, sizeof(*file->idx_derived_ptr->agg_level_end) * file->idx_ptr->variable[start_index]->patch_group_count);
      
      for (p = 0; p < file->idx_ptr->variable[start_index]->patch_group_count; p++)
      {
        file->idx_derived_ptr->agg_level_start[p] = malloc(sizeof(*file->idx_derived_ptr->agg_level_start[p]) * (end_index - start_index + 1));
        memset(file->idx_derived_ptr->agg_level_start[p], 0, sizeof(*file->idx_derived_ptr->agg_level_start[p]) * (end_index - start_index + 1));
        file->idx_derived_ptr->agg_level_end[p] = malloc(sizeof(*file->idx_derived_ptr->agg_level_end[p]) * (end_index - start_index + 1));
        memset(file->idx_derived_ptr->agg_level_end[p], 0, sizeof(*file->idx_derived_ptr->agg_level_end[p]) * (end_index - start_index + 1));
        
        for(var = start_index; var <= end_index; var++)
        {
          file->idx_derived_ptr->agg_level_start[p][var] = malloc(sizeof(*file->idx_derived_ptr->agg_level_start[p][var]) * (file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_to - file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_from));
          memset(file->idx_derived_ptr->agg_level_start[p][var], 0, sizeof(*file->idx_derived_ptr->agg_level_start[p][var]) * (file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_to - file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_from));
          file->idx_derived_ptr->agg_level_end[p][var] = malloc(sizeof(*file->idx_derived_ptr->agg_level_end[p][var]) * (file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_to - file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_from));
          memset(file->idx_derived_ptr->agg_level_end[p][var], 0, sizeof(*file->idx_derived_ptr->agg_level_end[p][var]) * (file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_to - file->idx_ptr->variable[start_index]->HZ_patch[p]->HZ_level_from));
        }
      }
      
      agg_2[vp] = PIDX_get_time();
      PIDX_agg_buf_create(file->agg_id);
      agg_3[vp] = PIDX_get_time();
      if (file->perform_agg == 1)
        PIDX_agg_write(file->agg_id);
      agg_4[vp] = PIDX_get_time();
      PIDX_hz_encode_buf_destroy(file->hz_id);
      agg_5[vp] = PIDX_get_time();
      
            
      /// Initialization ONLY ONCE for all TIME STEPS (caching across time)
      if (caching_state == 1 && file->idx_derived_ptr->agg_buffer->var_number == 0 && file->idx_derived_ptr->agg_buffer->sample_number == 0)
      {
        if (file->idx_ptr->enable_compression == 0)
        {
          PIDX_cache_headers(file);
          caching_state = 0;
        }
      }
      agg_6[vp] = PIDX_get_time();
    }
    agg_end[vp] = PIDX_get_time();
    ///---------------------------------------Agg end time---------------------------------------------------///
    

    ///---------------------------------------IO start time--------------------------------------------------///
    io_start[vp] = PIDX_get_time();
    if(do_agg == 1)
    {
      if (file->perform_io == 1)
      {
        if (time_step_caching == 1)
        {
          if (file->idx_ptr->enable_compression == 0)
            PIDX_io_cached_data(cached_header_copy);
        }
      
        PIDX_io_aggregated_write(file->io_id);
      }
      PIDX_agg_buf_destroy(file->agg_id);
      free(file->idx_derived_ptr->agg_buffer);
    }
    else
    {
      if (file->perform_io == 1)
        PIDX_io_per_process_write(file->io_id);
    }
    io_end[vp] = PIDX_get_time();
    ///---------------------------------------IO end time---------------------------------------------------///
    
    
    ///--------------------------------------cleanup start time---------------------------------------------///
    cleanup_start[vp] = PIDX_get_time();
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        if(global_do_rst == 0)
        {
          for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->box_count; j++)
            free(file->idx_ptr->variable[var]->patch_group_ptr[p]->box[j]);
          free(file->idx_ptr->variable[var]->patch_group_ptr[p]->box);
        }
        free(file->idx_ptr->variable[var]->patch_group_ptr[p]);
      }
      free(file->idx_ptr->variable[var]->patch_group_ptr);
    }
    
#if PIDX_HAVE_MPI
    free(file->idx_ptr->variable[start_index]->rank_r_offset);
    free(file->idx_ptr->variable[start_index]->rank_r_count);
#endif
    
    cleanup_end[vp] = PIDX_get_time();
    ///-------------------------------------cleanup end time------------------------------------------------///
    
    
    ///-------------------------------------finalize start time---------------------------------------------///
    finalize_start[vp] = PIDX_get_time();
    PIDX_io_finalize(file->io_id);
    if(do_agg == 1)
      PIDX_agg_finalize(file->agg_id);
    PIDX_hz_encode_finalize(file->hz_id);
    PIDX_block_rst_finalize(file->block_rst_id);
#if PIDX_HAVE_MPI
    if(global_do_rst == 1)
      PIDX_rst_finalize(file->rst_id);
#endif
    
    finalize_end[vp] = PIDX_get_time();
    ///------------------------------------finalize end time------------------------------------------------///
    
    vp++;
  }
#endif
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_flush(PIDX_file file)
{
  if (file->flags == PIDX_file_trunc)
    PIDX_write(file);
  else if (file->flags == PIDX_file_rdonly)
    PIDX_read(file);
  else
    return PIDX_err_unsupported_flags; //unsupoorted flags
    
  PIDX_cleanup(file);
  return PIDX_success;
}

/////////////////////////////////////////////////
static PIDX_return_code PIDX_cleanup(PIDX_file file)
{
  int i, p;
  for (i = file->local_variable_index; i < file->local_variable_index + file->local_variable_count; i++) 
  {
#ifdef PIDX_VAR_SLOW_LOOP
    free(file->idx_ptr->variable[i]->VAR_existing_file_index);
    file->idx_ptr->variable[i]->VAR_existing_file_index = 0;
#else
    free(file->idx_derived_ptr->existing_file_index);
    file->idx_derived_ptr->existing_file_index = 0;
#endif
    
    for(p = 0; p < file->idx_ptr->variable[i]->patch_count; p++)
    {
      //free(file->idx_ptr->variable[i]->HZ_patch[p]);
      //file->idx_ptr->variable[i]->HZ_patch[p] = 0;
    
      free(file->idx_ptr->variable[i]->patch[p]);
      file->idx_ptr->variable[i]->patch[p] = 0;
    }
  }
  
  file->local_variable_index = file->idx_ptr->variable_index_tracker;//file->idx_ptr->variable_count;
  file->local_variable_count = 0;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_close(PIDX_file file) 
{
  file->write_on_close = 1;
  PIDX_flush(file);
  
  sim_end = PIDX_get_time();
  
  int p = 0, i = 0;
  double total_time = sim_end - sim_start;
  double max_time = total_time;
  int sample_sum = 0, var = 0, rank = 0, nprocs = 1;
  
#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
#endif

  if (file->idx_count[0] *  file->idx_count[1] * file->idx_count[2] == 1)
  {
    if (max_time == total_time)
    {
      if (file->IDX_WRITE == 1 && file->IDX_READ == 0)
        fprintf(stdout, "\n------------- WRITE -------------\n");

      if (file->IDX_WRITE == 0 && file->IDX_READ == 1)
        fprintf(stdout, "\n------------- READ -------------\n");
        
      for (var = 0; var < file->idx_ptr->variable_count; var++)
        sample_sum = sample_sum + file->idx_ptr->variable[var]->values_per_sample;
      
      int64_t total_data = file->idx_ptr->global_bounds[0] * file->idx_ptr->global_bounds[1] * file->idx_ptr->global_bounds[2] * file->idx_ptr->global_bounds[3] * file->idx_ptr->global_bounds[4] * sample_sum * 8;
      fprintf(stdout, "\n==========================================================================================================\n");
      fprintf(stdout, "[%d] Time step %d File name %s\n", rank, file->idx_ptr->current_time_step, file->idx_ptr->filename);
      fprintf(stdout, "Cores %d Global Data %lld %lld %lld Variables %d IDX count %d = %d x %d x %d\n", nprocs, (long long) file->idx_ptr->global_bounds[0], (long long) file->idx_ptr->global_bounds[1], (long long) file->idx_ptr->global_bounds[2], file->idx_ptr->variable_count, file->idx_count[0] * file->idx_count[1] * file->idx_count[2], file->idx_count[0], file->idx_count[1], file->idx_count[2]);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d Aggregation Factor %d Aggregator Count %d\n", file->idx_ptr->blocks_per_file, file->idx_ptr->bits_per_block, file->idx_derived_ptr->existing_file_count, file->idx_derived_ptr->aggregation_factor, file->idx_ptr->variable_count * file->idx_derived_ptr->existing_file_count * file->idx_derived_ptr->aggregation_factor);
      fprintf(stdout, "Blocks Restructuring Size %d %d %d %d %d\n", (int)file->idx_ptr->compression_block_size[0], (int)file->idx_ptr->compression_block_size[1], (int)file->idx_ptr->compression_block_size[2], (int)file->idx_ptr->compression_block_size[3], (int)file->idx_ptr->compression_block_size[4]);
      fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", max_time, (float) total_data / (1000 * 1000 * max_time));
      fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
      //printf("File creation time %f\n", write_init_end - write_init_start);
      
      for (var = 0; var < hp; var++)
        fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
      
      for (var = 0; var < vp; var++)
      {
        fprintf(stdout, "------------------------------------------------VG %d (START)----------------------------------------------\n", var);
        fprintf(stdout, "Buffer init time %f\n", (var_init_end[var] - var_init_start[var]));
        fprintf(stdout, "Init time  [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_init_end[var] - rst_init_start[var]), (block_rst_init_end[var] - block_rst_init_start[var]), (hz_init_end[var] - hz_init_start[var]), (agg_init_end[var] - agg_init_start[var]), (io_init_end[var] - io_init_start[var]), (rst_init_end[var] - rst_init_start[var]) + (block_rst_init_end[var] - block_rst_init_start[var]) + (hz_init_end[var] - hz_init_start[var]) + (agg_init_end[var] - agg_init_start[var]) + (io_init_end[var] - io_init_start[var]));
        
        fprintf(stdout, "Write time [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_end[var] - rst_start[var]), (block_rst_end[var] - block_rst_start[var]), (hz_end[var] - hz_start[var]), (agg_end[var] - agg_start[var]), (io_end[var] - io_start[var]), (rst_end[var] - rst_start[var]) + (block_rst_end[var] - block_rst_start[var]) + (hz_end[var] - hz_start[var]) + (agg_end[var] - agg_start[var]) + (io_end[var] - io_start[var]));
        
        fprintf(stdout, "Block Restructuring time %f = %f + %f\n", (block_rst_end[var] - block_rst_start[var]), (block_2[var] - block_1[var]), (block_3[var] - block_2[var]));
        
        fprintf(stdout, "Agg time %f = %f + %f + %f + %f + %f\n", (agg_end[var] - agg_start[var]), (agg_2[var] - agg_1[var]), (agg_3[var] - agg_2[var]), (agg_4[var] - agg_3[var]), (agg_5[var] - agg_4[var]), (agg_6[var] - agg_5[var]));
        
        fprintf(stdout, "Cleanup time %f\n", cleanup_end[var] - cleanup_start[var]);
        fprintf(stdout, "-------------------------------------------------VG %d (END)-----------------------------------------------\n", var);
      }
      
      /*
      double total_agg_time = 0, all_time = 0;
      for (p = 0; p < file->idx_ptr->variable[0]->patch_group_count; p++)
        for (var = 0; var < file->idx_ptr->variable_count; var++)
          for (i = file->idx_ptr->variable[var]->HZ_patch[p]->HZ_level_from; i < file->idx_ptr->variable[var]->HZ_patch[p]->HZ_level_to; i++)
          {
            printf("Agg Time [Patch %d Var %d Level %d] = %f\n", p, var, i, (file->idx_derived_ptr->agg_level_end[p][var][i] - file->idx_derived_ptr->agg_level_start[p][var][i]));
            total_agg_time = total_agg_time + (file->idx_derived_ptr->agg_level_end[p][var][i] - file->idx_derived_ptr->agg_level_start[p][var][i]);
          }
      
      all_time = total_agg_time + (file->idx_derived_ptr->win_time_end - file->idx_derived_ptr->win_time_start) + (file->idx_derived_ptr->win_free_time_end - file->idx_derived_ptr->win_free_time_start);
      printf("Total Agg Time %f = [Network + Win_Create + Win_free] %f + %f + %f\n", all_time, total_agg_time, (file->idx_derived_ptr->win_time_end - file->idx_derived_ptr->win_time_start), (file->idx_derived_ptr->win_free_time_end - file->idx_derived_ptr->win_free_time_start));
      */
      fprintf(stdout, "==========================================================================================================\n");
    }
  }
  else
  {
    
    for (var = 0; var < file->idx_ptr->variable_count; var++)
      sample_sum = sample_sum + file->idx_ptr->variable[var]->values_per_sample;
    
    int64_t total_data = file->idx_ptr->global_bounds[0] * file->idx_ptr->global_bounds[1] * file->idx_ptr->global_bounds[2] * file->idx_ptr->global_bounds[3] * file->idx_ptr->global_bounds[4] * sample_sum * 8 * file->idx_count[0] * file->idx_count[1] * file->idx_count[2];
    if (max_time == total_time)
      fprintf(stdout, "[EXTRA INFO] %d %s Time %f Seconds Throughput %f MB/sec\n", file->idx_ptr->current_time_step, file->idx_ptr->filename, max_time, (float) total_data / (1000 * 1000 * max_time));
    
    int global_rank = 0, global_nprocs = 1;
    double global_max_time = 0;
    
#if PIDX_HAVE_MPI
    MPI_Allreduce(&max_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, file->global_comm);
    MPI_Comm_rank(file->global_comm, &global_rank);
    MPI_Comm_size(file->global_comm, &global_nprocs);
#endif
    
    if (global_max_time == total_time)
    {
      if (file->IDX_WRITE == 1 && file->IDX_READ == 0)
        fprintf(stdout, "\n------------- WRITE -------------\n");

      if (file->IDX_WRITE == 0 && file->IDX_READ == 1)
        fprintf(stdout, "\n------------- READ -------------\n");
      
      fprintf(stdout, "\n==========================================================================================================\n");
      fprintf(stdout, "[%d] Combined Time step %d File name %s\n", global_rank, file->idx_ptr->current_time_step, file->idx_ptr->filename);
      fprintf(stdout, "Cores %d Global Data %lld %lld %lld Variables %d IDX Count %d = %d x %d x %d\n", global_nprocs, (long long) file->idx_ptr->global_bounds[0] * file->idx_count[0], (long long) file->idx_ptr->global_bounds[1] * file->idx_count[1], (long long) file->idx_ptr->global_bounds[2] * file->idx_count[2], file->idx_ptr->variable_count, file->idx_count[0] * file->idx_count[1] * file->idx_count[2], file->idx_count[0], file->idx_count[1], file->idx_count[2]);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d Aggregation Factor %d Aggregator Count %d\n", file->idx_ptr->blocks_per_file, file->idx_ptr->bits_per_block, file->idx_derived_ptr->existing_file_count, file->idx_derived_ptr->aggregation_factor, file->idx_ptr->variable_count * file->idx_derived_ptr->existing_file_count * file->idx_derived_ptr->aggregation_factor );
      fprintf(stdout, "Blocks Restructuring Size %d %d %d %d %d\n", (int)file->idx_ptr->compression_block_size[0], (int)file->idx_ptr->compression_block_size[1], (int)file->idx_ptr->compression_block_size[2], (int)file->idx_ptr->compression_block_size[3], (int)file->idx_ptr->compression_block_size[4]);
      fprintf(stdout, "Time Taken: %f Seconds Throughput %f MB/sec\n", global_max_time, (float) total_data / (1000 * 1000 * global_max_time));
      fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
      //printf("File creation time %f\n", write_init_end - write_init_start);
      
      for (var = 0; var < hp; var++)
        fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
      
      for (var = 0; var < vp; var++)
      {
        fprintf(stdout, "------------------------------------------------VG %d (START)----------------------------------------------\n", var);
        fprintf(stdout, "Buffer init time %f\n", (var_init_end[var] - var_init_start[var]));
        
        fprintf(stdout, "Init time  [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_init_end[var] - rst_init_start[var]), (block_rst_init_end[var] - block_rst_init_start[var]), (hz_init_end[var] - hz_init_start[var]), (agg_init_end[var] - agg_init_start[var]), (io_init_end[var] - io_init_start[var]), (rst_init_end[var] - rst_init_start[var]) + (block_rst_init_end[var] - block_rst_init_start[var]) + (hz_init_end[var] - hz_init_start[var]) + (agg_init_end[var] - agg_init_start[var]) + (io_init_end[var] - io_init_start[var]));
        
        fprintf(stdout, "Write time [RST + BRST + HZ + AGG + IO] %f + %f + %f + %f + %f = %f\n", (rst_end[var] - rst_start[var]), (block_rst_end[var] - block_rst_start[var]), (hz_end[var] - hz_start[var]), (agg_end[var] - agg_start[var]), (io_end[var] - io_start[var]), (rst_end[var] - rst_start[var]) + (block_rst_end[var] - block_rst_start[var]) + (hz_end[var] - hz_start[var]) + (agg_end[var] - agg_start[var]) + (io_end[var] - io_start[var]));
        
        fprintf(stdout, "Block Restructuring time %f = %f + %f\n", (block_rst_end[var] - block_rst_start[var]), (block_2[var] - block_1[var]), (block_3[var] - block_2[var]));
        
        fprintf(stdout, "Agg time %f = %f + %f + %f + %f = %f\n", (agg_end[var] - agg_start[var]), (agg_2[var] - agg_1[var]), (agg_3[var] - agg_2[var]), (agg_4[var] - agg_3[var]), (agg_5[var] - agg_4[var]), (agg_6[var] - agg_5[var]));
        
        fprintf(stdout, "Cleanup time %f\n", cleanup_end[var] - cleanup_start[var]);
        fprintf(stdout, "--------------------------------------------------VG %d (END)-----------------------------------------------\n", var);
      }
      fprintf(stdout, "==========================================================================================================\n");
    }
    
  }
  
#if 1
  
  for (p = 0; p < file->idx_ptr->variable[0]->patch_group_count; p++)
  {
    for(var = 0; var < file->idx_ptr->variable_count; var++)
    {
      free(file->idx_derived_ptr->agg_level_start[p][var]);
      free(file->idx_derived_ptr->agg_level_end[p][var]);
    }
    free(file->idx_derived_ptr->agg_level_start[p]);
    free(file->idx_derived_ptr->agg_level_end[p]);

  }
  free(file->idx_derived_ptr->agg_level_start);
  free(file->idx_derived_ptr->agg_level_end);
  
#ifdef PIDX_VAR_SLOW_LOOP
  for (i = 0; i < file->idx_ptr->variable_count; i++) 
  {
    free(file->idx_ptr->variable[i]->VAR_blocks_per_file);
    file->idx_ptr->variable[i]->VAR_blocks_per_file = 0;
    
    destroyBlockBitmap(file->idx_ptr->variable[i]->VAR_global_block_layout);
    free(file->idx_ptr->variable[i]->VAR_global_block_layout);
    file->idx_ptr->variable[i]->VAR_global_block_layout = 0;
    
    for (p = 0; p < file->idx_ptr->variable[id->start_var_index]->patch_group_count; p++)
    {
      free(file->idx_ptr->variable[i]->HZ_patch[p]);
      file->idx_ptr->variable[i]->HZ_patch[p] = 0;
    }
  }
#else

  for (i = 0; i < file->idx_ptr->variable_count; i++)
  {
    for (p = 0; p < file->idx_ptr->variable[i]->patch_group_count; p++)
    {
      free(file->idx_ptr->variable[i]->HZ_patch[p]);
      file->idx_ptr->variable[i]->HZ_patch[p] = 0;
    }
  }
  
  free(file->idx_derived_ptr->existing_blocks_index_per_file);
  file->idx_derived_ptr->existing_blocks_index_per_file = 0;
  
  PIDX_blocks_free_layout(file->idx_derived_ptr->global_block_layout);
  free(file->idx_derived_ptr->global_block_layout);
  file->idx_derived_ptr->global_block_layout = 0;
#endif
  for (i = 0; i < 1024; i++)
  {
    free(file->idx_ptr->variable[i]);  
    file->idx_ptr->variable[i] = 0; 
  }
  
  file->idx_ptr->variable_count = 0;
  
  //free(file->idx_ptr->global_bounds);         file->idx_ptr->global_bounds = 0;
  free(file->idx_ptr);                        file->idx_ptr = 0;
  free(file->idx_derived_ptr->file_bitmap);   file->idx_derived_ptr->file_bitmap = 0;
  free(file->idx_derived_ptr);                file->idx_derived_ptr = 0;

#endif
  
#if PIDX_HAVE_MPI
  MPI_Comm_free(&(file->comm));
  if (file->idx_count[0] * file->idx_count[1] * file->idx_count[2] != 1)
    MPI_Comm_free(&(file->global_comm));
#endif
  
  free(file);
  
  vp = 0;
  hp = 0;
  
  free(write_init_start);               write_init_start        = 0;
  free(write_init_end);                 write_init_end          = 0;
  free(rst_init_start);                 rst_init_start          = 0;
  free(rst_init_end);                   rst_init_end            = 0;
  free(hz_init_start);                  hz_init_start           = 0;
  free(hz_init_end);                    hz_init_end             = 0;
  free(agg_init_start);                 agg_init_start          = 0;
  free(rst_start);                      rst_start               = 0;
  free(rst_end);                        rst_end                 = 0;
  free(hz_start);                       hz_start                = 0;
  free(hz_end);                         hz_end                  = 0;
  free(agg_init_end);                   agg_init_end            = 0;
  free(agg_start);                      agg_start               = 0;
  free(agg_end);                        agg_end                 = 0;
  free(io_start);                       io_start                = 0;
  free(io_end);                         io_end                  = 0;
  free(io_init_start);                  io_init_start           = 0;
  free(io_init_end);                    io_init_end             = 0;
  free(var_init_start);                 var_init_start          = 0;
  free(var_init_end);                   var_init_end            = 0;
  free(cleanup_start);                  cleanup_start           = 0;
  free(cleanup_end);                    cleanup_end             = 0;
  free(finalize_start);                 finalize_start          = 0;
  free(finalize_end);                   finalize_end            = 0;
  free(buffer_start);                   buffer_start            = 0;
  free(buffer_end);                     buffer_end              = 0;
  free(block_rst_init_start);           block_rst_init_start    = 0;
  free(block_rst_init_end);             block_rst_init_end      = 0;
  free(block_rst_start);                block_rst_start         = 0;
  free(block_rst_end);                  block_rst_end           = 0;
  free(compression_init_start);         compression_init_start  = 0;
  free(compression_init_end);           compression_init_end    = 0;
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

/////////////////////////////////////////////////
PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable)
{
  if(!variable)
    return PIDX_err_variable;
  
  variable->dump_meta_data_ON = 1;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable)
{
  if(!variable)
    return PIDX_err_variable;
  
  variable->dump_meta_data_ON = 0;
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_bits_per_sample(PIDX_type type_name, unsigned int bits_per_sample)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_box_count(PIDX_file file, int* box_count)
{
  return PIDX_err_not_implemented;
}
/////////////////////////////////////////////////
PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point dims)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_box_count_with_rank(PIDX_file file, int MPI_rank, int* box_count)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, int MPI_rank, PIDX_point offset, PIDX_point dims)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int variable_index)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_index <= 0)
    return PIDX_err_size;
  
  if(file->idx_ptr->variable_index_tracker >= file->idx_ptr->variable_count)
    return PIDX_err_count;
  
  file->idx_ptr->variable_index_tracker = variable_index;
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable* variable)
{
  if(!file)
    return PIDX_err_file;
  
  if(file->idx_ptr->variable_index_tracker >= file->idx_ptr->variable_count)
    return PIDX_err_count;
  
  (*variable) = file->idx_ptr->variable[file->idx_ptr->variable_index_tracker];
  
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_read_variable(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_write_variable(PIDX_variable variable, PIDX_point offset, PIDX_point dims, const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}

/////////////////////////////////////////////////
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
  MPI_Allgather(variable->patch[0]->Ndim_box_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
  MPI_Allgather(variable->patch[0]->Ndim_box_size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
#else
  memcpy(rank_r_offset, variable->patch[0]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
  memcpy(rank_r_count, variable->patch[0]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
#endif
  
  meta_data_file = fopen("Box_info.txt", "w");
  if (!meta_data_file) 
    return PIDX_err_name;
  
  //for(i = 0; i < nprocs; i++)
  //  fprintf(meta_data_file, "[%d]: %lld %lld %lld %lld %lld - %lld %lld %lld %lld %lld\n", i, rank_r_offset[PIDX_MAX_DIMENSIONS * i + 0], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 1], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 2], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 3], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 4], rank_r_count[PIDX_MAX_DIMENSIONS * i + 0], rank_r_count[PIDX_MAX_DIMENSIONS * i + 1], rank_r_count[PIDX_MAX_DIMENSIONS * i + 2], rank_r_count[PIDX_MAX_DIMENSIONS * i + 3], rank_r_count[PIDX_MAX_DIMENSIONS * i + 4]);
  
  fclose(meta_data_file);
  
  return PIDX_success;
}
