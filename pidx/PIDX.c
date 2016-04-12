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

#include "PIDX.h"

#define PIDX_DEBUG_OUTPUT 0

static int caching_state = 0;
static int time_step_caching = 0;
static uint32_t* cached_header_copy;

//static int static_var_counter = 0;

//static void PIDX_init_timming_buffers1();
//static void PIDX_init_timming_buffers2(PIDX_file file);
static PIDX_return_code PIDX_validate(PIDX_file file);


/// PIDX File descriptor (equivalent to the descriptor returned by)
/// POSIX or any other IO framework
struct PIDX_file_descriptor
{
  int flags;

  int io_type;
  PIDX_io io;
  PIDX_partitioned_idx_io partitioned_idx_io;
  PIDX_idx_io idx_io;
  PIDX_raw_io raw_io;

#if PIDX_HAVE_MPI
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
  MPI_Comm global_comm;                        ///< MPI super-communicator (includes all processes)
#endif

  PIDX_access access;                          ///< serial or parallel access


  int local_variable_index;                    ///<
  int local_variable_count;                    ///<

  int flush_used;
  int write_on_close;                          ///< HPC Writes
  int one_time_initializations;                ///<

  int ROI_writes;

  int small_agg_comm;

  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;

  int debug_output;

  //int layout_count;
  //int reduced_res_from;
  //int reduced_res_to;

  //PIDX_time time;
};




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
  if (flags != PIDX_MODE_CREATE && flags != PIDX_MODE_EXCL)
    return PIDX_err_unsupported_flags;
    
  if (flags == PIDX_MODE_EXCL)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }
  
  uint64_t i = 0;
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
  (*file)->idx_d = malloc(sizeof (*((*file)->idx_d)));
  memset((*file)->idx_d, 0, sizeof (*((*file)->idx_d)));

  (*file)->idx_d->time = malloc(sizeof (*((*file)->idx_d->time)));
  memset((*file)->idx_d->time, 0, sizeof (*((*file)->idx_d->time)));

  (*file)->idx_d->time->sim_start = PIDX_get_time();
  
  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));

  (*file)->flags = flags;

  (*file)->access = access_type;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    (*file)->idx_d->idx_count[i] = 1;//access_type->idx_count[i];

  //memcpy ((*file)->idx_d->idx_count, access_type->idx_count, sizeof(int) * PIDX_MAX_DIMENSIONS);

  (*file)->idx_dbg->debug_do_rst = 1;
  (*file)->idx_dbg->debug_do_chunk = 1;
  (*file)->idx_dbg->debug_do_compress = 1;
  (*file)->idx_dbg->debug_do_hz = 1;
  (*file)->idx_dbg->debug_do_agg = 1;
  (*file)->idx_dbg->debug_do_io = 1;
  (*file)->idx_dbg->debug_rst = 0;
  (*file)->idx_dbg->debug_hz = 0;

  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;

  (*file)->idx_d->color = 0;
  (*file)->small_agg_comm = 0;
  (*file)->io_type = PIDX_IDX_IO;
  //(*file)->enable_raw_dump = 0;
  (*file)->idx_d->data_core_count = -1;
  (*file)->debug_output = 0;

  (*file)->idx_d->reduced_res_from = 0;
  (*file)->idx_d->reduced_res_to = 0;

  (*file)->idx_d->parallel_mode = access_type->parallel;
  (*file)->idx_d->file_zero_optimization = 0;

#if PIDX_HAVE_MPI
  if (access_type->parallel)
  {
    MPI_Comm_rank(access_type->comm, &rank);
    (*file)->comm = access_type->comm;
    (*file)->idx->enable_rst = 1;
    (*file)->idx->enable_agg = 1;
  }
  else
  {
    (*file)->idx->enable_rst = 1;
    (*file)->idx->enable_agg = 1;
  }
#else
  (*file)->idx->enable_rst = 0;
  (*file)->idx->enable_agg = 0;
#endif

  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx->variable_index_tracker = 0;

  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';
  sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);

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

  (*file)->idx_d->agg_type = 1;
  (*file)->idx_d->layout_count = 0;
  (*file)->idx_d->reduced_res_from = 0;
  (*file)->idx_d->reduced_res_to = 0;

  if (rank == 0)
  {
    //TODO: close and delete the file (there is a way to do this automatically by fopen...)
    struct stat stat_buf;
    FILE *dummy = fopen(".dummy.txt", "w");
    fclose(dummy);
    ret = stat(".dummy.txt", &stat_buf);
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

  (*file)->idx_d->time->file_create_time = PIDX_get_time();

  return PIDX_success;
}


/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  int i;
  //int ret;
  char file_name_skeleton[1024];
  int rank = 0;

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
  
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  //(*file)->time = malloc(sizeof((*file)->time));
  //memset((*file)->time, 0, sizeof((*file)->time));

  (*file)->flags = flags;
    
  (*file)->idx = (idx_dataset)malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));
  
  (*file)->idx_d = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_d)));
  memset((*file)->idx_d, 0, sizeof (*((*file)->idx_d)));

  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));

  (*file)->idx_d->time = malloc(sizeof (*((*file)->idx_d->time)));
  memset((*file)->idx_d->time, 0, sizeof (*((*file)->idx_d->time)));
  (*file)->idx_d->time->sim_start = PIDX_get_time();


  (*file)->access = access_type;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    (*file)->idx_d->idx_count[i] = 1;

  (*file)->idx_dbg->debug_do_rst = 1;
  (*file)->idx_dbg->debug_do_chunk = 1;
  (*file)->idx_dbg->debug_do_compress = 1;
  (*file)->idx_dbg->debug_do_hz = 1;
  (*file)->idx_dbg->debug_do_agg = 1;
  (*file)->idx_dbg->debug_do_io = 1;
  (*file)->idx_dbg->debug_rst = 0;
  (*file)->idx_dbg->debug_hz = 0;

  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;
  (*file)->one_time_initializations = 0;

  (*file)->idx_d->reduced_res_from = 0;
  (*file)->idx_d->reduced_res_to = 0;

  (*file)->idx_d->file_zero_optimization = 0;

  (*file)->small_agg_comm = 0;
  (*file)->debug_output = 0;

  (*file)->idx_d->color = 0;

  (*file)->io_type = PIDX_IDX_IO;
  //(*file)->enable_raw_dump = 0;

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
      memcpy ((*file)->idx_d->idx_count, access_type->idx_count, sizeof(int) * PIDX_MAX_DIMENSIONS);

      colors = malloc(sizeof(*colors) * (*file)->idx_d->idx_count[0] * (*file)->idx_d->idx_count[1] * (*file)->idx_d->idx_count[2]);
      memset(colors, 0, sizeof(*colors) * (*file)->idx_d->idx_count[0] * (*file)->idx_d->idx_count[1] * (*file)->idx_d->idx_count[2]);

      for (k = 0; k < (*file)->idx_d->idx_count[2]; k++)
        for (j = 0; j < (*file)->idx_d->idx_count[1]; j++)
          for (i = 0; i < (*file)->idx_d->idx_count[0]; i++)
            colors[((*file)->idx_d->idx_count[0] * (*file)->idx_d->idx_count[1] * k) + ((*file)->idx_d->idx_count[0] * j) + i] = ((*file)->idx_d->idx_count[0] * (*file)->idx_d->idx_count[1] * k) + ((*file)->idx_d->idx_count[0] * j) + i;

      int index_x = 0, index_y = 0, index_z = 0;
      for (i = 0; i < access_type->sub_div[0]; i = i + (access_type->sub_div[0] / (*file)->idx_d->idx_count[0]))
      {
        if (access_type->rank_component[0] >= i && access_type->rank_component[0] < i + (access_type->sub_div[0] / (*file)->idx_d->idx_count[0]))
        {
          index_x = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[1]; i = i + (access_type->sub_div[1] / (*file)->idx_d->idx_count[1]))
      {
        if (access_type->rank_component[1] >= i && access_type->rank_component[1] < i + (access_type->sub_div[1] / (*file)->idx_d->idx_count[1]))
        {
          index_y = i;
          break;
        }
      }
      for (i = 0; i < access_type->sub_div[2]; i = i + (access_type->sub_div[2] / (*file)->idx_d->idx_count[2]))
      {
        if (access_type->rank_component[2] >= i && access_type->rank_component[2] < i + (access_type->sub_div[2] / (*file)->idx_d->idx_count[2]))
        {
          index_z = i;
          break;
        }
      }

      (*file)->idx_d->color = colors[((*file)->idx_d->idx_count[0] * (*file)->idx_d->idx_count[1] * (index_z/(access_type->sub_div[2] / (*file)->idx_d->idx_count[2]))) + ((*file)->idx_d->idx_count[0] * (index_y/ (access_type->sub_div[1] / (*file)->idx_d->idx_count[1]))) + (index_x / (access_type->sub_div[0] / (*file)->idx_d->idx_count[0]))];

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

  if ((*file)->idx_d->idx_count[0] == 1 && (*file)->idx_d->idx_count[1] == 1 && (*file)->idx_d->idx_count[2] == 1)
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
  (*file)->idx_d->data_core_count = -1;
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

      if (strcmp(line, "(raw_dump)") == 0)
      {
        (*file)->io_type = PIDX_RAW_IO;
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          (*file)->idx->reg_patch_size[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
        (*file)->idx->reg_patch_size[3] = 1;
        (*file)->idx->reg_patch_size[4] = 1;
      }
      if (strcmp(line, "(cores)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        (*file)->idx_d->data_core_count = atoi(line);
      }

      if (strcmp(line, "(fields)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        count = 0;
        variable_counter = 0;
        
        while (line[0] != '(')
        {

          (*file)->idx->variable[variable_counter] = malloc(sizeof (*((*file)->idx->variable[variable_counter])));
          if ((*file)->idx->variable[variable_counter] == NULL)
            return PIDX_err_file;

          memset((*file)->idx->variable[variable_counter], 0, sizeof (*((*file)->idx->variable[variable_counter])));
          
          pch1 = strtok(line, " +");
          while (pch1 != NULL)
          {
            if (count == 0)
            {
              char* temp_name = strdup(pch1);
              strcpy((*file)->idx->variable[variable_counter]->var_name, /*strdup(pch1)*/temp_name);
              free(temp_name);
            }

            if (count == 1)
            {
              len = strlen(pch1) - 1;
              if (pch1[len] == '\n')
                pch1[len] = 0;

              strcpy((*file)->idx->variable[variable_counter]->type_name, pch1);
              int ret;
              int bits_per_sample = 0;
              ret = PIDX_default_bits_per_datatype((*file)->idx->variable[variable_counter]->type_name, &bits_per_sample);
              if (ret != PIDX_success)  return PIDX_err_file;

              (*file)->idx->variable[variable_counter]->bits_per_value = bits_per_sample;
              (*file)->idx->variable[variable_counter]->values_per_sample = 1;
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
        if ((*file)->idx->compression_type != PIDX_NO_COMPRESSION)
        {
          int i1 = 0;
          for (i1 = 0; i1 < PIDX_MAX_DIMENSIONS; i1++)
            (*file)->idx->chunk_size[i1] = 4;
        }
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
    MPI_Bcast((*file)->idx->reg_patch_size, 5, MPI_LONG_LONG, 0, (*file)->comm);
    MPI_Bcast((*file)->idx->chunk_size, 5, MPI_LONG_LONG, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->blocks_per_file), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->bits_per_block), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->variable_count), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx_d->data_core_count), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->compression_bit_rate), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->idx->compression_type), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast(&((*file)->io_type), 1, MPI_INT, 0, (*file)->comm);

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
  
  (*file)->idx_d->time->file_create_time = PIDX_get_time();
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

  return PIDX_success;
}



PIDX_return_code PIDX_set_partition_size(PIDX_file file, int count_x, int count_y, int count_z)
{
  if(count_x < 0 || count_y < 0 || count_z < 0)
    return PIDX_err_box;
  
  if(file == NULL)
    return PIDX_err_file;
  
  file->idx_d->idx_size[0] = count_x;
  file->idx_d->idx_size[1] = count_y;
  file->idx_d->idx_size[2] = count_z;
  file->idx_d->idx_size[3] = 1;
  file->idx_d->idx_size[4] = 1;

  //printf("{S %d %d %d\n", file->idx_d->idx_size[0], file->idx_d->idx_size[1], file->idx_d->idx_size[2]);
  return PIDX_validate(file);
}


PIDX_return_code PIDX_set_aggregator_multiplier(PIDX_file file, int count_aggregator_multiplier)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->aggregator_multiplier = count_aggregator_multiplier;

  return PIDX_success;
}



PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims)
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(dims, file->idx->bounds, (sizeof(int64_t) * PIDX_MAX_DIMENSIONS));
  
  return PIDX_success;
}



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



PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_dbg->debug_do_chunk = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_compression(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_compress = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_hz(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_dbg->debug_do_hz = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_dbg->debug_do_agg = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_io(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_dbg->debug_do_io = 0;
  
  return PIDX_success;
}



PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst)
{
  if(!file)
    return PIDX_err_file;
  
  file->idx_dbg->debug_rst = debug_rst;
  
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
  
  file->idx_dbg->debug_hz = debug_hz;
  
  return PIDX_success;
}




PIDX_return_code PIDX_set_resolution(PIDX_file file, int hz_from, int hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  file->idx_d->reduced_res_from = hz_from;
  file->idx_d->reduced_res_to = hz_to;
  
  return PIDX_success;
}


PIDX_return_code PIDX_get_resolution(PIDX_file file, int *hz_from, int *hz_to)
{
  if(file == NULL)
    return PIDX_err_file;
  
  *hz_from = file->idx_d->reduced_res_from;
  *hz_to = file->idx_d->reduced_res_to;
  
  return PIDX_success;
}




PIDX_return_code PIDX_enable_raw_io(PIDX_file file)
{
  if(file == NULL)
    return PIDX_err_file;

  file->io_type = PIDX_RAW_IO;

  return PIDX_success;
}




PIDX_return_code PIDX_optimize_for_file_zero(PIDX_file file)
{
  if(file == NULL)
    return PIDX_err_file;

  file->idx_d->file_zero_optimization = 1;

  return PIDX_success;
}




PIDX_return_code PIDX_enable_partitioned_io(PIDX_file file)
{
  if(file == NULL)
    return PIDX_err_file;

  file->io_type = PIDX_PARTITIONED_IDX_IO;

  return PIDX_success;
}




PIDX_return_code PIDX_enable_partition_merge_io(PIDX_file file)
{
  if(file == NULL)
    return PIDX_err_file;

  file->io_type = PIDX_PARTITION_MERGE_IDX_IO;

  return PIDX_success;
}




PIDX_return_code PIDX_flush(PIDX_file file)
{
  int i, p;
  int ret;
  if (file->idx->variable_count <= 0)
    return PIDX_err_variable;

  file->io = PIDX_io_init(file->idx, file->idx_d, file->idx_dbg);
  if (file->io == NULL)
    return PIDX_err_flush;

  ret = PIDX_io_set_communicator(file->io, file->comm);
  if (ret != PIDX_success)
    return PIDX_err_io;

  ret = PIDX_io_io(file->io, file->flags, file->io_type, file->local_variable_index, (file->local_variable_index + file->local_variable_count));
  if (ret != PIDX_success)
    return PIDX_err_io;

  ret = PIDX_io_finalize(file->io);
  if (ret != PIDX_success)
    return PIDX_err_io;

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

  PIDX_time time = file->idx_d->time;
  time->sim_end = PIDX_get_time();

  int rank;
  MPI_Comm_rank(file->comm, &rank);

#if 1
  if (file->debug_output == 1)
  {
    if (rank == 0)
    {
      fprintf(stdout, "\n==========================================================================================================\n");
      fprintf(stdout, "Time step %d File name %s\n", file->idx->current_time_step, file->idx->filename);
      fprintf(stdout, "Global Data %lld %lld %lld Variables %d\n", (long long) file->idx->bounds[0], (long long) file->idx->bounds[1], (long long) file->idx->bounds[2], file->idx->variable_count);
      fprintf(stdout, "Partition count %d = %d x %d x %d\n", file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2], file->idx_d->idx_count[0], file->idx_d->idx_count[1], file->idx_d->idx_count[2]);
      fprintf(stdout, "Rst = %d Comp = %d\n", file->idx->enable_rst, file->idx->compression_type);
      fprintf(stdout, "Blocks Per File %d Bits per block %d File Count %d\n", file->idx->blocks_per_file, file->idx->bits_per_block, file->idx_d->max_file_count);
      fprintf(stdout, "Chunk Size %d %d %d %d %d\n", (int)file->idx->chunk_size[0], (int)file->idx->chunk_size[1], (int)file->idx->chunk_size[2], (int)file->idx->chunk_size[3], (int)file->idx->chunk_size[4]);
      fprintf(stdout, "Restructuring Box Size %d %d %d %d %d\n", (int)file->idx->reg_patch_size[0], (int)file->idx->reg_patch_size[1], (int)file->idx->reg_patch_size[2], (int)file->idx->reg_patch_size[3], (int)file->idx->reg_patch_size[4]);
      fprintf(stdout, "Aggregation Type = %d\n", file->idx_d->agg_type);
      fprintf(stdout, "File zero optimization = %d\n", file->idx_d->file_zero_optimization);
    }

    if (file->io_type == PIDX_IDX_IO)
      PIDX_print_idx_io_timing(file->comm, time, file->idx->variable_count, file->idx_d->layout_count);

    else if (file->io_type == PIDX_RAW_IO)
      PIDX_print_raw_io_timing(file->comm, time, time->variable_counter, file->idx_d->perm_layout_count);

    else if (file->io_type == PIDX_PARTITIONED_IDX_IO)
      PIDX_print_partition_timing(file->comm, time, file->idx->variable_count, file->idx_d->perm_layout_count);

    else if (file->io_type == PIDX_PARTITION_MERGE_IDX_IO)
      PIDX_print_partition_merge_timing(file->comm, time, file->idx->variable_count, file->idx_d->perm_layout_count);

    if (rank == 0)
      fprintf(stdout, "==========================================================================================================\n");
  }

  PIDX_delete_timming_buffers1(file->idx_d->time);
  PIDX_delete_timming_buffers2(file->idx_d->time, file->idx->variable_count);
#endif

  for (i = 0; i < file->idx->variable_count; i++)
    free(file->idx->variable[i]);

  file->idx->variable_count = 0;

  free(file->idx);                  file->idx = 0;
  free(file->idx_d->time);          file->idx_d->time = 0;
  free(file->idx_d);                file->idx_d = 0;
  free(file->idx_dbg);              file->idx_dbg = 0;

  free(file);
  
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

  file->idx_d->var_pipe_length = var_pipe_length;

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



PIDX_return_code PIDX_set_ROI_writes(PIDX_file file)
{
  if (file == NULL)
    return PIDX_err_file;

  file->ROI_writes = 1;
  file->idx->enable_rst = 0;
  //file->idx->enable_agg = 0;

  return PIDX_success;
}
