/******************************************************
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
 
#include "PIDX.h"

static double sim_start = 0, sim_end = 0;

const int PIDX_reccommended_bits_per_block       = 15;
const int PIDX_reccommended_blocks_per_file      = 256;

const unsigned PIDX_MODE_EXCL               = 1;     // Error creating a file that already exists.
const unsigned PIDX_MODE_CREATE             = 2;     // Create the file if it does not exist.

const unsigned PIDX_MODE_RDONLY             = 4;     // Read only.
const unsigned PIDX_MODE_WRONLY             = 16;    // Write only.

const unsigned PIDX_MODE_RDWR               = 8;     // Reading and writing.


const unsigned PIDX_MODE_APPEND             = 32;
const unsigned PIDX_MODE_DELETE_ON_CLOSE    = 64;
const unsigned PIDX_MODE_SEQUENTIAL         = 128;

int start_ts = 0;
int end_ts = 1;
int current_ts = 0;
char* ts_template;

static PIDX_return_code PIDX_cleanup(PIDX_file pidx);
static PIDX_return_code PIDX_write(PIDX_file pidx);

struct PIDX_file_descriptor 
{
  MPI_Comm comm; //MPI Communicator
  
  int IDX_WRITE; // 1 for write	0 for read
  int IDX_READ; // 0 for write	1 for read

  int flags;
  
  int current_time_step;
  
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  PIDX_rst_id rst_id;
  PIDX_hz_encode_id hz_id;
  PIDX_agg_id agg_id;
  PIDX_io_id io_id;
  
  Agg_buffer agg_buffer;
  
  int start_variable_index;
  int end_variable_index;
};

PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_file file)
{
  sim_start = MPI_Wtime();
  
  if(flags != PIDX_file_excl || flags != PIDX_file_trunc)
    return PIDX_err_unsopperted_flags;
    
  if(flags == PIDX_file_excl)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }
  
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 || !filename)
    return PIDX_err_name;
    
  file = malloc(sizeof (*file));
  memset(file, 0, sizeof (*file));

  file->flags = flags;
  
  file->idx_ptr = (idx_dataset)malloc(sizeof (*file->idx_ptr));
  memset(file->idx_ptr, 0, sizeof (*file->idx_ptr));
  
  file->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof (*file->idx_derived_ptr));
  memset(file->idx_derived_ptr, 0, sizeof (*file->idx_derived_ptr));

  file->idx_ptr->filename = strdup(filename);
  file->idx_ptr->global_bounds = malloc(sizeof(int) * PIDX_MAX_DIMENSIONS);
  
  file->idx_ptr->variable_count = 0;
  file->idx_ptr->variable_index_tracker = 0;
  
  file->start_variable_index = 0;
  file->end_variable_index = 0;
   
  file->idx_ptr->bits_per_block = PIDX_reccommended_bits_per_block;
  file->idx_derived_ptr->samples_per_block = pow(2, PIDX_reccommended_bits_per_block);
  file->idx_ptr->blocks_per_file = PIDX_reccommended_blocks_per_file;
  
  return PIDX_success;
}

PIDX_return_code PIDX_file_open(const char* data_set_path, PIDX_flags  flags, PIDX_file   file)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_set_bounding_box_size(PIDX_file file, PIDX_point box_size)
{
  if(!box_size || box_size[0] < 0 || box_size[1] < 0 || box_size[2] < 0 || box_size[3] < 0 || box_size[4] < 0)
    return PIDX_err_box;
  if(!file)
    return PIDX_err_file;
   
  memcpy(file->idx_ptr->global_bounds, box_size, (sizeof (int) * PIDX_MAX_DIMENSIONS));
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_bounding_box_size(PIDX_file file, PIDX_point box_size)
{
  if(!file)
    return PIDX_err_file;
  
  if(!box_size)
    return PIDX_err_box;
   
  memcpy(box_size, file->idx_ptr->global_bounds, (sizeof (int) * PIDX_MAX_DIMENSIONS));
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_current_time_step(PIDX_file file, const int current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  if(current_time_step < 0)
    return PIDX_err_time;
   
  file->idx_ptr->current_time_step = current_time_step;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_current_time_step(PIDX_file file, int* current_time_step)
{
  if(!file)
    return PIDX_err_file;
  
  *current_time_step = file->idx_ptr->current_time_step;
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_block_size(PIDX_file file, const int bits_per_block)
{
  if(!file)
    return PIDX_err_file;
  
  if(bits_per_block <= 0)
    return PIDX_err_block;
   
  file->idx_ptr->bits_per_block = bits_per_block;
  file->idx_derived_ptr->samples_per_block = pow(2, bits_per_block);
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_block_size(PIDX_file file, int* bits_per_block)
{ 
  if(!file)
    return PIDX_err_file;
  
  *bits_per_block = file->idx_ptr->bits_per_block;
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_block_count(PIDX_file file, const int blocks_per_file)
{
  if(!file)
    return PIDX_err_file;
  
  if(blocks_per_file <= 0)
    return PIDX_err_block;
   
  file->idx_ptr->blocks_per_file = blocks_per_file;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_block_count(PIDX_file file, int* blocks_per_file)
{ 
  if(!file)
    return PIDX_err_file;
  
  *blocks_per_file = file->idx_ptr->blocks_per_file;
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_communicator(PIDX_file file, MPI_Comm comm)
{
  if(!file)
    return PIDX_err_file;
  
  if(!comm)
    return PIDX_err_comm;
   
  MPI_Comm_dup(comm, &file->comm);
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_communicator(PIDX_file file, MPI_Comm* comm)
{
  if(!file)
    return PIDX_err_file;
  
  if(!file->comm)
    return PIDX_err_comm;
   
  MPI_Comm_dup(file->comm, comm);
  
  return PIDX_success;
}

PIDX_return_code PIDX_set_variable_count(PIDX_file file, int  variable_count)
{
  if(!file)
    return PIDX_err_file;
  
  if(variable_count <= 0)
    return PIDX_err_count;
   
  file->idx_ptr->variable_count = variable_count;
  
  return PIDX_success;
}

PIDX_return_code PIDX_get_variable_count(PIDX_file file, int* variable_count)
{ 
  if(!file)
    return PIDX_err_file;
  
  *variable_count = file->idx_ptr->variable_count;
  
  return PIDX_success;
}

PIDX_return_code PIDX_variable_create(PIDX_file file, char* variable_name, unsigned int bits_per_sample, PIDX_type type_name, PIDX_variable variable)
{
  if(!file)
    return PIDX_err_file;
  
  if(!variable_name)
    return PIDX_err_name;

  if(bits_per_sample <= 0)
    return PIDX_err_size;
  
  if(!type_name)
    return PIDX_err_type;
  
  if(file->idx_ptr->variable_index_tracker > file->idx_ptr->variable_count)
    return PIDX_err_count;
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker] = malloc(sizeof (*file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]));
  memset(file->idx_ptr->variable[file->idx_ptr->variable_index_tracker], 0, sizeof (*file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]));
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->bits_per_sample = bits_per_sample;
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->values_per_sample = 1;
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->bits_per_value = (bits_per_sample);
  
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->type_name = strdup(type_name);    
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->var_name = strdup(variable_name);
  file->idx_ptr->variable[file->idx_ptr->variable_index_tracker]->patch_count = 0;
 
  variable = file->idx_ptr->variable[file->idx_ptr->variable_index_tracker];
  
  file->idx_ptr->variable_index_tracker++;
  file->end_variable_index++;
  
  return PIDX_success;
}

PIDX_return_code PIDX_append_and_write_variable(PIDX_variable variable, PIDX_point offset, PIDX_point box_size, const void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  if(!variable)
    return PIDX_err_variable;
  
  if(!offset || offset[0] < 0 || offset[1] < 0 || offset[2] < 0 || offset[3] < 0 || offset[4] < 0)
    return PIDX_err_offset;
  
  if(!box_size || box_size[0] < 0 || box_size[1] < 0 || box_size[2] < 0 || box_size[3] < 0 || box_size[4] < 0)
    return PIDX_err_count;
  
  const void *temp_buffer;
  variable->patch[variable->patch_count] = malloc(sizeof(*(variable->patch[variable->patch_count])));
  memcpy(variable->patch[variable->patch_count]->offset, offset, PIDX_MAX_DIMENSIONS * sizeof(int));
  memcpy(variable->patch[variable->patch_count]->count, box_size, PIDX_MAX_DIMENSIONS * sizeof(int));
  
  temp_buffer = read_from_this_buffer;
  variable->patch[variable->patch_count]->buffer = (unsigned char*)temp_buffer;
  
  variable->data_layout = data_layout;
  variable->patch_count = variable->patch_count + 1;
  
  return PIDX_success; 
}

PIDX_return_code populate_idx_dataset(PIDX_file file)
{
  int var, i, j, p, ctr, counter = 0, file_number = 0, level_count = 1;
    int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };
  
  PointND global_bounds_point;
  global_bounds_point.x = file->idx_ptr->global_bounds[0];
  global_bounds_point.y = file->idx_ptr->global_bounds[1];
  global_bounds_point.z = file->idx_ptr->global_bounds[2];
  global_bounds_point.u = file->idx_ptr->global_bounds[3];
  global_bounds_point.v = file->idx_ptr->global_bounds[4];
  GuessBitmaskPattern(file->idx_ptr->bitSequence, global_bounds_point);
  file->idx_derived_ptr->maxh = strlen(file->idx_ptr->bitSequence);
  
  for (i = 0; i <= file->idx_derived_ptr->maxh; i++)
    file->idx_ptr->bitPattern[i] = RegExBitmaskBit(file->idx_ptr->bitSequence, i);
    
  for (var = file->start_variable_index; var < file->end_variable_index; var++)
  {
    block_layout* all_patch_local_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(all_patch_local_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    file->idx_ptr->variable[var]->global_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(file->idx_ptr->variable[var]->global_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    for(p = 0 ; p < file->idx_ptr->variable[var]->patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
	bounding_box[0][i] = file->idx_ptr->variable[var]->patch[p]->offset[i];
	bounding_box[1][i] = file->idx_ptr->variable[var]->patch[p]->count[i] + file->idx_ptr->variable[var]->patch[p]->offset[i];
      }
      
      block_layout* per_patch_local_block_layout = (block_layout*) malloc(sizeof (block_layout));
      createBlockBitmap(bounding_box, file->idx_ptr->blocks_per_file, file->idx_ptr->bits_per_block, file->idx_derived_ptr->maxh, file->idx_ptr->bitPattern, per_patch_local_block_layout);
      
      ctr = 1;
      for (i = 1; i < (all_patch_local_block_layout->levels); i++)
      {
	for(j = 0 ; j < ctr ; j++)
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
    for (i = 1; i < (file->idx_ptr->variable[var]->global_block_layout->levels); i++) 
    {
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i], level_count,
		MPI_INT, MPI_BOR, file->comm);
      level_count = level_count * 2;
    }
    
    destroyBlockBitmap(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
    
    ctr = 1;
    for (i = 1; i < (file->idx_ptr->variable[var]->global_block_layout->levels); i++)
    {
      for(j = 0 ; j < ctr ; j++)
      {
	if(file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i][j] != 0)
	  file->idx_ptr->variable[var]->global_block_layout->hz_block_count_array[i]++;
      }    
      
      counter = 0;
      for(j = 0 ; j < ctr ; j++)
      {
	if(file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i][j] != 0)
	{
	  file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i][counter] = file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i][j];
	  counter++;
	}
      }
      ctr = ctr * 2;
    }
  }
  
  file->idx_derived_ptr->max_file_count = (getPowerOf2(file->idx_ptr->global_bounds[0]) * getPowerOf2(file->idx_ptr->global_bounds[1]) * getPowerOf2(file->idx_ptr->global_bounds[2]) * getPowerOf2(file->idx_ptr->global_bounds[3]) * getPowerOf2(file->idx_ptr->global_bounds[4])) / ((unsigned long long) file->idx_derived_ptr->samples_per_block * (unsigned long long) file->idx_ptr->blocks_per_file);
  if ((getPowerOf2(file->idx_ptr->global_bounds[0]) * getPowerOf2(file->idx_ptr->global_bounds[1]) * getPowerOf2(file->idx_ptr->global_bounds[2]) * getPowerOf2(file->idx_ptr->global_bounds[3]) * getPowerOf2(file->idx_ptr->global_bounds[4])) % ((unsigned long long) file->idx_derived_ptr->samples_per_block * (unsigned long long) file->idx_ptr->blocks_per_file))
    file->idx_derived_ptr->max_file_count++;
  
  file->idx_derived_ptr->file_bitmap = (int*) malloc(file->idx_derived_ptr->max_file_count * sizeof (int));
  memset(file->idx_derived_ptr->file_bitmap, 0, file->idx_derived_ptr->max_file_count * sizeof (int));
  
  for (var = file->start_variable_index; var < file->end_variable_index; var++)
  { 
    int *temp_file_index = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
    memset(temp_file_index, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));
    
    file->idx_ptr->variable[var]->blocks_per_file = malloc(sizeof(int) * (file->idx_derived_ptr->max_file_count));
    memset(file->idx_ptr->variable[var]->blocks_per_file, 0, sizeof(int) * (file->idx_derived_ptr->max_file_count));
  
    temp_file_index[0] = 1;
    file->idx_derived_ptr->file_bitmap[0] = 1;
    file->idx_ptr->variable[var]->blocks_per_file[0] = 1;
    
    for (i = 1; i < file->idx_ptr->variable[var]->global_block_layout->levels; i++) 
    {
      for (j = 0; j < file->idx_ptr->variable[var]->global_block_layout->hz_block_count_array[i]; j++) 
      {
	file_number = file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i][j] / file->idx_ptr->blocks_per_file;
	file->idx_derived_ptr->file_bitmap[file_number] = 1;
	temp_file_index[file_number] = 1;
	file->idx_ptr->variable[var]->blocks_per_file[file_number]++;
      }
    }
    
    file->idx_ptr->variable[var]->existing_file_count = 0;
    for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
      if (temp_file_index[i] == 1)
	file->idx_ptr->variable[var]->existing_file_count++;
    
    file->idx_ptr->variable[var]->existing_file_index = (int*) malloc(file->idx_ptr->variable[var]->existing_file_count * sizeof (int));
    memset(file->idx_ptr->variable[var]->existing_file_index, 0, file->idx_ptr->variable[var]->existing_file_count * sizeof (int));
    
    int count = 0;
    for (i = 0; i < file->idx_derived_ptr->max_file_count; i++)
      if (temp_file_index[i] == 1)
      {
	file->idx_ptr->variable[var]->existing_file_index[count] = i;
	count++;
      }
  }
  
  return PIDX_success;
}

PIDX_return_code PIDX_read(PIDX_file file)
{

  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_write(PIDX_file file)
{ 
  int i = 0, p, var = 0;
  
  populate_idx_dataset(file);
  
  for(i = file->start_variable_index; i < file->end_variable_index/*file->idx_ptr->variable_count*/; i++)
  {
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, i, i);
    file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, file->comm, i, i);
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, file->comm, i, i);
    
    for (var = i; var < (i+1); var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
    PIDX_io_file_create(file->io_id, current_ts, file->idx_ptr->filename, PIDX_WRITE);
    PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable, PIDX_WRITE);
    PIDX_hz_encode_write_var(file->hz_id, file->idx_ptr->variable);
    
    file->agg_buffer = malloc(sizeof(*file->agg_buffer));
    PIDX_agg_aggregate(file->agg_id, file->agg_buffer);
    PIDX_agg_aggregate_write_read(file->agg_id, file->agg_buffer);
    
    
    PIDX_io_aggregated_IO(file->io_id, file->agg_buffer, PIDX_WRITE);
    //PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_WRITE);
    
    PIDX_agg_buf_destroy(file->agg_buffer);
    PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
    
    PIDX_io_finalize(file->io_id);
    PIDX_agg_finalize(file->agg_id);
    PIDX_hz_encode_finalize(file->hz_id);
  }
  
  
#if 0
  file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, file->start_variable_index, (file->end_variable_index - 1));
  file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, file->comm, file->start_variable_index, (file->end_variable_index - 1));
  
  for (var = file->start_variable_index; var < file->end_variable_index; var++)
  {
    for (p = 0; p < file->idx_ptr->variable[var]->patch_count; p++)
    {
      file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
      memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
    }
  }
  PIDX_io_file_create(file->io_id, current_ts, file->idx_ptr->data_set_path, PIDX_WRITE);
  PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable, PIDX_WRITE);
  PIDX_hz_encode_write_var(file->hz_id, file->idx_ptr->variable);
  
  
  PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_WRITE);
  PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
  PIDX_io_finalize(file->io_id);
  PIDX_hz_encode_finalize(file->hz_id);
#endif
  
  return PIDX_success;
}

PIDX_return_code PIDX_sync(PIDX_file pidx)
{
  PIDX_write(pidx);
  PIDX_cleanup(pidx);
  return PIDX_success;
}

PIDX_return_code PIDX_cleanup(PIDX_file pidx) 
{
  int i, p;
  for (i = pidx->start_variable_index; i < pidx->end_variable_index; i++) 
  {
    for(p = 0; p < pidx->idx_ptr->variable[i]->patch_count; p++)
    {
      free(pidx->idx_ptr->variable[i]->HZ_patch[p]);
      pidx->idx_ptr->variable[i]->HZ_patch[p] = 0;
    
      free(pidx->idx_ptr->variable[i]->patch[p]);
      pidx->idx_ptr->variable[i]->patch[p] = 0;
    }
  }
  
  pidx->start_variable_index = pidx->idx_ptr->variable_count;
  pidx->end_variable_index = pidx->idx_ptr->variable_count;
  
  return PIDX_success;
}

PIDX_return_code PIDX_close(PIDX_file pidx) 
{
  PIDX_sync(pidx);
  
  sim_end = MPI_Wtime();
  
  double total_time = sim_end - sim_start;
  double max_time;
  int sample_sum = 0, var = 0, i = 0;
  
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, pidx->comm);

  if (max_time == total_time) 
  {
    if (pidx->IDX_WRITE == 1 && pidx->IDX_READ == 0)
      printf("\n------------- WRITE -------------\n");

    if (pidx->IDX_WRITE == 0 && pidx->IDX_READ == 1)
      printf("\n------------- READ -------------\n");
      
    for (var = 0; var < pidx->idx_ptr->variable_count; var++)
      sample_sum = sample_sum + pidx->idx_ptr->variable[var]->values_per_sample;
    
    long long total_data = (long long) pidx->idx_ptr->global_bounds[0] * pidx->idx_ptr->global_bounds[1] * pidx->idx_ptr->global_bounds[2] * pidx->idx_ptr->global_bounds[3] * pidx->idx_ptr->global_bounds[4] * sample_sum * 8;
    
    printf("Global Data [%d %d %d] Variables [%d]\nTime Taken: %f Seconds Throughput %f MiB/sec\n", pidx->idx_ptr->global_bounds[0], pidx->idx_ptr->global_bounds[1], pidx->idx_ptr->global_bounds[2], pidx->idx_ptr->variable_count, max_time, (float) total_data / (1024 * 1024 * max_time));
  }
  
  for (i = 0; i < pidx->idx_ptr->variable_count; i++) 
  {
    free(pidx->idx_ptr->variable[i]->var_name);
    pidx->idx_ptr->variable[i]->var_name = 0;
    
    free(pidx->idx_ptr->variable[i]->type_name);
    pidx->idx_ptr->variable[i]->type_name = 0;
    
    destroyBlockBitmap(pidx->idx_ptr->variable[i]->global_block_layout);
    free(pidx->idx_ptr->variable[i]->global_block_layout);
    pidx->idx_ptr->variable[i]->global_block_layout = 0;
  }
  
  for (i = 0; i < 1024; i++)
  {
    free(pidx->idx_ptr->variable[i]);
    pidx->idx_ptr->variable[i] = 0; 
  }
  pidx->idx_ptr->variable_count = 0;
  
  free(pidx->idx_ptr->filename);
  pidx->idx_ptr->filename = 0;

  free(pidx->idx_ptr->global_bounds);
  pidx->idx_ptr->global_bounds = 0;
  
  free(pidx->idx_ptr);
  pidx->idx_ptr = 0;
  
  free(pidx->idx_derived_ptr);
  pidx->idx_derived_ptr = 0;
  
  MPI_Comm_free(&pidx->comm);
  
  free(pidx);
  pidx = 0;
  
  return PIDX_success;
}

PIDX_return_code PIDX_append_and_read_variable(PIDX_variable variable, PIDX_point offset, PIDX_point box_size, void* read_from_this_buffer, PIDX_data_layout data_layout)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_variable_set_box_metadata_on (PIDX_variable variable)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_variable_set_box_metadata_off(PIDX_variable variable)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_variable_get_box_metadata(PIDX_variable variable, int* on_off_bool)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_bits_per_sample(PIDX_type type_name, unsigned int bits_per_sample)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_box_count(PIDX_file file, int* box_count)
{
  return PIDX_err_not_implemented;
}
PIDX_return_code PIDX_get_box(PIDX_file file, int box_index, PIDX_point offset, PIDX_point box_size)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_box_count__with_rank(PIDX_file file, int* box_count, int MPI_rank)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_box_with_rank(PIDX_file file, int box_index, PIDX_point offset, PIDX_point box_size,  int MPI_rank)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_current_variable_index(PIDX_file file, int* variable_index)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_set_current_variable_index(PIDX_file file, int  variable_index)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_get_current_variable(PIDX_file file, PIDX_variable variable)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_set_current_variable(PIDX_file file, PIDX_variable variable)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_read_variable(PIDX_variable variable, PIDX_point offset, PIDX_point box_size, const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_write_variable(PIDX_variable variable, PIDX_point offset, PIDX_point box_size, 	const void* read_from_this_buffer, PIDX_data_layout layout)
{
  return PIDX_err_not_implemented;
}
