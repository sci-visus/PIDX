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

#include "PIDX_inc.h"
#define PIDX_ACTIVE_TARGET 1
//#define RANK_ORDER 1

struct PIDX_agg_struct 
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Win win;
#endif
  
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_var_index;
  int end_var_index;
  
  int aggregator_interval;
};

PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{  
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  /*
  agg_id->idx_ptr = (idx_dataset)malloc(sizeof(*(agg_id->idx_ptr)));
  memcpy(agg_id->idx_ptr, idx_meta_data, sizeof(*(agg_id->idx_ptr)));
  
  agg_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(agg_id->idx_derived_ptr)));
  memcpy(agg_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(agg_id->idx_derived_ptr)));
  */
  
  agg_id->idx_ptr = idx_meta_data;
  agg_id->idx_derived_ptr = idx_derived_ptr;
  agg_id->start_var_index = start_var_index;
  agg_id->end_var_index = end_var_index;
  
  return agg_id;
}

#if PIDX_HAVE_MPI
int PIDX_agg_set_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  agg_id->comm = comm;
  //MPI_Comm_dup(comm, &agg_id->comm);
  return 0;
}
#endif

int aggregate_write_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, int buffer_offset, int MODE)
{
  int ret;
  int rank = 0, itr;
  int bytes_per_datatype;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, values_per_sample;
  int target_rank = 0;
  long long start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, hz_start = 0, samples_in_file = 0;
  long long samples_per_file = (long long) agg_id->idx_derived_ptr->samples_per_block * agg_id->idx_ptr->blocks_per_file;
  //MPI_Aint target_disp_address;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

  values_per_sample = agg_id->idx_ptr->variable[variable_index]->values_per_sample; //number of samples for variable j

  //starting HZ index for the data buffer at level "level" and for regular box number "box"
  hz_start = hz_start_index;
  
  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / agg_id->idx_derived_ptr->samples_per_block;

  //number of empty blocks befor block "block_no" in the file "file_no"
#ifdef PIDX_VAR_SLOW_LOOP
  negative_block_offset = find_block_negative_offset(agg_id->idx_ptr->blocks_per_file, block_no, agg_id->idx_ptr->variable[variable_index]->VAR_global_block_layout);
  assert(negative_block_offset >= 0);
  
  //number of samples in file "file_no"
  samples_in_file = agg_id->idx_ptr->variable[variable_index]->VAR_blocks_per_file[file_no] * agg_id->idx_derived_ptr->samples_per_block;
  assert(samples_in_file <= samples_per_file);
#else
  negative_block_offset = find_block_negative_offset(agg_id->idx_ptr->blocks_per_file, block_no, agg_id->idx_derived_ptr->global_block_layout);
  assert(negative_block_offset >= 0);
  
  //number of samples in file "file_no"
  samples_in_file = agg_id->idx_derived_ptr->existing_blocks_index_per_file[file_no] * agg_id->idx_derived_ptr->samples_per_block;
  assert(samples_in_file <= samples_per_file);
#endif

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  assert(hz_start >= (samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block));
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block))) * values_per_sample)
    %
    (samples_in_file * values_per_sample);
  assert(target_disp >= 0);

  sample_index = target_disp / (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);
  assert(sample_index < agg_id->idx_ptr->variable[variable_index]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor);
  
  target_disp = target_disp % (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);

#if RANK_ORDER
  target_rank = agg_buffer->rank_holder[variable_index - agg_id->start_var_index][sample_index][file_no];
#else
  target_rank = agg_buffer->rank_holder[file_no][variable_index - agg_id->start_var_index][sample_index];
#endif
  target_count = hz_count * values_per_sample;
  
  bytes_per_datatype = agg_id->idx_ptr->variable[variable_index]->bits_per_value / 8;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;
  
  start_agg_index = target_disp / (long long) (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor);
  end_agg_index = ((target_disp + target_count - 1) / (long long) (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor));
  assert(start_agg_index >= 0 && end_agg_index >= 0 && end_agg_index >= start_agg_index);
  
  if (start_agg_index != end_agg_index)
  {
    if (target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      //target_disp_address = target_disp;
      if (MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[A] Count %lld Local Disp %d Target Disp %lld\n", ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp), 0, target_disp);
#endif
        ret = MPI_Put(hz_buffer, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }

#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    } 
    else
      if (MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MA] Count %lld Local Disp %d Target Disp %lld\n", target_disp, 0, ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp));
#endif
        
        memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, ( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) * bytes_per_datatype);
      
    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++) 
    {
      if (target_rank != rank)
      {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
        if (MODE == PIDX_WRITE)
        {
#if PIDX_PRINT_AGG
          if (rank == 0)
            printf("[B] Count %lld Local Dis %lld Target Disp %d\n", (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif
          
          ret = MPI_Put(hz_buffer + (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
        else
        {
          ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype, MPI_BYTE, agg_id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
      }
      else
      {
        if (MODE == PIDX_WRITE)
        {
#if PIDX_PRINT_AGG
          if (rank == 0)
            printf("[MB] Count %lld Local Dis %lld Target Disp %d\n", (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor), (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif
          
          memcpy( agg_buffer->buffer, hz_buffer + (( (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, ( samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype);
        }
        else
          memcpy( hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + (itr * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) * bytes_per_datatype);
      }
    }
      
    if (target_rank + agg_id->aggregator_interval != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
      if (MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[C] Count %lld Local Dis %lld Target Disp %d\n", (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif
        ret = MPI_Put(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype, 
                      MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype, 
                      MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
    }
    else
      if(MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MC] Count %lld Local Dis %lld Target Disp %d\n", (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor)) - target_disp))), (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))), 0);
#endif
        
        memcpy( agg_buffer->buffer, hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype);    
      }
      else
        memcpy( hz_buffer + (((samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor))) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / agg_id->idx_derived_ptr->aggregation_factor) - target_disp)) * bytes_per_datatype);
  }
  else 
  {
    if(target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      //target_disp_address = target_disp;
      if(MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[D] Count %lld Local Dis %d Target Disp %lld\n", hz_count, 0, target_disp);
#endif
        
        ret = MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }
      else
      {
        ret = MPI_Get(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return (-1);
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    }
    else
    {
      if(MODE == PIDX_WRITE)
      {
#if PIDX_PRINT_AGG
        if (rank == 0)
          printf("[MD] Count %lld Local Dis %d Target Disp %lld\n", hz_count, 0, target_disp);
#endif
        memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
      }
      else
        memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
    }
  }
  
  return PIDX_success;
}

int PIDX_agg_aggregate(PIDX_agg_id agg_id, Agg_buffer agg_buffer) 
{
  int i, j, k, var;
  int rank_counter = 0, no_of_aggregators = 0, nprocs = 1, rank = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif
  
  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;
  
#ifdef PIDX_VAR_SLOW_LOOP
  for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    no_of_aggregators = no_of_aggregators + agg_id->idx_ptr->variable[var]->values_per_sample * agg_id->idx_ptr->variable[var]->VAR_existing_file_count;
#else
  for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    no_of_aggregators = no_of_aggregators + agg_id->idx_ptr->variable[var]->values_per_sample * agg_id->idx_derived_ptr->existing_file_count;
#endif
  
  agg_id->aggregator_interval = nprocs/ (no_of_aggregators * agg_id->idx_derived_ptr->aggregation_factor);
  assert(agg_id->aggregator_interval != 0);
    
#if RANK_ORDER
  agg_buffer->rank_holder = malloc((agg_id->end_var_index - agg_id->start_var_index + 1) * sizeof (int**));
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++) 
  {
    agg_buffer->rank_holder[i - agg_id->start_var_index] = malloc( agg_id->idx_ptr->variable[i]->values_per_sample  * sizeof (int*));
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      agg_buffer->rank_holder[i - agg_id->start_var_index][j] = malloc(/*agg_id->idx_ptr->variable[i]->existing_file_count*/ agg_id->idx_derived_ptr->max_file_count * sizeof (int));
      memset(agg_buffer->rank_holder[i - agg_id->start_var_index][j], 0, agg_id->idx_derived_ptr->max_file_count * sizeof (int));
    }
  }
#else
  agg_buffer->rank_holder = malloc(agg_id->idx_derived_ptr->max_file_count * sizeof (int**));
  for (i = 0; i < agg_id->idx_derived_ptr->max_file_count; i++) 
  {
    agg_buffer->rank_holder[i] = malloc( (agg_id->end_var_index - agg_id->start_var_index + 1)  * sizeof (int*));
    for (j = agg_id->start_var_index; j <= agg_id->end_var_index; j++)
    {
      agg_buffer->rank_holder[i][j - agg_id->start_var_index] = malloc( agg_id->idx_ptr->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_derived_ptr->aggregation_factor);
      memset(agg_buffer->rank_holder[i][j - agg_id->start_var_index], 0, agg_id->idx_ptr->variable[j]->values_per_sample * sizeof (int) * agg_id->idx_derived_ptr->aggregation_factor);
    }
  }
#endif
  
  rank_counter = 0;
#if RANK_ORDER

#ifdef PIDX_VAR_SLOW_LOOP
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      for (k = 0; k < agg_id->idx_ptr->variable[i]->VAR_existing_file_count; k++)
      {
        agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;
        
        if(rank == agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]])
        {
          agg_buffer->file_number = agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;
          
          agg_buffer->buffer_size = agg_id->idx_ptr->variable[agg_buffer->var_number]->VAR_blocks_per_file[agg_buffer->file_number] * agg_id->idx_derived_ptr->samples_per_block * (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
          //printf("Aggregator Rank %d Buffer Size %d (Var no: %d) (Sample no: %d) (File no: %d) (%d x %d x %d)\n", rank, agg_buffer->buffer_size, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->file_number, agg_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number], agg_id->idx_derived_ptr->samples_per_block, (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8));
        }
      }
    }
  }
#else
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      for (k = 0; k < agg_id->idx_derived_ptr->existing_file_count; k++)
      {
        agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_derived_ptr->existing_file_index[k]] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;
        
        if(rank == agg_buffer->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_derived_ptr->existing_file_index[k]])
        {
          agg_buffer->file_number = agg_id->idx_derived_ptr->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;
          
          agg_buffer->buffer_size = agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_buffer->file_number] * agg_id->idx_derived_ptr->samples_per_block * (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
          //printf("Aggregator Rank %d Buffer Size %d (Var no: %d) (Sample no: %d) (File no: %d) (%d x %d x %d)\n", rank, agg_buffer->buffer_size, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->file_number, agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_buffer->file_number], agg_id->idx_derived_ptr->samples_per_block, (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8));
        }
      }
    }
  }
#endif

#else

#ifdef PIDX_VAR_SLOW_LOOP
  for (k = 0; k < agg_id->idx_ptr->variable[i]->VAR_existing_file_count; k++)
  {
    for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
    {
      for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
      {
        agg_buffer->rank_holder[agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]][i - agg_id->start_var_index][j] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;
        
        if(rank == agg_buffer->rank_holder[agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k]][i - agg_id->start_var_index][j])
        {
          agg_buffer->file_number = agg_id->idx_ptr->variable[i]->VAR_existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;
          
          agg_buffer->buffer_size = agg_id->idx_ptr->variable[agg_buffer->var_number]->VAR_blocks_per_file[agg_buffer->file_number] * agg_id->idx_derived_ptr->samples_per_block * (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8) / agg_id->idx_derived_ptr->aggregation_factor;
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
        }
      }
    }
  }
#else
  for (k = 0; k < agg_id->idx_derived_ptr->existing_file_count; k++)
  {
    for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
    {
      for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample * agg_id->idx_derived_ptr->aggregation_factor; j++)
      {
        agg_buffer->rank_holder[agg_id->idx_derived_ptr->existing_file_index[k]][i - agg_id->start_var_index][j] = rank_counter;
        rank_counter = rank_counter + agg_id->aggregator_interval;
                  
        if(rank == agg_buffer->rank_holder[agg_id->idx_derived_ptr->existing_file_index[k]][i - agg_id->start_var_index][j])
        {
          agg_buffer->file_number = agg_id->idx_derived_ptr->existing_file_index[k];
          agg_buffer->var_number = i;
          agg_buffer->sample_number = j;
          
          agg_buffer->buffer_size = agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_buffer->file_number] * (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor) * (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
          
          agg_buffer->buffer = malloc(agg_buffer->buffer_size);
          if (agg_buffer->buffer == NULL)
          {
            printf("[%d] [%d %d %d] : %lld (%d %d (%d/%d) %d)\n", rank, i, j, k, agg_buffer->buffer_size, agg_id->idx_derived_ptr->existing_blocks_index_per_file[agg_buffer->file_number], (agg_id->idx_derived_ptr->samples_per_block / agg_id->idx_derived_ptr->aggregation_factor), agg_id->idx_derived_ptr->samples_per_block, agg_id->idx_derived_ptr->aggregation_factor, (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8));
            
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", agg_buffer->buffer_size, __LINE__, __FILE__);
            return (-1);
          }
          //memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
        }
      }
    }
  }
#endif

#endif
  
  return PIDX_success;
}

int PIDX_agg_aggregate_write_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int MODE)
{
  int i, p, e1, var, ret = 0;
  int send_index = 0;
  long long index = 0, count = 0, hz_index = 0;
  int variable_order = 1;
  int aggregate_lower_levels = 0;
  int existing_levels = 0;
  int element_count = 0;
  int ***send_offset, ***send_count;
  unsigned char ***data_buffer;
  int counter = 0, dest_counter = 0;
  MPI_Datatype **chunk_data_type;
  
  int rank = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
#endif

#if PIDX_HAVE_MPI
  if(agg_buffer->buffer_size != 0)
    MPI_Win_create(agg_buffer->buffer, agg_buffer->buffer_size, agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  else
    MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));    
        
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
#endif
  
  send_offset = malloc(sizeof(*send_offset) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count);
  memset(send_offset, 0, (sizeof(*send_offset) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count));
  send_count = malloc(sizeof(*send_count) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count);
  memset(send_count, 0, (sizeof(*send_count) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count));
  data_buffer = malloc(sizeof(*data_buffer) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count);
  memset(data_buffer, 0, (sizeof(*data_buffer) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count));
  chunk_data_type = malloc(sizeof(*chunk_data_type) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count);
  memset(chunk_data_type, 0, (sizeof(*chunk_data_type) * agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count));
  
  for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_ptr[p]->box_group_type == 0)
    {
      for (i = 0; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i++) 
        hz_index = hz_index + agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i];
      
      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
        if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
        {
          for(e1 = 0; e1 < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
          {
            if(e1 == 0)
            {
              index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
              send_index = e1;
              count = 1;
              
              if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
              {
                for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                {
                  //printf("[A] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
                  ret = aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
                  if (ret == -1)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return (-1);
                  }
                }
              }
            }
            else
            {
              if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
              {
                count++;
                if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
                {
                  for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                  {
                    //printf("[B] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
                    aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
                    if (ret == -1)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return (-1);
                    }
                  }
                }
              }
              else
              {
                for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                {
                  //printf("[C] Size %lld Offset %lld\n", count, index);
                  aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
                  if (ret == -1)
                  {
                    fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                    return (-1);
                  }
                }

                if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
                {
                  for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
                  {
                    //printf("[D] Size %lld Offset %lld\n", count, index);
                    aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, MODE);
                    if (ret == -1)
                    {
                      fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                      return (-1);
                    }
                  }
                }
                index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
                count = 1;
                send_index = e1;
              }
            }
            hz_index++;
          }
        }
      }
    }
    else
    {
      if (aggregate_lower_levels == 1)
      {
        int intermediate_hz_level = ((agg_id->idx_ptr->bits_per_block + 1) >= agg_id->idx_derived_ptr->maxh) ? agg_id->idx_ptr->bits_per_block + 1 : agg_id->idx_ptr->bits_per_block + 2;
        //intermediate_hz_level = 24;
        
        for (i = 0; i < intermediate_hz_level; i++)
        {
          if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
          {
            existing_levels++;
            element_count = element_count + agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->start_hz_index[i] + 1;
          }
        }
        data_buffer[p] = malloc(sizeof(*data_buffer[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        memset(data_buffer[p], 0, sizeof(*data_buffer[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        
        send_offset[p] = malloc(sizeof(*send_offset[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        send_count[p] = malloc(sizeof(*send_count[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        memset(send_offset[p], 0, sizeof(*send_offset[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        memset(send_count[p], 0, sizeof(*send_count[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        for (var = 0; var <= (agg_id->end_var_index - agg_id->start_var_index); var++)
        {
          send_offset[p][var] = malloc(existing_levels * sizeof(*send_offset[p][var]));
          send_count[p][var] = malloc(existing_levels * sizeof(*send_count[p][var]));
          memset(send_offset[p][var], 0, existing_levels * sizeof(*send_offset[p][var]));
          memset(send_count[p][var], 0, existing_levels * sizeof(*send_count[p][var]));
        }
        
        int bytes_per_datatype;
        chunk_data_type[p] = malloc( sizeof(*chunk_data_type[p]) * (agg_id->end_var_index - agg_id->start_var_index + 1));
        
        for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
        {
          bytes_per_datatype = agg_id->idx_ptr->variable[var]->bits_per_value / 8;
          data_buffer[p][var] = malloc( element_count * sizeof(*data_buffer[p][var]) * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype );
          memset(data_buffer[p][var], 0, element_count *  sizeof(*data_buffer[p][var]) * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype);
          counter = 0;
          dest_counter = 0;
          for (i = 0; i < intermediate_hz_level; i++)
          {
            if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
            {
              memcpy(data_buffer[p][var] + dest_counter, 
                    agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 
                    (agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1) * agg_id->idx_ptr->variable[var]->values_per_sample  * bytes_per_datatype );
              send_count[p][var][counter] = (agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1) * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype;
              dest_counter = dest_counter + send_count[p][var][counter];

              send_offset[p][var][counter] = agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype;
              
              counter++;
            }
          }
          
          MPI_Type_indexed(existing_levels, send_count[p][var], send_offset[p][var], MPI_BYTE, &(chunk_data_type[p][var]));
          MPI_Type_commit(&(chunk_data_type[p][var]));

#if RANK_ORDER
          ret = MPI_Put(data_buffer[p][var], element_count * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype, MPI_BYTE, agg_buffer->rank_holder[var - agg_id->start_var_index][0][0], 0, 1, chunk_data_type[p][var], agg_id->win);
#else
          ret = MPI_Put(data_buffer[p][var], element_count * agg_id->idx_ptr->variable[var]->values_per_sample * bytes_per_datatype, MPI_BYTE, agg_buffer->rank_holder[0][var - agg_id->start_var_index][0], 0, 1, chunk_data_type[p][var], agg_id->win);
#endif
          if(ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return (-1);
          }
        }
        
        if (variable_order == 0)
        {
          for (i = intermediate_hz_level; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
          {
            if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
            {
              for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
              {
                index = 0;
                count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);
                
                aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
              }
            }
          }
        }
        else
        {
          for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
          {
            for (i = intermediate_hz_level; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
            {
              if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
              {
                index = 0;
                count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);
                
                aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
              }
            }
          }
        }
      }
      else
      {
        if (variable_order == 0)
        {
          for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
          {
            if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
            {
              for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
              {
                index = 0;
                count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);
                
                aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
                if (ret == -1)
                {
                  fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                  return (-1);
                }
              }
            }
          }
        }
        else
        {
          for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
          {
            for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
            {
              if (agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] != 0)
              {
                index = 0;
                count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1 - (agg_id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level[i] * agg_id->idx_derived_ptr->samples_per_block);
                
#if PIDX_PRINT_AGG
                if (rank == 0)
                  printf("[AGG] [Color %d] [VAR %d] [HZ %d] Size %lld Send Offset %lld\n", agg_id->idx_derived_ptr->color, var, i, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
#endif
                
                aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
                if (ret == -1)
                {
                  fprintf(stderr, " Error in aggregate_write_read Line %d File %s\n", __LINE__, __FILE__);
                  return (-1);
                }
              }
            }
          }
        }
      }
    }
  }

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);		//First Fence
#else
  //MPI_Win_create has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
  MPI_Win_free(&(agg_id->win));
#endif
  
  if (aggregate_lower_levels == 1)
  {
    for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
    {
      for (var = 0; var <= (agg_id->end_var_index - agg_id->start_var_index); var++)
      {
        
        free(send_offset[p][var]);
        send_offset[p][var] = 0;
        free(send_count[p][var]);
        send_count[p][var] = 0;
        free(data_buffer[p][var]);
        data_buffer[p][var] = 0;
        
      }
      free(send_offset[p]);
      send_offset[p] = 0;
      free(send_count[p]);
      send_count[p] = 0;
      free(data_buffer[p]);
      data_buffer[p] = 0;
      free(chunk_data_type[p]);
    }
    free(send_offset);
    send_offset = 0;
    free(send_count);
    send_count = 0;
    free(chunk_data_type);
    free(data_buffer);
    data_buffer = 0;
  }
  
  return PIDX_success;
}

int PIDX_agg_buf_destroy(PIDX_agg_id agg_id, Agg_buffer agg_buffer) 
{
  if (agg_buffer->buffer_size != 0) 
  {
    free(agg_buffer->buffer);
    agg_buffer->buffer = 0;
  }
  
  int i = 0, j = 0;
#if RANK_ORDER
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++) 
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      free(agg_buffer->rank_holder[i - agg_id->start_var_index][j]);
      agg_buffer->rank_holder[i - agg_id->start_var_index][j] = 0;
    }
    free(agg_buffer->rank_holder[i - agg_id->start_var_index]);
    agg_buffer->rank_holder[i - agg_id->start_var_index] = 0;
  }
#else
  for (i = 0; i < agg_id->idx_derived_ptr->max_file_count; i++) 
  {
    for (j = agg_id->start_var_index; j <= agg_id->end_var_index; j++)
    {
      free(agg_buffer->rank_holder[i][j - agg_id->start_var_index]);
      agg_buffer->rank_holder[i][j - agg_id->start_var_index] = 0;
    }
    free(agg_buffer->rank_holder[i]);
  }
#endif
  
  free(agg_buffer->rank_holder);
  agg_buffer->rank_holder = 0;
  
  return PIDX_success;
}

int PIDX_agg_finalize(PIDX_agg_id agg_id) 
{
/*
  free(agg_id->idx_ptr);
  agg_id->idx_ptr = 0;
  
  free(agg_id->idx_derived_ptr);
  agg_id->idx_derived_ptr = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_free(&agg_id->comm);
#endif
*/
  free(agg_id);
  agg_id = 0;

  return 0;
}