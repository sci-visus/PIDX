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

#include "PIDX_io.h"
#include "PIDX_utils.h"

struct PIDX_io_struct 
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
    
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_var_index;
  int end_var_index;
};

static int write_read_samples(PIDX_io_id io_id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, long long buffer_offset, int MODE);

//static int generate_file_name(PIDX_io_id io_id, int file_number, char* filename, int maxlen);

static int populate_meta_data(PIDX_io_id io_id, int file_number, char* bin_file)
{
  int total_header_size, ret = 0, block_negative_offset = 0, block_limit = 0;
  int bytes_per_sample, bytes_per_sample_previous;
  
#if PIDX_HAVE_MPI
  MPI_File fh;
  MPI_Status status;
#else
  int fh = 0;
#endif  
  
  int i = 0, m = 0, n = 0, b = 0;
  uint32_t* headers;
  long long initial_offset = 0;
  off_t data_offset = 0, max_offset = 0;  
  uint32_t little_data_offset;
  
  total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * io_id->idx_ptr->variable_count;
  headers = (uint32_t*)malloc(total_header_size);
  
#if PIDX_HAVE_MPI      
  ret = MPI_File_open(MPI_COMM_SELF, bin_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_open() failed on %s\n", __FILE__, __LINE__, bin_file);
    return 1;
  }
#else
  fh = open(bin_file, O_CREAT | O_WRONLY, 0664);
#endif
  
  //printf("FS block: %d %d\n", io_id->start_fs_block, io_id->idx_derived_ptr->fs_block_size);
  for (m = 0; m < 10; m++)
    headers[m] = htonl(0);
  for (n = 0; n < io_id->idx_ptr->variable_count; n++)
  {
    bytes_per_sample = io_id->idx_ptr->variable[n]->bits_per_value / 8;
    initial_offset = 0;
    for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++) 
    {
      if (is_block_present((i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_ptr->variable[n]->global_block_layout))
      {
	block_negative_offset = find_block_negative_offset(io_id->idx_ptr->blocks_per_file, (i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_ptr->variable[n]->global_block_layout);
	if (n == 0) 
	{
	  block_limit = i - block_negative_offset;
	  data_offset = ((i) - block_negative_offset) * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample;
	  data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;
	}
	else 
	{
	  if (i == 0)
	    for (b = 0; b < n; b++)
	    {
	      bytes_per_sample_previous = io_id->idx_ptr->variable[b]->bits_per_value / 8;
	      initial_offset = initial_offset + ((block_limit + 1) * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample_previous * io_id->idx_ptr->variable[b]->values_per_sample);
	    }

	  data_offset = initial_offset + ((i) - block_negative_offset) * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample;
	  data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;
	}
	max_offset = data_offset + io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample;
	  
	little_data_offset = 0;
	little_data_offset += data_offset;

	headers[10 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[11 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[12 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(little_data_offset);
	headers[13 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[14 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample);
	
	//printf("[%d] [%d] :: [%d %d %d] Offste and Count %d %d\n", file_number, i, io_id->idx_derived_ptr->samples_per_block, bytes_per_sample, io_id->idx_ptr->variable[n]->values_per_sample, little_data_offset, io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample);
	
	for (m = 15; m < 20; m++)
	  headers[m + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
      } 
      else 
      {
	for (m = 10; m < 20; m++)
	  headers[m + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
      }
    }
  }
  
#if PIDX_HAVE_MPI
  ret = MPI_File_write_at(fh, 0, headers, total_header_size, MPI_BYTE, &status);
  if (ret != MPI_SUCCESS) 
  {
    fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
    return 1;
  }
#else
  pwrite(fh, headers, total_header_size, 0);
#endif
  
  free(headers);

#if 1 //sid: some problem detecting file size, temporarily disable this bit (or get rid of it)
  if (max_offset != 0) 
  {
#if PIDX_HAVE_MPI
    ret = MPI_File_set_size(fh, max_offset);
    if (ret != MPI_SUCCESS) 
    {
      fprintf(stderr, "[%s] [%d] MPI_File_set_size() failed.\n", __FILE__, __LINE__);
      return 1;
    }
#else
    ftruncate(fh, max_offset);
#endif
  }
#endif
  
#if PIDX_HAVE_MPI  
  ret = MPI_File_close(&fh);
  if (ret != MPI_SUCCESS) 
  {
    fprintf(stderr, "[%s] [%d] MPI_File_open() failed on %s\n", __FILE__, __LINE__, bin_file);
    return 1;
  }
#else
  close(fh);
#endif
  
  return 0;
}

static int write_read_samples(PIDX_io_id io_id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, long long buffer_offset, int MODE)
{
  int samples_per_file, block_number, file_index, file_count, ret = 0, block_negative_offset = 0, file_number;
  int bytes_per_sample, bytes_per_datatype;
  int i = 0, l = 0;
  char file_name[PATH_MAX];
  off_t data_offset = 0;

#if PIDX_HAVE_MPI
  MPI_File fh;
  MPI_Status status;
#else
  int fh;
#endif
    
  samples_per_file = io_id->idx_derived_ptr->samples_per_block * io_id->idx_ptr->blocks_per_file;    

  bytes_per_datatype = io_id->idx_ptr->variable[variable_index]->bits_per_value / 8;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * io_id->idx_ptr->variable[variable_index]->values_per_sample;
  
  while (hz_count) 
  {
    block_number = hz_start_index / io_id->idx_derived_ptr->samples_per_block;
    file_number = hz_start_index / samples_per_file;
    file_index = hz_start_index % samples_per_file;
    file_count = samples_per_file - file_index;
    
    if (file_count > hz_count)
      file_count = hz_count;

    // build file name
    ret = generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, file_number, file_name, PATH_MAX);
    if (ret == 1)
    {
      fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
      return 1;
    }
    //printf("File name = %s\n", file_name);
#if PIDX_HAVE_MPI
    if(MODE == PIDX_WRITE)
    {
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS) 
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
	return 1;
      }
    }
    else
    {
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS) 
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
	return 1;
      }
    }
#else
    fh = open(file_name, O_WRONLY);
#endif
    
    data_offset = 0;
    bytes_per_sample = io_id->idx_ptr->variable[variable_index]->bits_per_value / 8;
    data_offset = file_index * bytes_per_sample * io_id->idx_ptr->variable[variable_index]->values_per_sample;
    data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;
    block_negative_offset = find_block_negative_offset(io_id->idx_ptr->blocks_per_file, block_number, io_id->idx_ptr->variable[variable_index]->global_block_layout);
    data_offset -= block_negative_offset * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[variable_index]->values_per_sample;
      
    //printf("NEGATIVE = %d\n", block_negative_offset);
    for (l = 0; l < variable_index; l++) 
    {
      bytes_per_sample = io_id->idx_ptr->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
	if (is_block_present((i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_ptr->variable[l]->global_block_layout))
	  data_offset = data_offset + (io_id->idx_ptr->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_derived_ptr->samples_per_block);
    }
    
#if 0
    int u = 0;
    printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->idx_derived_ptr->start_fs_block, io_id->idx_derived_ptr->fs_block_size, data_offset, file_count, file_number, hz_start_index, samples_per_file);
    if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(double) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "float64") == 0)
    {
      double value;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&value, hz_buffer + (u * sizeof(double)), sizeof(double));
	printf("[DOUBLE] Value at %d = %f\n", u,  value);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(int) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "float32") == 0)
    {
      float fvalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&fvalue, hz_buffer + (u * sizeof(float)), sizeof(float));
	printf("[FLOAT] Value at %d = %f\n", u,  fvalue);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(float) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "int32") == 0)
    {
      int ivalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&ivalue, hz_buffer + (u * sizeof(int)), sizeof(int));
	printf("[INT] Value at %d = %f\n", u,  ivalue);
      }
    }
    else if(io_id->idx_ptr->variable[variable_index]->bits_per_value / 8 == (int)sizeof(unsigned char) && strcmp(io_id->idx_ptr->variable[variable_index]->type_name, "uint8") == 0)
    {
      unsigned char uivalue;
      printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	memcpy(&uivalue, hz_buffer + (u * sizeof(unsigned char)), sizeof(unsigned char));
	printf("[U INT] Value at %d = %f\n", u,  uivalue);
      }
    }
#endif

    if(MODE == PIDX_WRITE)
    {
#if PIDX_HAVE_MPI
      
      //printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->idx_derived_ptr->start_fs_block, io_id->idx_derived_ptr->fs_block_size, data_offset, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), file_number, hz_start_index, samples_per_file);
      
      MPI_File_write_at(fh, data_offset, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
#else
      pwrite(fh, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), data_offset);
#endif
    }
    if(MODE == PIDX_READ)
    {
#if PIDX_HAVE_MPI
      ret = MPI_File_read_at(fh, data_offset, hz_buffer, file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * (io_id->idx_ptr->variable[variable_index]->bits_per_value/8), MPI_BYTE, &status);
      if (ret != MPI_SUCCESS) 
      {
	fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
	return 1;
      }
#endif
    }
        
    //file_count = file_count / io_id->idx_ptr->variable[variable_index]->values_per_sample;
    
    hz_count -= file_count;
    hz_start_index += file_count;
    hz_buffer += file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample * bytes_per_datatype;

#if PIDX_HAVE_MPI
    MPI_File_close(&fh);
#else
    close(fh);
#endif

  }
  return (0);
}

PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index) 
{
  PIDX_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_io_id)malloc(sizeof (*io_id));
  if (!io_id) 
  memset(io_id, 0, sizeof (*io_id));
  
  io_id->idx_ptr = idx_meta_data;
  io_id->idx_derived_ptr = idx_derived_ptr;
  
  io_id->idx_ptr = (idx_dataset)malloc(sizeof(*(io_id->idx_ptr)));
  memcpy(io_id->idx_ptr, idx_meta_data, sizeof(*(io_id->idx_ptr)));
  
  io_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(io_id->idx_derived_ptr)));
  memcpy(io_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(io_id->idx_derived_ptr)));
    
  io_id->start_var_index = start_var_index;
  io_id->end_var_index = end_var_index;
      
  return io_id;
}

#if PIDX_HAVE_MPI
int PIDX_io_set_communicator(PIDX_io_id io_id, MPI_Comm comm)
{
  MPI_Comm_dup(comm, &io_id->comm);
  return 0;
}
#endif

void print_buffer(unsigned char *buffer, int offset, int count, int sample_size, int datatype_size, char* name)
{
  int i;
  if(datatype_size/8 == sizeof(int) && strcmp(name, "int32") == 0)
  {
    int ival = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&ival, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("Int value at %d = %d\n", i, ival);
    }
  }
  else if(datatype_size/8 == sizeof(double) && strcmp(name, "float64") == 0)
  {
    double dval = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&dval, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Double value at %d = %f\n", count, i, dval);
    }
  }
  else if (datatype_size/8 == sizeof(float) && strcmp(name, "float32") == 0)
  {
    float fval = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&fval, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Float value at %d = %f\n", count, i, fval);
    }
  }
  else if (datatype_size/8 == sizeof(unsigned char) && strcmp(name, "uint8") == 0)
  {
    unsigned int uival = 0;
    for(i = 0; i < count * sample_size; i++)
    {
      memcpy(&uival, buffer + (offset * sample_size + i) * (datatype_size/8), (datatype_size/8));
      printf("[%d] Unsigned int value at %d = %d\n", count, i, uival);
    }
  }
}

int PIDX_io_aggregated_IO(PIDX_io_id io_id, Agg_buffer agg_buffer, int MODE)
{
  int i, k;
#if PIDX_HAVE_MPI
  MPI_File fh;
  MPI_Status status;
#else
  int fh;
#endif
  
  char file_name[PATH_MAX];
  off_t data_offset;
  
  if (agg_buffer->var_number != -1 && agg_buffer->sample_number != -1 && agg_buffer->file_number != -1) 
  {
    generate_file_name(io_id->idx_ptr->blocks_per_file, io_id->idx_ptr->filename_template, (unsigned int) agg_buffer->file_number, file_name, PATH_MAX);

#if PIDX_HAVE_MPI
    MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
#else
    fh = open(file_name, O_WRONLY);
#endif
      
    data_offset = 0;
    data_offset += io_id->idx_derived_ptr->start_fs_block * io_id->idx_derived_ptr->fs_block_size;
    
    for (k = 0; k < agg_buffer->var_number; k++) 
      for (i = 0; i < io_id->idx_ptr->variable[k]->values_per_sample; i++)
	data_offset = data_offset + io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[k]->bits_per_value/8);
    
    for (i = 0; i < agg_buffer->sample_number; i++)
      data_offset = data_offset + io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
  
    //printf("[%d] :: [FS %d %d] [%d %d %d]\n", rank, io_id->idx_derived_ptr->start_fs_block, io_id->idx_derived_ptr->fs_block_size, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->file_number);

#if PIDX_HAVE_MPI
    if(MODE == PIDX_WRITE)
      MPI_File_write_at(fh, data_offset, agg_buffer->buffer, ((io_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number]) * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8)), MPI_BYTE, &status);
    else
      MPI_File_read_at(fh, data_offset, agg_buffer->buffer, ((io_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number]) * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8)), MPI_BYTE, &status);
#else
      pwrite(fh, agg_buffer->buffer, ((io_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number]) * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8)), data_offset);
#endif

#if PIDX_HAVE_MPI
    MPI_File_close(&fh);
#else
    close(fh);
#endif
    
  }
  return 0;
}

int PIDX_io_independent_IO_var(PIDX_io_id io_id, PIDX_variable* variable, int MODE)
{
  int e1 = 0, i = 0, p = 0, var = 0;
  int send_index = 0;
  long long hz_index = 0;
  long long index = 0, count = 0;
   
  for (p = 0; p < io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->type == 0)
    {
      for (i = 0; i < variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i++) 
	hz_index = hz_index + variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i];
      
      for (i = variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for(e1 = 0; e1 < variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
	{
	  if(e1 == 0)
	  {
	    index = variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
	    send_index = e1;
	    count = 1;
	    
	    if(variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
	    {
	      for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
	      {
		write_read_samples(io_id, var, index, count, variable[var]->HZ_patch[p]->buffer[i], (send_index), MODE);
		//printf("X");
		//print_buffer(variable[var]->HZ_patch[p]->buffer[i], send_index, count, variable[var]->values_per_sample, variable[var]->bits_per_value, variable[var]->type_name);
	      }
	    }
	  }
	  else
	  {
	    if(variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
	    {
	      count++;
	      if(e1 == variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
	      {
		for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
		{
		  write_read_samples(io_id, var, index, count, variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
		  //printf("Y");
		  //print_buffer(variable[var]->HZ_patch[p]->buffer[i], send_index, count, variable[var]->values_per_sample, variable[var]->bits_per_value, variable[var]->type_name);
		}
	      }
	    }
	    else
	    {
	      for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
	      {
		write_read_samples(io_id, var, index, count, variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
		//printf("Z");
		//print_buffer(variable[var]->HZ_patch[p]->buffer[i], send_index, count, variable[var]->values_per_sample, variable[var]->bits_per_value, variable[var]->type_name);
	      }
	      
	      if(e1 == variable[io_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
	      {
		for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
		{
		  write_read_samples(io_id, var, variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, variable[var]->HZ_patch[p]->buffer[i], e1, MODE);
		  //printf("M");
		  //print_buffer(variable[var]->HZ_patch[p]->buffer[i], e1, 1, variable[var]->values_per_sample, variable[var]->bits_per_value, variable[var]->type_name);
		}
	      }
	      index = variable[io_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
	      count = 1;
	      send_index = e1;
	    }
	  }
	  hz_index++;
	}
      }
    }
    
    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->type == 1)
    {
      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for(var = io_id->start_var_index; var <= io_id->end_var_index; var++)
	{
	  index = 0;
	  count =  io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1;
	  write_read_samples(io_id, var, io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
	}
      }
    }
    
    else if(io_id->idx_ptr->variable[io_id->start_var_index]->patch_group_ptr[p]->type == 2)
    {
      int start_block_index, end_block_index, bl;
      for (i = io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < io_id->idx_ptr->variable[io_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for (var = io_id->start_var_index; var <= io_id->end_var_index; var++)
	{
	  start_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
	  end_block_index = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] / io_id->idx_derived_ptr->samples_per_block;
	  assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);
	  
	  send_index = 0;
	  for (bl = start_block_index; bl <= end_block_index; bl++) 
	  {
	    if (end_block_index == start_block_index) 
	    {
	      index = 0;
	      count = (io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i] + 1);
	    } 
	    else
	    {
	      if (bl == start_block_index) 
	      {
		index = 0;
		count = ((start_block_index + 1) * io_id->idx_derived_ptr->samples_per_block) - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i];
	      } 
	      else if (bl == end_block_index) 
	      {
		index = (end_block_index * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
		count = io_id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index[i] - ((end_block_index) * io_id->idx_derived_ptr->samples_per_block) + 1;
	      }
	      else 
	      {
		index = (bl * io_id->idx_derived_ptr->samples_per_block - io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i]);
		count = io_id->idx_derived_ptr->samples_per_block;
	      }
	    }
	    //printf("[%d] Offset %d Count %d\n", i, index + io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count);
	    //printf(" %d [B%d] INDEX and COUNT and SEND INDEX %d %d %d\n", bl, i, index, count, send_index); 
	    if(bl == 14)
	    {
	      double check;
	      memcpy(&check, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i] + send_index * (io_id->idx_ptr->variable[var]->bits_per_value / 8) * io_id->idx_ptr->variable[var]->values_per_sample, sizeof(double));
	      printf("%d %d %d (%d) VALUE = %f\n", var, p, i, (send_index  * io_id->idx_ptr->variable[var]->values_per_sample), check);
	    }
	    
	    write_read_samples(io_id, var, index + io_id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index[i], count, io_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
	    
	    send_index = send_index + count;
	  }
	}
      }
    }
  }
  return 0;
}

int PIDX_io_finalize(PIDX_io_id io_id) 
{
#if PIDX_HAVE_MPI
  MPI_Comm_free(&io_id->comm);
#endif
  
  free(io_id->idx_derived_ptr);
  io_id->idx_derived_ptr = 0;
  
  free(io_id->idx_ptr);
  io_id->idx_ptr = 0;
  
  free(io_id);
  io_id = 0;

  return 0;
}