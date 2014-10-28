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
#define MAX_TEMPLATE_DEPTH 6

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
  
  char filename_template[1024];
  int fs_block_size;
  off_t start_fs_block;
};

static int write_read_samples(PIDX_io_id io_id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, long long buffer_offset, int MODE);

static int generate_file_name(PIDX_io_id io_id, int file_number, char* filename, int maxlen);

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
  
  //printf("FS block: %d %d\n", io_id->start_fs_block, io_id->fs_block_size);
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
	  data_offset += io_id->start_fs_block * io_id->fs_block_size;
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
	  data_offset += io_id->start_fs_block * io_id->fs_block_size;
	}
	max_offset = data_offset + io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample;
	  
	little_data_offset = 0;
	little_data_offset += data_offset;

	headers[10 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[11 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[12 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(little_data_offset);
	headers[13 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[14 + ((i + (io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample);
	
	//printf("[%d] :: [%d %d %d] Offste and Count %d %d\n", file_number, io_id->idx_derived_ptr->samples_per_block, bytes_per_sample, io_id->idx_ptr->variable[n]->values_per_sample, little_data_offset, io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[n]->values_per_sample);
	
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

#if 0 //sid: some problem detecting file size, temporarily disable this bit (or get rid of it)
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
    ret = generate_file_name(io_id, file_number, file_name, PATH_MAX);
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
    data_offset += io_id->start_fs_block * io_id->fs_block_size;
    block_negative_offset = find_block_negative_offset(io_id->idx_ptr->blocks_per_file, block_number, io_id->idx_ptr->variable[variable_index]->global_block_layout);
    //file_count = file_count * io_id->idx_ptr->variable[variable_index]->values_per_sample;

    data_offset -= block_negative_offset * io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * io_id->idx_ptr->variable[variable_index]->values_per_sample;
      
    for (l = 0; l < variable_index; l++) 
    {
      bytes_per_sample = io_id->idx_ptr->variable[l]->bits_per_value / 8;
      for (i = 0; i < io_id->idx_ptr->blocks_per_file; i++)
	if (is_block_present((i + (io_id->idx_ptr->blocks_per_file * file_number)), io_id->idx_ptr->variable[l]->global_block_layout))
	{
	  data_offset = data_offset + (io_id->idx_ptr->variable[l]->values_per_sample * bytes_per_sample * io_id->idx_derived_ptr->samples_per_block);
	}
    }
    
#if 0
    int u = 0;
    printf("[Offset Count] [%d %d] :: [%d %d] File Number %d HZ Start %lld Total Samples %d\n", io_id->start_fs_block, io_id->fs_block_size, data_offset, file_count, file_number, hz_start_index, samples_per_file);
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
      
      double value;
      int u;
      //printf("[D] Offset %d\n", data_offset);
      for(u = 0 ; u < file_count ; u++)
      {
	//memcpy(&value, hz_buffer + (u * sizeof(double)), sizeof(double));
	//printf("[DOUBLE] Value at %d = %f\n", u,  value);
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

PIDX_io_id PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, 
#if PIDX_HAVE_MPI
			MPI_Comm in_comm,
#endif
			int start_var_index, int end_var_index) 
{
  int ret;
  PIDX_io_id io_id;

  //Creating the IO ID
  io_id = (PIDX_io_id)malloc(sizeof (*io_id));
  if (!io_id) 
  memset(io_id, 0, sizeof (*io_id));

#if PIDX_HAVE_MPI
  MPI_Comm_dup(in_comm, &io_id->comm);
#endif
  
  io_id->idx_ptr = idx_meta_data;
  io_id->idx_derived_ptr = idx_derived_ptr;
  
  io_id->idx_ptr = (idx_dataset)malloc(sizeof(*(io_id->idx_ptr)));
  memcpy(io_id->idx_ptr, idx_meta_data, sizeof(*(io_id->idx_ptr)));
  
  io_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(io_id->idx_derived_ptr)));
  memcpy(io_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(io_id->idx_derived_ptr)));
    
  io_id->start_var_index = start_var_index;
  io_id->end_var_index = end_var_index;
    
  io_id->fs_block_size == 0;

  return io_id;
}

int PIDX_io_file_create(PIDX_io_id io_id, int time_step, char* data_set_path, int MODE) 
{
  int i = 0, l = 0, rank = 0, nprocs = 1, j, ret, N;
  FILE* idx_file_p;
  char bin_file[PATH_MAX];
  char last_path[PATH_MAX] = {0};
  char this_path[PATH_MAX] = {0};
  char tmp_path[PATH_MAX] = {0};
  char* pos;
  char dirname[1024], basename[1024];
  char* idx_file;
  struct stat stat_buf;
  int total_header_size;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(io_id->comm, &rank);
  MPI_Comm_size(io_id->comm, &nprocs);
#endif
  
  int nbits_blocknumber = (io_id->idx_derived_ptr->maxh - io_id->idx_ptr->bits_per_block - 1);
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
  strcpy(io_id->filename_template, data_set_path);
  for (N = strlen(io_id->filename_template) - 1; N >= 0; N--) 
  {
    int ch = io_id->filename_template[N];
    io_id->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
  {
    strcat(io_id->filename_template, "/%01x.bin");
  } 
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
    {
      strcat(io_id->filename_template, "/%02x.bin"); //no directories, 256 files
    } 
    else if (nbits_blocknumber <= 12) 
    {
      strcat(io_id->filename_template, "/%03x.bin"); //no directories, 4096 files
    } 
    else if (nbits_blocknumber <= 16) 
    {
      strcat(io_id->filename_template, "/%04x.bin"); //no directories, 65536  files
    } 
    else 
    {
      while (nbits_blocknumber > 16) 
      {
	strcat(io_id->filename_template, "/%02x"); //256 subdirectories
	nbits_blocknumber -= 8;
      }
      strcat(io_id->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      assert(nbits_blocknumber <= 0);
    }
  }
  
  if (strncmp(".idx", &data_set_path[strlen(data_set_path) - 4], 4) != 0) 
  {
    fprintf(stderr, "[%s] [%d] Bad file name extension.\n", __FILE__, __LINE__);
    return 1;
  }
  
  idx_file = strdup(data_set_path);
  if (!idx_file) 
  {
    fprintf(stderr, "[%s] [%d] idx_file is corrupt.\n", __FILE__, __LINE__);
    return 1;
  }
  
  // only write idx file once (from first process, timestep 0, once all the variables have been added)
  if (rank == 0)
  {
    fprintf(stderr, "current_time_step = %d, end_var_index = %d, variable_count = %d\n",io_id->idx_ptr->current_time_step, io_id->end_var_index, io_id->idx_ptr->variable_count);
  }
  if (rank == 0 && io_id->idx_ptr->current_time_step == 2 && io_id->end_var_index == 9/*godhelpus io_id->idx_ptr->variable_count - 1*/)
  {
    fprintf(stderr, "writing IDX file...\n", __FILE__, __LINE__);

    idx_file_p = fopen(idx_file, "w");
    if (!idx_file_p) 
    {
      fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
      return 1;
    }

    fprintf(idx_file_p, "(version)\n6\n");
    fprintf(idx_file_p, "(box)\n0 %d 0 %d 0 %d 0 %d 0 %d\n", (io_id->idx_ptr->global_bounds[0] - 1), (io_id->idx_ptr->global_bounds[1] - 1), (io_id->idx_ptr->global_bounds[2] - 1), (io_id->idx_ptr->global_bounds[3] - 1), (io_id->idx_ptr->global_bounds[4] - 1));
    fprintf(idx_file_p, "(fields)\n");  
    
    for (l = 0; l < io_id->idx_ptr->variable_count; l++) 
    {
      fprintf(idx_file_p, "%s %s", io_id->idx_ptr->variable[l]->var_name, io_id->idx_ptr->variable[l]->type_name);
      
      /*
      if(io_id->idx_ptr->variable[l]->bits_per_value/8 == sizeof(double) && strcmp(io_id->idx_ptr->variable[l]->type_name, "float64") == 0)
	fprintf(idx_file_p, "float64 ");
      else if (io_id->idx_ptr->variable[l]->bits_per_value/8 == sizeof(int) && strcmp(io_id->idx_ptr->variable[l]->type_name, "int32") == 0)
	fprintf(idx_file_p, "int32 ");
      else if (io_id->idx_ptr->variable[l]->bits_per_value/8 == sizeof(unsigned char) && strcmp(io_id->idx_ptr->variable[l]->type_name, "uint8") == 0)
	fprintf(idx_file_p, "uint8 ");
      else if (io_id->idx_ptr->variable[l]->bits_per_value/8 == sizeof(float) && strcmp(io_id->idx_ptr->variable[l]->type_name, "float32") == 0)
	fprintf(idx_file_p, "float32 ");
      */
      if (l != io_id->idx_ptr->variable_count - 1)
	fprintf(idx_file_p, " + \n");
    }
    
    fprintf(idx_file_p, "\n(bits)\n%s\n", io_id->idx_ptr->bitSequence);
    fprintf(idx_file_p, "(bitsperblock)\n%d\n(blocksperfile)\n%d\n", io_id->idx_ptr->bits_per_block, io_id->idx_ptr->blocks_per_file);
    fprintf(idx_file_p, "(filename_template)\n./%s\n", io_id->filename_template);
    fprintf(idx_file_p, "(time)\n1 %d time%%06d/"/*note: uintah starts at timestep 1, but we shouldn't assume...*/, 99999 /*io_id->idx_ptr->current_time_step*/); //fix #1: need to add * notation to idx reader (because we can't write the .idx every time)
    fclose(idx_file_p);
  }

  //determine hdd block size
  if (rank == 0 && io_id->fs_block_size == 0)
  {
    FILE *dummy = fopen("dummy.txt", "w"); //TODO: close and delete the file (there is a way to do this automatically by fopen...)
    fclose(dummy);
    ret = stat("dummy.txt", &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] MPI_File_open() failed on %s\n", __FILE__, __LINE__, "dummy.txt");
      return 1;
    }
    io_id->fs_block_size = stat_buf.st_blksize;
  }
  
#if PIDX_HAVE_MPI
  MPI_Bcast(&io_id->fs_block_size, 1, MPI_INT, 0, io_id->comm);
#endif
  
  //printf("[%d] FS block: %d %d\n", (int)total_header_size, io_id->start_fs_block, io_id->fs_block_size);
  total_header_size = (10 + (10 * io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * /*io_id->idx_ptr->variable_count*/16;
  io_id->start_fs_block = total_header_size / io_id->fs_block_size;
    
  if (total_header_size % io_id->fs_block_size)
    io_id->start_fs_block++;
  
  //io_id->start_fs_block = 16;
  
  free(idx_file);
  
  char* directory_path;
  char* data_set_path_copy;
  directory_path = (char*) malloc(sizeof (char) * 1024);
  assert(directory_path);
  memset(directory_path, 0, sizeof (char) * 1024);

  data_set_path_copy = (char*) malloc(sizeof (char) * 1024);
  assert(data_set_path_copy);
  memset(data_set_path_copy, 0, sizeof (char) * 1024);

  assert(data_set_path);
  strncpy(directory_path, data_set_path, strlen(data_set_path) - 4);

  sprintf(data_set_path_copy, "%s/time%06d.idx", directory_path, io_id->idx_ptr->current_time_step);

  nbits_blocknumber = (io_id->idx_derived_ptr->maxh - io_id->idx_ptr->bits_per_block - 1);
  VisusSplitFilename(data_set_path_copy, dirname, basename);

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
  strcpy(io_id->filename_template, data_set_path_copy);
  for (N = strlen(io_id->filename_template) - 1; N >= 0; N--) 
  {
    int ch = io_id->filename_template[N];
    io_id->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
  {
    strcat(io_id->filename_template, "/%01x.bin");
  } 
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
    {
      strcat(io_id->filename_template, "/%02x.bin"); //no directories, 256 files
    } 
    else if (nbits_blocknumber <= 12) 
    {
      strcat(io_id->filename_template, "/%03x.bin"); //no directories, 4096 files
    } 
    else if (nbits_blocknumber <= 16) 
    {
      strcat(io_id->filename_template, "/%04x.bin"); //no directories, 65536  files
    } 
    else 
    {
      while (nbits_blocknumber > 16) 
      {
	strcat(io_id->filename_template, "/%02x"); //256 subdirectories
	nbits_blocknumber -= 8;
      }
      strcat(io_id->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      assert(nbits_blocknumber <= 0);
    }
  }
  if (strncmp(".idx", &data_set_path_copy[strlen(data_set_path_copy) - 4], 4) != 0) 
  {
    fprintf(stderr, "Error: bad file name extension.\n");
    return 1;
  }
  
  for (i = 0; i < io_id->idx_derived_ptr->max_file_count; i++) 
  {
#if PIDX_HAVE_MPI
    if (i % nprocs == rank && io_id->idx_derived_ptr->file_bitmap[i] == 1) 
#else
    if (rank == 0 && io_id->idx_derived_ptr->file_bitmap[i] == 1) 
#endif
    {
      ret = generate_file_name(io_id, i, bin_file, PATH_MAX);
      if (ret == 1) 
      {
	fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
	return 1;
      }

      //TODO: the logic for creating the subdirectory hierarchy could
      //be made to be more efficient than this. This implementation
      //walks up the tree attempting to mkdir() each segment every
      //time we switch to a new directory when creating binary files.

      // see if we need to make parent directory 
      strcpy(this_path, bin_file);
      if ((pos = rindex(this_path, '/'))) 
      {
	pos[1] = '\0';
	if (!strcmp(this_path, last_path) == 0) 
	{
	  //this file is in a previous directory than the last
	  //one; we need to make sure that it exists and create
	  //it if not.
	  strcpy(last_path, this_path);
	  memset(tmp_path, 0, PATH_MAX * sizeof (char));
	  //walk up path and mkdir each segment 
	  for (j = 0; j < strlen(this_path); j++) 
	  {
	    if (j > 0 && this_path[j] == '/') 
	    {
	      ret = mkdir(tmp_path, 0770);
	      if (ret != 0 && errno != EEXIST) 
	      {
		  perror("mkdir");
		  fprintf(stderr, "Error: failed to mkdir %s\n", tmp_path);
		  return 1;
	      }
	    }
	    tmp_path[j] = this_path[j];
	  }
	}
      }
      populate_meta_data(io_id, i, bin_file);
    }
  }
  free(directory_path);
  free(data_set_path_copy);
  
#if PIDX_HAVE_MPI
  MPI_Barrier(io_id->comm);
#endif
  
  return 0;
}

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

#if PIDX_HAVE_MPI
int PIDX_io_aggregated_IO(PIDX_io_id io_id, Agg_buffer agg_buffer, int MODE)
{
  int i, k;
  int rank, nprocs;
  MPI_File fh;
  MPI_Status status;
  char file_name[PATH_MAX];
  off_t data_offset;

  MPI_Comm_rank(io_id->comm, &rank);
  MPI_Comm_size(io_id->comm, &nprocs);

  if (agg_buffer->var_number != -1 && agg_buffer->sample_number != -1 && agg_buffer->file_number != -1) 
  {
    generate_file_name(io_id, (unsigned int) agg_buffer->file_number, file_name, PATH_MAX);
      
    MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      
    data_offset = 0;
    data_offset += io_id->start_fs_block * io_id->fs_block_size;
    
    for (k = 0; k < agg_buffer->var_number; k++) 
      for (i = 0; i < io_id->idx_ptr->variable[k]->values_per_sample; i++)
	data_offset = data_offset + io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[k]->bits_per_value/8);
    
    for (i = 0; i < agg_buffer->sample_number; i++)
      data_offset = data_offset + io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number] * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
  
    //printf("[%d] :: [FS %d %d] [%d %d %d]\n", rank, io_id->start_fs_block, io_id->fs_block_size, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->file_number);
    if(MODE == PIDX_WRITE)
      MPI_File_write_at(fh, data_offset, agg_buffer->buffer, ((io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number]) * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8)), MPI_BYTE, &status);
    else
      MPI_File_read_at(fh, data_offset, agg_buffer->buffer, ((io_id->idx_ptr->variable[k]->blocks_per_file[agg_buffer->file_number]) * io_id->idx_derived_ptr->samples_per_block * (io_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8)), MPI_BYTE, &status);
     
    MPI_File_close(&fh);
  }
  return 0;
}
#endif

int PIDX_io_independent_IO_var(PIDX_io_id io_id, PIDX_variable* variable, int MODE)
{
  int e1 = 0, i = 0, p = 0, var = 0;
  int send_index = 0;
  long long hz_index = 0;
  long long index = 0, count = 0;
  
  for (p = 0; p < io_id->idx_ptr->variable[io_id->start_var_index]->patch_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    
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
  return 0;
}

static int generate_file_name(PIDX_io_id io_id, int file_number, char* filename, int maxlen) 
{
  long long address = 0;
  unsigned int segs[MAX_TEMPLATE_DEPTH] = {0};
  int seg_count = 0;
  char* pos;
  int ret;

  //printf("[generate_file_name]: %d %s %d :: %s\n", file_number, filename, maxlen, io_id->filename_template);
  // determine the first HZ address for the file in question 
  address = file_number * io_id->idx_ptr->blocks_per_file;

  // walk backwards through the file name template to find places where we need to substitute strings
  for (pos = &io_id->filename_template[strlen(io_id->filename_template) - 1];
	  pos != io_id->filename_template;
	  pos--) 
  {
    // be careful not to lo0 past the end of the array 
    if (pos - io_id->filename_template > (strlen(io_id->filename_template) - 3))
      continue;

    if (pos[0] == '%' && pos[1] == '0' && pos[3] == 'x') 
    {
      // TODO: for now we have a hard coded max depth 
      if (seg_count >= MAX_TEMPLATE_DEPTH) 
      {
	fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", io_id->filename_template);
	return 1;
      }

      // found an occurance of %0 in the template; check the next character to see how many bits to use here 

      switch (pos[2]) 
      {
	case '1':
	    segs[seg_count] += address & 0xf;
	    address = address >> 4;
	    break;
	case '2':
	    segs[seg_count] += address & 0xff;
	    address = address >> 8;
	    break;
	case '3':
	    segs[seg_count] += address & 0xfff;
	    address = address >> 12;
	    break;
	case '4':
	    segs[seg_count] += address & 0xffff;
	    address = address >> 16;
	    break;
	case '5':
	    segs[seg_count] += address & 0xfffff;
	    address = address >> 20;
	    break;
	default:
	    // TODO: generalize this to work for any value 
	    fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", io_id->filename_template);
	    return 1;
      }
      seg_count++;
    }
  }
  switch (seg_count) 
  {
    case 0:
	ret = strlen(io_id->filename_template);
	if (ret < maxlen) {
	    strcpy(filename, io_id->filename_template);
	}
	break;
    case 1:
	ret = snprintf(filename, maxlen, io_id->filename_template, segs[0]);
	break;
    case 2:
	ret = snprintf(filename, maxlen, io_id->filename_template,
		segs[1], segs[0]);
	break;
    case 3:
	ret = snprintf(filename, maxlen, io_id->filename_template,
		segs[2], segs[1], segs[0]);
	break;
    case 4:
	ret = snprintf(filename, maxlen, io_id->filename_template,
		segs[3], segs[2], segs[1], segs[0]);
	break;
    case 5:
	ret = snprintf(filename, maxlen, io_id->filename_template,
		segs[4], segs[3], segs[2], segs[1], segs[0]);
	break;
    case 6:
	ret = snprintf(filename, maxlen, io_id->filename_template,
		segs[5], segs[4], segs[3], segs[2], segs[1], segs[0]);
	break;
    default:
	// TODO: generalize this 
	fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", io_id->filename_template);
	return 1;
	break;
  }
  // make sure that the resulting string fit into the buffer ok 
  if (ret >= maxlen - 1) 
  {
      fprintf(stderr, "Error: filename too short in generate_filename()\n");
      return 1;
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
