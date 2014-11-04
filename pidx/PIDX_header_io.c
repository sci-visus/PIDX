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

#include "PIDX_header_io.h"
#include "PIDX_utils.h"
#define MAX_TEMPLATE_DEPTH 6

static int populate_meta_data(PIDX_header_io_id header_io_id, int file_number, char* bin_file);

struct PIDX_header_io_struct 
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
};

PIDX_header_io_id PIDX_header_io_init(idx_dataset idx_meta_data,
			idx_dataset_derived_metadata idx_derived_ptr,
			int start_var_index, int end_var_index )
{
  PIDX_header_io_id header_io_id;

  //Creating the IO ID
  header_io_id = malloc(sizeof (*header_io_id));
  if (!header_io_id) 
  memset(header_io_id, 0, sizeof (*header_io_id));
  
  header_io_id->idx_ptr = idx_meta_data;
  header_io_id->idx_derived_ptr = idx_derived_ptr;
  
  header_io_id->idx_ptr = (idx_dataset)malloc(sizeof(*(header_io_id->idx_ptr)));
  memcpy(header_io_id->idx_ptr, idx_meta_data, sizeof(*(header_io_id->idx_ptr)));
  
  header_io_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(header_io_id->idx_derived_ptr)));
  memcpy(header_io_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(header_io_id->idx_derived_ptr)));
    
  header_io_id->start_var_index = start_var_index;
  header_io_id->end_var_index = end_var_index;
      
  return header_io_id;
}

#if PIDX_HAVE_MPI
int PIDX_header_io_set_communicator(PIDX_header_io_id header_io, MPI_Comm comm)
{
  MPI_Comm_dup(comm, &header_io->comm);
  return 0;
}
#endif

int PIDX_header_io_write_idx (PIDX_header_io_id header_io, char* data_set_path, int current_time_step)
{
  int l = 0, rank = 0, N;
  FILE* idx_file_p;
  char dirname[1024], basename[1024];
  
#if PIDX_HAVE_MPI
  MPI_Comm_rank(header_io->comm, &rank);
#endif
  
  int nbits_blocknumber = (header_io->idx_derived_ptr->maxh - header_io->idx_ptr->bits_per_block - 1);
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
  strcpy(header_io->filename_template, data_set_path);
  for (N = strlen(header_io->filename_template) - 1; N >= 0; N--) 
  {
    int ch = header_io->filename_template[N];
    header_io->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
    strcat(header_io->filename_template, "/%01x.bin");
   
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
      strcat(header_io->filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12) 
      strcat(header_io->filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16) 
      strcat(header_io->filename_template, "/%04x.bin"); //no directories, 65536  files
    else 
    {
      while (nbits_blocknumber > 16) 
      {
	strcat(header_io->filename_template, "/%02x"); //256 subdirectories
	nbits_blocknumber -= 8;
      }
      strcat(header_io->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      assert(nbits_blocknumber <= 0);
    }
  }
  
  if (strncmp(".idx", &data_set_path[strlen(data_set_path) - 4], 4) != 0) 
  {
    fprintf(stderr, "[%s] [%d] Bad file name extension.\n", __FILE__, __LINE__);
    return 1;
  }
    
  if (rank == 0)
  {
    //fprintf(stderr, "writing IDX file...\n", __FILE__, __LINE__);

    idx_file_p = fopen(data_set_path, "w");
    if (!idx_file_p) 
    {
      fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
      return 1;
    }

    fprintf(idx_file_p, "(version)\n6\n");
    fprintf(idx_file_p, "(box)\n0 %lld 0 %lld 0 %lld 0 %lld 0 %lld\n", (header_io->idx_ptr->global_bounds[0] - 1), (header_io->idx_ptr->global_bounds[1] - 1), (header_io->idx_ptr->global_bounds[2] - 1), (header_io->idx_ptr->global_bounds[3] - 1), (header_io->idx_ptr->global_bounds[4] - 1));
    fprintf(idx_file_p, "(fields)\n");  
    
    for (l = 0; l < header_io->end_var_index; l++) 
    {
      fprintf(idx_file_p, "%s %s", header_io->idx_ptr->variable[l]->var_name, header_io->idx_ptr->variable[l]->type_name);
      if (l != header_io->end_var_index - 1)
	fprintf(idx_file_p, " + \n");
    }
    
    fprintf(idx_file_p, "\n(bits)\n%s\n", header_io->idx_ptr->bitSequence);
    fprintf(idx_file_p, "(bitsperblock)\n%d\n(blocksperfile)\n%d\n", header_io->idx_ptr->bits_per_block, header_io->idx_ptr->blocks_per_file);
    fprintf(idx_file_p, "(filename_template)\n./%s\n", header_io->filename_template);
    fprintf(idx_file_p, "(time)\n0 %d time%%06d/"/*note: uintah starts at timestep 1, but we shouldn't assume...*/, header_io->idx_ptr->current_time_step); //fix #1: need to add * notation to idx reader (because we can't write the .idx every time)
    fclose(idx_file_p);
  }
  
  return 0;
}

int PIDX_header_io_file_create(PIDX_header_io_id header_io_id)
{
  int i = 0, rank = 0, nprocs = 1, j, ret;
  char bin_file[PATH_MAX];
  char last_path[PATH_MAX] = {0};
  char this_path[PATH_MAX] = {0};
  char tmp_path[PATH_MAX] = {0};
  char* pos;
  
#if PIDX_HAVE_MPI
  MPI_Comm_rank(header_io_id->comm, &rank);
  MPI_Comm_size(header_io_id->comm, &nprocs);
#endif
  
  for (i = 0; i < header_io_id->idx_derived_ptr->max_file_count; i++) 
  {
#if PIDX_HAVE_MPI
    if (i % nprocs == rank && header_io_id->idx_derived_ptr->file_bitmap[i] == 1) 
#else
    if (rank == 0 && header_io_id->idx_derived_ptr->file_bitmap[i] == 1) 
#endif
    {
      ret = generate_file_name(header_io_id->idx_ptr->blocks_per_file, header_io_id->idx_ptr->filename_template, i, bin_file, PATH_MAX);
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
      populate_meta_data(header_io_id, i, bin_file);
    }
  }
  
#if PIDX_HAVE_MPI
  MPI_Barrier(header_io_id->comm);
#endif
  
  return 0;
}

static int populate_meta_data(PIDX_header_io_id header_io_id, int file_number, char* bin_file)
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
  
  total_header_size = (10 + (10 * header_io_id->idx_ptr->blocks_per_file)) * sizeof (uint32_t) * header_io_id->end_var_index;
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
  
  //printf("FS block: %d %d\n", header_io_id->start_fs_block, header_io_id->idx_derived_ptr->fs_block_size);
  for (m = 0; m < 10; m++)
    headers[m] = htonl(0);
  for (n = 0; n < header_io_id->end_var_index; n++)
  {
    bytes_per_sample = header_io_id->idx_ptr->variable[n]->bits_per_value / 8;
    initial_offset = 0;
    for (i = 0; i < header_io_id->idx_ptr->blocks_per_file; i++) 
    {
      if (is_block_present((i + (header_io_id->idx_ptr->blocks_per_file * file_number)), header_io_id->idx_ptr->variable[n]->global_block_layout))
      {
	block_negative_offset = find_block_negative_offset(header_io_id->idx_ptr->blocks_per_file, (i + (header_io_id->idx_ptr->blocks_per_file * file_number)), header_io_id->idx_ptr->variable[n]->global_block_layout);
	if (n == 0) 
	{
	  block_limit = i - block_negative_offset;
	  data_offset = ((i) - block_negative_offset) * header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * header_io_id->idx_ptr->variable[n]->values_per_sample;
	  data_offset += header_io_id->idx_derived_ptr->start_fs_block * header_io_id->idx_derived_ptr->fs_block_size;
	}
	else 
	{
	  if (i == 0)
	    for (b = 0; b < n; b++)
	    {
	      bytes_per_sample_previous = header_io_id->idx_ptr->variable[b]->bits_per_value / 8;
	      initial_offset = initial_offset + ((block_limit + 1) * header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample_previous * header_io_id->idx_ptr->variable[b]->values_per_sample);
	    }

	  data_offset = initial_offset + ((i) - block_negative_offset) * header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * header_io_id->idx_ptr->variable[n]->values_per_sample;
	  data_offset += header_io_id->idx_derived_ptr->start_fs_block * header_io_id->idx_derived_ptr->fs_block_size;
	}
	max_offset = data_offset + header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * header_io_id->idx_ptr->variable[n]->values_per_sample;
	  
	little_data_offset = 0;
	little_data_offset += data_offset;

	headers[10 + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[11 + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[12 + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(little_data_offset);
	headers[13 + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
	headers[14 + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * header_io_id->idx_ptr->variable[n]->values_per_sample);
	
	//printf("[%d] [%d] (%d) :: [%d %d %d] Offste and Count %d %d\n", file_number, i, header_io_id->idx_derived_ptr->start_fs_block, header_io_id->idx_derived_ptr->samples_per_block, bytes_per_sample, header_io_id->idx_ptr->variable[n]->values_per_sample, little_data_offset, header_io_id->idx_derived_ptr->samples_per_block * bytes_per_sample * header_io_id->idx_ptr->variable[n]->values_per_sample);
	
	for (m = 15; m < 20; m++)
	  headers[m + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
      } 
      else 
      {
	for (m = 10; m < 20; m++)
	  headers[m + ((i + (header_io_id->idx_ptr->blocks_per_file * n))*10)] = htonl(0);
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

int PIDX_header_io_finalize(PIDX_header_io_id header_io_id)
{
  #if PIDX_HAVE_MPI
  MPI_Comm_free(&header_io_id->comm);
#endif
  
  free(header_io_id->idx_derived_ptr);
  header_io_id->idx_derived_ptr = 0;
  
  free(header_io_id->idx_ptr);
  header_io_id->idx_ptr = 0;
  
  free(header_io_id);
  header_io_id = 0;

  return 0;
}