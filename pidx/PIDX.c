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

static PIDX_return_code PIDX_cleanup(PIDX_file file);
static PIDX_return_code PIDX_write(PIDX_file file);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* file_name, int current_time_step);

///
/// PIDX File descriptor (equivalent to the descriptor returned by)
/// POSIX or any other IO framework
///
struct PIDX_file_descriptor 
{
#if PIDX_HAVE_MPI
  MPI_Comm comm; //MPI Communicator
#endif
  
  int IDX_WRITE; // 1 for write 0 for read
  int IDX_READ;  // 0 for write 1 for read

  int flags;

  PIDX_access access;
  
  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  PIDX_header_io_id header_io_id;
  
#if PIDX_HAVE_MPI
  PIDX_rst_id rst_id;
#endif
  
  PIDX_hz_encode_id hz_id;
  PIDX_agg_id agg_id;
  PIDX_io_id io_id;
  
  Agg_buffer agg_buffer;
  
  int local_variable_index;
  int local_variable_count;
  int variable_pipelining_factor;
  
  /// HPC Writes
  int write_on_close;
  
  ///
  int one_time_initializations;
};

///
/// Returns elapsed time
///
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

///
/// Function to create IDX file descriptor (based on flags and access)
///
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  sim_start = PIDX_get_time();

  if(flags != PIDX_file_excl && flags != PIDX_file_trunc)
    return PIDX_err_unsupported_flags;
    
  if(flags == PIDX_file_excl)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }
  
  int ret, rank;
  
  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;
    
  *file = malloc(sizeof (*(*file)) );
  memset(*file, 0, sizeof (*(*file)) );

  (*file)->flags = flags;
  
  (*file)->idx_ptr = (idx_dataset)malloc(sizeof (*((*file)->idx_ptr)));
  memset((*file)->idx_ptr, 0, sizeof (*((*file)->idx_ptr)));
  
  (*file)->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof (*((*file)->idx_derived_ptr)));
  memset((*file)->idx_derived_ptr, 0, sizeof (*((*file)->idx_derived_ptr)));

  (*file)->idx_ptr->filename = strdup(filename);
  (*file)->idx_ptr->global_bounds = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
  { int i; for (i=0;i<PIDX_MAX_DIMENSIONS;i++) (*file)->idx_ptr->global_bounds[i]=65535; }

  //initialize logic_to_physic transform to identity
  (*file)->idx_ptr->transform[0]  = 1.0;
  (*file)->idx_ptr->transform[5]  = 1.0;
  (*file)->idx_ptr->transform[10] = 1.0;
  (*file)->idx_ptr->transform[15] = 1.0;
  
  (*file)->idx_ptr->variable_count = -1;
  (*file)->idx_ptr->variable_index_tracker = 0;

  (*file)->access = access_type;
  (*file)->idx_ptr->current_time_step = 0;
    
#if PIDX_HAVE_MPI
  if (access_type->parallel)
  {
    MPI_Comm_dup( access_type->comm , &((*file)->comm));
    MPI_Comm_rank((*file)->comm, &rank);
  }
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
  MPI_Bcast(&((*file)->idx_derived_ptr->fs_block_size), 1, MPI_INT, 0, (*file)->comm);
#endif
    
  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
  (*file)->write_on_close = 0; 
  (*file)->one_time_initializations = 0;
  
  (*file)->idx_ptr->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx_derived_ptr->samples_per_block = pow(2, PIDX_default_bits_per_block);
  (*file)->idx_ptr->blocks_per_file = PIDX_default_blocks_per_file;
    
  write_init_start = malloc (sizeof(double) * 64); memset(write_init_start, 0, sizeof(double) * 64);
  write_init_end = malloc (sizeof(double) * 64);   memset(write_init_end, 0, sizeof(double) * 64);
  rst_init_start = malloc (sizeof(double) * 64);   memset(rst_init_start, 0, sizeof(double) * 64);
  rst_init_end = malloc (sizeof(double) * 64);     memset(rst_init_end, 0, sizeof(double) * 64);
  hz_init_start = malloc (sizeof(double) * 64);    memset(hz_init_start, 0, sizeof(double) * 64);
  hz_init_end = malloc (sizeof(double) * 64);      memset(hz_init_end, 0, sizeof(double) * 64);
  agg_init_start = malloc (sizeof(double) * 64);   memset(agg_init_start, 0, sizeof(double) * 64);
  agg_init_end = malloc (sizeof(double) * 64);     memset(agg_init_end, 0, sizeof(double) * 64);
  io_init_start = malloc (sizeof(double) * 64);    memset(io_init_start, 0, sizeof(double) * 64);
  io_init_end = malloc (sizeof(double) * 64);      memset(io_init_end, 0, sizeof(double) * 64);
  var_init_start = malloc (sizeof(double) * 64);   memset(var_init_start, 0, sizeof(double) * 64);
  var_init_end = malloc (sizeof(double) * 64);     memset(var_init_end, 0, sizeof(double) * 64);
  rst_start = malloc (sizeof(double) * 64);        memset(rst_start, 0, sizeof(double) * 64);
  rst_end = malloc (sizeof(double) * 64);          memset(rst_end, 0, sizeof(double) * 64);
  hz_start = malloc (sizeof(double) * 64);         memset(hz_start, 0, sizeof(double) * 64);
  hz_end = malloc (sizeof(double) * 64);           memset(hz_end, 0, sizeof(double) * 64);
  agg_start = malloc (sizeof(double) * 64);        memset(agg_start, 0, sizeof(double) * 64);
  agg_end = malloc (sizeof(double) * 64);          memset(agg_end, 0, sizeof(double) * 64);
  io_start = malloc (sizeof(double) * 64);         memset(io_start, 0, sizeof(double) * 64);
  io_end = malloc (sizeof(double) * 64);           memset(io_end, 0, sizeof(double) * 64);
  cleanup_start = malloc (sizeof(double) * 64);    memset(cleanup_start, 0, sizeof(double) * 64);
  cleanup_end = malloc (sizeof(double) * 64);      memset(cleanup_end, 0, sizeof(double) * 64);
  finalize_start = malloc (sizeof(double) * 64);   memset(finalize_start, 0, sizeof(double) * 64);
  finalize_end = malloc (sizeof(double) * 64);     memset(finalize_end, 0, sizeof(double) * 64);
    
  return PIDX_success;
}

/// Function to get file descriptor when opening an existing IDX file
PIDX_return_code PIDX_file_open(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_file* file)
{
  if(flags != PIDX_file_rdonly)
    return PIDX_err_unsupported_flags;
    
  if(flags == PIDX_file_rdonly)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
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

  (*file)->idx_ptr->filename = strdup(filename);
  (*file)->idx_ptr->global_bounds = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
  
  (*file)->idx_ptr->variable_count = 0;
  (*file)->idx_ptr->variable_index_tracker = 0;
  
  (*file)->local_variable_index = 0;
  (*file)->local_variable_count = 0;
 
  int var = 0, variable_counter = 0, count = 0, len = 0;
  int rank;
  char *pch, *pch1;
  char line [ 512 ];
  
#if PIDX_HAVE_MPI
  if ((*file)->access->parallel)
  {
    MPI_Comm_dup(access_type->comm, &((*file)->comm));
    MPI_Comm_rank((*file)->comm, &rank);
  }
#endif
  
  if (rank == 0) 
  {
    FILE *fp = fopen(filename, "r");
    while (fgets(line, sizeof (line), fp) != NULL) 
    {
      printf("%s", line);
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
  
  MPI_Bcast((*file)->idx_ptr->global_bounds, 5, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->blocks_per_file), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->bits_per_block), 1, MPI_INT, 0, (*file)->comm);
  MPI_Bcast(&((*file)->idx_ptr->variable_count), 1, MPI_INT, 0, (*file)->comm);

  if(rank != 0)
  {
    for (var = 0; var < (*file)->idx_ptr->variable_count; var++) 
      (*file)->idx_ptr->variable[var] = malloc(sizeof (*((*file)->idx_ptr->variable[var])));
  }
  
  for (var = 0; var < (*file)->idx_ptr->variable_count; var++) 
  {
    MPI_Bcast(&((*file)->idx_ptr->variable[var]->values_per_sample), 1, MPI_INT, 0, (*file)->comm);
    MPI_Bcast((*file)->idx_ptr->variable[var]->var_name, 512, MPI_CHAR, 0, (*file)->comm);
    MPI_Bcast((*file)->idx_ptr->variable[var]->type_name, 512, MPI_CHAR, 0, (*file)->comm); 
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
  long long dims;
  if (PIDX_inner_product(file->idx_ptr->global_bounds, &dims))
    return PIDX_err_size;
  if (dims < file->idx_derived_ptr->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx_derived_ptr->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx_ptr->bits_per_block = getNumBits(file->idx_derived_ptr->samples_per_block) - 1;
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
  
  memcpy(file->idx_ptr->global_bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(long long));

  return PIDX_validate(file);
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_get_dims(PIDX_file file, PIDX_point dims)
{
  if(!file)
    return PIDX_err_file;
  
  memcpy(dims, file->idx_ptr->global_bounds, (sizeof(long long) * PIDX_MAX_DIMENSIONS));
  
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

/////////////////////////////////////////////////
PIDX_return_code PIDX_file_initialize_time_step(PIDX_file file, char* filename, int current_time_step)
{
  int N;
  char dirname[1024], basename[1024];
  int nbits_blocknumber;
  char* directory_path;
  char* data_set_path;
  
  directory_path = (char*) malloc(sizeof (char) * 1024);
  assert(directory_path);
  memset(directory_path, 0, sizeof (char) * 1024);

  data_set_path = (char*) malloc(sizeof (char) * 1024);
  assert(data_set_path);
  memset(data_set_path, 0, sizeof (char) * 1024);

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
  memcpy(variable->patch[variable->patch_count]->Ndim_box_offset, offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_size, dims, PIDX_MAX_DIMENSIONS * sizeof(long long));
  
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
  memcpy(variable->patch[variable->patch_count]->Ndim_box_offset, offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
  memcpy(variable->patch[variable->patch_count]->Ndim_box_size, dims, PIDX_MAX_DIMENSIONS * sizeof(long long));
  
  temp_buffer = read_from_this_buffer;
  variable->patch[variable->patch_count]->Ndim_box_buffer = (unsigned char*)temp_buffer;
  
  variable->data_layout = data_layout;
  variable->patch_count = variable->patch_count + 1;
  return PIDX_success; 
}

/////////////////////////////////////////////////
PIDX_return_code populate_idx_dataset(PIDX_file file)
{
  int var, i, j, p, ctr, counter = 0, file_number = 0, level_count = 1;
    int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };
  
  PointND global_bounds_point;
  global_bounds_point.x = (int) file->idx_ptr->global_bounds[0];
  global_bounds_point.y = (int) file->idx_ptr->global_bounds[1];
  global_bounds_point.z = (int) file->idx_ptr->global_bounds[2];
  global_bounds_point.u = (int) file->idx_ptr->global_bounds[3];
  global_bounds_point.v = (int) file->idx_ptr->global_bounds[4];
  GuessBitmaskPattern(file->idx_ptr->bitSequence, global_bounds_point);
  file->idx_derived_ptr->maxh = strlen(file->idx_ptr->bitSequence);
  
  for (i = 0; i <= file->idx_derived_ptr->maxh; i++)
    file->idx_ptr->bitPattern[i] = RegExBitmaskBit(file->idx_ptr->bitSequence, i);
    
  for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
  {
    block_layout* all_patch_local_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(all_patch_local_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    file->idx_ptr->variable[var]->global_block_layout = (block_layout*) malloc(sizeof (block_layout));
    initialize_block_layout(file->idx_ptr->variable[var]->global_block_layout, file->idx_derived_ptr->maxh, file->idx_ptr->bits_per_block);
    
    for(p = 0 ; p < file->idx_ptr->variable[var]->patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
	bounding_box[0][i] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[i];
	bounding_box[1][i] = file->idx_ptr->variable[var]->patch[p]->Ndim_box_size[i] + file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset[i];
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
#if PIDX_HAVE_MPI
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i], level_count,
		MPI_INT, MPI_BOR, file->comm);
#else
      memcpy(file->idx_ptr->variable[var]->global_block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
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
  
  for (var = file->local_variable_index; var < file->local_variable_index + file->local_variable_count; var++)
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
    free(temp_file_index);
  }
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_read(PIDX_file file)
{
  int j = 0, p, var = 0;
  int do_agg = 1;
  int rank = 0;
  int local_do_rst = 0, global_do_rst = 0;
  int var_used_in_binary_file, total_header_size;
  
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
  /// All variables are written in one go 
  /// (basically no intermittent PIDX_flush involved).
  if(file->local_variable_count == file->idx_ptr->variable_index_tracker && file->write_on_close == 1)
  {
    if (file->idx_ptr->variable_count == -1)
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    else if (file->idx_ptr->variable_count < file->idx_ptr->variable_index_tracker)
    {
      if(rank == 0)
	fprintf(stderr, "(Warning !!!!) Variable count set at %d while attempting to write %d variables\n", file->idx_ptr->variable_count, file->idx_ptr->variable_index_tracker);
      
      file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    }
    file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_count);
    PIDX_header_io_set_communicator(file->header_io_id, file->comm);
    PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
    PIDX_header_io_file_create(file->header_io_id);
    PIDX_header_io_finalize(file->header_io_id);
  }
  else
  {
    file->idx_ptr->variable_count = file->idx_ptr->variable_index_tracker;
    file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_index_tracker);
    PIDX_header_io_set_communicator(file->header_io_id, file->comm);
    PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
    PIDX_header_io_file_create(file->header_io_id);
    PIDX_header_io_finalize(file->header_io_id);
  }
    
#if PIDX_HAVE_MPI

  //if (file->idx_ptr->variable[var]->dump_meta_data_ON == 1)
  //    dump_meta_data(file->idx_ptr->variable[var], file->comm);
  int start_index = 0, end_index = 0;
  file->variable_pipelining_factor = 0;
  
  for (start_index = file->local_variable_index;
      start_index < file->local_variable_index + file->local_variable_count; 
      start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    if (file->idx_ptr->variable[var]->patch_count == 1)
      local_do_rst = 1;
    
    MPI_Allreduce(&local_do_rst, &global_do_rst, 1, MPI_INT, MPI_LOR, file->comm);
    
    global_do_rst = 1;
    if(global_do_rst == 1)
      file->rst_id = PIDX_rst_init(file->comm, file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_agg_set_communicator(file->agg_id, file->comm);
    }
    
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    PIDX_io_set_communicator(file->io_id, file->comm);
    //PIDX_io_file_create(file->io_id, file->idx_ptr->current_time_step, file->idx_ptr->filename, PIDX_WRITE);
    
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = 0;
      if(global_do_rst == 1)
	file->idx_ptr->variable[var]->patch_group_count = PIDX_rst_set_restructuring_box(file->rst_id, 0, NULL);
      else
	file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
    
    if(global_do_rst == 0)
    {
      for (var = start_index; var <= end_index; var++)
      {
	for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	{
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->count = 1;
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->type = 0;
	  
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->count);
	  for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->count; j++)
	  {
	    file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j])));
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(int));
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(int));
	    file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
	  }
	}
      }
    }
    else
    {
      PIDX_rst_restructure(file->rst_id, file->idx_ptr->variable);
      //PIDX_rst_buf_destroy(file->rst_id, file->idx_ptr->variable);
    }
    
    PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable);
    //if(global_do_rst == 1)
    //HELPER_Hz_encode(file->hz_id, file->idx_ptr->variable);

    if(do_agg == 1)
    {
      file->agg_buffer = malloc(sizeof(*file->agg_buffer));
      PIDX_agg_aggregate(file->agg_id, file->agg_buffer);
      
      PIDX_io_aggregated_IO(file->io_id, file->agg_buffer, PIDX_READ);
      PIDX_agg_aggregate_write_read(file->agg_id, file->agg_buffer, PIDX_READ);
      
      PIDX_agg_buf_destroy(file->agg_buffer);
    }
    else
      PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_READ);
    
    PIDX_hz_encode_read_var(file->hz_id, file->idx_ptr->variable);
    PIDX_rst_restructure_IO(file->rst_id, file->idx_ptr->variable, PIDX_READ);
    
    PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
    
    PIDX_io_finalize(file->io_id);
    
    if(do_agg == 1)
      PIDX_agg_finalize(file->agg_id);
    
    PIDX_hz_encode_finalize(file->hz_id);
  }
#else

  int start_index = 0, end_index = 0;
  file->variable_pipelining_factor = 0;
  
  for (start_index = file->local_variable_index;
      start_index < file->local_variable_index + file->local_variable_count; 
      start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_agg_set_communicator(file->agg_id, file->comm);
    }
    
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    PIDX_io_set_communicator(file->io_id, file->comm);
    
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
    
    
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
	file->idx_ptr->variable[var]->patch_group_ptr[p]->count = 1;
	file->idx_ptr->variable[var]->patch_group_ptr[p]->type = 0;
	
	file->idx_ptr->variable[var]->patch_group_ptr[p]->block = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->count);
	for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->count; j++)
	{
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j])));
	  memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->offset, file->idx_ptr->variable[var]->patch[p]->offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
	  memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->count, file->idx_ptr->variable[var]->patch[p]->count, PIDX_MAX_DIMENSIONS * sizeof(long long));
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->buffer = file->idx_ptr->variable[var]->patch[p]->buffer;
	}
      }
    }
    
    
    
    PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable);
    //if(global_do_rst == 1)
    //HELPER_Hz_encode(file->hz_id, file->idx_ptr->variable);

    if(do_agg == 1)
    {
      file->agg_buffer = malloc(sizeof(*file->agg_buffer));
      PIDX_agg_aggregate(file->agg_id, file->agg_buffer);
      
      PIDX_io_aggregated_IO(file->io_id, file->agg_buffer, PIDX_READ);
      PIDX_agg_aggregate_write_read(file->agg_id, file->agg_buffer, PIDX_READ);
      
      PIDX_agg_buf_destroy(file->agg_buffer);
    }
    else
      PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_READ);
    
    PIDX_hz_encode_read_var(file->hz_id, file->idx_ptr->variable);
    
    PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
    
    PIDX_io_finalize(file->io_id);
    
    if(do_agg == 1)
      PIDX_agg_finalize(file->agg_id);
    
    PIDX_hz_encode_finalize(file->hz_id);
  }

  
#endif

  return PIDX_success;
}

/////////////////////////////////////////////////
int dump_meta_data(PIDX_variable variable 
#if PIDX_HAVE_MPI
		   , MPI_Comm comm
#endif
		  )
{
  int i, nprocs = 1;
  FILE* meta_data_file;
  long long *rank_r_offset, *rank_r_count;
  
  rank_r_offset = malloc(sizeof (long long) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_offset, 0, (sizeof (long long) * nprocs * PIDX_MAX_DIMENSIONS));

  rank_r_count = malloc(sizeof (long long) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_count, 0, (sizeof (long long) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  MPI_Comm_size(comm, &nprocs);
  MPI_Allgather(variable->patch[0]->Ndim_box_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
  MPI_Allgather(variable->patch[0]->Ndim_box_size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, MPI_COMM_WORLD);
#else
  memcpy(rank_r_offset, variable->patch[0]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
  memcpy(rank_r_count, variable->patch[0]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(long long));
#endif
  
  meta_data_file = fopen("Box_info.txt", "w");
  if (!meta_data_file) 
    return PIDX_err_name;
  
  for(i = 0; i < nprocs; i++)
    fprintf(meta_data_file, "[%d]: %lld %lld %lld %lld %lld - %lld %lld %lld %lld %lld\n", i, rank_r_offset[PIDX_MAX_DIMENSIONS * i + 0], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 1], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 2], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 3], rank_r_offset[PIDX_MAX_DIMENSIONS * i + 4], rank_r_count[PIDX_MAX_DIMENSIONS * i + 0], rank_r_count[PIDX_MAX_DIMENSIONS * i + 1], rank_r_count[PIDX_MAX_DIMENSIONS * i + 2], rank_r_count[PIDX_MAX_DIMENSIONS * i + 3], rank_r_count[PIDX_MAX_DIMENSIONS * i + 4]);
  
  fclose(meta_data_file);
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_write(PIDX_file file)
{
  if (file->local_variable_index == file->idx_ptr->variable_count)
    return PIDX_success;
    
  int j = 0, p, var = 0;
  int do_agg = 1;
  int debug_rst = 0, debug_hz = 0;
  int rank = 0;
  int local_do_rst = 0, global_do_rst = 0;
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
    file->header_io_id = PIDX_header_io_init(file->idx_ptr, file->idx_derived_ptr, 0, file->idx_ptr->variable_count);
    PIDX_header_io_set_communicator(file->header_io_id, file->comm);
    PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
    PIDX_header_io_file_create(file->header_io_id);
    PIDX_header_io_file_write(file->header_io_id);
    PIDX_header_io_finalize(file->header_io_id);
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
    PIDX_header_io_set_communicator(file->header_io_id, file->comm);
    
    if(hp == 0)
      PIDX_header_io_file_create(file->header_io_id);
    
    if (file->idx_ptr->variable_count == -1 || (file->idx_ptr->variable_count == file->idx_ptr->variable_index_tracker))
    {
      PIDX_header_io_write_idx (file->header_io_id, file->idx_ptr->filename, file->idx_ptr->current_time_step);
      PIDX_header_io_file_write(file->header_io_id);
    }
        
    PIDX_header_io_finalize(file->header_io_id);
    
    write_init_end[hp] = PIDX_get_time();
    hp++;
  }
  
    
#if PIDX_HAVE_MPI

  //if (file->idx_ptr->variable[var]->dump_meta_data_ON == 1)
  //    dump_meta_data(file->idx_ptr->variable[var], file->comm);
  int start_index = 0, end_index = 0;
  file->variable_pipelining_factor = 0;
  for (start_index = file->local_variable_index; start_index < file->local_variable_index + file->local_variable_count; start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    ///----------------------------------- RST init start----------------------------------------------///
    rst_init_start[vp] = PIDX_get_time();                                            
    if (file->idx_ptr->variable[var]->patch_count == 1)
      local_do_rst = 1;  
    
    MPI_Allreduce(&local_do_rst, &global_do_rst, 1, MPI_INT, MPI_LOR, file->comm);
    if(global_do_rst == 1)
      file->rst_id = PIDX_rst_init(file->comm, file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    rst_init_end[vp] = PIDX_get_time();
    ///----------------------------------- RST init end------------------------------------------------///
    
    
    ///------------------------------------HZ init start-----------------------------------------------///
    hz_init_start[vp] = PIDX_get_time();                                             
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    hz_init_end[vp] = PIDX_get_time();
    ///------------------------------------HZ init end-------------------------------------------------///
    
    
    ///-----------------------------------AGG init start-----------------------------------------------///
    agg_init_start[vp] = PIDX_get_time();                                            
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_agg_set_communicator(file->agg_id, file->comm);
    }
    agg_init_end[vp] = PIDX_get_time();
    ///-----------------------------------AGG init end-------------------------------------------------///
    
    
    ///----------------------------------IO init start-------------------------------------------------///
    io_init_start[vp] = PIDX_get_time();
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    PIDX_io_set_communicator(file->io_id, file->comm);
    io_init_end[vp] = PIDX_get_time();
    ///----------------------------------IO init end---------------------------------------------------///
    
    
    ///------------------------------Var buffer init start---------------------------------------------///
    var_init_start[vp] = PIDX_get_time();
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = 0;
      if(global_do_rst == 1)
	file->idx_ptr->variable[var]->patch_group_count = PIDX_rst_set_restructuring_box(file->rst_id, 0, NULL);
      else
	file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
      
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
    var_init_end[vp] = PIDX_get_time();
    ///------------------------------Var buffer init end--------------------------------------------------///
    
    
    ///------------------------------------RST start time-------------------------------------------------///
    rst_start[vp] = PIDX_get_time();
    if(global_do_rst == 0)
    {
      for (var = start_index; var <= end_index; var++)
      {
	for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	{
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->count = 1;
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->type = 0;
	  
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->count);
	  for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->count; j++)
	  {
	    file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j])));
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_size, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(long long));
	    
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->power_two_offset, file->idx_ptr->variable[var]->patch[p]->Ndim_box_offset, PIDX_MAX_DIMENSIONS * sizeof(long long));
	    memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->power_two_count, file->idx_ptr->variable[var]->patch[p]->Ndim_box_size, PIDX_MAX_DIMENSIONS * sizeof(long long));
	    file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->Ndim_box_buffer = file->idx_ptr->variable[var]->patch[p]->Ndim_box_buffer;
	  }
	}
      }
    }
    else
    {
      PIDX_rst_restructure(file->rst_id, file->idx_ptr->variable);
      PIDX_rst_restructure_IO(file->rst_id, file->idx_ptr->variable, PIDX_WRITE);
      if(global_do_rst == 1 && debug_rst == 1)
	HELPER_rst(file->rst_id, file->idx_ptr->variable);
    }
    rst_end[vp] = PIDX_get_time();
    ///--------------------------------------RST end time---------------------------------------------------///
    
    
    ///-------------------------------------HZ start time---------------------------------------------------///
    hz_start[vp] = PIDX_get_time();
    PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable);
    PIDX_hz_encode_write_var(file->hz_id, file->idx_ptr->variable);
    if(global_do_rst == 1 && debug_hz == 1)
      HELPER_Hz_encode(file->hz_id, file->idx_ptr->variable);
    if(global_do_rst == 1)
      PIDX_rst_buf_destroy(file->rst_id);
    hz_end[vp] = PIDX_get_time();
    ///------------------------------------HZ end time------------------------------------------------------///
    
    
    ///------------------------------------Agg start time---------------------------------------------------///
    agg_start[vp] = PIDX_get_time();
    if(do_agg == 1)
    {
      file->agg_buffer = malloc(sizeof(*file->agg_buffer));
      PIDX_agg_aggregate(file->agg_id, file->agg_buffer);
      PIDX_agg_aggregate_write_read(file->agg_id, file->agg_buffer, PIDX_WRITE);
      PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
    }
    agg_end[vp] = PIDX_get_time();
    ///---------------------------------------Agg end time---------------------------------------------------///
     
    
    ///---------------------------------------IO start time--------------------------------------------------///
    io_start[vp] = PIDX_get_time();
    if(do_agg == 1)
    {
      PIDX_io_aggregated_IO(file->io_id, file->agg_buffer, PIDX_WRITE);
      PIDX_agg_buf_destroy(file->agg_buffer);
      free(file->agg_buffer);
    }
    else
      PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_WRITE);
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
	  for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->count; j++)
	    free(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]);
	  free(file->idx_ptr->variable[var]->patch_group_ptr[p]->block);
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
    if(global_do_rst == 1)
      PIDX_rst_finalize(file->rst_id);
    finalize_end[vp] = PIDX_get_time();
    ///------------------------------------finalize end time------------------------------------------------///
    
    vp++;
  }
  

#else

  //if (file->idx_ptr->variable[var]->dump_meta_data_ON == 1)
  //    dump_meta_data(file->idx_ptr->variable[var], file->comm);
  int start_index = 0, end_index = 0;
  file->variable_pipelining_factor = 0;
  
  for (start_index = file->local_variable_index;
      start_index < file->local_variable_index + file->local_variable_count; 
      start_index = start_index + (file->variable_pipelining_factor + 1))
  {
    end_index = ((start_index + file->variable_pipelining_factor) >= (file->local_variable_index + file->local_variable_count)) ? ((file->local_variable_index + file->local_variable_count) - 1) : (start_index + file->variable_pipelining_factor);
    
    file->hz_id = PIDX_hz_encode_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    
    if(do_agg == 1)
    {
      file->agg_id = PIDX_agg_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
      PIDX_agg_set_communicator(file->agg_id, file->comm);
    }
    
    file->io_id = PIDX_io_init(file->idx_ptr, file->idx_derived_ptr, start_index, end_index);
    PIDX_io_set_communicator(file->io_id, file->comm);
    
    for (var = start_index; var <= end_index; var++)
    {
      file->idx_ptr->variable[var]->patch_group_count = file->idx_ptr->variable[var]->patch_count;
      
      file->idx_ptr->variable[var]->patch_group_ptr = malloc(file->idx_ptr->variable[var]->patch_group_count * sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr)));
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
	file->idx_ptr->variable[var]->patch_group_ptr[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p])));
      
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
        file->idx_ptr->variable[var]->HZ_patch[p] = malloc(sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
        assert(file->idx_ptr->variable[var]->HZ_patch[p] != NULL);
        memset(file->idx_ptr->variable[var]->HZ_patch[p], 0, sizeof(*(file->idx_ptr->variable[var]->HZ_patch[p])));
      }
    }
    
    for (var = start_index; var <= end_index; var++)
    {
      for (p = 0; p < file->idx_ptr->variable[var]->patch_group_count; p++)
      {
	file->idx_ptr->variable[var]->patch_group_ptr[p]->count = 1;
	file->idx_ptr->variable[var]->patch_group_ptr[p]->type = 0;
	
	file->idx_ptr->variable[var]->patch_group_ptr[p]->block = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block)) * file->idx_ptr->variable[var]->patch_group_ptr[p]->count);
	for(j = 0; j < file->idx_ptr->variable[var]->patch_group_ptr[p]->count; j++)
	{
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j] = malloc(sizeof(*(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j])));
	  memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->offset, file->idx_ptr->variable[var]->patch[p]->offset, PIDX_MAX_DIMENSIONS * sizeof(int));
	  memcpy(file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->count, file->idx_ptr->variable[var]->patch[p]->count, PIDX_MAX_DIMENSIONS * sizeof(int));
	  file->idx_ptr->variable[var]->patch_group_ptr[p]->block[j]->buffer = file->idx_ptr->variable[var]->patch[p]->buffer;
	}
      }
    }
    
    PIDX_hz_encode_var(file->hz_id, file->idx_ptr->variable);
    PIDX_hz_encode_write_var(file->hz_id, file->idx_ptr->variable);
    
    if(do_agg == 1)
    {
      file->agg_buffer = malloc(sizeof(*file->agg_buffer));
      PIDX_agg_aggregate(file->agg_id, file->agg_buffer);
      PIDX_agg_aggregate_write_read(file->agg_id, file->agg_buffer, PIDX_WRITE);
      
      PIDX_io_aggregated_IO(file->io_id, file->agg_buffer, PIDX_WRITE);
      PIDX_agg_buf_destroy(file->agg_buffer);
    }
    else
      PIDX_io_independent_IO_var(file->io_id, file->idx_ptr->variable, PIDX_WRITE);
    
    PIDX_hz_encode_buf_destroy_var(file->hz_id, file->idx_ptr->variable);
    PIDX_io_finalize(file->io_id);
    
    if(do_agg == 1)
      PIDX_agg_finalize(file->agg_id);
    
    PIDX_hz_encode_finalize(file->hz_id);
  }
  
#endif
  
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_flush(PIDX_file file)
{
  PIDX_write(file);
  PIDX_cleanup(file);
  return PIDX_success;
}

/////////////////////////////////////////////////
PIDX_return_code PIDX_cleanup(PIDX_file file)
{
  int i, p;
  for (i = file->local_variable_index; i < file->local_variable_index + file->local_variable_count; i++) 
  {
    free(file->idx_ptr->variable[i]->existing_file_index);
    file->idx_ptr->variable[i]->existing_file_index = 0;
      
    for(p = 0; p < file->idx_ptr->variable[i]->patch_count; p++)
    {   
      free(file->idx_ptr->variable[i]->HZ_patch[p]);
      file->idx_ptr->variable[i]->HZ_patch[p] = 0;
    
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
  
  double total_time = sim_end - sim_start;
  double max_time = total_time;
  int sample_sum = 0, var = 0, i = 0;
  
#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, file->comm);
#endif

  if (max_time == total_time) 
  {
    if (file->IDX_WRITE == 1 && file->IDX_READ == 0)
      fprintf(stdout, "\n------------- WRITE -------------\n");

    if (file->IDX_WRITE == 0 && file->IDX_READ == 1)
      fprintf(stdout, "\n------------- READ -------------\n");
      
    for (var = 0; var < file->idx_ptr->variable_count; var++)
      sample_sum = sample_sum + file->idx_ptr->variable[var]->values_per_sample;
    
    long long total_data = file->idx_ptr->global_bounds[0] * file->idx_ptr->global_bounds[1] * file->idx_ptr->global_bounds[2] * file->idx_ptr->global_bounds[3] * file->idx_ptr->global_bounds[4] * sample_sum * 8;
    fprintf(stdout, "\n=======================================================================================\n");
    fprintf(stdout, "Global Data [%lld %lld %lld] Variables [%d]\nTime Taken: %f Seconds Throughput %f MiB/sec\n", file->idx_ptr->global_bounds[0], file->idx_ptr->global_bounds[1], file->idx_ptr->global_bounds[2], file->idx_ptr->variable_count, max_time, (float) total_data / (1024 * 1024 * max_time));
    fprintf(stdout, "---------------------------------------------------------------------------------------\n");
    //printf("File creation time %f\n", write_init_end - write_init_start);
    
    for (var = 0; var < hp; var++)
      fprintf(stdout, "File Create time (+ header IO) %f\n", (write_init_end[var] - write_init_start[var]));
    
    for (var = 0; var < vp; var++)
    {
      fprintf(stdout, "---------------------------------------VG %d (START)------------------------------------\n", var);
      fprintf(stdout, "Buffer init time %f\n", (var_init_end[var] - var_init_start[var]));
      fprintf(stdout, "Init time [RST + HZ + AGG + IO] [%f + %f + %f + %f] = %f\n", (rst_init_end[var] - rst_init_start[var]), (hz_init_end[var] - hz_init_start[var]), (agg_init_end[var] - agg_init_start[var]), (io_init_end[var] - io_init_start[var]), (rst_init_end[var] - rst_init_start[var]) + (hz_init_end[var] - hz_init_start[var]) + (agg_init_end[var] - agg_init_start[var]) + (io_init_end[var] - io_init_start[var]));
      
      fprintf(stdout, "Write time [RST + HZ + AGG + IO] [%f + %f + %f + %f] = %f\n", (rst_end[var] - rst_start[var]), (hz_end[var] - hz_start[var]), (agg_end[var] - agg_start[var]), (io_end[var] - io_start[var]), (rst_end[var] - rst_start[var]) + (hz_end[var] - hz_start[var]) + (agg_end[var] - agg_start[var]) + (io_end[var] - io_start[var]));
      
      fprintf(stdout, "Cleanup time %f\n", cleanup_end[var] - cleanup_start[var]);
      fprintf(stdout, "----------------------------------------VG %d (END)-------------------------------------\n", var);
    }
    fprintf(stdout, "=======================================================================================\n");
  }
  
  vp = 0;
  hp = 0;
  free(write_init_start);  write_init_start = 0;
  free(write_init_end);    write_init_end   = 0;
  free(rst_init_start);    rst_init_start   = 0;
  free(rst_init_end);      rst_init_end     = 0;
  free(hz_init_start);     hz_init_start    = 0;
  free(hz_init_end);       hz_init_end      = 0;
  free(agg_init_start);    agg_init_start   = 0;
  free(agg_init_end);      agg_init_end     = 0;
  free(io_init_start);     io_init_start    = 0;
  free(io_init_end);       io_init_end      = 0;
  free(var_init_start);    var_init_start   = 0;
  free(var_init_end);      var_init_end     = 0;
  free(cleanup_start);     cleanup_start    = 0;
  free(cleanup_end);       cleanup_end      = 0;
  free(finalize_start);    finalize_start   = 0;
  free(finalize_end);      finalize_end     = 0;
  
  for (i = 0; i < file->idx_ptr->variable_count; i++) 
  {
    free(file->idx_ptr->variable[i]->blocks_per_file);
    file->idx_ptr->variable[i]->blocks_per_file = 0;
    
    destroyBlockBitmap(file->idx_ptr->variable[i]->global_block_layout);
    free(file->idx_ptr->variable[i]->global_block_layout);
    file->idx_ptr->variable[i]->global_block_layout = 0;
  }
  
  for (i = 0; i < 1024; i++)
  {
    free(file->idx_ptr->variable[i]);  
    file->idx_ptr->variable[i] = 0; 
  }
  
  file->idx_ptr->variable_count = 0;
  
  free(file->idx_ptr->filename);              file->idx_ptr->filename = 0;
  free(file->idx_ptr->global_bounds);         file->idx_ptr->global_bounds = 0;
  free(file->idx_ptr);                        file->idx_ptr = 0;
  free(file->idx_derived_ptr->file_bitmap);   file->idx_derived_ptr->file_bitmap = 0;
  free(file->idx_derived_ptr);                file->idx_derived_ptr = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_free(&(file->comm));
#endif
  
  free(file);
  
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
