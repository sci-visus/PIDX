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

#include "PIDX_rst.h"
#define PIDX_MAX_NEIGHBOR_PROC 1024

struct NDim_chunk_bound 
{
  int lower_bound[PIDX_MAX_DIMENSIONS];
  int upper_bound[PIDX_MAX_DIMENSIONS];
};
typedef struct NDim_chunk_bound* NDim_chunk;

struct Ndim_power_two_buffer_struct
{
  int type;
  int count;
  Ndim_buffer block[512];
  int rank[512];
  int recieve[512];
  int max_rank;
  int power_two_offset[PIDX_MAX_DIMENSIONS];
  int power_two_count[PIDX_MAX_DIMENSIONS];
};
typedef struct Ndim_power_two_buffer_struct* Ndim_power_two_buffer;

//Struct for restructuring ID
struct PIDX_rst_struct 
{
  //Passed by PIDX API
  MPI_Comm comm; //Communicator

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  //dimension of the power-two volume imposed box
  int regular_box_dim[PIDX_MAX_DIMENSIONS];
  
  int owned_regular_box_count;
  
  Ndim_power_two_buffer *regular_box_buffer;
  
  int start_variable_index;
  int end_variable_index;
  
};

//Function to check if NDimensional data chunks A and B intersects
int intersectNDChunk(NDim_chunk A, NDim_chunk B) 
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
    check_bit = check_bit || A->upper_bound[d] < B->lower_bound[d] || B->upper_bound[d] < A->lower_bound[d];
  
  return !(check_bit);
}

//Function to find the power of 2 of an integer value (example 5->8)
int getPowerOftwo(int x) 
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}

//Function to find the dimension of the imposing regular box
void set_default_box_size(PIDX_rst_id rst_id, int* process_bounds, int nprocs) 
{
  int i = 0, average_count = 0, j = 0;
  int check_bit = 0;
  int* max_dim_length;
  int equal_partiton = 1;

  max_dim_length = (int*) malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
  assert(max_dim_length);
  memset(max_dim_length, 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++) 
  {
    max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * 0 + i];
    for (j = 0; j < nprocs; j++) 
    {
      if (max_dim_length[i] <= process_bounds[PIDX_MAX_DIMENSIONS * j + i])
	max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * j + i];
    }
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    average_count = average_count + max_dim_length[i];
  }
  average_count = average_count / PIDX_MAX_DIMENSIONS;
  average_count = getPowerOftwo(average_count);
  
  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    check_bit = check_bit || ((double) rst_id->idx_ptr->global_bounds[i] / average_count > (double) rst_id->idx_ptr->global_bounds[i] / max_dim_length[i]);

  while (check_bit) 
  {
    average_count = average_count * 2;
    check_bit = 0;
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      check_bit = check_bit || ((double) rst_id->idx_ptr->global_bounds[i] / average_count > (double) rst_id->idx_ptr->global_bounds[i] / max_dim_length[i]);
  }
  //regular_box_dim =  average_count;
  if (equal_partiton == 1) 
  {
    rst_id->regular_box_dim[0] = average_count * 1;
    rst_id->regular_box_dim[1] = average_count * 1;
    rst_id->regular_box_dim[2] = average_count * 1;
    rst_id->regular_box_dim[3] = average_count * 1;
    rst_id->regular_box_dim[4] = average_count * 1;
  } 
  else 
  {
    rst_id->regular_box_dim[0] = getPowerOftwo(process_bounds[0]) * 1;
    rst_id->regular_box_dim[1] = getPowerOftwo(process_bounds[1]) * 1;
    rst_id->regular_box_dim[2] = getPowerOftwo(process_bounds[2]) * 1;
    rst_id->regular_box_dim[3] = getPowerOftwo(process_bounds[3]) * 1;
    rst_id->regular_box_dim[4] = getPowerOftwo(process_bounds[4]) * 1;
  }
  free(max_dim_length);
  max_dim_length = 0;
  //regular_box_dim = regular_box_dim * 4;
}

int* PIDX_rst_get_box_dimension(PIDX_rst_id id) 
{
  return id->regular_box_dim;
}

/* output value: num_output_buffers (number of buffers this process will hold after restructuring given the above parameters) */
PIDX_rst_id PIDX_rst_init(MPI_Comm comm, idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int var_start_index, int var_end_index)
{
  int ret;
  //Creating the restructuring ID
  PIDX_rst_id rst_id;
  rst_id = (PIDX_rst_id)malloc(sizeof (*rst_id));
  if (!rst_id) PIDX_rst_print_error("Error Creating Rst ID", __FILE__, __LINE__);
  memset(rst_id, 0, sizeof (*rst_id));

  rst_id->idx_ptr = (idx_dataset)malloc(sizeof(*(rst_id->idx_ptr)));
  memcpy(rst_id->idx_ptr, idx_meta_data, sizeof(*(rst_id->idx_ptr)));

  rst_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(rst_id->idx_derived_ptr)));
  memcpy(rst_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(rst_id->idx_derived_ptr)));

  rst_id->start_variable_index = var_start_index;
  rst_id->end_variable_index = var_end_index;
  
  ret = MPI_Comm_dup(comm, &rst_id->comm);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("Communicator Duplication", __FILE__, __LINE__);

  return (rst_id);
}

int PIDX_rst_set_restructuring_box(PIDX_rst_id rst_id, int set_box_dim, int* box_dim)
{
  int num_output_buffers;
  int r, d, i, j, k, l, m, c, nprocs, rank, ret;
  int *rank_r_offset, *rank_r_count;
  int max_rank, max_vol, regular_box_count;
  
  ret = MPI_Comm_rank(rst_id->comm, &rank);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("Rank ", __FILE__, __LINE__);

  ret = MPI_Comm_size(rst_id->comm, &nprocs);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("nprocs ", __FILE__, __LINE__);

  //creating rank_r_count and rank_r_offset to hold the offset and count of every process
  rst_id->owned_regular_box_count = 0;

  rank_r_offset = (int*) malloc(sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS);
  if (!rank_r_offset) PIDX_rst_print_error("Memory : rank_r_offset", __FILE__, __LINE__);
  memset(rank_r_offset, 0, (sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS));

  rank_r_count = (int*) malloc(sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS);
  if (!rank_r_count) PIDX_rst_print_error("Memory : rank_r_count", __FILE__, __LINE__);
  memset(rank_r_count, 0, (sizeof (int) * nprocs * PIDX_MAX_DIMENSIONS));

  //STEP 1 : Doing an all to all Communication to get extents of all processes.
  ret = MPI_Allgather(rst_id->idx_ptr->variable[0]->patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_INT, rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_INT, MPI_COMM_WORLD);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Allgather : rank_r_offset", __FILE__, __LINE__);

  ret = MPI_Allgather(rst_id->idx_ptr->variable[0]->patch[0]->count, PIDX_MAX_DIMENSIONS, MPI_INT, rank_r_count, PIDX_MAX_DIMENSIONS, MPI_INT, MPI_COMM_WORLD);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Allgather : rank_r_count", __FILE__, __LINE__);

  //STEP 2 : Compute the dimension of the regular BOX
  if(set_box_dim == 0)
   set_default_box_size(rst_id, rank_r_count, nprocs);
  else
    memcpy(rst_id->regular_box_dim, box_dim, PIDX_MAX_DIMENSIONS * sizeof(int));
    
  if(rank == 0)
    printf("[%d] Imposed Box Dimension : %d %d %d %d %d\n", rank, rst_id->regular_box_dim[0], rst_id->regular_box_dim[1], rst_id->regular_box_dim[2],
	  rst_id->regular_box_dim[3], rst_id->regular_box_dim[4]);
  
  //extents for the local process(rank)
  NDim_chunk local_proc_bound = malloc(sizeof (*local_proc_bound));
  if (!local_proc_bound) 
  {
    fprintf(stderr, "[Rank : %d] [File : %s] [Line : %d] local_proc_bound\n", rank, __FILE__, __LINE__);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  memset(local_proc_bound, 0, sizeof (*local_proc_bound));
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
  {
    local_proc_bound->lower_bound[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * rank + d];
    local_proc_bound->upper_bound[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * rank + d] + rank_r_count[PIDX_MAX_DIMENSIONS * rank + d] - 1;
  }

  rst_id->owned_regular_box_count = 0;
  for (i = 0; i < rst_id->idx_ptr->global_bounds[0]; i = i + rst_id->regular_box_dim[0])
    for (j = 0; j < rst_id->idx_ptr->global_bounds[1]; j = j + rst_id->regular_box_dim[1])
      for (k = 0; k < rst_id->idx_ptr->global_bounds[2]; k = k + rst_id->regular_box_dim[2])
	for (l = 0; l < rst_id->idx_ptr->global_bounds[3]; l = l + rst_id->regular_box_dim[3])
	  for (m = 0; m < rst_id->idx_ptr->global_bounds[4]; m = m + rst_id->regular_box_dim[4]) 
	  {
	    NDim_chunk regular_box_bound = malloc(sizeof (*regular_box_bound));
	    if (!regular_box_bound) PIDX_rst_print_error("Memory : regular_box_bound", __FILE__, __LINE__);
	    memset(regular_box_bound, 0, sizeof (*regular_box_bound));

	    //Interior regular boxes
	    regular_box_bound->lower_bound[0] = i;
	    regular_box_bound->lower_bound[1] = j;
	    regular_box_bound->lower_bound[2] = k;
	    regular_box_bound->lower_bound[3] = l;
	    regular_box_bound->lower_bound[4] = m;
	    regular_box_bound->upper_bound[0] = i + rst_id->regular_box_dim[0] - 1;
	    regular_box_bound->upper_bound[1] = j + rst_id->regular_box_dim[1] - 1;
	    regular_box_bound->upper_bound[2] = k + rst_id->regular_box_dim[2] - 1;
	    regular_box_bound->upper_bound[3] = l + rst_id->regular_box_dim[3] - 1;
	    regular_box_bound->upper_bound[4] = m + rst_id->regular_box_dim[4] - 1;

	    //Edge regular boxes
	    if ((i + rst_id->regular_box_dim[0]) > rst_id->idx_ptr->global_bounds[0])
		regular_box_bound->upper_bound[0] = rst_id->idx_ptr->global_bounds[0] - 1;
	    if ((j + rst_id->regular_box_dim[1]) > rst_id->idx_ptr->global_bounds[1])
		regular_box_bound->upper_bound[1] = rst_id->idx_ptr->global_bounds[1] - 1;
	    if ((k + rst_id->regular_box_dim[2]) > rst_id->idx_ptr->global_bounds[2])
		regular_box_bound->upper_bound[2] = rst_id->idx_ptr->global_bounds[2] - 1;
	    if ((l + rst_id->regular_box_dim[3]) > rst_id->idx_ptr->global_bounds[3])
		regular_box_bound->upper_bound[3] = rst_id->idx_ptr->global_bounds[3] - 1;
	    if ((m + rst_id->regular_box_dim[4]) > rst_id->idx_ptr->global_bounds[4])
		regular_box_bound->upper_bound[4] = rst_id->idx_ptr->global_bounds[4] - 1;

	    //STEP 4: If local process intersects with regular box, then find all other process that intersects with the regular box.
	    if (intersectNDChunk(regular_box_bound, local_proc_bound))
	      rst_id->owned_regular_box_count++;
	  }
  
  rst_id->regular_box_buffer = malloc(sizeof(*rst_id->regular_box_buffer) * rst_id->owned_regular_box_count);
  
  regular_box_count = 0;
  //STEP 3 : iterate through extents of all imposed regular boxes, and find all the regular boxes a process (local_proc_bound) intersects with
  for (i = 0; i < rst_id->idx_ptr->global_bounds[0]; i = i + rst_id->regular_box_dim[0])
    for (j = 0; j < rst_id->idx_ptr->global_bounds[1]; j = j + rst_id->regular_box_dim[1])
      for (k = 0; k < rst_id->idx_ptr->global_bounds[2]; k = k + rst_id->regular_box_dim[2])
	for (l = 0; l < rst_id->idx_ptr->global_bounds[3]; l = l + rst_id->regular_box_dim[3])
	  for (m = 0; m < rst_id->idx_ptr->global_bounds[4]; m = m + rst_id->regular_box_dim[4]) 
	  {
	    NDim_chunk regular_box_bound = malloc(sizeof (*regular_box_bound));
	    if (!regular_box_bound) PIDX_rst_print_error("Memory : regular_box_bound", __FILE__, __LINE__);
	    memset(regular_box_bound, 0, sizeof (*regular_box_bound));

	    //Interior regular boxes
	    regular_box_bound->lower_bound[0] = i;
	    regular_box_bound->lower_bound[1] = j;
	    regular_box_bound->lower_bound[2] = k;
	    regular_box_bound->lower_bound[3] = l;
	    regular_box_bound->lower_bound[4] = m;
	    regular_box_bound->upper_bound[0] = i + rst_id->regular_box_dim[0] - 1;
	    regular_box_bound->upper_bound[1] = j + rst_id->regular_box_dim[1] - 1;
	    regular_box_bound->upper_bound[2] = k + rst_id->regular_box_dim[2] - 1;
	    regular_box_bound->upper_bound[3] = l + rst_id->regular_box_dim[3] - 1;
	    regular_box_bound->upper_bound[4] = m + rst_id->regular_box_dim[4] - 1;

	    //Edge regular boxes
	    if ((i + rst_id->regular_box_dim[0]) > rst_id->idx_ptr->global_bounds[0])
		regular_box_bound->upper_bound[0] = rst_id->idx_ptr->global_bounds[0] - 1;
	    if ((j + rst_id->regular_box_dim[1]) > rst_id->idx_ptr->global_bounds[1])
		regular_box_bound->upper_bound[1] = rst_id->idx_ptr->global_bounds[1] - 1;
	    if ((k + rst_id->regular_box_dim[2]) > rst_id->idx_ptr->global_bounds[2])
		regular_box_bound->upper_bound[2] = rst_id->idx_ptr->global_bounds[2] - 1;
	    if ((l + rst_id->regular_box_dim[3]) > rst_id->idx_ptr->global_bounds[3])
		regular_box_bound->upper_bound[3] = rst_id->idx_ptr->global_bounds[3] - 1;
	    if ((m + rst_id->regular_box_dim[4]) > rst_id->idx_ptr->global_bounds[4])
		regular_box_bound->upper_bound[4] = rst_id->idx_ptr->global_bounds[4] - 1;

	    //STEP 4: If local process intersects with regular box, then find all other process that intersects with the regular box.
	    if (intersectNDChunk(regular_box_bound, local_proc_bound))
	    {      
	      rst_id->regular_box_buffer[regular_box_count] = malloc(sizeof(*(rst_id->regular_box_buffer[regular_box_count])));
	      rst_id->regular_box_buffer[regular_box_count]->count = 0;
	      
	      //Iterate through all processes
	      for (r = 0; r < nprocs; r++)
	      {
		//Extent of process with rank r
		NDim_chunk rank_r_bound = malloc(sizeof (*rank_r_bound));
		if (!rank_r_bound) PIDX_rst_print_error("Memory : rank_r_bound", __FILE__, __LINE__);
		memset(rank_r_bound, 0, sizeof (*rank_r_bound));

		for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
		{
		  rank_r_bound->lower_bound[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * r + d];
		  rank_r_bound->upper_bound[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * r + d] + rank_r_count[PIDX_MAX_DIMENSIONS * r + d] - 1;
		}

		//If process with rank r intersects with the regular box, then calculate the offset, count and volume of the intersecting volume
		if (intersectNDChunk(regular_box_bound, rank_r_bound)) 
		{
		  //if(rank == 0)
		  //printf("[%d]: (%d %d %d %d %d :: %d %d %d %d %d) -- (%d %d %d %d %d :: %d %d %d %d %d)\n", rank, regular_box_bound->lower_bound[0], regular_box_bound->lower_bound[1], regular_box_bound->lower_bound[2], regular_box_bound->lower_bound[3], regular_box_bound->lower_bound[4], regular_box_bound->upper_bound[0], regular_box_bound->upper_bound[1], regular_box_bound->upper_bound[2], regular_box_bound->upper_bound[3], regular_box_bound->upper_bound[4], rank_r_bound->lower_bound[0], rank_r_bound->lower_bound[1], rank_r_bound->lower_bound[2], rank_r_bound->lower_bound[3], rank_r_bound->lower_bound[4], rank_r_bound->upper_bound[0], rank_r_bound->upper_bound[1], rank_r_bound->upper_bound[2], rank_r_bound->upper_bound[3], rank_r_bound->upper_bound[4]);
		  
		  rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count] = malloc(sizeof(*(rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count])));
		  
		  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
		  {
		    //STEP 5 : offset and count of intersecting chunk of process with rank r and regular box
		    if (rank_r_bound->lower_bound[d] <= regular_box_bound->lower_bound[d] && rank_r_bound->upper_bound[d] <= regular_box_bound->upper_bound[d]) 
		    {
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->offset[d] = regular_box_bound->lower_bound[d];
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->count[d] = rank_r_bound->upper_bound[d] - regular_box_bound->lower_bound[d] + 1;    
		    } 
		    else if (regular_box_bound->lower_bound[d] <= rank_r_bound->lower_bound[d] && rank_r_bound->upper_bound[d] >= regular_box_bound->upper_bound[d]) 
		    {
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->offset[d] = rank_r_bound->lower_bound[d];
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->count[d] = regular_box_bound->upper_bound[d] - rank_r_bound->lower_bound[d] + 1;
		    } 
		    else if (regular_box_bound->upper_bound[d] <= rank_r_bound->upper_bound[d] && regular_box_bound->lower_bound[d] >= rank_r_bound->lower_bound[d]) 
		    {
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->offset[d] = regular_box_bound->lower_bound[d];
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->count[d] = regular_box_bound->upper_bound[d] - regular_box_bound->lower_bound[d] + 1;
		    } 
		    else if (rank_r_bound->upper_bound[d] <= regular_box_bound->upper_bound[d] && rank_r_bound->lower_bound[d] >= regular_box_bound->lower_bound[d]) 
		    {
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->offset[d] = rank_r_bound->lower_bound[d];
		      rst_id->regular_box_buffer[regular_box_count]->block[rst_id->regular_box_buffer[regular_box_count]->count]->count[d] = rank_r_bound->upper_bound[d] - rank_r_bound->lower_bound[d] + 1;
		    }
		    //offset and count of intersecting regular box
		    
		    rst_id->regular_box_buffer[regular_box_count]->power_two_offset[d] = regular_box_bound->lower_bound[d];
		    rst_id->regular_box_buffer[regular_box_count]->power_two_count[d] = regular_box_bound->upper_bound[d] - regular_box_bound->lower_bound[d] + 1;		    
		  }

		  rst_id->regular_box_buffer[regular_box_count]->rank[rst_id->regular_box_buffer[regular_box_count]->count] = r;		  
		  rst_id->regular_box_buffer[regular_box_count]->count++;
		}
		free(rank_r_bound);
	      }
	      
	      max_rank = rst_id->regular_box_buffer[regular_box_count]->rank[0];
	      max_vol = 1;
	      for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
		max_vol = max_vol * rst_id->regular_box_buffer[regular_box_count]->block[0]->count[d];
	      int c_vol = 1;
	      for(c = 1; c < rst_id->regular_box_buffer[regular_box_count]->count ; c++)
	      {
		c_vol = 1;
		for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
		  c_vol = c_vol * rst_id->regular_box_buffer[regular_box_count]->block[c]->count[d];
		if(c_vol > max_vol)
		{
		  max_vol = c_vol;
		  max_rank = rst_id->regular_box_buffer[regular_box_count]->rank[c];
		}
	      }

	      if(rank == max_rank)
		num_output_buffers = num_output_buffers + 1;
	      
	      for(c = 0; c < rst_id->regular_box_buffer[regular_box_count]->count ; c++)
	      {
		if (rst_id->regular_box_buffer[regular_box_count]->rank[c] == max_rank)
		  rst_id->regular_box_buffer[regular_box_count]->recieve[c] = 1;
		else
		  rst_id->regular_box_buffer[regular_box_count]->recieve[c] = 0;
	      }
	      rst_id->regular_box_buffer[regular_box_count]->max_rank = max_rank;
	      
	      regular_box_count++;
	    }
	    free(regular_box_bound);
	  }

  free(local_proc_bound);
  free(rank_r_offset);
  free(rank_r_count);
  
  return num_output_buffers;
}

/* actually do the restructuring, using pre-calculated data associated with the rst_id */
int PIDX_rst_restructure(PIDX_rst_id rst_id, int samples_per_variable, MPI_Datatype datatype, Ndim_buffer* in_buf, Ndim_buffer_group* out_buf_array, int num_output_buffers) 
{
  int j = 0, i, cnt = 0, ret = 0;
  int rank, nprocs, bytes_per_sample;

  //rank and nprocs
  ret = MPI_Comm_rank(rst_id->comm, &rank);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("Rank ", __FILE__, __LINE__);

  MPI_Comm_size(rst_id->comm, &nprocs);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("nprocs ", __FILE__, __LINE__);

  //Bytes per sample for this datatype
  MPI_Type_size(datatype, &bytes_per_sample);

  for (i = 0; i < rst_id->owned_regular_box_count; i++)
  {
    if (rank == rst_id->regular_box_buffer[i]->max_rank)
    {
      //printf("[Cnt] %d: %d\n", rank, rst_id->regular_box_buffer[i]->count);
      out_buf_array[cnt]->count = rst_id->regular_box_buffer[i]->count;
      out_buf_array[cnt]->block = malloc(sizeof(*(out_buf_array[cnt]->block)) * rst_id->regular_box_buffer[i]->count);
      for(j = 0; j < rst_id->regular_box_buffer[i]->count; j++)
      {
	out_buf_array[cnt]->block[j] = malloc(sizeof(*(out_buf_array[cnt]->block[j])));
	
	//printf("i = %d j = %d cnt = %d\n", i, j, cnt);
	//printf("O: %d %d %d %d %d\n", rst_id->regular_box_buffer[i]->block[j]->offset[0], rst_id->regular_box_buffer[i]->block[j]->offset[1], rst_id->regular_box_buffer[i]->block[j]->offset[2], rst_id->regular_box_buffer[i]->block[j]->offset[3], rst_id->regular_box_buffer[i]->block[j]->offset[4]);
	//printf("C: %d %d %d %d %d\n", rst_id->regular_box_buffer[i]->block[j]->count[0], rst_id->regular_box_buffer[i]->block[j]->count[1], rst_id->regular_box_buffer[i]->block[j]->count[2], rst_id->regular_box_buffer[i]->block[j]->count[3], rst_id->regular_box_buffer[i]->block[j]->count[4]);
	
	memcpy(out_buf_array[cnt]->block[j]->offset, rst_id->regular_box_buffer[i]->block[j]->offset, PIDX_MAX_DIMENSIONS * sizeof(int));
	memcpy(out_buf_array[cnt]->block[j]->count, rst_id->regular_box_buffer[i]->block[j]->count, PIDX_MAX_DIMENSIONS * sizeof(int));
	out_buf_array[cnt]->block[j]->buffer = malloc(out_buf_array[cnt]->block[j]->count[0] * out_buf_array[cnt]->block[j]->count[1] * out_buf_array[cnt]->block[j]->count[2] * out_buf_array[cnt]->block[j]->count[3] * out_buf_array[cnt]->block[j]->count[4] * bytes_per_sample * samples_per_variable);
      }
      memcpy(out_buf_array[cnt]->power_two_offset, rst_id->regular_box_buffer[i]->power_two_offset, sizeof(int) * PIDX_MAX_DIMENSIONS);
      memcpy(out_buf_array[cnt]->power_two_count, rst_id->regular_box_buffer[i]->power_two_count, sizeof(int) * PIDX_MAX_DIMENSIONS);
      cnt++;
    }
  } 
  assert(cnt == num_output_buffers);
  return 0;
}

int PIDX_rst_restructure_IO(PIDX_rst_id rst_id, int samples_per_variable, MPI_Datatype datatype, Ndim_buffer* in_buf, Ndim_buffer_group* out_buf_array, int num_output_buffers)
{  
  int i, j, a1 = 0, b1 = 0, k1 = 0, i1 = 0, j1 = 0, index, count1 = 0, ret = 0, req_count = 0;
  int *send_count, *send_offset;
  int rank, nprocs, send_c = 0, send_o = 0, counter = 0, req_counter = 0, bytes_per_sample;

  MPI_Request *req;
  MPI_Status *status;

  //rank and nprocs
  ret = MPI_Comm_rank(rst_id->comm, &rank);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("Rank", __FILE__, __LINE__);

  MPI_Comm_size(rst_id->comm, &nprocs);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("nprocs", __FILE__, __LINE__);

  MPI_Type_size(datatype, &bytes_per_sample);
  
  for (i = 0; i < rst_id->owned_regular_box_count; i++)
    for(j = 0; j < rst_id->regular_box_buffer[i]->count; j++)
      req_count++;
    
  //creating ample requests and statuses
  req = (MPI_Request*) malloc(sizeof (*req) * req_count * 2);
  if (!req) PIDX_rst_print_error("Memory Error : req", __FILE__, __LINE__);

  status = (MPI_Status*) malloc(sizeof (*status) * req_count * 2);
  if (!status) PIDX_rst_print_error("Memory Error : status", __FILE__, __LINE__);

  for (i = 0; i < rst_id->owned_regular_box_count; i++)
  {
    if (rank == rst_id->regular_box_buffer[i]->max_rank)
    {
      for(j = 0; j < rst_id->regular_box_buffer[i]->count; j++)
      {
	if(rank == rst_id->regular_box_buffer[i]->rank[j])
	{
	  count1 = 0;
	  for (a1 = rst_id->regular_box_buffer[i]->block[j]->offset[4]; a1 < rst_id->regular_box_buffer[i]->block[j]->offset[4] + rst_id->regular_box_buffer[i]->block[j]->count[4]; a1++)
	    for (b1 = rst_id->regular_box_buffer[i]->block[j]->offset[3]; b1 < rst_id->regular_box_buffer[i]->block[j]->offset[3] + rst_id->regular_box_buffer[i]->block[j]->count[3]; b1++)
	      for (k1 = rst_id->regular_box_buffer[i]->block[j]->offset[2]; k1 < rst_id->regular_box_buffer[i]->block[j]->offset[2] + rst_id->regular_box_buffer[i]->block[j]->count[2]; k1++)
		for (j1 = rst_id->regular_box_buffer[i]->block[j]->offset[1]; j1 < rst_id->regular_box_buffer[i]->block[j]->offset[1] + rst_id->regular_box_buffer[i]->block[j]->count[1]; j1++)
		  for (i1 = rst_id->regular_box_buffer[i]->block[j]->offset[0]; i1 < rst_id->regular_box_buffer[i]->block[j]->offset[0] + rst_id->regular_box_buffer[i]->block[j]->count[0]; i1 = i1 + rst_id->regular_box_buffer[i]->block[j]->count[0]) 
		  {
		    index = (in_buf[0]->count[0] * in_buf[0]->count[1] * in_buf[0]->count[2] * in_buf[0]->count[3] * (a1 - in_buf[0]->offset[4])) +
			    (in_buf[0]->count[0] * in_buf[0]->count[1] * in_buf[0]->count[2] * (b1 - in_buf[0]->offset[3])) +
			    (in_buf[0]->count[0] * in_buf[0]->count[1] * (k1 - in_buf[0]->offset[2])) +
			    (in_buf[0]->count[0] * (j1 - in_buf[0]->offset[1])) +
			    (i1 - in_buf[0]->offset[0]);

		    send_o = index * samples_per_variable;
		    send_c = rst_id->regular_box_buffer[i]->block[j]->count[0] * samples_per_variable;
		    memcpy(out_buf_array[counter]->block[j]->buffer + (count1 * send_c * bytes_per_sample), in_buf[0]->buffer + send_o * bytes_per_sample, send_c * bytes_per_sample);
		    count1++;
		  }
	}
	else
	{
	  ret = MPI_Irecv(out_buf_array[counter]->block[j]->buffer, (rst_id->regular_box_buffer[i]->block[j]->count[0] * rst_id->regular_box_buffer[i]->block[j]->count[1] * rst_id->regular_box_buffer[i]->block[j]->count[2] * rst_id->regular_box_buffer[i]->block[j]->count[3] * rst_id->regular_box_buffer[i]->block[j]->count[4]) * samples_per_variable, datatype, rst_id->regular_box_buffer[i]->rank[j], 123, rst_id->comm, &req[req_counter]);
	  if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Irecv", __FILE__, __LINE__);
	  
	  req_counter++;
	}
      }
      counter++;
    }
    else
    {
      for(j = 0; j < rst_id->regular_box_buffer[i]->count; j++)
      {
	if(rank == rst_id->regular_box_buffer[i]->rank[j])
	{
	  send_offset = (int*) malloc(sizeof (int) * (rst_id->regular_box_buffer[i]->block[j]->count[1] * rst_id->regular_box_buffer[i]->block[j]->count[2] * rst_id->regular_box_buffer[i]->block[j]->count[3] * rst_id->regular_box_buffer[i]->block[j]->count[4]));
	  if (!send_offset) PIDX_rst_print_error("Memory Error : send_offset", __FILE__, __LINE__);
	  memset(send_offset, 0, sizeof (int) * (rst_id->regular_box_buffer[i]->block[j]->count[1] * rst_id->regular_box_buffer[i]->block[j]->count[2] * rst_id->regular_box_buffer[i]->block[j]->count[3] * rst_id->regular_box_buffer[i]->block[j]->count[4]));

	  send_count = (int*) malloc(sizeof (int) * (rst_id->regular_box_buffer[i]->block[j]->count[1] * rst_id->regular_box_buffer[i]->block[j]->count[2] * rst_id->regular_box_buffer[i]->block[j]->count[3] * rst_id->regular_box_buffer[i]->block[j]->count[4]));
	  if (!send_count) PIDX_rst_print_error("Memory Error : send_count", __FILE__, __LINE__);
	  memset(send_count, 0, sizeof (int) * (rst_id->regular_box_buffer[i]->block[j]->count[1] * rst_id->regular_box_buffer[i]->block[j]->count[2] * rst_id->regular_box_buffer[i]->block[j]->count[3] * rst_id->regular_box_buffer[i]->block[j]->count[4]));
	  
	  count1 = 0;
	  int tot_cnt = 0;
	  for (a1 = rst_id->regular_box_buffer[i]->block[j]->offset[4]; a1 < rst_id->regular_box_buffer[i]->block[j]->offset[4] + rst_id->regular_box_buffer[i]->block[j]->count[4]; a1++)
	    for (b1 = rst_id->regular_box_buffer[i]->block[j]->offset[3]; b1 < rst_id->regular_box_buffer[i]->block[j]->offset[3] + rst_id->regular_box_buffer[i]->block[j]->count[3]; b1++)
	      for (k1 = rst_id->regular_box_buffer[i]->block[j]->offset[2]; k1 < rst_id->regular_box_buffer[i]->block[j]->offset[2] + rst_id->regular_box_buffer[i]->block[j]->count[2]; k1++)
		for (j1 = rst_id->regular_box_buffer[i]->block[j]->offset[1]; j1 < rst_id->regular_box_buffer[i]->block[j]->offset[1] + rst_id->regular_box_buffer[i]->block[j]->count[1]; j1++)
		  for (i1 = rst_id->regular_box_buffer[i]->block[j]->offset[0]; i1 < rst_id->regular_box_buffer[i]->block[j]->offset[0] + rst_id->regular_box_buffer[i]->block[j]->count[0]; i1 = i1 + rst_id->regular_box_buffer[i]->block[j]->count[0]) 
		  {
		    index = (in_buf[0]->count[0] * in_buf[0]->count[1] * in_buf[0]->count[2] * in_buf[0]->count[3] * (a1 - in_buf[0]->offset[4])) +
			    (in_buf[0]->count[0] * in_buf[0]->count[1] * in_buf[0]->count[2] * (b1 - in_buf[0]->offset[3])) +
			    (in_buf[0]->count[0] * in_buf[0]->count[1] * (k1 - in_buf[0]->offset[2])) +
			    (in_buf[0]->count[0] * (j1 - in_buf[0]->offset[1])) +
			    (i1 - in_buf[0]->offset[0]);
		    send_offset[count1] = index * samples_per_variable;
		    send_count[count1] = rst_id->regular_box_buffer[i]->block[j]->count[0] * samples_per_variable;
		    tot_cnt = tot_cnt + send_count[count1];
		    count1++;
		  }

	  MPI_Datatype chunk_data_type;
	  MPI_Type_indexed(count1, send_count, send_offset, MPI_DOUBLE, &chunk_data_type);
	  MPI_Type_commit(&chunk_data_type);

	  
	  ret = MPI_Isend(in_buf[0]->buffer, 1, chunk_data_type, rst_id->regular_box_buffer[i]->max_rank, 123, rst_id->comm, &req[req_counter]);
	  if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Isend", __FILE__, __LINE__);
	  req_counter++;
	  
	  MPI_Type_free(&chunk_data_type);
	  free(send_offset);
	  free(send_count);
	}
      }
    }
  }

  ret = MPI_Waitall(req_counter, req, status);
  if (ret != MPI_SUCCESS) PIDX_rst_print_error("MPI_Waitall", __FILE__, __LINE__);

  free(req);
  req = 0;
  free(status);
  status = 0;
  
  return 0;
}

/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
int PIDX_rst_buf_destroy(int count, Ndim_buffer_group* out_buf_array)
{
  int i, j;
  for(i = 0; i < count; i++)
  {
    for(j = 0; j < out_buf_array[i]->count; j++)
    {
      free(out_buf_array[i]->block[j]->buffer);
      out_buf_array[i]->block[j]->buffer = 0;
      
      free(out_buf_array[i]->block[j]);
      out_buf_array[i]->block[j] = 0;
    }
    free(out_buf_array[i]->block);
    out_buf_array[i]->block = 0;
  }
  
  return 0;
}

/* tear down whatever was calculated for this particular combination of dimensions and bounds */
int PIDX_rst_finalize(PIDX_rst_id id) 
{
  free(id->idx_ptr);
  id->idx_ptr = 0;
  
  free(id->idx_derived_ptr);
  id->idx_derived_ptr = 0;
  
  free(id);
  id = 0;
  
  return 0;
}

int HELPER_rst(Ndim_buffer_group* out_buf_array1, PIDX_rst_id rst_id, int num_output_buffers, int spv)
{
  int i, j, k, rank = 0, v = 0, u = 0, s = 0, cnt = 0, m, n;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  long long element_count = 0;
  long long lost_element_count = 0;
  long long per_process_exact = 0;
  double **temp_buffer;
  int buffer_count = 0;
  
  temp_buffer = (double**) malloc(sizeof (*temp_buffer) * buffer_count);
  
  for (m = 0; m < num_output_buffers; m++)
    buffer_count = buffer_count + rst_id->regular_box_buffer[m]->count;
  
  for (m = 0; m < num_output_buffers; m++)
  {
    //printf("[%d] counts %d\n", rank, out_buf_array1[m]->count);
    for(n = 0; n < out_buf_array1[m]->count; n++)
    {
      //temp_buffer[cnt] = (double*) out_buf_array1[m]->block[n]->buffer;
      //printf("[A %d] [%d] Count: %d %d %d %d %d\n", rank, rst_id->regular_box_buffer[m]->count, out_buf_array1[m]->block[n]->count[0], out_buf_array1[m]->block[n]->count[1], out_buf_array1[m]->block[n]->count[2], out_buf_array1[m]->block[n]->count[3], out_buf_array1[m]->block[n]->count[4]);
      
      temp_buffer[cnt] = (double*) malloc(out_buf_array1[m]->block[n]->count[4] * out_buf_array1[m]->block[n]->count[3] * out_buf_array1[m]->block[n]->count[2] * out_buf_array1[m]->block[n]->count[1] * out_buf_array1[m]->block[n]->count[0] * sizeof(double));
      
      memcpy(temp_buffer[cnt], out_buf_array1[m]->block[n]->buffer, (out_buf_array1[m]->block[n]->count[4] * out_buf_array1[m]->block[n]->count[3] * out_buf_array1[m]->block[n]->count[2] * out_buf_array1[m]->block[n]->count[1] * out_buf_array1[m]->block[n]->count[0] * sizeof(double)));
      
      
      for (v = 0; v < out_buf_array1[m]->block[n]->count[4]; v++) 
	for (u = 0; u < out_buf_array1[m]->block[n]->count[3]; u++)
	  for (k = 0; k < out_buf_array1[m]->block[n]->count[2]; k++) 
	    for (j = 0; j < out_buf_array1[m]->block[n]->count[1]; j++) 
	      for (i = 0; i < out_buf_array1[m]->block[n]->count[0]; i++) 
	      {
		int index = (out_buf_array1[m]->block[n]->count[0] * out_buf_array1[m]->block[n]->count[1] * out_buf_array1[m]->block[n]->count[2] * out_buf_array1[m]->block[n]->count[3] * v) +
			(out_buf_array1[m]->block[n]->count[0] * out_buf_array1[m]->block[n]->count[1] * out_buf_array1[m]->block[n]->count[2] * u) +
			(out_buf_array1[m]->block[n]->count[0] * out_buf_array1[m]->block[n]->count[1] * k) +
			(out_buf_array1[m]->block[n]->count[0] * j) +
			i;
		int check_bit = 1;
		for (s = 0; s < spv; s++)
		    check_bit = check_bit && ((int) temp_buffer[cnt][spv * index + s] == s + 100 + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2] * rst_id->idx_ptr->global_bounds[3]*(out_buf_array1[m]->block[n]->offset[4] + v)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]*(out_buf_array1[m]->block[n]->offset[3] + u)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * (out_buf_array1[m]->block[n]->offset[2] + k)) + (rst_id->idx_ptr->global_bounds[0] * (out_buf_array1[m]->block[n]->offset[1] + j)) + out_buf_array1[m]->block[n]->offset[0] + i);

		if (check_bit == 0) 
		{
		  lost_element_count++;
		  printf("LOST Element : %f %d\n", temp_buffer[cnt][1 * index + 0], (s + 100 + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2] * rst_id->idx_ptr->global_bounds[3]*(out_buf_array1[m]->block[n]->offset[4] + v)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]*(out_buf_array1[m]->block[n]->offset[3] + u)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * (out_buf_array1[m]->block[n]->offset[2] + k)) + (rst_id->idx_ptr->global_bounds[0] * (out_buf_array1[m]->block[n]->offset[1] + j)) + out_buf_array1[m]->block[n]->offset[0] + i));
		} 
		else 
		{
		  element_count++;
		  //printf("Element : %d %d\n", (int)temp_buffer[cnt][1 * index + 0], (s + 100 + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2] * rst_id->idx_ptr->global_bounds[3]*(out_buf_array1[m]->block[n]->offset[4] + v)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]*(out_buf_array1[m]->block[n]->offset[3] + u)) + (rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * (out_buf_array1[m]->block[n]->offset[2] + k)) + (rst_id->idx_ptr->global_bounds[0] * (out_buf_array1[m]->block[n]->offset[1] + j)) + out_buf_array1[m]->block[n]->offset[0] + i));
		}
	      }
      
    }
    cnt++;
  }
  
  long long global_volume;
  MPI_Allreduce(&element_count, &global_volume, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (global_volume != (long long) rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]) 
  {
    fprintf(stderr, "[%d] RST Volume Error %lld %lld\n", rank, global_volume, (long long) rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]);
      MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  if (rank == 0)
    if (global_volume == (long long) rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2])
      printf("[%d] RST Volume %lld %lld\n", rank, global_volume, (long long) rst_id->idx_ptr->global_bounds[0] * rst_id->idx_ptr->global_bounds[1] * rst_id->idx_ptr->global_bounds[2]);
  
  return 1;
}

void PIDX_rst_print_error(char *error_message, char* file, int line) 
{
  fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);
  MPI_Abort(MPI_COMM_WORLD, -1);
}
