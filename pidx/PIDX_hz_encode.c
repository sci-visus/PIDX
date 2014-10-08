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

#include "PIDX_hz_encode.h"

struct hz_tupple_double
{
  //pointer to variables -> pointer to samples per variable -> actual data
  unsigned char*** value;
  int index;
};
typedef struct hz_tupple_double hz_tupple;

struct PIDX_hz_encode_struct 
{
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_var_index;
  int end_var_index;
  
};

int compare( const hz_tupple* a, const hz_tupple* b)
{
  long long int_a = a->index;
  long long int_b = b->index;

  if ( int_a == int_b ) return 0;
  else if ( int_a < int_b ) return -1;
  else return 1;
}

int compare_long (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}


PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{
  PIDX_hz_encode_id hz_id;
  hz_id = (PIDX_hz_encode_id)malloc(sizeof (*hz_id));
  memset(hz_id, 0, sizeof (*hz_id));

  hz_id->idx_ptr = (idx_dataset)malloc(sizeof(*(hz_id->idx_ptr)));
  memcpy(hz_id->idx_ptr, idx_meta_data, sizeof(*(hz_id->idx_ptr)));
  
  hz_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(hz_id->idx_derived_ptr)));
  memcpy(hz_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(hz_id->idx_derived_ptr)));
  
  hz_id->start_var_index = start_var_index;
  hz_id->end_var_index = end_var_index;
  
  return hz_id;
}

int PIDX_hz_encode_var(PIDX_hz_encode_id id, PIDX_variable* variable, int MODE)
{
  int** userBox;
  int i = 0, j = 0, k = 0, d = 0;
  
  userBox = (int**) malloc(2 * sizeof (int*));
  userBox[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  userBox[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  assert(userBox);

  if(variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return 1;
  }
  //printf("Start and End Var Index = %d and %d\n", id->start_var_index, id->end_var_index);
  for(i = id->start_var_index; i <= id->end_var_index; i++)
  {
    for (k = 0; k < variable[i]->patch_count; k++)
    {
      variable[i]->HZ_patch[k]->HZ_level_from = 0;
      variable[i]->HZ_patch[k]->HZ_level_to = id->idx_derived_ptr->maxh;
      variable[i]->HZ_patch[k]->allign_offset = (int**) malloc(sizeof (int*) * (id->idx_derived_ptr->maxh));
      variable[i]->HZ_patch[k]->allign_count = (int**) malloc(sizeof (int*) * (id->idx_derived_ptr->maxh));
      variable[i]->HZ_patch[k]->allign_start_hz = (long long*) malloc(sizeof (long long) * (id->idx_derived_ptr->maxh));
      variable[i]->HZ_patch[k]->allign_end_hz = (long long*) malloc(sizeof (long long) * (id->idx_derived_ptr->maxh));
      
      for (j = 0; j < id->idx_derived_ptr->maxh; j++) 
      {
	variable[i]->HZ_patch[k]->allign_offset[j] = (int*) malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
	variable[i]->HZ_patch[k]->allign_count[j] = (int*) malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
	memset(variable[i]->HZ_patch[k]->allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
	memset(variable[i]->HZ_patch[k]->allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
      }
      
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
	userBox[0][d] = variable[i]->patch[k]->offset[d];
	userBox[1][d] = variable[i]->patch[k]->offset[d] + variable[i]->patch[k]->count[d] - 1;
      }
      
      for (j = 0; j < id->idx_derived_ptr->maxh; j++) 
      {
	Align((id->idx_derived_ptr->maxh - 1), j, id->idx_ptr->bitPattern, userBox, variable[i]->HZ_patch[k]->allign_offset, variable[i]->HZ_patch[k]->allign_count);
	
	PointND startXYZ;
	startXYZ.x = variable[i]->HZ_patch[k]->allign_offset[j][0];
	startXYZ.y = variable[i]->HZ_patch[k]->allign_offset[j][1];
	startXYZ.z = variable[i]->HZ_patch[k]->allign_offset[j][2];
	startXYZ.u = variable[i]->HZ_patch[k]->allign_offset[j][3];
	startXYZ.v = variable[i]->HZ_patch[k]->allign_offset[j][4];
	variable[i]->HZ_patch[k]->allign_start_hz[j] = xyz_to_HZ(id->idx_ptr->bitPattern, id->idx_derived_ptr->maxh - 1, startXYZ);

	PointND endXYZ;
	endXYZ.x = variable[i]->HZ_patch[k]->allign_count[j][0];
	endXYZ.y = variable[i]->HZ_patch[k]->allign_count[j][1];
	endXYZ.z = variable[i]->HZ_patch[k]->allign_count[j][2];
	endXYZ.u = variable[i]->HZ_patch[k]->allign_count[j][3];
	endXYZ.v = variable[i]->HZ_patch[k]->allign_count[j][4];
	variable[i]->HZ_patch[k]->allign_end_hz[j] = xyz_to_HZ(id->idx_ptr->bitPattern, id->idx_derived_ptr->maxh - 1, endXYZ); 
      }
    }
  }
  free(userBox[0]);
  userBox[0] = 0;
  free(userBox[1]);
  userBox[1] = 0;
  free(userBox);
  userBox = 0;
  
  for (k = 0; k < variable[id->start_var_index]->patch_count; k++) 
  {
    variable[id->start_var_index]->HZ_patch[k]->samples_per_level = (int*) malloc( id->idx_derived_ptr->maxh * sizeof (int));
    memset(variable[id->start_var_index]->HZ_patch[k]->samples_per_level, 0, id->idx_derived_ptr->maxh * sizeof (int));
  }
  
  for(i = id->start_var_index; i <= id->end_var_index; i++)
  {
    for (k = 0; k < variable[id->start_var_index]->patch_count; k++)
    {
      variable[i]->HZ_patch[k]->buffer = (unsigned char**)malloc( id->idx_derived_ptr->maxh * sizeof (unsigned char*));
      memset(variable[i]->HZ_patch[k]->buffer, 0,  id->idx_derived_ptr->maxh * sizeof (unsigned char*));
    }
    
    for (k = 0; k < variable[id->start_var_index]->patch_count; k++) 
    {
      variable[i]->HZ_patch[k]->buffer_index = (long long*) malloc(((variable[i]->patch[k]->count[0]) * (variable[i]->patch[k]->count[1]) * (variable[i]->patch[k]->count[2])) * sizeof (long long));
      
      memset(variable[i]->HZ_patch[k]->buffer_index, 0, ((variable[i]->patch[k]->count[0]) * (variable[i]->patch[k]->count[1]) * (variable[i]->patch[k]->count[2])) * sizeof (long long));
    }
  }
  
  /*
  if(MODE == PIDX_READ)
  {
    long long z_order = 0, hz_order = 0, index = 0;
    long long l_x, l_y, l_z, l_u, l_v;
    int i = 0, j = 0, k = 0, d = 0, level = 0, u = 0, v = 0, cnt = 0, y = 0, number_levels = 0;
    int index_count = 0;
  
    id->index = malloc(variable[id->start_var_index]->patch_count * sizeof(long long*));
    for (y = 0; y < variable[id->start_var_index]->patch_count; y++)
    { 
      id->index[y] = malloc((variable[id->start_var_index]->patch[y]->count[0] * variable[id->start_var_index]->patch[y]->count[1] * variable[id->start_var_index]->patch[y]->count[2] * variable[id->start_var_index]->patch[y]->count[3] * variable[id->start_var_index]->patch[y]->count[4]) * sizeof (long long));
      
      index_count = 0;
      
      l_x = variable[id->start_var_index]->patch[y]->count[0];
      l_y = variable[id->start_var_index]->patch[y]->count[1];
      l_z = variable[id->start_var_index]->patch[y]->count[2];
      l_u = variable[id->start_var_index]->patch[y]->count[3];
      l_v = variable[id->start_var_index]->patch[y]->count[4];
      
      number_levels = id->idx_derived_ptr->maxh - 1;

      PointND xyzuv_Index;
      for (v = variable[id->start_var_index]->patch[y]->offset[4]; v < variable[id->start_var_index]->patch[y]->offset[4] + variable[id->start_var_index]->patch[y]->count[4]; v++)
      {
	for (u = variable[id->start_var_index]->patch[y]->offset[3]; u < variable[id->start_var_index]->patch[y]->offset[3] + variable[id->start_var_index]->patch[y]->count[3]; u++)
	{
	  for (k = variable[id->start_var_index]->patch[y]->offset[2]; k < variable[id->start_var_index]->patch[y]->offset[2] + variable[id->start_var_index]->patch[y]->count[2]; k++)
	  {
	    for (j = variable[id->start_var_index]->patch[y]->offset[1]; j < variable[id->start_var_index]->patch[y]->offset[1] + variable[id->start_var_index]->patch[y]->count[1]; j++)
	    {
	      for (i = variable[id->start_var_index]->patch[y]->offset[0]; i < variable[id->start_var_index]->patch[y]->offset[0] + variable[id->start_var_index]->patch[y]->count[0]; i++) 
	      {
		xyzuv_Index.x = i;
		xyzuv_Index.y = j;
		xyzuv_Index.z = k;
		xyzuv_Index.u = u;
		xyzuv_Index.v = v;

		z_order = 0;
		PointND zero;
		zero.x = 0;
		zero.y = 0;
		zero.z = 0;
		zero.u = 0;
		zero.v = 0;
		memset(&zero, 0, sizeof (PointND));

		for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--) 
		{
		  int bit = id->idx_ptr->bitPattern[number_levels];
		  z_order |= ((long long) PGET(xyzuv_Index, bit) & 1) << cnt;
		  PGET(xyzuv_Index, bit) >>= 1;
		}

		number_levels = id->idx_derived_ptr->maxh - 1;
		long long lastbitmask = ((long long) 1) << number_levels;
		z_order |= lastbitmask;
		while (!(1 & z_order)) z_order >>= 1;
		z_order >>= 1;

		hz_order = z_order;
		level = getLeveL(hz_order);
		
		index = (l_x * l_y * l_z * l_u * (v - variable[id->start_var_index]->patch[y]->offset[4])) + (l_x * l_y * l_z * (u - variable[id->start_var_index]->patch[y]->offset[3])) + (l_x * l_y * (k - variable[id->start_var_index]->patch[y]->offset[2])) + (l_x * (j - variable[id->start_var_index]->patch[y]->offset[1])) + (i - variable[id->start_var_index]->patch[y]->offset[0]);
		
		variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
		
		id->index[y][index_count] = index;
		variable[id->start_var_index]->HZ_patch[y]->buffer_index[index_count] = hz_order;
		
		printf("[%d][%d : %d] : IJK %lld HZ %lld\n", index_count, level, variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level], index, hz_order);
		
		index_count++;
	      }
	    }
	  }
	}
      }
      qsort( variable[id->start_var_index]->HZ_patch[y]->buffer_index, ( variable[id->start_var_index]->patch[y]->count[0] * variable[id->start_var_index]->patch[y]->count[1] * variable[id->start_var_index]->patch[y]->count[2] * variable[id->start_var_index]->patch[y]->count[3] * variable[id->start_var_index]->patch[y]->count[4] ), sizeof(long long), compare_long );
    }
  }
  */ 
  return 0;
}

int PIDX_hz_encode_write_var(PIDX_hz_encode_id id, PIDX_variable* variable)
{
  long long z_order = 0, hz_order = 0, index = 0;
  long long l_x, l_y, l_z, l_u, l_v;
  int i = 0, j = 0, k = 0, level = 0, u = 0, v = 0, cnt = 0, c = 0, s = 0, y = 0, number_levels = 0, var = 0;
  int index_count = 0;
  int bytes_for_datatype;
  
  if(variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return 1;
  }
  
  for (y = 0; y < variable[id->start_var_index]->patch_count; y++)
  {
    hz_tupple *tupple = malloc((variable[id->start_var_index]->patch[y]->count[0] * variable[id->start_var_index]->patch[y]->count[1] * variable[id->start_var_index]->patch[y]->count[2] * variable[id->start_var_index]->patch[y]->count[3] * variable[id->start_var_index]->patch[y]->count[4]) /*values_per_sample*/ * sizeof (hz_tupple));
    
    index_count = 0;
    l_x = variable[id->start_var_index]->patch[y]->count[0];
    l_y = variable[id->start_var_index]->patch[y]->count[1];
    l_z = variable[id->start_var_index]->patch[y]->count[2];
    l_u = variable[id->start_var_index]->patch[y]->count[3];
    l_v = variable[id->start_var_index]->patch[y]->count[4];
	
    //printf("Patch %d [%d]: %lld %lld %lld %lld %lld :: %d %d %d %d %d\n", y, id->idx_derived_ptr->patch_count, l_x, l_y, l_z, l_u, l_v, variable[id->start_var_index]->patch[y]->offset[0], variable[id->start_var_index]->patch[y]->offset[1], variable[id->start_var_index]->patch[y]->offset[2], variable[id->start_var_index]->patch[y]->offset[3], variable[id->start_var_index]->patch[y]->offset[4]);
    
    number_levels = id->idx_derived_ptr->maxh - 1;

    PointND xyzuv_Index;
    for (v = variable[id->start_var_index]->patch[y]->offset[4]; v < variable[id->start_var_index]->patch[y]->offset[4] + variable[id->start_var_index]->patch[y]->count[4]; v++)
    {
      for (u = variable[id->start_var_index]->patch[y]->offset[3]; u < variable[id->start_var_index]->patch[y]->offset[3] + variable[id->start_var_index]->patch[y]->count[3]; u++)
      {
	for (k = variable[id->start_var_index]->patch[y]->offset[2]; k < variable[id->start_var_index]->patch[y]->offset[2] + variable[id->start_var_index]->patch[y]->count[2]; k++)
	{
	  for (j = variable[id->start_var_index]->patch[y]->offset[1]; j < variable[id->start_var_index]->patch[y]->offset[1] + variable[id->start_var_index]->patch[y]->count[1]; j++)
	  {
	    for (i = variable[id->start_var_index]->patch[y]->offset[0]; i < variable[id->start_var_index]->patch[y]->offset[0] + variable[id->start_var_index]->patch[y]->count[0]; i++) 
	    {   
	      xyzuv_Index.x = i;
	      xyzuv_Index.y = j;
	      xyzuv_Index.z = k;
	      xyzuv_Index.u = u;
	      xyzuv_Index.v = v;

	      z_order = 0;
	      PointND zero;
	      zero.x = 0;
	      zero.y = 0;
	      zero.z = 0;
	      zero.u = 0;
	      zero.v = 0;
	      memset(&zero, 0, sizeof (PointND));

	      for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--) 
	      {
		int bit = id->idx_ptr->bitPattern[number_levels];
		z_order |= ((long long) PGET(xyzuv_Index, bit) & 1) << cnt;
		PGET(xyzuv_Index, bit) >>= 1;
	      }

	      number_levels = id->idx_derived_ptr->maxh - 1;
	      long long lastbitmask = ((long long) 1) << number_levels;
	      z_order |= lastbitmask;
	      while (!(1 & z_order)) z_order >>= 1;
	      z_order >>= 1;

	      hz_order = z_order;
	      level = getLeveL(hz_order);
	      
	      index = (l_x * l_y * l_z * l_u * (v - variable[id->start_var_index]->patch[y]->offset[4])) + (l_x * l_y * l_z * (u - variable[id->start_var_index]->patch[y]->offset[3])) + (l_x * l_y * (k - variable[id->start_var_index]->patch[y]->offset[2])) + (l_x * (j - variable[id->start_var_index]->patch[y]->offset[1])) + (i - variable[id->start_var_index]->patch[y]->offset[0]);
	      
	      tupple[index_count].value = malloc((id->end_var_index - id->start_var_index + 1)* sizeof(unsigned char**));
	      variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
	      
	      for(var = id->start_var_index; var <= id->end_var_index; var++)
	      {
		tupple[index_count].value[var - id->start_var_index] = malloc(variable[var]->values_per_sample * sizeof(unsigned char*));
		bytes_for_datatype = variable[var]->bits_per_value / 8;
		
		for (s = 0; s < variable[var]->values_per_sample; s++) 
		{
		  tupple[index_count].value[var - id->start_var_index][s] = malloc(bytes_for_datatype);
		  memcpy(tupple[index_count].value[var - id->start_var_index][s], variable[var]->patch[y]->buffer + ((index * variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
		  
		  
		  //if(bytes_for_datatype == sizeof(double))
		  //{
		    //double dvalue;
		    //memcpy(&value, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
		    //memcpy(&value, variable[var]->patch[y]->buffer + ((index * variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
		    //printf("HZD: [%d %d %d %d] %d %f\n", y, index_count, var, s, bytes_for_datatype, value);
		  //}
		  //else
		  //{
		    //int ivalue;
		    //memcpy(&ivalue, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
		    //memcpy(&ivalue, variable[var]->patch[y]->buffer + ((index * variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
		    //printf("HZI: [%d %d %d %d] %d %d\n", y, index_count, var, s, bytes_for_datatype, ivalue);
		  //}
		  
		  tupple[index_count].index = hz_order;
		}
	      }
	      index_count++;
	    }
	  }
	}
      }
    }
    qsort( tupple, (l_x * l_y * l_z * l_u * l_v), sizeof(hz_tupple), compare );
    
    for(var = id->start_var_index; var <= id->end_var_index; var++)
      for(c = 0 ; c < id->idx_derived_ptr->maxh ; c++)
      {
	bytes_for_datatype = variable[var]->bits_per_value / 8;
	variable[var]->HZ_patch[y]->buffer[c] = malloc(bytes_for_datatype * variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c] * variable[var]->values_per_sample);
      }
    
    cnt = 0;
    
    for(c = 0; c < id->idx_derived_ptr->maxh; c++)
    {
      for(s = 0; s < variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c]; s++)
      {
	for(var = id->start_var_index; var <= id->end_var_index; var++)
	{
	  for (i = 0; i < variable[var]->values_per_sample; i++) 
	  {
	    bytes_for_datatype = variable[var]->bits_per_value / 8;
	    memcpy(variable[var]->HZ_patch[y]->buffer[c] + ((s * variable[var]->values_per_sample + i) * bytes_for_datatype), tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype);
	    variable[var]->HZ_patch[y]->buffer_index[cnt] = tupple[cnt].index;
	    
	    
	    //if(variable[var]->bits_per_value / 8 == sizeof(double))
	    //{
	    //  double check;
	    //  memcpy(&check, tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype);
	      //memcpy(&check, variable[var]->HZ_patch[y]->buffer[c] + ((s * variable[var]->values_per_sample + i) * bytes_for_datatype), bytes_for_datatype);
	      //printf("[HD%d] %d %d %d %d %d %f\n", var, cnt, y, c, (s * variable[var]->values_per_sample + i), bytes_for_datatype, check );
	    //}
	    //else
	    //{
	    //  int icheck;
	    //  memcpy(&icheck, tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype);
	      //memcpy(&icheck, variable[var]->HZ_patch[y]->buffer[c] + ((s * variable[var]->values_per_sample + i) * bytes_for_datatype), bytes_for_datatype);
	      //printf("[HI%d] %d %d %d %d %d %d\n", var, cnt, y, c, (s * variable[var]->values_per_sample + i), bytes_for_datatype, icheck );
	    //}
	    
	    free(tupple[cnt].value[var - id->start_var_index][i]);
	  }
	  free(tupple[cnt].value[var - id->start_var_index]);
	}
	free(tupple[cnt].value);
	cnt++;
      }
    }
    
    free(tupple);
    tupple = 0;
  }
  return 0;
}

int PIDX_hz_encode_read_var(PIDX_hz_encode_id id, PIDX_variable* variable)
{
  /*
  int i = 0, cnt = 0, c = 0, s = 0, y = 0, var = 0;
  int bytes_for_datatype;
  
  if(variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return 1;
  }
  
  for (y = 0; y < variable[id->start_var_index]->patch_count; y++)
  { 
    hz_tupple *tupple = malloc((variable[id->start_var_index]->patch[y]->count[0] * variable[id->start_var_index]->patch[y]->count[1] * variable[id->start_var_index]->patch[y]->count[2] * variable[id->start_var_index]->patch[y]->count[3] * variable[id->start_var_index]->patch[y]->count[4])  * sizeof (hz_tupple));
   
    for(var = id->start_var_index; var <= id->end_var_index; var++)
    {
      for(c = 0 ; c < id->idx_derived_ptr->maxh ; c++)
      {
	bytes_for_datatype = variable[var]->bits_per_value / 8;
	variable[var]->HZ_patch[y]->buffer[c] = malloc(bytes_for_datatype * variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c] * variable[var]->values_per_sample);
      }
    }
    
    cnt = 0;
    for(c = 0; c < id->idx_derived_ptr->maxh; c++)
    {
      for(s = 0; s < variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c]; s++)
      {
	tupple[cnt].value = malloc(id->end_var_index - id->start_var_index + 1);
	for(var = id->start_var_index; var <= id->end_var_index; var++)
	{
	  tupple[cnt].value[var] = malloc(variable[var]->values_per_sample);
	  for (i = 0; i < variable[var]->values_per_sample; i++) 
	  {
	    bytes_for_datatype = variable[var]->bits_per_value / 8;
	    tupple[cnt].value[var][i] = malloc(bytes_for_datatype);
	    memcpy(tupple[cnt].value[var][i], variable[var]->HZ_patch[y]->buffer[c] + ((s * variable[var]->values_per_sample + i) * bytes_for_datatype), bytes_for_datatype);
	    tupple[cnt].index = variable[var]->HZ_patch[y]->buffer_index[cnt];
	  }
	}
	cnt++;
      }
    }
    
    for (v = variable[id->start_var_index]->patch[y]->offset[4]; v < variable[id->start_var_index]->patch[y]->offset[4] + variable[id->start_var_index]->patch[y]->count[4]; v++)
    
      for (u = variable[id->start_var_index]->patch[y]->offset[3]; u < variable[id->start_var_index]->patch[y]->offset[3] + variable[id->start_var_index]->patch[y]->count[3]; u++)
      
	for (k = variable[id->start_var_index]->patch[y]->offset[2]; k < variable[id->start_var_index]->patch[y]->offset[2] + variable[id->start_var_index]->patch[y]->count[2]; k++)
	
	  for (j = variable[id->start_var_index]->patch[y]->offset[1]; j < variable[id->start_var_index]->patch[y]->offset[1] + variable[id->start_var_index]->patch[y]->count[1]; j++)
	  
	    for (i = variable[id->start_var_index]->patch[y]->offset[0]; i < variable[id->start_var_index]->patch[y]->offset[0] + variable[id->start_var_index]->patch[y]->count[0]; i++) 
	    
	      for(var = id->start_var_index; var <= id->end_var_index; var++)
	      {
		
		index = (l_x * l_y * l_z * l_u * (v - variable[id->start_var_index]->patch[y]->offset[4])) + (l_x * l_y * l_z * (u - variable[id->start_var_index]->patch[y]->offset[3])) + (l_x * l_y * (k - variable[id->start_var_index]->patch[y]->offset[2])) + (l_x * (j - variable[id->start_var_index]->patch[y]->offset[1])) + (i - variable[id->start_var_index]->patch[y]->offset[0]);
		
		printf("HZZZZ %d\n", tupple[i].index);
		bytes_for_datatype = variable[var]->bits_per_value / 8;
		for (s = 0; s < variable[var]->values_per_sample; s++) 
		{
		  memcpy(variable[var]->patch[y]->buffer + ((index * variable[var]->values_per_sample) + s) * bytes_for_datatype, tupple[i].value[var][s], bytes_for_datatype);
		  free(tupple[i].value[var][s]);
		}
		free(tupple[i].value[var]);
	      }
	      free(tupple[i].value);
	    }
    
  free(tupple);
  tupple = 0;
	    
  free(id->index[y]);
  id->index[y] = 0;
  
  return 0;
  */
}

/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
int PIDX_hz_encode_buf_destroy_var(PIDX_hz_encode_id id, PIDX_variable* variable) 
{
  int itr = 0, p = 0, var = 0;
  
  if(variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh (%d) not set.\n", __FILE__, __LINE__, id->idx_derived_ptr->maxh);
    return 1;
  }
  
  for (p = 0; p < variable[id->start_var_index]->patch_count; p++)
    free(variable[id->start_var_index]->HZ_patch[p]->samples_per_level);
  
  for(var = id->start_var_index; var <= id->end_var_index; var++)
  {
    for (p = 0; p < variable[id->start_var_index]->patch_count; p++)
    {
      free(variable[var]->HZ_patch[p]->allign_start_hz);
      free(variable[var]->HZ_patch[p]->allign_end_hz);
      free(variable[var]->HZ_patch[p]->buffer_index);
      
      for(itr = 0 ; itr < id->idx_derived_ptr->maxh ; itr++)
      {
	free(variable[var]->HZ_patch[p]->buffer[itr]);
	variable[var]->HZ_patch[p]->buffer[itr] = 0;
      }
      free(variable[var]->HZ_patch[p]->buffer);
      variable[var]->HZ_patch[p]->buffer = 0;
    
      for (itr = 0; itr < variable[var]->HZ_patch[p]->HZ_level_to - variable[var]->HZ_patch[p]->HZ_level_from; itr++) 
      {
	free(variable[var]->HZ_patch[p]->allign_offset[itr]);
	variable[var]->HZ_patch[p]->allign_offset[itr] = 0;
	free(variable[var]->HZ_patch[p]->allign_count[itr]);
	variable[var]->HZ_patch[p]->allign_count[itr] = 0;
      }
      free(variable[var]->HZ_patch[p]->allign_offset);
      variable[var]->HZ_patch[p]->allign_offset = 0;
      free(variable[var]->HZ_patch[p]->allign_count);
      variable[var]->HZ_patch[p]->allign_count = 0;
    
      free(variable[var]->HZ_patch[p]);
      variable[var]->HZ_patch[p] = 0;
    }
  }
  
  return 0;
}

int PIDX_hz_encode_finalize(PIDX_hz_encode_id id) 
{  
  if(id == NULL)
  {
    fprintf(stderr, "[%s] [%d] hz id is null.\n", __FILE__, __LINE__);
    return 1;
  }
  
  free(id->idx_ptr);
  id->idx_ptr = 0;
  
  free(id->idx_derived_ptr);
  id->idx_derived_ptr = 0;
  
  free(id);
  id = 0;
  
  return 0;
}

#if PIDX_HAVE_MPI
int HELPER_Hz_encode(PIDX_hz_encode_id hz_id, HZ_buffer* out_hz_array, MPI_Datatype datatype, int values_per_sample)
{
  int i = 0, k = 0, b = 0, rank;
  long long global_hz, element_counts = 0, lost_element_count = 0;
  long long* ZYX;
  int check_bit = 1, s = 0;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Patch Count: %d\n", hz_id->idx_ptr->variable[hz_id->start_var_index]->patch_count);
  for (b = 0; b < hz_id->idx_ptr->variable[hz_id->start_var_index]->patch_count; b++)
  {
    printf("maxh = %d\n", (out_hz_array[b]->HZ_level_to - out_hz_array[b]->HZ_level_from));
    for (i = 0; i < out_hz_array[b]->HZ_level_to - out_hz_array[b]->HZ_level_from; i++) 
    {
      printf("samples at %d = %lld %lld\n", i, out_hz_array[b]->allign_start_hz[i], out_hz_array[b]->allign_end_hz[i]);
      for (k = 0; k <= (out_hz_array[b]->allign_end_hz[i] - out_hz_array[b]->allign_start_hz[i]) * 1; k++) 
      {
	global_hz = out_hz_array[b]->allign_start_hz[i] + k;
	
	ZYX = malloc(PIDX_MAX_DIMENSIONS * sizeof(long long));
	Hz_to_xyz(hz_id->idx_ptr->bitPattern, hz_id->idx_derived_ptr->maxh - 1, global_hz, ZYX);

	if ((ZYX[0] < hz_id->idx_ptr->global_bounds[0] && ZYX[1] < hz_id->idx_ptr->global_bounds[1] && ZYX[2] < hz_id->idx_ptr->global_bounds[2] && ZYX[3] < hz_id->idx_ptr->global_bounds[3] && ZYX[4] < hz_id->idx_ptr->global_bounds[4])) 
	{
	  check_bit = 1, s = 0;
	  
	  for (s = 0; s < values_per_sample; s++)
	  {
	    if(datatype == MPI_DOUBLE)
	    {
	      check_bit = check_bit && (*(*((double**)out_hz_array[b]->buffer + i) + ((k * values_per_sample) + s))  == (double)s + 100 + (hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1]*(ZYX[2]))+(hz_id->idx_ptr->global_bounds[0]*(ZYX[1])) + ZYX[0]);

	      if (check_bit == 0)
	      {
		printf("Elements: %f\n", *(*((double**)out_hz_array[b]->buffer + i) + ((k * values_per_sample) + s)));
		lost_element_count++;
	      }
	      else
		element_counts++;
	    }
	    else
	    {
	      check_bit = check_bit && (*(*((int**)out_hz_array[b]->buffer + i) + ((k * values_per_sample) + s))  == (int)s + 100 + (hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1]*(ZYX[2]))+(hz_id->idx_ptr->global_bounds[0]*(ZYX[1])) + ZYX[0]);

	      if (check_bit == 0)
	      {
		//printf("Elements:; %f %f\n", out_hz_array[b]->buffer[i][(k * 1) + 0]);
		lost_element_count++;
	      }
	      else
		element_counts++;
	    }
	  }
	}
	free(ZYX);
	ZYX = 0;
      }    
    }
  }

  long long global_volume = 0;
  MPI_Allreduce(&element_counts, &global_volume, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (global_volume != (long long) hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1] * hz_id->idx_ptr->global_bounds[2]) {
    //if(rank == 0)
    fprintf(stderr, "Volume Error %lld %lld\n", global_volume, (long long) hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1] * hz_id->idx_ptr->global_bounds[2]);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  assert(global_volume == (long long) hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1] * hz_id->idx_ptr->global_bounds[2]);
  if (rank == 0)
    if (global_volume == (long long) hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1] * hz_id->idx_ptr->global_bounds[2])
      printf("HZ : Volume %lld %lld\n", global_volume, (long long) hz_id->idx_ptr->global_bounds[0] * hz_id->idx_ptr->global_bounds[1] * hz_id->idx_ptr->global_bounds[2]);
    
  return 0;
}
#endif
