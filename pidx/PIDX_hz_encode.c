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
  
  int** index;
  
  int start_var_index;
  int end_var_index;
  
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
};


#if PIDX_HAVE_MPI
int PIDX_hz_encode_set_communicator(PIDX_hz_encode_id id, MPI_Comm comm)
{
  id->comm = comm;
  return 0;
}
#endif

int compare( const void* a, const void* b)
{
  int64_t int_a = ((const hz_tupple*)a)->index;
  int64_t int_b = ((const hz_tupple*)b)->index;

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

  hz_id->idx_ptr = idx_meta_data;
  hz_id->idx_derived_ptr = idx_derived_ptr;
  hz_id->start_var_index = start_var_index;
  hz_id->end_var_index = end_var_index;
  
  return hz_id;
}

int PIDX_hz_encode_buf_create(PIDX_hz_encode_id id)
{
  int** userBox;
  int **allign_offset;
  int **allign_count;
  int start_block_no, end_block_no, b;
  int i = 0, j = 0, k = 0, d = 0, c = 0, bytes_for_datatype = 0, count = 0;

  userBox = (int**) malloc(2 * sizeof (int*));
  userBox[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  userBox[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  assert(userBox);

  if(id->idx_ptr->variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return 1;
  }

  for(i = id->start_var_index; i <= id->end_var_index; i++)
  {
    for (k = 0; k < id->idx_ptr->variable[i]->patch_group_count; k++)
    {
      id->idx_ptr->variable[i]->HZ_patch[k]->HZ_level_from = 0;
      id->idx_ptr->variable[i]->HZ_patch[k]->HZ_level_to = id->idx_derived_ptr->maxh;

      id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index = malloc(sizeof (int64_t) * (id->idx_derived_ptr->maxh));
      id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index = malloc(sizeof (int64_t) * (id->idx_derived_ptr->maxh));
      memset(id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index, 0, sizeof (int64_t) * (id->idx_derived_ptr->maxh));
      memset(id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index, 0, sizeof (int64_t) * (id->idx_derived_ptr->maxh));

      id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level = malloc(sizeof (int) * (id->idx_derived_ptr->maxh));
      id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level = malloc(sizeof (int*) * (id->idx_derived_ptr->maxh));
      memset(id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level, 0, sizeof (int*) * (id->idx_derived_ptr->maxh));
      memset(id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level, 0, sizeof (int) * (id->idx_derived_ptr->maxh));

      allign_offset = malloc(sizeof (int*) * (id->idx_derived_ptr->maxh));
      allign_count = malloc(sizeof (int*) * (id->idx_derived_ptr->maxh));
      memset(allign_offset, 0, sizeof (int*) * id->idx_derived_ptr->maxh);
      memset(allign_count, 0, sizeof (int*) * id->idx_derived_ptr->maxh);

      for (j = 0; j < id->idx_derived_ptr->maxh; j++)
      {
    allign_offset[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    allign_count[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    memset(allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
    memset(allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

        id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j] = malloc(sizeof (int) * (id->idx_ptr->blocks_per_file));
        memset(id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j], 0, (sizeof (int) * (id->idx_ptr->blocks_per_file)));
      }

      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
    userBox[0][d] = (int) id->idx_ptr->variable[i]->patch_group_ptr[k]->enclosing_box_offset[d] / id->idx_ptr->compression_block_size[d];
    userBox[1][d] = (int) (id->idx_ptr->variable[i]->patch_group_ptr[k]->enclosing_box_offset[d] / id->idx_ptr->compression_block_size[d]) + (id->idx_ptr->variable[i]->patch_group_ptr[k]->enclosing_box_size[d] / id->idx_ptr->compression_block_size[d]) - 1;
      }

      for (j = 0; j < id->idx_derived_ptr->maxh; j++)
      {
        Align((id->idx_derived_ptr->maxh - 1), j, id->idx_ptr->bitPattern, userBox, allign_offset, allign_count);

        PointND startXYZ;
        startXYZ.x = allign_offset[j][0];
        startXYZ.y = allign_offset[j][1];
        startXYZ.z = allign_offset[j][2];
        startXYZ.u = allign_offset[j][3];
        startXYZ.v = allign_offset[j][4];
        id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[j] = xyz_to_HZ(id->idx_ptr->bitPattern, id->idx_derived_ptr->maxh - 1, startXYZ);

        PointND endXYZ;
        endXYZ.x = allign_count[j][0];
        endXYZ.y = allign_count[j][1];
        endXYZ.z = allign_count[j][2];
        endXYZ.u = allign_count[j][3];
        endXYZ.v = allign_count[j][4];
        id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[j] = xyz_to_HZ(id->idx_ptr->bitPattern, id->idx_derived_ptr->maxh - 1, endXYZ);

        if (id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[k]->box_group_type == 2)
        {
          start_block_no = id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[j] / id->idx_derived_ptr->samples_per_block;
          end_block_no = id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[j] / id->idx_derived_ptr->samples_per_block;

          count = 0;
          for (b = start_block_no; b <= end_block_no; b++)
          {
#ifdef PIDX_VAR_SLOW_LOOP
            if (PIDX_blocks_is_block_present(b, id->idx_ptr->id->idx_ptr->variable[id->start_var_index]->VAR_global_block_layout) == 0)
            {
              //printf("[%d] missing blocks == %d\n", j, b);
              id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level[j]++;
              id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j][count] = b;
              count++;
            }
#else
            if (PIDX_blocks_is_block_present(b, id->idx_derived_ptr->global_block_layout) == 0)
            {
              id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level[j]++;
              id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j][count] = b;
              //printf("XXXXXXXXX [%d] [%d] (%d %d) (%d %d) (%lld %lld) missing blocks == %d\n", j, count,  id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level[j], id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j][count], start_block_no, end_block_no, (long long)id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[j], (long long)id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[j], b);
              count++;
            }
#endif
          }
          //if(j == 19)
          //{
            //id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level[j]++;
            //id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j][count] = 15;
            //printf("[%d] [%d] (%d %d) (%d %d) (%lld %lld) missing blocks == %d\n", j, count, id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_count_per_level[j], id->idx_ptr->variable[i]->HZ_patch[k]->missing_block_index_per_level[j][count], start_block_no, end_block_no, id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[j], id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[j], b);
          //}
        }

        free(allign_offset[j]);
        free(allign_count[j]);
      }
      free(allign_offset);
      free(allign_count);
    }
  }
  free(userBox[0]);
  userBox[0] = 0;
  free(userBox[1]);
  userBox[1] = 0;
  free(userBox);
  userBox = 0;

  for (k = 0; k < id->idx_ptr->variable[id->start_var_index]->patch_group_count; k++)
  {
    id->idx_ptr->variable[id->start_var_index]->HZ_patch[k]->samples_per_level = malloc( id->idx_derived_ptr->maxh * sizeof (int64_t));
    memset(id->idx_ptr->variable[id->start_var_index]->HZ_patch[k]->samples_per_level, 0, id->idx_derived_ptr->maxh * sizeof (int64_t));
  }

  for (i = id->start_var_index; i <= id->end_var_index; i++)
  {
    //printf("p g count %d\n", id->idx_ptr->variable[id->start_var_index]->patch_group_count);
    for (k = 0; k < id->idx_ptr->variable[id->start_var_index]->patch_group_count; k++)
    {
      id->idx_ptr->variable[i]->HZ_patch[k]->buffer = (unsigned char**)malloc( id->idx_derived_ptr->maxh * sizeof (unsigned char*));
      memset(id->idx_ptr->variable[i]->HZ_patch[k]->buffer, 0,  id->idx_derived_ptr->maxh * sizeof (unsigned char*));
      //printf("p g type %d = %d\n", k, id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[k]->box_group_type);
      if(id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[k]->box_group_type == 1 || id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[k]->box_group_type == 2)
      {
        for (c = 0 ; c < id->idx_derived_ptr->maxh ; c++)
        {
          //printf("buffer at %d = %d\n", c, bytes_for_datatype * (id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[c] - id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[c] + 1) * id->idx_ptr->variable[i]->values_per_sample * id->idx_ptr->compression_block_size[0] * id->idx_ptr->compression_block_size[1] * id->idx_ptr->compression_block_size[2] * id->idx_ptr->compression_block_size[3] * id->idx_ptr->compression_block_size[4]);
          bytes_for_datatype = id->idx_ptr->variable[i]->bits_per_value / 8;
          id->idx_ptr->variable[i]->HZ_patch[k]->buffer[c] = malloc(bytes_for_datatype * (id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[c] - id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[c] + 1) * id->idx_ptr->variable[i]->values_per_sample * id->idx_ptr->compression_block_size[0] * id->idx_ptr->compression_block_size[1] * id->idx_ptr->compression_block_size[2] * id->idx_ptr->compression_block_size[3] * id->idx_ptr->compression_block_size[4]);
          memset(id->idx_ptr->variable[i]->HZ_patch[k]->buffer[c], 0, bytes_for_datatype * (id->idx_ptr->variable[i]->HZ_patch[k]->end_hz_index[c] - id->idx_ptr->variable[i]->HZ_patch[k]->start_hz_index[c] + 1) * id->idx_ptr->variable[i]->values_per_sample * id->idx_ptr->compression_block_size[0] * id->idx_ptr->compression_block_size[1] * id->idx_ptr->compression_block_size[2] * id->idx_ptr->compression_block_size[3] * id->idx_ptr->compression_block_size[4]);
        }
      }
    }

    for (k = 0; k < id->idx_ptr->variable[id->start_var_index]->patch_group_count; k++)
    {
      if(id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[k]->box_group_type == 0)
      {
    id->idx_ptr->variable[i]->HZ_patch[k]->buffer_index = (int64_t*) malloc(((id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[0]/id->idx_ptr->compression_block_size[0]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[1]/id->idx_ptr->compression_block_size[1]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[2]/id->idx_ptr->compression_block_size[2]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[3]/id->idx_ptr->compression_block_size[3]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[4]/id->idx_ptr->compression_block_size[4])) * sizeof (int64_t));

    memset(id->idx_ptr->variable[i]->HZ_patch[k]->buffer_index, 0, ((id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[0]/id->idx_ptr->compression_block_size[0]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[1]/id->idx_ptr->compression_block_size[1]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[2]/id->idx_ptr->compression_block_size[2]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[3]/id->idx_ptr->compression_block_size[3]) * (id->idx_ptr->variable[i]->patch_group_ptr[k]->box[0]->Ndim_box_size[4]/id->idx_ptr->compression_block_size[4])) * sizeof (int64_t));
      }
    }
  }

  return 0;
}

int PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  int64_t z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, c = 0, s = 0, y = 0, n = 0, m = 0, number_levels = 0, var = 0;
  int64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0;
  int index_count = 0;
  int bytes_for_datatype;
  int64_t hz_index;
  int64_t total_compressed_patch_size;
  
  int64_t total_compression_block_size = id->idx_ptr->compression_block_size[0] * id->idx_ptr->compression_block_size[1] * id->idx_ptr->compression_block_size[2] * id->idx_ptr->compression_block_size[3] * id->idx_ptr->compression_block_size[4];
  
  int compressed_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int compressed_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  
  if(id->idx_ptr->variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return 1;
  }
  
  for (y = 0; y < id->idx_ptr->variable[id->start_var_index]->patch_group_count; y++)
  { 
    if(id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box_group_type == 0)
    {
      for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
      {
        compressed_patch_offset[l] = id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[l] / id->idx_ptr->compression_block_size[l];
        compressed_patch_size[l] = id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[l] / id->idx_ptr->compression_block_size[l];
      }
      
      total_compressed_patch_size = (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[0] / id->idx_ptr->compression_block_size[0]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[1] / id->idx_ptr->compression_block_size[1]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[2] / id->idx_ptr->compression_block_size[2]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[3] / id->idx_ptr->compression_block_size[3]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[4] / id->idx_ptr->compression_block_size[4]);
    
      index_count = 0;
      hz_tupple *tupple = malloc(total_compressed_patch_size * sizeof(hz_tupple));
      
      //printf("Patch %d [%d]: %lld %lld %lld %lld %lld :: %d %d %d %d %d\n", y, id->idx_derived_ptr->patch_count, l_x, l_y, l_z, l_u, l_v, id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[0], id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[1], id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[2], id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[3], id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[4]);
      
      number_levels = id->idx_derived_ptr->maxh - 1;
      PointND xyzuv_Index;
      
      if(id->idx_ptr->variable[id->start_var_index]->data_layout == PIDX_row_major)
      {
        for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
        {
          for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
          {
            for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
            {
              for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
              {
                for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++) 
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
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_derived_ptr->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  level = getLeveL(hz_order);
                  
                  index = (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * compressed_patch_size[3] * (v - compressed_patch_offset[4])) 
                        + (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * (u - compressed_patch_offset[3])) 
                        + (compressed_patch_size[0] * compressed_patch_size[1] * (k - compressed_patch_offset[2])) 
                        + (compressed_patch_size[0] * (j - compressed_patch_offset[1])) 
                        + (i - compressed_patch_offset[0]);
                  
                  tupple[index_count].value = malloc((id->end_var_index - id->start_var_index + 1)* sizeof(unsigned char**));
                  id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                  
                  for(var = id->start_var_index; var <= id->end_var_index; var++)
                  {
                    tupple[index_count].value[var - id->start_var_index] = malloc(id->idx_ptr->variable[var]->values_per_sample * sizeof(unsigned char*));
                    bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                    
                    for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++) 
                    {
                      tupple[index_count].value[var - id->start_var_index][s] = malloc(bytes_for_datatype * total_compression_block_size);
                      memcpy(tupple[index_count].value[var - id->start_var_index][s], id->idx_ptr->variable[var]->patch[y]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size, bytes_for_datatype * total_compression_block_size);
                      
#if 0
                      if(bytes_for_datatype == sizeof(double))
                      {
                        double dvalue;
                        memcpy(&value, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
                        memcpy(&value, id->idx_ptr->variable[var]->patch[y]->buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                        printf("HZD: [%d %d %d %d] %d %f\n", y, index_count, var, s, bytes_for_datatype, value);
                      }
                      else
                      {
                        int ivalue;
                        memcpy(&ivalue, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
                        memcpy(&ivalue, id->idx_ptr->variable[var]->patch[y]->buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                        printf("HZI: [%d %d %d %d] %d %d\n", y, index_count, var, s, bytes_for_datatype, ivalue);
                      }
#endif
                      tupple[index_count].index = hz_order;
                    }
                  }
                  index_count++;
                }
              }
            }
          }
        }
      }
      else
      {
        for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
        {
          for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
          {
            for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
            {
              for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
              {
                for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++)
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
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_derived_ptr->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  level = getLeveL(hz_order);
                                  
                  index = (compressed_patch_size[2] * compressed_patch_size[1] * (i - compressed_patch_offset[0])) 
                        + (compressed_patch_size[2] * (j - compressed_patch_offset[1])) 
                        + (k - compressed_patch_offset[2]);
                  
                  tupple[index_count].value = malloc((id->end_var_index - id->start_var_index + 1)* sizeof(unsigned char**));
                  id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                  
                  for(var = id->start_var_index; var <= id->end_var_index; var++)
                  {
                    tupple[index_count].value[var - id->start_var_index] = malloc(id->idx_ptr->variable[var]->values_per_sample * sizeof(unsigned char*));
                    bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                    
                    for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++) 
                    {
                      tupple[index_count].value[var - id->start_var_index][s] = malloc(bytes_for_datatype * total_compression_block_size);
                      memcpy(tupple[index_count].value[var - id->start_var_index][s], id->idx_ptr->variable[var]->patch[y]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size, bytes_for_datatype * total_compression_block_size);

#if 0
                      if(bytes_for_datatype == sizeof(double))
                      {
                        double dvalue;
                        memcpy(&value, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
                        memcpy(&value, id->idx_ptr->variable[var]->patch[y]->buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                        printf("HZD: [%d %d %d %d] %d %f\n", y, index_count, var, s, bytes_for_datatype, value);
                      }
                      else
                      {
                        int ivalue;
                        memcpy(&ivalue, tupple[index_count].value[var - id->start_var_index][s], bytes_for_datatype);
                        memcpy(&ivalue, id->idx_ptr->variable[var]->patch[y]->buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                        printf("HZI: [%d %d %d %d] %d %d\n", y, index_count, var, s, bytes_for_datatype, ivalue);
                      }
#endif
                      
                      tupple[index_count].index = hz_order;
                    }
                  }
                  index_count++;
                }
              }
            }
          }
        }
      }
      qsort( tupple, (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * compressed_patch_size[3] * compressed_patch_size[4]), sizeof(hz_tupple), compare );
      
      for(var = id->start_var_index; var <= id->end_var_index; var++)
        for(c = 0 ; c < id->idx_derived_ptr->maxh ; c++)
        {
          bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
          id->idx_ptr->variable[var]->HZ_patch[y]->buffer[c] = malloc(bytes_for_datatype * id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c] * id->idx_ptr->variable[var]->values_per_sample * total_compression_block_size);
          memset(id->idx_ptr->variable[var]->HZ_patch[y]->buffer[c], 0, bytes_for_datatype * id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c] * id->idx_ptr->variable[var]->values_per_sample * total_compression_block_size);
        }
      
      cnt = 0;
      
      for(c = 0; c < id->idx_derived_ptr->maxh; c++)
      {
        for(s = 0; s < id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[c]; s++)
        {
          for(var = id->start_var_index; var <= id->end_var_index; var++)
          {
            for (i = 0; i < id->idx_ptr->variable[var]->values_per_sample; i++) 
            {
              bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
              memcpy(id->idx_ptr->variable[var]->HZ_patch[y]->buffer[c] + ((s * id->idx_ptr->variable[var]->values_per_sample + i) * bytes_for_datatype * total_compression_block_size), tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype * total_compression_block_size);
              id->idx_ptr->variable[var]->HZ_patch[y]->buffer_index[cnt] = tupple[cnt].index;
      
#if 0
              if(id->idx_ptr->variable[var]->bits_per_value / 8 == sizeof(double))
              {
                double check;
                memcpy(&check, tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype);
                memcpy(&check, id->idx_ptr->variable[var]->HZ_patch[y]->buffer[c] + ((s * id->idx_ptr->variable[var]->values_per_sample + i) * bytes_for_datatype), bytes_for_datatype);
                printf("[HD%d] %d %d %d %d %d %f\n", var, cnt, y, c, (s * id->idx_ptr->variable[var]->values_per_sample + i), bytes_for_datatype, check );
              }
              else
              {
                int icheck;
                memcpy(&icheck, tupple[cnt].value[var - id->start_var_index][i], bytes_for_datatype);
                memcpy(&icheck, id->idx_ptr->variable[var]->HZ_patch[y]->buffer[c] + ((s * id->idx_ptr->variable[var]->values_per_sample + i) * bytes_for_datatype), bytes_for_datatype);
                printf("[HI%d] %d %d %d %d %d %d\n", var, cnt, y, c, (s * id->idx_ptr->variable[var]->values_per_sample + i), bytes_for_datatype, icheck );
              }
#endif
              
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
    else
    {
      for (b = 0; b < id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box_count; b++) 
      {
        index_count = 0;
        total_compressed_patch_size = 0;
        
        for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
        {
          compressed_patch_offset[l] = id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box[b]->Ndim_box_offset[l] / id->idx_ptr->compression_block_size[l];
          compressed_patch_size[l] = id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box[b]->Ndim_box_size[l] / id->idx_ptr->compression_block_size[l];  
          total_compressed_patch_size = total_compressed_patch_size * compressed_patch_size[l];
        }
        
        number_levels = id->idx_derived_ptr->maxh - 1;
        PointND xyzuv_Index;
      
        if(id->idx_ptr->variable[id->start_var_index]->data_layout == PIDX_row_major)
        {
          for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
          {
            for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
            {
              for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
              {
                for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
                {
                  for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++)
                  {
                    index = (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * compressed_patch_size[3] * (v - compressed_patch_offset[4])) 
                          + (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * (u - compressed_patch_offset[3])) 
                          + (compressed_patch_size[0] * compressed_patch_size[1] * (k - compressed_patch_offset[2])) 
                          + (compressed_patch_size[0] * (j - compressed_patch_offset[1])) 
                          + (i - compressed_patch_offset[0]);
                          
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
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = id->idx_derived_ptr->maxh - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    
                    level = getLeveL(hz_order);

                    //if (level >= id->idx_derived_ptr->maxh - id->idx_derived_ptr->resolution_to - 1)
                    //  continue;
                    id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                    hz_index = hz_order - id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->start_hz_index[level];
                    for(var = id->start_var_index; var <= id->end_var_index; var++)
                    {
                      bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                      for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++)
                      {
                        memcpy(id->idx_ptr->variable[var]->HZ_patch[y]->buffer[level] + ((hz_index * id->idx_ptr->variable[var]->values_per_sample + s) * bytes_for_datatype * total_compression_block_size), 
                                id->idx_ptr->variable[var]->post_rst_block[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size,
                                bytes_for_datatype * total_compression_block_size);
                        
#if 0
                        //if (hz_order == 985668 || hz_order == 985669 || hz_order == 1971392 || hz_order == 1971394)
                        //if (i == 31 && j == 127 && k == 155)
                        //int64_t dvalue;
                        //memcpy(&dvalue, id->idx_ptr->variable[var]->patch_group_ptr[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                        //if(dvalue == 638051)
                        //{
                          double dvalue;
                          double dvalue2;
                          memcpy(&dvalue, id->idx_ptr->variable[var]->patch_group_ptr[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);
                          memcpy(&dvalue2, id->idx_ptr->variable[var]->HZ_patch[y]->buffer[level] + ((hz_index * id->idx_ptr->variable[var]->values_per_sample + s) * bytes_for_datatype), bytes_for_datatype);
                          printf("HZD: [%lld%lld%lld] (%lld %lld) : (%f %f)\n", k, j, i, index, hz_order, dvalue, dvalue2);
                        //}
#endif
                      }
                    }
                  }
                }
              }
            }
          }
          
        }
        else
        {
          for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
          {
            for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
            {
              for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
              {
                for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
                {
                  for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++)
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
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = id->idx_derived_ptr->maxh - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    level = getLeveL(hz_order);
                    
                    index = (compressed_patch_size[2] * compressed_patch_size[1] * (i - compressed_patch_offset[0])) 
                          + (compressed_patch_size[2] * (j - compressed_patch_offset[1])) 
                          + (k - compressed_patch_offset[2]);

                    id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                    for(var = id->start_var_index; var <= id->end_var_index; var++)
                    {
                      hz_index = hz_order - id->idx_ptr->variable[var]->HZ_patch[y]->start_hz_index[level];
                      bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                      for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++)
                      {
                        memcpy(id->idx_ptr->variable[var]->HZ_patch[y]->buffer[level] + ((hz_index * id->idx_ptr->variable[var]->values_per_sample + s) * bytes_for_datatype * total_compression_block_size), 
                               id->idx_ptr->variable[var]->post_rst_block[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size, 
                               bytes_for_datatype * total_compression_block_size);                        
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box_group_type == 2)
    {
      for (i = id->start_var_index; i <= id->end_var_index; i++)
      {
        bytes_for_datatype = id->idx_ptr->variable[i]->bits_per_value / 8;  
        for (j = 0; j < id->idx_derived_ptr->maxh; j++) 
        {
          int level_blocks = 0;
          if (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] != 0)
          {
            int adjusted_buffer_size = (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1 - (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample;
            
            //printf("Initial Size %lld Adjusted size %lld (%lld - %d)\n", (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1), adjusted_buffer_size/(bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample), (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1), (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] * id->idx_derived_ptr->samples_per_block));
            
            int initial_m = 0;
            for (n = 0; n < id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j]; n++)
            {
              //printf("[%d] [L %d] start %lld end %lld\n", n, j, id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j], id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j]);
              for (m = initial_m; m < (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1); m++)
              {
                if (m + id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] == id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n] * id->idx_derived_ptr->samples_per_block)
                {
                  //printf("[MISS] (%d + %d) %d = %d (%d x %d)\n", m, id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j], (m + id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j]), id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n] * id->idx_derived_ptr->samples_per_block, id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n], id->idx_derived_ptr->samples_per_block);                  
                  //printf("[SOURCE] %d [DEST] %d [COUNT] %d\n", (m - (level_blocks * id->idx_derived_ptr->samples_per_block)), ((m + id->idx_derived_ptr->samples_per_block) - (level_blocks * id->idx_derived_ptr->samples_per_block)), (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1 - ((m + level_blocks) + id->idx_derived_ptr->samples_per_block) ));
                  
                  memmove(id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + (m - (level_blocks * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample, id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + ((m + id->idx_derived_ptr->samples_per_block) - (level_blocks * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample, (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - (id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + m + id->idx_derived_ptr->samples_per_block) + 1) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample /*((id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1 - ((m + level_blocks) + id->idx_derived_ptr->samples_per_block) ) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample)*/);
                  level_blocks++;
                  initial_m = m;
                }
              }
            }
            
            unsigned char* temp_buffer = realloc(id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j], adjusted_buffer_size);
            if (temp_buffer == NULL) 
            {
              
            }
            else
            {
              id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] = temp_buffer;
            }
            //int64_t dvalue;
            //if (j == 20)
            //{
            //  memcpy(&dvalue, id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + 887364 * bytes_for_datatype, bytes_for_datatype);
            //  printf("[XXX] Value at 887364 = %lld\n", dvalue);
            //}
          }
        }
      }
    }
  }
  return 0;
}

int PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  int64_t z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, s = 0, y = 0, n = 0, m = 0, number_levels = 0, var = 0;
  int64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0;
  int bytes_for_datatype;
  int64_t hz_index;
  int64_t total_compressed_patch_size;
  
  int64_t total_compression_block_size = id->idx_ptr->compression_block_size[0] * id->idx_ptr->compression_block_size[1] * id->idx_ptr->compression_block_size[2] * id->idx_ptr->compression_block_size[3] * id->idx_ptr->compression_block_size[4];
  
  int compressed_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int compressed_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  
  if(id->idx_ptr->variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return -1;
  }
  
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh not set.\n", __FILE__, __LINE__);
    return -1;
  }
  
  
  for (y = 0; y < id->idx_ptr->variable[id->start_var_index]->patch_group_count; y++)
  {
    for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
    {
      compressed_patch_offset[l] = id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_offset[l] / id->idx_ptr->compression_block_size[l];
      compressed_patch_size[l] = id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[l] / id->idx_ptr->compression_block_size[l];
    }
    
    total_compressed_patch_size = (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[0] / id->idx_ptr->compression_block_size[0]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[1] / id->idx_ptr->compression_block_size[1]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[2] / id->idx_ptr->compression_block_size[2]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[3] / id->idx_ptr->compression_block_size[3]) * (id->idx_ptr->variable[id->start_var_index]->patch[y]->Ndim_box_size[4] / id->idx_ptr->compression_block_size[4]);
  
    for (b = 0; b < id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box_count; b++) 
    {
      total_compressed_patch_size = 0;
      
      for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
      {
        compressed_patch_offset[l] = id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box[b]->Ndim_box_offset[l] / id->idx_ptr->compression_block_size[l];
        compressed_patch_size[l] = id->idx_ptr->variable[id->start_var_index]->post_rst_block[y]->box[b]->Ndim_box_size[l] / id->idx_ptr->compression_block_size[l];  
        total_compressed_patch_size = total_compressed_patch_size * compressed_patch_size[l];
      }
      
      number_levels = id->idx_derived_ptr->maxh - 1;
      PointND xyzuv_Index;
    
      if(id->idx_ptr->variable[id->start_var_index]->data_layout == PIDX_row_major)
      {
        for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
          for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
            for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
              for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
                for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++)
                {
                  index = (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * compressed_patch_size[3] * (v - compressed_patch_offset[4])) 
                        + (compressed_patch_size[0] * compressed_patch_size[1] * compressed_patch_size[2] * (u - compressed_patch_offset[3])) 
                        + (compressed_patch_size[0] * compressed_patch_size[1] * (k - compressed_patch_offset[2])) 
                        + (compressed_patch_size[0] * (j - compressed_patch_offset[1])) 
                        + (i - compressed_patch_offset[0]);
                        
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
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_derived_ptr->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  
                  level = getLeveL(hz_order);
                  id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                  hz_index = hz_order - id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->start_hz_index[level];
                  for(var = id->start_var_index; var <= id->end_var_index; var++)
                  {
                    bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                    for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++)
                    {                        
                      memcpy(id->idx_ptr->variable[var]->post_rst_block[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size,
                             id->idx_ptr->variable[var]->HZ_patch[y]->buffer[level] + ((hz_index * id->idx_ptr->variable[var]->values_per_sample + s) * bytes_for_datatype * total_compression_block_size),
                             bytes_for_datatype * total_compression_block_size);
                    }
                  }
                }
      }
      else
      {
        for (v = compressed_patch_offset[4]; v < compressed_patch_offset[4] + compressed_patch_size[4]; v++)
          for (u = compressed_patch_offset[3]; u < compressed_patch_offset[3] + compressed_patch_size[3]; u++)
            for (k = compressed_patch_offset[2]; k < compressed_patch_offset[2] + compressed_patch_size[2]; k++)
              for (j = compressed_patch_offset[1]; j < compressed_patch_offset[1] + compressed_patch_size[1]; j++)
                for (i = compressed_patch_offset[0]; i < compressed_patch_offset[0] + compressed_patch_size[0]; i++)
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
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_derived_ptr->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  level = getLeveL(hz_order);
                  
                  index = (compressed_patch_size[2] * compressed_patch_size[1] * (i - compressed_patch_offset[0])) 
                        + (compressed_patch_size[2] * (j - compressed_patch_offset[1])) 
                        + (k - compressed_patch_offset[2]);

                  id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] = id->idx_ptr->variable[id->start_var_index]->HZ_patch[y]->samples_per_level[level] + 1;
                  for(var = id->start_var_index; var <= id->end_var_index; var++)
                  {
                    hz_index = hz_order - id->idx_ptr->variable[var]->HZ_patch[y]->start_hz_index[level];
                    bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                    for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++)
                    {
                      bytes_for_datatype = id->idx_ptr->variable[var]->bits_per_value / 8;
                      memcpy(id->idx_ptr->variable[var]->post_rst_block[y]->box[b]->Ndim_box_buffer + ((index * id->idx_ptr->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_compression_block_size, id->idx_ptr->variable[var]->HZ_patch[y]->buffer[level] + ((hz_index * id->idx_ptr->variable[var]->values_per_sample + s) * bytes_for_datatype * total_compression_block_size), bytes_for_datatype * total_compression_block_size);                        
                    }
                  }
                }
      }
    }
    
    if (id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box_group_type == 2)
    {
      for (i = id->start_var_index; i <= id->end_var_index; i++)
      {
        bytes_for_datatype = id->idx_ptr->variable[i]->bits_per_value / 8;  
        for (j = 0; j < id->idx_derived_ptr->maxh; j++) 
        {
          int level_blocks = 0;
          if (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] != 0)
          {
            int adjusted_buffer_size = (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1 - (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample;
            
            //printf("Initial Size %lld Adjusted size %lld (%lld - %d)\n", (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1), adjusted_buffer_size/(bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample), (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1), (id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j] * id->idx_derived_ptr->samples_per_block));
            
            int initial_m = 0;
            for (n = 0; n < id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_count_per_level[j]; n++)
            {
              //printf("[%d] [L %d] start %lld end %lld\n", n, j, id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j], id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j]);
              for (m = initial_m; m < (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1); m++)
              {
                if (m + id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] == id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n] * id->idx_derived_ptr->samples_per_block)
                {
                  //printf("[MISS] (%d + %d) %d = %d (%d x %d)\n", m, id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j], (m + id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j]), id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n] * id->idx_derived_ptr->samples_per_block, id->idx_ptr->variable[i]->HZ_patch[y]->missing_block_index_per_level[j][n], id->idx_derived_ptr->samples_per_block);                  
                  //printf("[SOURCE] %d [DEST] %d [COUNT] %d\n", (m - (level_blocks * id->idx_derived_ptr->samples_per_block)), ((m + id->idx_derived_ptr->samples_per_block) - (level_blocks * id->idx_derived_ptr->samples_per_block)), (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + 1 - ((m + level_blocks) + id->idx_derived_ptr->samples_per_block) ));
                  
                  memmove(id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + (m - (level_blocks * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample, 
                          id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + ((m + id->idx_derived_ptr->samples_per_block) - (level_blocks * id->idx_derived_ptr->samples_per_block)) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample, 
                          (id->idx_ptr->variable[i]->HZ_patch[y]->end_hz_index[j] - (id->idx_ptr->variable[i]->HZ_patch[y]->start_hz_index[j] + m + id->idx_derived_ptr->samples_per_block) + 1) * bytes_for_datatype * id->idx_ptr->variable[i]->values_per_sample );
                  level_blocks++;
                  initial_m = m;
                }
              }
            }
            
            unsigned char* temp_buffer = realloc(id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j], adjusted_buffer_size);
            if (temp_buffer == NULL) 
            {
              
            }
            else
            {
              id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] = temp_buffer;
            }
            
            //int64_t dvalue;
            //if (j == 20)
            //{
            //  memcpy(&dvalue, id->idx_ptr->variable[i]->HZ_patch[y]->buffer[j] + 887364 * bytes_for_datatype, bytes_for_datatype);
            //  printf("[XXX] Value at 887364 = %lld\n", dvalue);
            //}
          }
        }
      }
    }
  }
  return 0;
  
  /*
  int64_t hz_index = 0, index = 0;
  int64_t mins[PIDX_MAX_DIMENSIONS], maxes[PIDX_MAX_DIMENSIONS], l_x, l_y, l_z, l_u, l_v, hzaddress, hz_mins, hz_maxes;
  int i =0, j = 0, k = 0, d = 0, s = 0, buffer_count = 0, rank = 0, y = 0, spv, number_levels = 0, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int64_t lastbitmask = ((int64_t) 1) << id->maxh - 1;
  int64_t n = 0, m = 0;
    
  if(rank == 0)
    printf("[read] Cache Status: %d\n", var_cache);
  
  if(var_cache == 1)
  {
    for (y = 0; y < id->number_of_buffers; y++) 
    {
      if (num_regular_blocks > 1)
	if (y == id->regular_block_index[buffer_count + 1]) 
	  buffer_count++;
    }
    
    index_cache = (int***)malloc(sizeof(int**) * id->number_of_buffers);
    memset(index_cache, -1, sizeof(int**) * id->number_of_buffers);
    for(i = 0; i < id->number_of_buffers; i++)
    {
      index_cache[i] = (int**)malloc(sizeof(int*) * (id->maxh-id->resh));
      memset(index_cache[i], -1, sizeof(int*) * (id->maxh-id->resh));
      for (j = 0; j < (id->maxh-id->resh); j++) 
      {
	if ((out_buf_array[buffer_count]->n_samples[j][0] * out_buf_array[buffer_count]->n_samples[j][1] * out_buf_array[buffer_count]->n_samples[j][2] * out_buf_array[buffer_count]->n_samples[j][3] * out_buf_array[buffer_count]->n_samples[j][4]) != 0) 
	{
	  hz_mins = out_buf_array[buffer_count]->start_hz_index[j];
	  hz_maxes = out_buf_array[buffer_count]->end_hz_index[j] + 1;
      
	  index_cache[i][j] = (int*)malloc(sizeof(int) * (hz_maxes - hz_mins));
	  memset(index_cache[i][j], -1, sizeof(int) * (hz_maxes - hz_mins));
	}
      }
    }
      
    for(i = 0; i < id->number_of_buffers; i++)
    {
      for (j = 0; j < (id->maxh-id->resh); j++) 
      {
	if ((out_buf_array[buffer_count]->n_samples[j][0] * out_buf_array[buffer_count]->n_samples[j][1] * out_buf_array[buffer_count]->n_samples[j][2] * out_buf_array[buffer_count]->n_samples[j][3] * out_buf_array[buffer_count]->n_samples[j][4]) != 0) 
	{
	  hz_mins = out_buf_array[buffer_count]->start_hz_index[j];
	  hz_maxes = out_buf_array[buffer_count]->end_hz_index[j] + 1;
	    
	  for(k = 0; k < hz_maxes-hz_mins ; k++)
	    index_cache[i][j][k] = -1;//(int*)malloc(sizeof(int) * (hz_maxes - hz_mins));
	}
      }
    }
  }
  
  for (y = 0; y < id->idx_ptr->variable[id->start_var_index]->patch_group_count; y++)
  {
    for (b = 0; b < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->count; b++)
    {
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++) 
      {
	mins[d] = in_buf[y]->lower_bounds[d];
	maxes[d] = in_buf[y]->upper_bounds[d];
      }
      l_x = maxes[0] - mins[0];
      l_y = maxes[1] - mins[1];
      l_z = maxes[2] - mins[2];
      l_u = maxes[3] - mins[3];
      l_v = maxes[4] - mins[4];
      
      for (j = 0; j < id->maxh-id->resh; j++) 
      {
	hz_mins = id->idx_ptr->variable[var]->HZ_patch[y]->start_hz_index[j];
	hz_maxes = id->idx_ptr->variable[var]->HZ_patch[y]->end_hz_index[j] + 1;

	number_levels = id->idx_derived_ptr->maxh - 1;

	for (m = hz_mins; m < hz_maxes; m++) 
	{
	  if(var_cache == 1)
	  {
	    hzaddress = m;
	    hzaddress <<= 1;
	    hzaddress |= 1;
	    while ((lastbitmask & hzaddress) == 0) hzaddress <<= 1;
	    hzaddress &= lastbitmask - 1;

	    PointND cnt;
	    PointND p;
	    n = 0;

	    memset(&cnt, 0, sizeof (PointND));
	    memset(&p, 0, sizeof (PointND));

	    for (; hzaddress; hzaddress >>= 1, ++n, number_levels--) 
	    {
	      int bit = id->bitPattern[number_levels];
	      PGET(p, bit) |= (hzaddress & 1) << PGET(cnt, bit);
	      ++PGET(cnt, bit);
	    }
	    number_levels = id->maxh - 1;
	    
	    if (p.x >= id->idx_ptr->global_bounds[0] || p.y >= id->idx_ptr->global_bounds[1] || p.z >= id->idx_ptr->global_bounds[2] || p.u >= id->idx_ptr->global_bounds[3] || p.v >= id->idx_ptr->global_bounds[4])
	      continue;
	    
	    if (p.x < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[0] || p.y < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[1] || p.z < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[2] || p.u < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[3] || p.v < id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[4]) 
	      continue;
	    
	    if (p.x >= id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[0] + id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->count[0] || p.y >= id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[1] + id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->count[1] || p.z >= id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[2] + id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->count[2] || p.u >= id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[3] + id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->count[3] || p.v >= id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->Ndim_box_offset[4] + id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[y]->box[b]->count[4]) 
	      continue;
	    

	    hz_index = m - hz_mins;
	    index = (l_x * l_y * (p.z - in_buf[y]->lower_bounds[2])) + (l_x * (p.y - in_buf[y]->lower_bounds[1])) + (p.x - in_buf[y]->lower_bounds[0]);
	    
	    if(m-hz_mins >= 0 && m-hz_mins < (hz_maxes - hz_mins))
	      index_cache[y][j][m-hz_mins] = index;
	    
	    spv = out_buf_array[buffer_count]->sample_per_id->idx_ptr->variable;
	    for (s = 0; s < spv; s++) 
	      in_buf[y]->buffer[(index * spv) + s] = out_buf_array[buffer_count]->buffer[j][(hz_index * spv) + s];
		      
		      
	    if(j == (id->maxh-id->resh -1) && m == (hz_maxes - 1) && y == (id->number_of_buffers - 1))
	      var_cache = 0;
	  }
	  else
	  {
	    if(index_cache[y][j][m-hz_mins] == -1)
		continue;
	    hz_index = m - hz_mins;
	    spv = out_buf_array[buffer_count]->sample_per_id->idx_ptr->variable;
	    memcpy(in_buf[y]->buffer + (index_cache[y][j][m-hz_mins] * spv), out_buf_array[buffer_count]->buffer[j] + (hz_index * spv), sizeof(double)*spv );
	  }
	}
	
      }
    }
  }
  */
  return 1;
}


/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
int PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id) 
{
  int itr = 0, p = 0, var = 0;
  
  if(id->idx_ptr->variable[id->start_var_index]->patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_derived_ptr->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_derived_ptr->maxh (%d) not set.\n", __FILE__, __LINE__, id->idx_derived_ptr->maxh);
    return 1;
  }
  
  for (p = 0; p < id->idx_ptr->variable[id->start_var_index]->patch_group_count; p++)
  {
    free(id->idx_ptr->variable[id->start_var_index]->HZ_patch[p]->samples_per_level);
    id->idx_ptr->variable[id->start_var_index]->HZ_patch[p]->samples_per_level = 0;
  }
  
  for(var = id->start_var_index; var <= id->end_var_index; var++)
  {
    for (p = 0; p < id->idx_ptr->variable[id->start_var_index]->patch_group_count; p++)
    {
      free(id->idx_ptr->variable[var]->HZ_patch[p]->start_hz_index);
      free(id->idx_ptr->variable[var]->HZ_patch[p]->end_hz_index);
      
      free(id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_count_per_level);     
      
      if(id->idx_ptr->variable[id->start_var_index]->patch_group_ptr[p]->box_group_type == 0)
	free(id->idx_ptr->variable[var]->HZ_patch[p]->buffer_index);
      
      for(itr = 0 ; itr < id->idx_derived_ptr->maxh ; itr++)
      {
	free(id->idx_ptr->variable[var]->HZ_patch[p]->buffer[itr]);
	id->idx_ptr->variable[var]->HZ_patch[p]->buffer[itr] = 0;
        
        free(id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_index_per_level[itr]);
        id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_index_per_level[itr] = 0;
      }
      
      free(id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_index_per_level);
      id->idx_ptr->variable[var]->HZ_patch[p]->missing_block_index_per_level = 0;
      
      free(id->idx_ptr->variable[var]->HZ_patch[p]->buffer);
      id->idx_ptr->variable[var]->HZ_patch[p]->buffer = 0;
      
      //free(id->idx_ptr->variable[var]->HZ_patch[p]);
      //id->idx_ptr->variable[var]->HZ_patch[p] = 0;
    }
  }
  
  return 0;
}

int PIDX_hz_encode_finalize(PIDX_hz_encode_id id) 
{
  
  free(id);
  id = 0;
  
  return 0;
}


int HELPER_Hz_encode(PIDX_hz_encode_id id)
{
  int i = 0, k = 0, b = 0, var = 0, rank;
  int64_t global_hz, element_count = 0, lost_element_count = 0;
  int64_t ZYX[PIDX_MAX_DIMENSIONS];
  int check_bit = 1, s = 0;
  
#if long_buffer
  uint64_t dvalue_1, dvalue_2;
#else
  double dvalue_1, dvalue_2;
#endif

#if PIDX_HAVE_MPI
  MPI_Comm_rank(id->comm, &rank);
#endif
  
  for(var = id->start_var_index; var <= id->end_var_index; var++)
  {
    for (b = 0; b < id->idx_ptr->variable[var]->patch_group_count; b++)
    {
      for (i = 0; i < id->idx_ptr->variable[var]->HZ_patch[b]->HZ_level_to - id->idx_ptr->variable[var]->HZ_patch[b]->HZ_level_from; i++) 
      {
        //printf("samples at level %d = %d\n", i, id->idx_ptr->variable[id->start_var_index]->HZ_patch[b]->samples_per_level[i]);
        if (id->idx_ptr->variable[id->start_var_index]->HZ_patch[b]->samples_per_level[i] != 0)
        {
          for (k = 0; k <= (id->idx_ptr->variable[var]->HZ_patch[b]->end_hz_index[i] - id->idx_ptr->variable[var]->HZ_patch[b]->start_hz_index[i]) * 1; k++) 
          {
            global_hz = id->idx_ptr->variable[var]->HZ_patch[b]->start_hz_index[i] + k;
            
            Hz_to_xyz(id->idx_ptr->bitPattern, id->idx_derived_ptr->maxh - 1, global_hz, ZYX);
            if (!(ZYX[0] >= id->idx_ptr->global_bounds[0] || ZYX[1] >= id->idx_ptr->global_bounds[1] || ZYX[2] >= id->idx_ptr->global_bounds[2])) 
            {
              check_bit = 1, s = 0;    
              for (s = 0; s < id->idx_ptr->variable[var]->values_per_sample; s++)
              {
                dvalue_1 = 100 + var + (id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1]*(ZYX[2]))+(id->idx_ptr->global_bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_derived_ptr->color * id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2]);
#if long_buffer
                dvalue_2 = *(*((uint64_t**)id->idx_ptr->variable[var]->HZ_patch[b]->buffer + i) + ((k * id->idx_ptr->variable[var]->values_per_sample) + s));
#else
                dvalue_2 = *(*((double**)id->idx_ptr->variable[var]->HZ_patch[b]->buffer + i) + ((k * id->idx_ptr->variable[var]->values_per_sample) + s));
#endif
                
                check_bit = check_bit && (dvalue_1  == dvalue_2);
                if (check_bit == 0)
                {
                  //printf("[HZ] %f %f (%lld :: %lld %lld %lld)\n", dvalue_1, dvalue_2, (long long)global_hz, (long long)ZYX[0], (long long)ZYX[1], (long long)ZYX[2]);
                  lost_element_count++;
                }
                else
                {
                  //printf("HZ [%d] %f %f\n", rank, dvalue_1, dvalue_2);
                  element_count++;
                }
              }
            }
          }
        }
      }
    }
  }
 
#if PIDX_HAVE_MPI
  int64_t global_volume = 0;
  MPI_Allreduce(&element_count, &global_volume, 1, MPI_LONG_LONG, MPI_SUM, id->comm);
  
  //printf("[HZ] Volume [%lld] and Volume [%lld]\n", global_volume, (int64_t)(id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2] * id->idx_ptr->global_bounds[3] * id->idx_ptr->global_bounds[4] * (id->end_var_index - id->start_var_index + 1)));
  
  if (global_volume != (int64_t) id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2] * (id->end_var_index - id->start_var_index + 1)) 
  {
    if (rank == 0)
      fprintf(stderr, "[HZ Debug FAILED!!!!] [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", id->idx_derived_ptr->color, (long long) global_volume, (long long) id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2] * (id->end_var_index - id->start_var_index + 1));
    
    printf("[HZ]  Rank %d Color %d [LOST ELEMENT COUNT %lld] [FOUND ELEMENT COUNT %lld] [TOTAL ELEMNTS %lld] \n", rank,  id->idx_derived_ptr->color, (long long) lost_element_count, (long long) element_count, (long long) (id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2] * id->idx_ptr->global_bounds[3] * id->idx_ptr->global_bounds[4]) * (id->end_var_index - id->start_var_index + 1));
    
    return (-1);
  }
  else
  {
    if (rank == 0)
      fprintf(stderr, "[HZ Debug PASSED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", id->idx_derived_ptr->color, (long long) global_volume, (long long) id->idx_ptr->global_bounds[0] * id->idx_ptr->global_bounds[1] * id->idx_ptr->global_bounds[2] * (id->end_var_index - id->start_var_index + 1));
  }
#endif
    
  return 0;
}
