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
  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;
  
  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;
  
  int** index;
  
  int init_index;
  int first_index;
  int last_index;
  
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
};

static int compare( const void* a, const void* b);

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_hz_encode_set_communicator(PIDX_hz_encode_id hz_id, MPI_Comm comm)
{
  if (hz_id == NULL)
    return PIDX_err_id;

  hz_id->comm = comm;

  return PIDX_success;
}
#endif


static int compare( const void* a, const void* b)
{
  int64_t int_a = ((const hz_tupple*)a)->index;
  int64_t int_b = ((const hz_tupple*)b)->index;

  if ( int_a == int_b ) return 0;
  else if ( int_a < int_b ) return -1;
  else return 1;
}



PIDX_hz_encode_id PIDX_hz_encode_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_hz_encode_id hz_id;
  hz_id = (PIDX_hz_encode_id)malloc(sizeof (*hz_id));
  memset(hz_id, 0, sizeof (*hz_id));

  hz_id->idx = idx_meta_data;
  hz_id->idx_d = idx_d;

  hz_id->init_index = init_index;
  hz_id->first_index = first_index;
  hz_id->last_index = last_index;
  
  return hz_id;
}


PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id)
{
  int **tpatch;
  int **allign_offset;
  int **allign_count;
  int start_block_no, end_block_no, b;
  int j = 0, p = 0, d = 0, c = 0, v = 0, bytes_for_datatype = 0, count = 0;

  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];
  PIDX_variable var0 = id->idx->variable[id->first_index];

  tpatch = (int**) malloc(2 * sizeof (int*));
  tpatch[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  tpatch[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  assert(tpatch);

  printf("id->idx->variable[id->first_index]->sim_patch_count %d\n", id->idx->variable[id->first_index]->sim_patch_count);
  if(id->idx->variable[id->first_index]->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->count not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  if(maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var0->patch_group_count; p++)
    {
      var->hz_buffer[p] = malloc(sizeof(*(var->hz_buffer[p])));
      memset(var->hz_buffer[p], 0, sizeof(*(var->hz_buffer[p])));
    }
  }

  for (p = 0; p < var0->patch_group_count; p++)
  {
    HZ_buffer hz_buf = var0->hz_buffer[p];
    if(var0->chunk_patch_group[p]->type == 0)
    {
      uint64_t buffer_size = 1;
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        buffer_size = buffer_size * (var0->chunk_patch_group[p]->patch[0]->size[d]/id->idx->chunk_size[d]);

      hz_buf->buffer_index = malloc(buffer_size * sizeof (int64_t));
      memset(hz_buf->buffer_index, 0, buffer_size);
    }

    hz_buf->type = var0->chunk_patch_group[p]->type;
    hz_buf->HZ_agg_from = 0;
    hz_buf->HZ_agg_to = maxH;
    hz_buf->HZ_io_from = 0;
    hz_buf->HZ_io_to = maxH;

    hz_buf->start_hz_index = malloc(sizeof (int64_t) * maxH);
    hz_buf->end_hz_index = malloc(sizeof (int64_t) * maxH);
    hz_buf->nsamples_per_level = malloc(sizeof (int*) * maxH);
    hz_buf->samples_per_level = malloc( maxH * sizeof (int64_t));
    memset(hz_buf->start_hz_index, 0, sizeof (int64_t) * maxH);
    memset(hz_buf->end_hz_index, 0, sizeof (int64_t) * maxH);
    memset(hz_buf->nsamples_per_level, 0, sizeof (int64_t) * maxH);
    memset(hz_buf->samples_per_level, 0, maxH * sizeof (int64_t));

    hz_buf->missing_block_count_per_level = malloc(sizeof (int) * maxH);
    hz_buf->missing_block_index_per_level = malloc(sizeof (int*) * maxH);
    memset(hz_buf->missing_block_index_per_level, 0, sizeof (int*) * maxH);
    memset(hz_buf->missing_block_count_per_level, 0, sizeof (int) * maxH);

    allign_offset = malloc(sizeof (int*) * maxH);
    allign_count = malloc(sizeof (int*) * maxH);
    memset(allign_offset, 0, sizeof (int*) * maxH);
    memset(allign_count, 0, sizeof (int*) * maxH);

    for (j = 0; j < maxH; j++)
    {
      allign_offset[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      allign_count[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      hz_buf->nsamples_per_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
      hz_buf->missing_block_index_per_level[j] = malloc(sizeof (int) * (id->idx->blocks_per_file));
      memset(hz_buf->missing_block_index_per_level[j], 0, (sizeof (int) * (id->idx->blocks_per_file)));
    }

    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      tpatch[0][d] = var0->chunk_patch_group[p]->reg_patch_offset[d] / id->idx->chunk_size[d];
      tpatch[1][d] = (var0->chunk_patch_group[p]->reg_patch_offset[d] / id->idx->chunk_size[d]) + (var0->chunk_patch_group[p]->reg_patch_size[d] / id->idx->chunk_size[d]) - 1;
    }

    for (j = 0; j < maxH; j++)
    {
      Align((maxH - 1), j, id->idx->bitPattern, tpatch, allign_offset, allign_count, hz_buf->nsamples_per_level);

      PointND startXYZ;
      startXYZ.x = allign_offset[j][0];
      startXYZ.y = allign_offset[j][1];
      startXYZ.z = allign_offset[j][2];
      startXYZ.u = allign_offset[j][3];
      startXYZ.v = allign_offset[j][4];
      hz_buf->start_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, startXYZ);

      PointND endXYZ;
      endXYZ.x = allign_count[j][0];
      endXYZ.y = allign_count[j][1];
      endXYZ.z = allign_count[j][2];
      endXYZ.u = allign_count[j][3];
      endXYZ.v = allign_count[j][4];

      hz_buf->end_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, endXYZ);

      if (id->idx->variable[id->first_index]->chunk_patch_group[p]->type == 2)
      {
        start_block_no = hz_buf->start_hz_index[j] / id->idx_d->samples_per_block;
        end_block_no = hz_buf->end_hz_index[j] / id->idx_d->samples_per_block;

        count = 0;
        for (b = start_block_no; b <= end_block_no; b++)
        {
          if (PIDX_blocks_is_block_present(b, id->idx->variable[id->init_index]->global_block_layout) == 0)
          {
            hz_buf->missing_block_count_per_level[j]++;
            hz_buf->missing_block_index_per_level[j][count] = b;
            count++;
          }
        }
      }
       free(allign_offset[j]);
       free(allign_count[j]);
    }
    free(allign_offset);
    free(allign_count);
  }

  free(tpatch[0]);
  tpatch[0] = 0;
  free(tpatch[1]);
  tpatch[1] = 0;
  free(tpatch);
  tpatch = 0;

  // Allocate actual HZ buffer for the variables
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var0->patch_group_count; p++)
    {
      var->hz_buffer[p]->buffer = (unsigned char**)malloc( maxH * sizeof (unsigned char*));
      memset(var->hz_buffer[p]->buffer, 0,  maxH * sizeof (unsigned char*));
      bytes_for_datatype = (var->bits_per_value / 8) * chunk_size * var->values_per_sample;

      if (var0->chunk_patch_group[p]->type == 1 || var0->chunk_patch_group[p]->type == 2)
      {
        for (c = 0 ; c < maxH ; c++)
        {
          int64_t samples_per_level = (var0->hz_buffer[p]->end_hz_index[c] - var0->hz_buffer[p]->start_hz_index[c] + 1);
          var->hz_buffer[p]->buffer[c] = malloc(bytes_for_datatype * samples_per_level);
          memset(var->hz_buffer[p]->buffer[c], 0, bytes_for_datatype * samples_per_level);
        }
      }
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  int64_t z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, c = 0, s = 0, y = 0, n = 0, m = 0, number_levels = 0, var = 0;
  int64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0, d = 0;
  int index_count = 0;
  int bytes_for_datatype;
  int64_t hz_index;
  int64_t total_chunked_patch_size = 1;
  
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];
  PIDX_variable var0 = id->idx->variable[id->first_index];
  
  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  
  if (var0->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->count not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  
  if (maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }
  
  for (y = 0; y < var0->patch_group_count; y++)
  {
    if (var0->chunk_patch_group[y]->type == 0)
    {
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
#if 1
        total_chunked_patch_size = 1;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          chunked_patch_offset[d] = var0->chunk_patch_group[y]->patch[b]->offset[d] / id->idx->chunk_size[d];
          chunked_patch_size[d] = var0->chunk_patch_group[y]->patch[b]->size[d] / id->idx->chunk_size[d];
          total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[d];
        }

        index_count = 0;
        hz_tupple *tupple = malloc(total_chunked_patch_size * sizeof(hz_tupple));
      
        number_levels = maxH - 1;
        PointND xyzuv_Index;
      
        if (var0->data_layout == PIDX_row_major)
        {
          for (m = chunked_patch_offset[4]; m < chunked_patch_offset[4] + chunked_patch_size[4]; m++)
            for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                  {
                    xyzuv_Index.x = i;
                    xyzuv_Index.y = j;
                    xyzuv_Index.z = k;
                    xyzuv_Index.u = u;
                    xyzuv_Index.v = m;

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
                      int bit = id->idx->bitPattern[number_levels];
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = maxH - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    level = getLeveL(hz_order);

                    index = (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * chunked_patch_size[3] * (m - chunked_patch_offset[4])) +
                            (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * (u - chunked_patch_offset[3])) +
                            (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                        + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                        + (i - chunked_patch_offset[0]);
                  
                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;

                    for(v = id->first_index; v <= id->last_index; v++)
                    {
                      PIDX_variable var = id->idx->variable[v];
                      tupple[index_count].value[v - id->first_index] = malloc(var->values_per_sample * sizeof(unsigned char*));
                      bytes_for_datatype = var->bits_per_value / 8;
                    
                      for (s = 0; s < var->values_per_sample; s++)
                      {
                        tupple[index_count].value[v - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memcpy(tupple[index_count].value[v- id->first_index][s], var->chunk_patch_group[y]->patch[b]->buffer + ((index * var->values_per_sample) + s) * bytes_for_datatype * chunk_size, bytes_for_datatype * chunk_size);

                        tupple[index_count].index = hz_order;
                      }
                    }
                    index_count++;
                  }
        }
        else
        {
          for (m = chunked_patch_offset[4]; m < chunked_patch_offset[4] + chunked_patch_size[4]; m++)
            for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                  {
                    xyzuv_Index.x = i;
                    xyzuv_Index.y = j;
                    xyzuv_Index.z = k;
                    xyzuv_Index.u = u;
                    xyzuv_Index.v = m;

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
                      int bit = id->idx->bitPattern[number_levels];
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = maxH - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    level = getLeveL(hz_order);
                                  
                    index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                        + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                        + (k - chunked_patch_offset[2]);
                  
                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;
                  
                    for(var = id->first_index; var <= id->last_index; var++)
                    {
                      tupple[index_count].value[var - id->first_index] = malloc(id->idx->variable[var]->values_per_sample * sizeof(unsigned char*));
                      bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                    
                      for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                      {
                        tupple[index_count].value[var - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memcpy(tupple[index_count].value[var - id->first_index][s], id->idx->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[var]->values_per_sample) + s) * bytes_for_datatype * chunk_size, bytes_for_datatype * chunk_size);

                        tupple[index_count].index = hz_order;
                      }
                    }
                    index_count++;
                  }

        }
        qsort( tupple, total_chunked_patch_size, sizeof(hz_tupple), compare );
      
        for(var = id->first_index; var <= id->last_index; var++)
        {
          for(c = 0 ; c < maxH ; c++)
          {
            bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
            id->idx->variable[var]->hz_buffer[y]->buffer[c] = malloc(bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * id->idx->variable[var]->values_per_sample * chunk_size);
            memset(id->idx->variable[var]->hz_buffer[y]->buffer[c], 0, bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * id->idx->variable[var]->values_per_sample * chunk_size);
          }
        }
      
        cnt = 0;
        for(c = 0; c < maxH; c++)
        {
          for(s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
          {
            for(var = id->first_index; var <= id->last_index; var++)
            {
              for (i = 0; i < id->idx->variable[var]->values_per_sample; i++)
              {
                bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                memcpy(id->idx->variable[var]->hz_buffer[y]->buffer[c] + ((s * id->idx->variable[var]->values_per_sample + i) * bytes_for_datatype * chunk_size), tupple[cnt].value[var - id->first_index][i], bytes_for_datatype * chunk_size);
                id->idx->variable[id->first_index]->hz_buffer[y]->buffer_index[cnt] = tupple[cnt].index;

                free(tupple[cnt].value[var - id->first_index][i]);
              }
              free(tupple[cnt].value[var - id->first_index]);
            }
            free(tupple[cnt].value);
            cnt++;
          }
        }
        free(tupple);
        tupple = 0;
#endif
      }
    }
    else
    {
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
        index_count = 0;
        total_chunked_patch_size = 0;
        
        for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
        {
          chunked_patch_offset[l] = var0->chunk_patch_group[y]->patch[b]->offset[l] / id->idx->chunk_size[l];
          chunked_patch_size[l] = var0->chunk_patch_group[y]->patch[b]->size[l] / id->idx->chunk_size[l];
          total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
        }
        
        number_levels = maxH - 1;
        PointND xyzuv_Index;
      
        if(var0->data_layout == PIDX_row_major)
        {
          for (v = chunked_patch_offset[4]; v < chunked_patch_offset[4] + chunked_patch_size[4]; v++)
            for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                  {
                    index = (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * chunked_patch_size[3] * (v - chunked_patch_offset[4]))
                          + (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * (u - chunked_patch_offset[3]))
                          + (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                          + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                          + (i - chunked_patch_offset[0]);
                          
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
                      int bit = id->idx->bitPattern[number_levels];
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = maxH - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    
                    level = getLeveL(hz_order);

                    //if (level >= maxH - id->idx_d->resolution_to - 1)
                    //  continue;
                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;
                    hz_index = hz_order - var0->hz_buffer[y]->start_hz_index[level];
                    for(var = id->first_index; var <= id->last_index; var++)
                    {
                      bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                      for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                      {
                        memcpy(id->idx->variable[var]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[var]->values_per_sample + s) * bytes_for_datatype * chunk_size),
                                id->idx->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[var]->values_per_sample) + s) * bytes_for_datatype * chunk_size,
                                bytes_for_datatype * chunk_size);
                        //printf("[%d] [%d %d] [%d %d %d %d %d] [%d %d %d %d %d :: %d %d %d %d %d] %d %d\n", var, b, y, i, j, k, u, v, chunked_patch_offset[0], chunked_patch_offset[1], chunked_patch_offset[2], chunked_patch_offset[3], chunked_patch_offset[4], chunked_patch_size[0], chunked_patch_size[1], chunked_patch_size[2], chunked_patch_size[3], chunked_patch_size[4],  (int)index, (int)hz_index);
                      }
                    }
                  }
        }
        else
        {
          for (v = chunked_patch_offset[4]; v < chunked_patch_offset[4] + chunked_patch_size[4]; v++)
            for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
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
                      int bit = id->idx->bitPattern[number_levels];
                      z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                      PGET(xyzuv_Index, bit) >>= 1;
                    }

                    number_levels = maxH - 1;
                    int64_t lastbitmask = ((int64_t) 1) << number_levels;
                    z_order |= lastbitmask;
                    while (!(1 & z_order)) z_order >>= 1;
                    z_order >>= 1;

                    hz_order = z_order;
                    level = getLeveL(hz_order);
                    
                    index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                          + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                          + (k - chunked_patch_offset[2]);

                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;
                    for(var = id->first_index; var <= id->last_index; var++)
                    {
                      hz_index = hz_order - id->idx->variable[var]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                      for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                      {
                        memcpy(id->idx->variable[var]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[var]->values_per_sample + s) * bytes_for_datatype * chunk_size),
                               id->idx->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[var]->values_per_sample) + s) * bytes_for_datatype * chunk_size,
                               bytes_for_datatype * chunk_size);
                      }
                    }
                  }
        }
      }
    }

    if (var0->chunk_patch_group[y]->type == 2)
    {
      for (i = id->first_index; i <= id->last_index; i++)
      {
        PIDX_variable var = id->idx->variable[i];
        bytes_for_datatype = (var->bits_per_value / 8) * var->values_per_sample;
        for (j = 0; j < maxH; j++)
        {
          int level_blocks = 0;
          if (var0->hz_buffer[y]->missing_block_count_per_level[j] != 0)
          {
            int samples_per_level = var0->hz_buffer[y]->end_hz_index[j] - var0->hz_buffer[y]->start_hz_index[j] + 1;
            int missing_sample_count = var0->hz_buffer[y]->missing_block_count_per_level[j] * id->idx_d->samples_per_block;
            int adjusted_buffer_size = (samples_per_level- missing_sample_count) * bytes_for_datatype;
            
            int initial_m = 0;
            for (n = 0; n < var0->hz_buffer[y]->missing_block_count_per_level[j]; n++)
            {
              for (m = initial_m; m < samples_per_level; m++)
              {
                if (var0->hz_buffer[y]->start_hz_index[j] + m == var0->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
                {
                  memmove(var->hz_buffer[y]->buffer[j] + (m - (level_blocks * id->idx_d->samples_per_block)) * bytes_for_datatype * var->values_per_sample, var->hz_buffer[y]->buffer[j] + ((m + id->idx_d->samples_per_block) - (level_blocks * id->idx_d->samples_per_block)) * bytes_for_datatype * var->values_per_sample, (var0->hz_buffer[y]->end_hz_index[j] - (var0->hz_buffer[y]->start_hz_index[j] + m + id->idx_d->samples_per_block) + 1) * bytes_for_datatype * var->values_per_sample);
                  level_blocks++;
                  initial_m = m;
                }
              }
            }
            
            unsigned char* temp_buffer = realloc(var->hz_buffer[y]->buffer[j], adjusted_buffer_size);
            if (temp_buffer == NULL)
            {
              fprintf(stderr, "realloc error File %s Line no %d\n", __FILE__, __LINE__);
              return PIDX_err_hz;
            }
            else
              var->hz_buffer[y]->buffer[j] = temp_buffer;
          }
        }
      }
    }
  }
  return PIDX_success;
}


PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  int64_t z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, s = 0, y = 0, n = 0, m = 0, number_levels = 0, var = 0;
  int64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0;
  int bytes_for_datatype;
  int64_t hz_index;
  int64_t total_chunked_patch_size;
  
  int64_t total_chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];
  
  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  
  if(id->idx->variable[id->first_index]->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->count not set.\n", __FILE__, __LINE__);
    return -1;
  }
  
  if(id->idx_d->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->maxh not set.\n", __FILE__, __LINE__);
    return -1;
  }
  
  
  for (y = 0; y < id->idx->variable[id->first_index]->patch_group_count; y++)
  {
    for (b = 0; b < id->idx->variable[id->first_index]->chunk_patch_group[y]->count; b++)
    {
      total_chunked_patch_size = 0;
      
      for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
      {
        chunked_patch_offset[l] = id->idx->variable[id->first_index]->chunk_patch_group[y]->patch[b]->offset[l] / id->idx->chunk_size[l];
        chunked_patch_size[l] = id->idx->variable[id->first_index]->chunk_patch_group[y]->patch[b]->size[l] / id->idx->chunk_size[l];
        total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
      }
      
      number_levels = id->idx_d->maxh - 1;
      PointND xyzuv_Index;
    
      if(id->idx->variable[id->first_index]->data_layout == PIDX_row_major)
      {
        for (v = chunked_patch_offset[4]; v < chunked_patch_offset[4] + chunked_patch_size[4]; v++)
          for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
            for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
              for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                {
                  index = (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * chunked_patch_size[3] * (v - chunked_patch_offset[4]))
                        + (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * (u - chunked_patch_offset[3]))
                        + (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                        + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                        + (i - chunked_patch_offset[0]);
                        
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
                    int bit = id->idx->bitPattern[number_levels];
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_d->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  
                  level = getLeveL(hz_order);
                  id->idx->variable[id->first_index]->hz_buffer[y]->samples_per_level[level] = id->idx->variable[id->first_index]->hz_buffer[y]->samples_per_level[level] + 1;
                  hz_index = hz_order - id->idx->variable[id->first_index]->hz_buffer[y]->start_hz_index[level];
                  for(var = id->first_index; var <= id->last_index; var++)
                  {
                    bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                    for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                    {                        
                      memcpy(id->idx->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_chunk_size,
                             id->idx->variable[var]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[var]->values_per_sample + s) * bytes_for_datatype * total_chunk_size),
                             bytes_for_datatype * total_chunk_size);
                    }
                  }
                }
      }
      else
      {
        for (v = chunked_patch_offset[4]; v < chunked_patch_offset[4] + chunked_patch_size[4]; v++)
          for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
            for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
              for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
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
                    int bit = id->idx->bitPattern[number_levels];
                    z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                    PGET(xyzuv_Index, bit) >>= 1;
                  }

                  number_levels = id->idx_d->maxh - 1;
                  int64_t lastbitmask = ((int64_t) 1) << number_levels;
                  z_order |= lastbitmask;
                  while (!(1 & z_order)) z_order >>= 1;
                  z_order >>= 1;

                  hz_order = z_order;
                  level = getLeveL(hz_order);
                  
                  index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                        + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                        + (k - chunked_patch_offset[2]);

                  id->idx->variable[id->first_index]->hz_buffer[y]->samples_per_level[level] = id->idx->variable[id->first_index]->hz_buffer[y]->samples_per_level[level] + 1;
                  for(var = id->first_index; var <= id->last_index; var++)
                  {
                    hz_index = hz_order - id->idx->variable[var]->hz_buffer[y]->start_hz_index[level];
                    bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                    for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                    {
                      bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;
                      memcpy(id->idx->variable[var]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[var]->values_per_sample) + s) * bytes_for_datatype * total_chunk_size, id->idx->variable[var]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[var]->values_per_sample + s) * bytes_for_datatype * total_chunk_size), bytes_for_datatype * total_chunk_size);
                    }
                  }
                }
      }
    }
    
    if (id->idx->variable[id->first_index]->chunk_patch_group[y]->type == 2)
    {
      for (i = id->first_index; i <= id->last_index; i++)
      {
        bytes_for_datatype = id->idx->variable[i]->bits_per_value / 8;
        for (j = 0; j < id->idx_d->maxh; j++)
        {
          int level_blocks = 0;
          if (id->idx->variable[i]->hz_buffer[y]->missing_block_count_per_level[j] != 0)
          {
            int adjusted_buffer_size = (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + 1 - (id->idx->variable[i]->hz_buffer[y]->missing_block_count_per_level[j] * id->idx_d->samples_per_block)) * bytes_for_datatype * id->idx->variable[i]->values_per_sample;
            
            //printf("Initial Size %lld Adjusted size %lld (%lld - %d)\n", (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + 1), adjusted_buffer_size/(bytes_for_datatype * id->idx->variable[i]->values_per_sample), (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + 1), (id->idx->variable[i]->hz_buffer[y]->missing_block_count_per_level[j] * id->idx_d->samples_per_block));
            
            int initial_m = 0;
            for (n = 0; n < id->idx->variable[i]->hz_buffer[y]->missing_block_count_per_level[j]; n++)
            {
              //printf("[%d] [L %d] start %lld end %lld\n", n, j, id->idx->variable[i]->hz_buffer[y]->start_hz_index[j], id->idx->variable[i]->hz_buffer[y]->end_hz_index[j]);
              for (m = initial_m; m < (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + 1); m++)
              {
                if (m + id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] == id->idx->variable[i]->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
                {
                  //printf("[MISS] (%d + %d) %d = %d (%d x %d)\n", m, id->idx->variable[i]->hz_buffer[y]->start_hz_index[j], (m + id->idx->variable[i]->hz_buffer[y]->start_hz_index[j]), id->idx->variable[i]->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block, id->idx->variable[i]->hz_buffer[y]->missing_block_index_per_level[j][n], id->idx_d->samples_per_block);
                  //printf("[SOURCE] %d [DEST] %d [COUNT] %d\n", (m - (level_blocks * id->idx_d->samples_per_block)), ((m + id->idx_d->samples_per_block) - (level_blocks * id->idx_d->samples_per_block)), (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + 1 - ((m + level_blocks) + id->idx_d->samples_per_block) ));
                  
                  memmove(id->idx->variable[i]->hz_buffer[y]->buffer[j] + (m - (level_blocks * id->idx_d->samples_per_block)) * bytes_for_datatype * id->idx->variable[i]->values_per_sample,
                          id->idx->variable[i]->hz_buffer[y]->buffer[j] + ((m + id->idx_d->samples_per_block) - (level_blocks * id->idx_d->samples_per_block)) * bytes_for_datatype * id->idx->variable[i]->values_per_sample,
                          (id->idx->variable[i]->hz_buffer[y]->end_hz_index[j] - (id->idx->variable[i]->hz_buffer[y]->start_hz_index[j] + m + id->idx_d->samples_per_block) + 1) * bytes_for_datatype * id->idx->variable[i]->values_per_sample );
                  level_blocks++;
                  initial_m = m;
                }
              }
            }
            
            unsigned char* temp_buffer = realloc(id->idx->variable[i]->hz_buffer[y]->buffer[j], adjusted_buffer_size);
            if (temp_buffer == NULL) 
            {
              
            }
            else
            {
              id->idx->variable[i]->hz_buffer[y]->buffer[j] = temp_buffer;
            }
            
            //int64_t dvalue;
            //if (j == 20)
            //{
            //  memcpy(&dvalue, id->idx->variable[i]->hz_buffer[y]->buffer[j] + 887364 * bytes_for_datatype, bytes_for_datatype);
            //  printf("[XXX] Value at 887364 = %lld\n", dvalue);
            //}
          }
        }
      }
    }
  }
  return 0;
  

  return 1;
}


/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id)
{
  int itr = 0, p = 0, v = 0;
  
  if(id->idx->variable[id->first_index]->sim_patch_count < 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->sim_patch_count not set.\n", __FILE__, __LINE__);
    return 1;
  }
  if(id->idx_d->maxh <= 0)
  {
    fprintf(stderr, "[%s] [%d] id->idx_d->maxh (%d) not set.\n", __FILE__, __LINE__, id->idx_d->maxh);
    return 1;
  }
  
  PIDX_variable var0 = id->idx->variable[id->first_index];
  for (p = 0; p < var0->patch_group_count; p++)
  {

    if (id->idx->enable_agg == 2)
    {
      free(var0->hz_buffer[p]->samples_per_level);
      var0->hz_buffer[p]->samples_per_level = 0;

      free(var0->hz_buffer[p]->start_hz_index);
      free(var0->hz_buffer[p]->end_hz_index);

      free(var0->hz_buffer[p]->missing_block_count_per_level);

      if (var0->hz_buffer[p]->type == 0)
        free(var0->hz_buffer[p]->buffer_index);

      for (itr = 0; itr < id->idx_d->maxh; itr++)
        free(var0->hz_buffer[p]->nsamples_per_level[itr]);
      free(var0->hz_buffer[p]->nsamples_per_level);
    }

    for (itr = 0; itr < id->idx_d->maxh; itr++)
    {
      free(var0->hz_buffer[p]->missing_block_index_per_level[itr]);
      var0->hz_buffer[p]->missing_block_index_per_level[itr] = 0;
    }
    free(var0->hz_buffer[p]->missing_block_index_per_level);
    var0->hz_buffer[p]->missing_block_index_per_level = 0;



  }
  
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var0->patch_group_count; p++)
    { 
      if (id->idx->enable_agg == 2 || id->idx->enable_agg == 1)
      {
        for (itr = var->hz_buffer[p]->HZ_agg_from; itr < var->hz_buffer[p]->HZ_agg_to; itr++)
        {
          free(var->hz_buffer[p]->buffer[itr]);
          var->hz_buffer[p]->buffer[itr] = 0;
        }
      
        /*
        if (id->idx->enable_agg == 2)
        {
          free(var->hz_buffer[p]->buffer);
          var->hz_buffer[p]->buffer = 0;
      
          free(var->hz_buffer[p]);
          var->hz_buffer[p] = 0;
        }
        */
      }
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_hz_encode_finalize(PIDX_hz_encode_id id)
{
  
  free(id);
  id = 0;
  
  return PIDX_success;
}


PIDX_return_code HELPER_Hz_encode(PIDX_hz_encode_id id)
{
  int i = 0, k = 0, b = 0, var = 0;
  int64_t global_hz, element_count = 0, lost_element_count = 0;
  int64_t ZYX[PIDX_MAX_DIMENSIONS];
  int check_bit = 1, s = 0;
  
#if long_buffer
  uint64_t dvalue_1, dvalue_2;
#else
  double dvalue_1, dvalue_2;
#endif

#if PIDX_HAVE_MPI
  int rank = 0;
  MPI_Comm_rank(id->comm, &rank);
#endif
  
  PIDX_variable var0 = id->idx->variable[id->first_index];
  for(var = id->first_index; var <= id->last_index; var++)
  {
    for (b = 0; b < var0->patch_group_count; b++)
    {
      for (i = 0; i < id->idx_d->maxh; i++)
      {
        //printf("samples at level %d = %d\n", i, id->idx->variable[id->first_index]->hz_buffer[b]->samples_per_level[i]);
        if (var0->hz_buffer[b]->samples_per_level[i] != 0)
        {
          for (k = 0; k <= (var0->hz_buffer[b]->end_hz_index[i] - var0->hz_buffer[b]->start_hz_index[i]) * 1; k++)
          {
            global_hz = var0->hz_buffer[b]->start_hz_index[i] + k;
            
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_hz, ZYX);
            if (!(ZYX[0] >= id->idx->bounds[0] || ZYX[1] >= id->idx->bounds[1] || ZYX[2] >= id->idx->bounds[2]))
            {
              check_bit = 1, s = 0;    
              for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
              {
                dvalue_1 = 100 + var + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
#if long_buffer
                dvalue_2 = *(*((uint64_t**)id->idx->variable[var]->hz_buffer[b]->buffer + i) + ((k * id->idx->variable[var]->values_per_sample) + s));
#else
                dvalue_2 = *(*((double**)id->idx->variable[var]->hz_buffer[b]->buffer + i) + ((k * id->idx->variable[var]->values_per_sample) + s));
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
  
  //printf("[HZ] Volume [%lld] and Volume [%lld]\n", global_volume, (int64_t)(id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * id->idx->bounds[3] * id->idx->bounds[4] * (id->last_index - id->first_index + 1)));
  
  if (global_volume != (int64_t) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1))
  {
    if (rank == 0)
      fprintf(stderr, "[HZ Debug FAILED!!!!] [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", id->idx_d->color, (long long) global_volume, (long long) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1));
    
    printf("[HZ]  Rank %d Color %d [LOST ELEMENT COUNT %lld] [FOUND ELEMENT COUNT %lld] [TOTAL ELEMNTS %lld] \n", rank,  id->idx_d->color, (long long) lost_element_count, (long long) element_count, (long long) (id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * id->idx->bounds[3] * id->idx->bounds[4]) * (id->last_index - id->first_index + 1));
    
    return (-1);
  }
  else
  {
    if (rank == 0)
      fprintf(stderr, "[HZ Debug PASSED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", id->idx_d->color, (long long) global_volume, (long long) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1));
  }
#endif
    
  return PIDX_success;
}
