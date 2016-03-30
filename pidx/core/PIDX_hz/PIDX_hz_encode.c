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

#include "../../PIDX_inc.h"

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

  int resolution_from;
  int resolution_to;
  
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Comm global_comm;
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

PIDX_return_code PIDX_hz_encode_set_global_communicator(PIDX_hz_encode_id hz_id, MPI_Comm comm)
{
  if (hz_id == NULL)
    return PIDX_err_id;

  hz_id->global_comm = comm;

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


PIDX_return_code PIDX_hz_encode_set_resolution(PIDX_hz_encode_id id, int resolution_from, int resolution_to)
{
  if (resolution_from < 0 || resolution_from < 0)
    return PIDX_err_hz;

  id->resolution_from = resolution_from;
  id->resolution_to = resolution_to;

  return PIDX_success;
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
  
  hz_id->resolution_from = 0;
  hz_id->resolution_to = 0;

  return hz_id;
}

PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id)
{
  int **tpatch;
  int **allign_offset;
  int **allign_count;
  int j = 0, p = 0, d = 0, v = 0;
  /*
  int rank, nprocs;
  MPI_Comm_rank(id->global_comm, &rank);
  MPI_Comm_size(id->global_comm, &nprocs);
  */

  int maxH = id->idx_d->maxh;

  tpatch = (int**) malloc(2 * sizeof (int*));
  memset(tpatch, 0, 2 * sizeof (int*));
  tpatch[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  tpatch[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[0], 0, PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[1], 0, PIDX_MAX_DIMENSIONS * sizeof (int));

  if(maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH not set.\n", __FILE__, __LINE__);
    return PIDX_err_hz;
  }

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    var->hz_buffer = malloc(sizeof(*(var->hz_buffer)) * var->patch_group_count);
    memset(var->hz_buffer, 0, sizeof(*(var->hz_buffer)) * var->patch_group_count);

    for (p = 0; p < var->patch_group_count; p++)
    {
      var->hz_buffer[p] = malloc(sizeof(*(var->hz_buffer[p])));
      memset(var->hz_buffer[p], 0, sizeof(*(var->hz_buffer[p])));

      HZ_buffer hz_buf = var->hz_buffer[p];
      if (var->chunk_patch_group[p]->type == 0)
      {
        uint64_t buffer_size = 1;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          buffer_size = buffer_size * (var->chunk_patch_group[p]->patch[0]->size[d]/id->idx->chunk_size[d]);

        hz_buf->buffer_index = malloc(buffer_size * sizeof (int64_t));
        memset(hz_buf->buffer_index, 0, buffer_size * sizeof (int64_t));
      }

      hz_buf->type = var->chunk_patch_group[p]->type;

      hz_buf->start_hz_index = malloc(sizeof (int64_t) * maxH);
      hz_buf->end_hz_index = malloc(sizeof (int64_t) * maxH);
      memset(hz_buf->start_hz_index, 0, sizeof (int64_t) * maxH);
      memset(hz_buf->end_hz_index, 0, sizeof (int64_t) * maxH);

      hz_buf->nsamples_per_level = malloc(sizeof (int*) * maxH);
      memset(hz_buf->nsamples_per_level, 0, sizeof (int*) * maxH);

      if (var->chunk_patch_group[p]->type == 0)
      {
        hz_buf->samples_per_level = malloc( maxH * sizeof (int64_t));
        memset(hz_buf->samples_per_level, 0, maxH * sizeof (int64_t));
      }

      allign_offset = malloc(sizeof (int*) * maxH);
      allign_count = malloc(sizeof (int*) * maxH);
      memset(allign_offset, 0, sizeof (int*) * maxH);
      memset(allign_count, 0, sizeof (int*) * maxH);

      for (j = 0; j < maxH; j++)
      //for (j = id->resolution_from; j < maxH - id->resolution_to; j++)
      {
        allign_offset[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
        memset(allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

        allign_count[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
        memset(allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

        hz_buf->nsamples_per_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
        memset(hz_buf->nsamples_per_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
      }

      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        //assert(var->chunk_patch_group[p]->reg_patch_offset[d] % id->idx->chunk_size[d] == 0);
        tpatch[0][d] = var->chunk_patch_group[p]->reg_patch->offset[d] / id->idx->chunk_size[d];

        if (var->chunk_patch_group[p]->reg_patch->size[d] % id->idx->chunk_size[d] == 0)
          tpatch[1][d] = (var->chunk_patch_group[p]->reg_patch->offset[d] / id->idx->chunk_size[d]) + (var->chunk_patch_group[p]->reg_patch->size[d] / id->idx->chunk_size[d]) - 1;
        else
          tpatch[1][d] = (var->chunk_patch_group[p]->reg_patch->offset[d] / id->idx->chunk_size[d]) + ((var->chunk_patch_group[p]->reg_patch->size[d] / id->idx->chunk_size[d]) + 1) - 1;
      }

      //for (j = 0; j < maxH; j++)
      for (j = id->resolution_from; j < maxH - id->resolution_to; j++)
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
        //if (rank == 1)
        //printf("XX%d [%d (%d): %d %d %d %d %d - %d %d %d %d %d] [%s] SE %d %d T %d %d %d\n", v, j, maxH, tpatch[0][0], tpatch[0][1], tpatch[0][2], tpatch[0][3], tpatch[0][4], tpatch[1][0], tpatch[1][1], tpatch[1][2], tpatch[1][3], tpatch[1][4], id->idx->bitPattern, hz_buf->start_hz_index[j], hz_buf->end_hz_index[j], hz_buf->nsamples_per_level[j][0], hz_buf->nsamples_per_level[j][1], hz_buf->nsamples_per_level[j][2]);
      }


      for (j = 0; j < maxH; j++)
      {
        free(allign_offset[j]);
        free(allign_count[j]);
      }
      free(allign_offset);
      free(allign_count);

    }
  }

  free(tpatch[0]);
  tpatch[0] = 0;
  free(tpatch[1]);
  tpatch[1] = 0;
  free(tpatch);
  tpatch = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_hz_encode_buf_create(PIDX_hz_encode_id id)
{

  int p = 0, c = 0, v = 0, bytes_for_datatype = 0;

  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];

  //int rank = 0;
  //MPI_Comm_rank(id->global_comm, &rank);

  // Allocate actual HZ buffer for the variables
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
#if !SIMULATE_IO
      var->hz_buffer[p]->buffer = (unsigned char**)malloc( maxH * sizeof (unsigned char*));
      memset(var->hz_buffer[p]->buffer, 0,  maxH * sizeof (unsigned char*));
#endif

      bytes_for_datatype = ((var->bits_per_value / 8) * chunk_size * var->values_per_sample) / id->idx->compression_factor;
      if (var->chunk_patch_group[p]->type == 1 || var->chunk_patch_group[p]->type == 2)
      {
        //for (c = 0 ; c < maxH ; c++)
        for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
        {
          int64_t samples_per_level = (var->hz_buffer[p]->end_hz_index[c] - var->hz_buffer[p]->start_hz_index[c] + 1);

#if !SIMULATE_IO
          var->hz_buffer[p]->buffer[c] = malloc(bytes_for_datatype * samples_per_level);
          memset(var->hz_buffer[p]->buffer[c], 0, bytes_for_datatype * samples_per_level);
          //if (rank == 1)
          //printf("[%d] Buffer size on level %d = %d %d\n", v, c, samples_per_level, bytes_for_datatype);
#endif
        }
      }
    }
  }

  return PIDX_success;
}

PIDX_return_code PIDX_hz_encode_write(PIDX_hz_encode_id id)
{
  uint64_t z_order = 0, hz_order = 0, index = 0;
  int b = 0, level = 0, cnt = 0, c = 0, s = 0, y = 0, m = 0, number_levels = 0;
  uint64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0, d = 0;
  int v1 = 0;
  int index_count = 0;
  int bytes_for_datatype;
  uint64_t hz_index;
  uint64_t total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];
  PIDX_variable var0 = id->idx->variable[id->first_index];
  
  //int rank = 0;
  //MPI_Comm_rank(id->global_comm, &rank);

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
        memset(tupple, 0, total_chunked_patch_size * sizeof(hz_tupple));
      
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

                    //if (level >= maxH - id->resolution_to)
                    //if (level >= maxH)
                    //  continue;

                    index = (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * chunked_patch_size[3] * (m - chunked_patch_offset[4])) +
                            (chunked_patch_size[0] * chunked_patch_size[1] * chunked_patch_size[2] * (u - chunked_patch_offset[3])) +
                            (chunked_patch_size[0] * chunked_patch_size[1] * (k - chunked_patch_offset[2]))
                        + (chunked_patch_size[0] * (j - chunked_patch_offset[1]))
                        + (i - chunked_patch_offset[0]);
                  
                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    for(v = id->first_index; v <= id->last_index; v++)
                    {
                      PIDX_variable var = id->idx->variable[v];

                      if (level < maxH - id->resolution_to)
                        var->hz_buffer[y]->samples_per_level[level] = var->hz_buffer[y]->samples_per_level[level] + 1;

                      tupple[index_count].value[v - id->first_index] = malloc(var->values_per_sample * sizeof(unsigned char*));
                      memset(tupple[index_count].value[v - id->first_index], 0, var->values_per_sample * sizeof(unsigned char*));

                      //bytes_for_datatype = var->bits_per_value / 8;
                      bytes_for_datatype = ((var->bits_per_value / 8) * chunk_size) / id->idx->compression_factor;
                      for (s = 0; s < var->values_per_sample; s++)
                      {
                        tupple[index_count].value[v - id->first_index][s] = malloc(bytes_for_datatype);
                        memset(tupple[index_count].value[v - id->first_index][s], 0, bytes_for_datatype);

                        memcpy(tupple[index_count].value[v- id->first_index][s], var->chunk_patch_group[y]->patch[b]->buffer + ((index * var->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);

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

                    //if (level >= maxH - id->resolution_to)
                    //if (level >= maxH)
                    //  continue;
                                  
                    index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                        + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                        + (k - chunked_patch_offset[2]);
                  
                    tupple[index_count].value = malloc((id->last_index - id->first_index + 1)* sizeof(unsigned char**));
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      if (level < maxH - id->resolution_to)
                        id->idx->variable[v1]->hz_buffer[y]->samples_per_level[level] = id->idx->variable[v1]->hz_buffer[y]->samples_per_level[level] + 1;

                      tupple[index_count].value[v1 - id->first_index] = malloc(id->idx->variable[v1]->values_per_sample * sizeof(unsigned char*));
                      memset(tupple[index_count].value[v1 - id->first_index], 0, id->idx->variable[v1]->values_per_sample * sizeof(unsigned char*));

                      //bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
                      bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size) / id->idx->compression_factor;
                      for (s = 0; s < id->idx->variable[v1]->values_per_sample; s++)
                      {
                        tupple[index_count].value[v1 - id->first_index][s] = malloc(bytes_for_datatype);
                        memset(tupple[index_count].value[v1 - id->first_index][s], 0, bytes_for_datatype);

                        memcpy(tupple[index_count].value[v1 - id->first_index][s], id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[v1]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);

                        tupple[index_count].index = hz_order;
                      }
                    }
                    index_count++;
                  }
        }
        qsort( tupple, total_chunked_patch_size, sizeof(hz_tupple), compare );
      
        for(v1 = id->first_index; v1 <= id->last_index; v1++)
        {
          for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
          //for(c = 0 ; c < maxH ; c++)
          {
            //bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
            bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size) / id->idx->compression_factor;

            id->idx->variable[v1]->hz_buffer[y]->buffer[c] = malloc(bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * id->idx->variable[v1]->values_per_sample);
            memset(id->idx->variable[v1]->hz_buffer[y]->buffer[c], 0, bytes_for_datatype * var0->hz_buffer[y]->samples_per_level[c] * id->idx->variable[v1]->values_per_sample);
          }
        }
      

        cnt = 0;
        for (c = id->resolution_from; c < maxH - id->resolution_to; c++)
        //for (c = 0; c < maxH; c++)
        {
          for (s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
          {
            for (v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              for (i = 0; i < id->idx->variable[v1]->values_per_sample; i++)
              {
                //bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
                bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size) / id->idx->compression_factor;
                memcpy(id->idx->variable[v1]->hz_buffer[y]->buffer[c] + ((s * id->idx->variable[v1]->values_per_sample + i) * bytes_for_datatype), tupple[cnt].value[v1 - id->first_index][i], bytes_for_datatype);
                id->idx->variable[/*id->first_index*/v1]->hz_buffer[y]->buffer_index[cnt] = tupple[cnt].index;
              }
            }
            cnt++;
          }
        }

        /*
        cnt = 0;
        for (c = 0; c < maxH; c++)
        {
          for (s = 0; s < var0->hz_buffer[y]->samples_per_level[c]; s++)
          {
            for (v1 = id->first_index; v1 <= id->last_index; v1++)
            {
              for (i = 0; i < id->idx->variable[v1]->values_per_sample; i++)
              {
                free(tupple[cnt].value[v1 - id->first_index][i]);
              }
              free(tupple[cnt].value[v1 - id->first_index]);
            }
            free(tupple[cnt].value);
            cnt++;
          }
        }
        free(tupple);
        tupple = 0;
        */

        if (var0->data_layout == PIDX_row_major)
        {
          index_count = 0;
          for (m = chunked_patch_offset[4]; m < chunked_patch_offset[4] + chunked_patch_size[4]; m++)
            for (u = chunked_patch_offset[3]; u < chunked_patch_offset[3] + chunked_patch_size[3]; u++)
              for (k = chunked_patch_offset[2]; k < chunked_patch_offset[2] + chunked_patch_size[2]; k++)
                for (j = chunked_patch_offset[1]; j < chunked_patch_offset[1] + chunked_patch_size[1]; j++)
                  for (i = chunked_patch_offset[0]; i < chunked_patch_offset[0] + chunked_patch_size[0]; i++)
                  {
                    for(v = id->first_index; v <= id->last_index; v++)
                    {
                      PIDX_variable var = id->idx->variable[v];
                      for (s = 0; s < var->values_per_sample; s++)
                      {
                        free(tupple[index_count].value[v - id->first_index][s]);
                      }
                      free(tupple[index_count].value[v - id->first_index]);
                    }
                    free(tupple[index_count].value);
                    index_count++;
                  }
        }
        free(tupple);



#endif
      }
    }
    else
    {
#if !SIMULATE_IO
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
                    //printf("%lld %lld -> %lld %d\n", (long long)i, (long long)j, (long long)hz_order, level);

                    //if (level >= maxH)
                    if (level >= maxH - id->resolution_to)
                      continue;

                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      hz_index = hz_order - id->idx->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size * id->idx->variable[v1]->values_per_sample) / id->idx->compression_factor;
                      //printf("bytes_for_datatype = %d\n", bytes_for_datatype);
                      //for (s = 0; s < id->idx->variable[v1]->values_per_sample; s++)
                      //{
                      /*
                      int nprocs;
                      MPI_Comm_size(id->comm, &nprocs);
                        double x;
                        memcpy(&x, id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), bytes_for_datatype);
                        if (nprocs == 1)
                        printf("Dest %d %d Src %d Count %d byte size %d V %f\n", level, (int)((hz_index * id->idx->variable[v1]->values_per_sample + s) * chunk_size), (int)((index * id->idx->variable[v1]->values_per_sample) + s) * chunk_size, (int)chunk_size, id->idx->variable[v1]->bits_per_value, x);
                        */

#if !SIMULATE_IO
                        //memcpy(id->idx->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[v1]->values_per_sample + s) * bytes_for_datatype), id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[v1]->values_per_sample) + s) * bytes_for_datatype, bytes_for_datatype);

                      //if (rank == 1)
                      //{
                      memcpy(id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype),
                              id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                              bytes_for_datatype);
                      //}
                      /*
                      double x;
                      memcpy(&x, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);
                      double y1;
                      memcpy(&y1, id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), bytes_for_datatype);
                      if (rank == 1)
                      printf("R %d ----> HZ %d (%d) [V%d Y%d B%d] [%f %f] [Bytes %d]\n", index, hz_index, level, v1, y, b, x, y1, bytes_for_datatype);
                      */
#endif
                      //}
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

                    if (level >= maxH - id->resolution_to)
                    //if (level >= maxH)
                      continue;
                    
                    index = (chunked_patch_size[2] * chunked_patch_size[1] * (i - chunked_patch_offset[0]))
                          + (chunked_patch_size[2] * (j - chunked_patch_offset[1]))
                          + (k - chunked_patch_offset[2]);

                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      hz_index = hz_order - id->idx->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      //bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
                      bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size) / (id->idx->variable[v1]->bits_per_value / id->idx->compression_bit_rate);
                      for (s = 0; s < id->idx->variable[v1]->values_per_sample; s++)
                      {
#if !SIMULATE_IO
                        memcpy(id->idx->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[v1]->values_per_sample + s) * bytes_for_datatype),
                               id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[v1]->values_per_sample) + s) * bytes_for_datatype,
                               bytes_for_datatype);
#endif
                      }
                    }
                  }
        }
      }
#endif
    }

    /*
    int rank = 0;
    MPI_Comm_rank(id->comm, &rank);
    int missing_sample_count = 0;
    if (var0->chunk_patch_group[y]->type == 2)
    {
      for (i = id->first_index; i <= id->last_index; i++)
      {
        PIDX_variable var = id->idx->variable[i];
        bytes_for_datatype = (var->bits_per_value / 8) * var->values_per_sample;
        for (j = id->idx_d->res_from; j < maxH - id->idx_d->res_to; j++)
        {
          int skipped_blocks = 0;
          int skipped_samples = 0;
          if (var0->hz_buffer[y]->missing_block_count_per_level[j] != 0)
          {
            int samples_per_level = var0->hz_buffer[y]->end_hz_index[j] - var0->hz_buffer[y]->start_hz_index[j] + 1;
            missing_sample_count = var0->hz_buffer[y]->missing_block_count_per_level[j] * id->idx_d->samples_per_block;
            int adjusted_buffer_size = (samples_per_level- missing_sample_count) * bytes_for_datatype;

            int sample_start_index = 0;
            for (n = 0; n < var0->hz_buffer[y]->missing_block_count_per_level[j]; n++)
            {
              int sample = 0;
              for (sample = sample_start_index; sample < samples_per_level; sample++)
              {
                if (var0->hz_buffer[y]->start_hz_index[j] + sample == var0->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
                {
                  skipped_samples = skipped_blocks * id->idx_d->samples_per_block;
                  int dest = (sample - skipped_samples) * bytes_for_datatype;
                  int src = (sample + id->idx_d->samples_per_block - skipped_samples) * bytes_for_datatype;
                  int count = (samples_per_level - sample - id->idx_d->samples_per_block) * bytes_for_datatype;
                  //printf("[%d] src dest count %d %d %d (%d - %d - %d)\n", rank, src, dest, count, samples_per_level, sample, id->idx_d->samples_per_block);
                  if (count < 0)
                    break;

#if !SIMULATE_IO
                  memmove(var->hz_buffer[y]->buffer[j] + dest, var->hz_buffer[y]->buffer[j] + src, count);
#endif
                  skipped_blocks++;
                  sample_start_index = sample;
                  break;
                }
                sample++;
              }
            }
#if !SIMULATE_IO
            //printf("adjusted_buffer_size = %d\n", adjusted_buffer_size);
            unsigned char* temp_buffer = realloc(var->hz_buffer[y]->buffer[j], adjusted_buffer_size);
            if (temp_buffer == NULL)
            {
              fprintf(stderr, "realloc error File %s Line no %d\n", __FILE__, __LINE__);
              return PIDX_err_hz;
            }
            else
              var->hz_buffer[y]->buffer[j] = temp_buffer;
#endif
          }
        }
      }
    }
    */
  }
  return PIDX_success;
}




PIDX_return_code PIDX_hz_encode_write_inverse(PIDX_hz_encode_id id, int start_hz_index, int end_hz_index)
{
  int rank, nprocs;
  uint64_t index = 0;
  int b = 0, y = 0, m = 0;
  uint64_t j = 0, l = 0;
  int v1 = 0;
  int bytes_for_datatype;
  uint64_t hz_index;
  uint64_t total_chunked_patch_size = 1;
  int maxH = id->idx_d->maxh;
  int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2] * id->idx->chunk_size[3] * id->idx->chunk_size[4];
  PIDX_variable var0 = id->idx->variable[id->first_index];

  int chunked_patch_offset[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int chunked_patch_size[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};

  MPI_Comm_rank(id->global_comm, &rank);
  MPI_Comm_size(id->global_comm, &nprocs);

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


  //
  for (y = 0; y < var0->patch_group_count; y++)
  {
#if !SIMULATE_IO
    for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
    {
      total_chunked_patch_size = 0;

      for (l = 0; l < PIDX_MAX_DIMENSIONS; l++)
      {
        chunked_patch_offset[l] = var0->chunk_patch_group[y]->patch[b]->offset[l] / id->idx->chunk_size[l];
        chunked_patch_size[l] = var0->chunk_patch_group[y]->patch[b]->size[l] / id->idx->chunk_size[l];
        total_chunked_patch_size = total_chunked_patch_size * chunked_patch_size[l];
      }

      if(var0->data_layout == PIDX_row_major)
      {
        uint64_t hz_mins = 0, hz_maxes = 0, mindex = 0;
        //for (j = 0; j < id->idx_d->maxh - id->resolution_to; j++)
        for (j = start_hz_index; j < end_hz_index; j++)
        {

          //if (rank == 1)
          //  printf("At level %d\n", j);
          if (id->idx->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][0] * id->idx->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][1] * id->idx->variable[id->first_index]->hz_buffer[y]->nsamples_per_level[j][2] != 0)
          {
            hz_mins = id->idx->variable[id->first_index]->hz_buffer[y]->start_hz_index[j];// out_buf_array[buffer_count]->allign_start_hz[j];
            hz_maxes = id->idx->variable[id->first_index]->hz_buffer[y]->end_hz_index[j] + 1;// out_buf_array[buffer_count]->allign_end_hz[j] + 1;

            //if (rank == 1)
            //  printf("[%d (%d - %d)] ----> %d %d\n", j, maxH, id->resolution_to, hz_mins, hz_maxes);

            for (m = hz_mins; m < hz_maxes; m++)
            {
              mindex = m;
              maxH = id->idx_d->maxh - 1;
              int64_t lastbitmask=((int64_t)1)<<maxH;

              mindex <<= 1;
              mindex  |= 1;
              while ((lastbitmask & mindex) == 0) mindex <<= 1;
                mindex &= lastbitmask - 1;

              PointND cnt;
              PointND p  ;
              int n = 0;

              memset(&cnt,0,sizeof(PointND));
              memset(&p  ,0,sizeof(PointND));

              for (;mindex; mindex >>= 1,++n, maxH--)
              {
                int bit= id->idx->bitPattern[maxH];
                PGET(p,bit) |= (mindex & 1) << PGET(cnt,bit);
                ++PGET(cnt,bit);
              }

              if (p.x >= id->idx->bounds[0] || p.y >= id->idx->bounds[1] || p.z >= id->idx->bounds[2])
              {
                //printf("A [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              if (p.x < chunked_patch_offset[0] || p.y < chunked_patch_offset[1] || p.z < chunked_patch_offset[2])
              {
                //printf("B [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              if (p.x >= chunked_patch_offset[0] + chunked_patch_size[0] || p.y >= chunked_patch_offset[1] + chunked_patch_size[1] || p.z >= chunked_patch_offset[2] + chunked_patch_size[2])
              {
                //printf("C [%d %d %d]\n", p.x, p.y, p.z);
                continue;
              }

              hz_index = m - hz_mins;
              index = (chunked_patch_size[0] * chunked_patch_size[1] * (p.z - chunked_patch_offset[2])) + (chunked_patch_size[0] * (p.y - chunked_patch_offset[1])) + (p.x - chunked_patch_offset[0]);

#if 1
              for(v1 = id->first_index; v1 <= id->last_index; v1++)
              {
                //hz_index = hz_order - id->idx->variable[v1]->hz_buffer[y]->start_hz_index[level];
                bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size * id->idx->variable[v1]->values_per_sample) / id->idx->compression_factor;
#if !SIMULATE_IO
                memcpy(id->idx->variable[v1]->hz_buffer[y]->buffer[j] + (hz_index * bytes_for_datatype),
                        id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype),
                        bytes_for_datatype);

                /*
                if (rank == 1)
                {
                  double xx;
                  memcpy(&xx, id->idx->variable[v1]->hz_buffer[y]->buffer[j] + (hz_index * bytes_for_datatype), sizeof(double));
                  printf("[%d] value %d %d - %f\n", j, index, hz_index, xx);
                }
                */
#endif
              }
#endif

            }
          }

        }
      }
    }
#endif
  }
  //
  return PIDX_success;
}


#if 0
int PIDX_hz_encode_read(PIDX_hz_encode_id id, PIDX_Ndim_buffer** in_buf, PIDX_HZ_buffer** out_buf_array, int num_regular_blocks)
{
    long long hz_index = 0, index = 0;
    long long mins[PIDX_MAX_DIMENSIONS], maxes[PIDX_MAX_DIMENSIONS], l_x, l_y, l_z, l_u, l_v, hzaddress, hz_mins, hz_maxes;
    int i =0, j = 0, k = 0, d = 0, s = 0, buffer_count = 0, rank = 0, y = 0, spv, number_levels = 0, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    long long lastbitmask = ((long long) 1) << id->maxh - 1;
    long long n = 0, m = 0;


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
    //printf("Buffer Count %d\n", buffer_count);
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
            hz_mins = out_buf_array[buffer_count]->allign_start_hz[j];
            hz_maxes = out_buf_array[buffer_count]->allign_end_hz[j] + 1;

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
            hz_mins = out_buf_array[buffer_count]->allign_start_hz[j];
            hz_maxes = out_buf_array[buffer_count]->allign_end_hz[j] + 1;

            for(k = 0; k < hz_maxes-hz_mins ; k++)
            index_cache[i][j][k] = -1;//(int*)malloc(sizeof(int) * (hz_maxes - hz_mins));
        }
        }
    }
    }

    for (y = 0; y < id->number_of_buffers; y++)
    {
    if (num_regular_blocks > 1)
        if (y == id->regular_block_index[buffer_count + 1])
          buffer_count++;

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
        //printf("[%d] Samples at Level [%d] = %d x %d x %d x %d x %d\n", rank, j, out_buf_array[buffer_count]->n_samples[j][0], out_buf_array[buffer_count]->n_samples[j][1], out_buf_array[buffer_count]->n_samples[j][2], out_buf_array[buffer_count]->n_samples[j][3], out_buf_array[buffer_count]->n_samples[j][4]);
        if ((out_buf_array[buffer_count]->n_samples[j][0] * out_buf_array[buffer_count]->n_samples[j][1] * out_buf_array[buffer_count]->n_samples[j][2] * out_buf_array[buffer_count]->n_samples[j][3] * out_buf_array[buffer_count]->n_samples[j][4]) != 0)
        {
        hz_mins = out_buf_array[buffer_count]->allign_start_hz[j];
        hz_maxes = out_buf_array[buffer_count]->allign_end_hz[j] + 1;

        //printf("[%d] Range at %d : %lld %lld : [%d %d %d %d %d]\n", rank, j, hz_mins, hz_maxes, l_x, l_y, l_z, l_u, l_v);
        number_levels = id->maxh - 1;

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

            for (; hzaddress; hzaddress >>= 1, ++n, number_levels--) {
              int bit = id->bitPattern[number_levels];
              PGET(p, bit) |= (hzaddress & 1) << PGET(cnt, bit);
              ++PGET(cnt, bit);
            }
            number_levels = id->maxh - 1;
            //printf("HZ [%lld] : [%lld %lld %lld] : %f\n", m, p.x, p.y, p.z, out_buf_array[buffer_count]->buffer[j][m - hz_mins]);
            //index = (l_x * l_y * (k - in_buf[y]->lower_bounds[2])) + (l_x * (j - in_buf[y]->lower_bounds[1])) + (i - in_buf[y]->lower_bounds[0]);

            if (p.x >= id->global_extents[0] || p.y >= id->global_extents[1] || p.z >= id->global_extents[2] || p.u >= id->global_extents[3] || p.v >= id->global_extents[4]){
              //printf("Case 1\n");
              continue;
            }
            if (p.x < in_buf[y]->lower_bounds[0] || p.y < in_buf[y]->lower_bounds[1] || p.z < in_buf[y]->lower_bounds[2] || p.u < in_buf[y]->lower_bounds[3] || p.v < in_buf[y]->lower_bounds[4]) {
              //printf("Case 2\n");
              continue;
            }
            if (p.x >= in_buf[y]->upper_bounds[0] || p.y >= in_buf[y]->upper_bounds[1] || p.z >= in_buf[y]->upper_bounds[2] || p.u >= in_buf[y]->upper_bounds[3] || p.v >= in_buf[y]->upper_bounds[4]) {
              //printf("Case 3\n");
              continue;
            }

            hz_index = m - hz_mins;
            index = (l_x * l_y * (p.z - in_buf[y]->lower_bounds[2])) + (l_x * (p.y - in_buf[y]->lower_bounds[1])) + (p.x - in_buf[y]->lower_bounds[0]);

            if(m-hz_mins >= 0 && m-hz_mins < (hz_maxes - hz_mins))
            {
                index_cache[y][j][m-hz_mins] = index;
            }

            //if(rank == 0)
            //printf("R [%d %d] [%d] [%d] [%lld] [%d %d %d %d %d] : HZ Index %lld : Index %lld : [L : %d %d %d] [B : %d %d %d]\n", y, buffer_count, rank, j, m, p.x, p.y, p.z, p.u, p.v, hz_index, index, l_x, l_y, l_z, in_buf[y]->lower_bounds[0], in_buf[y]->lower_bounds[1], in_buf[y]->lower_bounds[2]);
            spv = out_buf_array[buffer_count]->sample_per_variable;
            for (s = 0; s < spv; s++) {
              in_buf[y]->buffer[(index * spv) + s] = out_buf_array[buffer_count]->buffer[j][(hz_index * spv) + s];
            }

            if(j == (id->maxh-id->resh -1) && m == (hz_maxes - 1) && y == (id->number_of_buffers - 1))
                var_cache = 0;
            }
            else
            {
            if(index_cache[y][j][m-hz_mins] == -1)
                continue;
            hz_index = m - hz_mins;populate_idx_dataset
            spv = out_buf_array[buffer_count]->sample_per_variable;
            memcpy(in_buf[y]->buffer + (index_cache[y][j][m-hz_mins] * spv), out_buf_array[buffer_count]->buffer[j] + (hz_index * spv), sizeof(double)*spv );
            /*
            for (s = 0; s < spv; s++)
            {
                //TODO: Endian
                in_buf[y]->buffer[(index_cache[y][j][m-hz_mins] * spv) + s] = out_buf_array[buffer_count]->buffer[j][(hz_index * spv) + s];
            }
            */
            }
        }
        }
    }
    }
    //TODO: return ?
    return 1;
}
#endif

PIDX_return_code PIDX_hz_encode_read(PIDX_hz_encode_id id)
{
  int64_t z_order = 0, hz_order = 0, index = 0;
  //int n = 0, m = 0, d = 0, c = 0;
  //int index_count = 0;
  int b = 0, level = 0, cnt = 0, s = 0, y = 0, number_levels = 0;
  int64_t i = 0, j = 0, k = 0, u = 0, v = 0, l = 0;
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
    /*
    if (var0->chunk_patch_group[y]->type == 2)
    {
      for (i = id->first_index; i <= id->last_index; i++)
      {
        PIDX_variable var = id->idx->variable[i];
        bytes_for_datatype = (var->bits_per_value / 8) * var->values_per_sample;
        for (j = 0; j < maxH; j++)
        {
          if (var0->hz_buffer[y]->missing_block_count_per_level[j] != 0)
          {
            int samples_per_level = var0->hz_buffer[y]->end_hz_index[j] - var0->hz_buffer[y]->start_hz_index[j] + 1;

            int sample_start_index = 0;
            for (n = 0; n < var0->hz_buffer[y]->missing_block_count_per_level[j]; n++)
            {
              int sample = 0;
              for (sample = sample_start_index; sample < samples_per_level; sample++)
              {
                if (var0->hz_buffer[y]->start_hz_index[j] + sample == var0->hz_buffer[y]->missing_block_index_per_level[j][n] * id->idx_d->samples_per_block)
                {
                  int src = (sample) * bytes_for_datatype;
                  int dest = (sample + id->idx_d->samples_per_block) * bytes_for_datatype;
                  int count = (samples_per_level - sample - id->idx_d->samples_per_block) * bytes_for_datatype;
                  //printf("Dest %d Source %d Count %d (%d)\n", src, dest, count, (var0->hz_buffer[y]->end_hz_index[j] - (var0->hz_buffer[y]->start_hz_index[j] + sample + id->idx_d->samples_per_block) + 1) * bytes_for_datatype);

                  memmove(var->hz_buffer[y]->buffer[j] + dest, var->hz_buffer[y]->buffer[j] + src, count);
                  //memset(var->hz_buffer[y]->buffer[j] + dest, 0, count);
                  sample_start_index = sample;
                  break;
                }
                sample++;
              }
            }
          }
        }
      }
    }
    */

    if (var0->chunk_patch_group[y]->type == 0)
    {
#if 0
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
        memset(tupple, 0, total_chunked_patch_size * sizeof(hz_tupple));

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
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;

                    for(v = id->first_index; v <= id->last_index; v++)
                    {
                      PIDX_variable var = id->idx->variable[v];
                      tupple[index_count].value[v - id->first_index] = malloc(var->values_per_sample * sizeof(unsigned char*));
                      memset(tupple[index_count].value[v - id->first_index], 0, var->values_per_sample * sizeof(unsigned char*));

                      bytes_for_datatype = var->bits_per_value / 8;

                      for (s = 0; s < var->values_per_sample; s++)
                      {
                        tupple[index_count].value[v - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memset(tupple[index_count].value[v - id->first_index][s], 0, bytes_for_datatype * chunk_size);

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
                    memset(tupple[index_count].value, 0, (id->last_index - id->first_index + 1)* sizeof(unsigned char**));

                    var0->hz_buffer[y]->samples_per_level[level] = var0->hz_buffer[y]->samples_per_level[level] + 1;

                    for(var = id->first_index; var <= id->last_index; var++)
                    {
                      tupple[index_count].value[var - id->first_index] = malloc(id->idx->variable[var]->values_per_sample * sizeof(unsigned char*));
                      memset(tupple[index_count].value[var - id->first_index], 0, id->idx->variable[var]->values_per_sample * sizeof(unsigned char*));

                      bytes_for_datatype = id->idx->variable[var]->bits_per_value / 8;

                      for (s = 0; s < id->idx->variable[var]->values_per_sample; s++)
                      {
                        tupple[index_count].value[var - id->first_index][s] = malloc(bytes_for_datatype * chunk_size);
                        memset(tupple[index_count].value[var - id->first_index][s], 0, bytes_for_datatype * chunk_size);

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
#endif
    }
    else
    {
      for (b = 0; b < var0->chunk_patch_group[y]->count; b++)
      {
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

                    hz_index = hz_order - var0->hz_buffer[y]->start_hz_index[level];
                    int v1;
                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      hz_index = hz_order - id->idx->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = ((id->idx->variable[v1]->bits_per_value / 8) * chunk_size * id->idx->variable[v1]->values_per_sample) / id->idx->compression_factor;


  #if !SIMULATE_IO
                      /*
                      double x1, x2, x3, x4;
                      memcpy(&x1, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), sizeof(double));
                      memcpy(&x2, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + sizeof(double), sizeof(double));
                      memcpy(&x3, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + 2*sizeof(double), sizeof(double));
                      memcpy(&x4, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype) + 3*sizeof(double), sizeof(double));
                      printf("XXXXXXXXXXXXx x1 = %f %f %f %f\n", x1, x2, x3, x4);
                      */
                      memcpy(id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + (index * bytes_for_datatype), id->idx->variable[v1]->hz_buffer[y]->buffer[level] + (hz_index * bytes_for_datatype), bytes_for_datatype);
  #endif
                      /*
                      bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
                      for (s = 0; s < id->idx->variable[v1]->values_per_sample; s++)
                      {
                        memcpy(id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[v1]->values_per_sample) + s) * bytes_for_datatype * chunk_size, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[v1]->values_per_sample + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);
                      }
                      */
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

                    int v1 = 0;
                    for(v1 = id->first_index; v1 <= id->last_index; v1++)
                    {
                      hz_index = hz_order - id->idx->variable[v1]->hz_buffer[y]->start_hz_index[level];
                      bytes_for_datatype = id->idx->variable[v1]->bits_per_value / 8;
                      for (s = 0; s < id->idx->variable[v1]->values_per_sample; s++)
                      {
                        memcpy(id->idx->variable[v1]->chunk_patch_group[y]->patch[b]->buffer + ((index * id->idx->variable[v1]->values_per_sample) + s) * bytes_for_datatype * chunk_size, id->idx->variable[v1]->hz_buffer[y]->buffer[level] + ((hz_index * id->idx->variable[v1]->values_per_sample + s) * bytes_for_datatype * chunk_size), bytes_for_datatype * chunk_size);
                      }
                    }
                  }
        }
      }
    }
  }

  return PIDX_success;
}



/* tear down the various buffer structs. In the case of the output structs this function should also free the memory buffers as well */
PIDX_return_code PIDX_hz_encode_buf_destroy(PIDX_hz_encode_id id)
{
#if !SIMULATE_IO
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

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    { 
      for (itr = id->resolution_from; itr < id->idx_d->maxh - id->resolution_to; itr++)
      {
        free(var->hz_buffer[p]->buffer[itr]);
        var->hz_buffer[p]->buffer[itr] = 0;
      }
    }
  }
#endif
  return PIDX_success;
}

PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id)
{
  int itr = 0, p = 0, v = 0;
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      if (var->chunk_patch_group[p]->type == 0)
      {
        free(var->hz_buffer[p]->samples_per_level);
        var->hz_buffer[p]->samples_per_level = 0;
      }
      free(var->hz_buffer[p]->start_hz_index);
      free(var->hz_buffer[p]->end_hz_index);

      if (var->hz_buffer[p]->type == 0)
        free(var->hz_buffer[p]->buffer_index);

      for (itr = 0; itr < id->idx_d->maxh; itr++)
        free(var->hz_buffer[p]->nsamples_per_level[itr]);
      free(var->hz_buffer[p]->nsamples_per_level);

      free(var->hz_buffer[p]->buffer);
      var->hz_buffer[p]->buffer = 0;

      free(var->hz_buffer[p]);
      var->hz_buffer[p] = 0;
    }

    free(var->hz_buffer);
    var->hz_buffer = 0;
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
#if !SIMULATE_IO
  int i = 0, k = 0, b = 0, v = 0;
  int64_t global_hz, element_count = 0, lost_element_count = 0;
  int64_t ZYX[PIDX_MAX_DIMENSIONS];
  int check_bit = 1, s = 0;
  double dvalue_1, dvalue_2;
  float fvalue_1, fvalue_2;
  uint64_t uvalue_1, uvalue_2;

#if PIDX_HAVE_MPI
  int rank = 0;
  if (id->idx_d->parallel_mode == 1)
    MPI_Comm_rank(id->comm, &rank);
#endif
  
  for(v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];
    //printf("[%d %d] var->type_name = %s\n", v, (id->last_index - id->first_index + 1), var->type_name);
    for (b = 0; b < var->patch_group_count; b++)
    {
      for (i = 0; i < id->idx_d->maxh; i++)
      {
        if (var->hz_buffer[b]->nsamples_per_level[i][0] * var->hz_buffer[b]->nsamples_per_level[i][1] * var->hz_buffer[b]->nsamples_per_level[i][2] != 0)
        {
          //printf("[%d %d] - %d\n", b, i, (var->hz_buffer[b]->end_hz_index[i] - var->hz_buffer[b]->start_hz_index[i]));
          for (k = 0; k <= (var->hz_buffer[b]->end_hz_index[i] - var->hz_buffer[b]->start_hz_index[i]) * 1; k++)
          {
            global_hz = var->hz_buffer[b]->start_hz_index[i] + k;
            
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_hz, ZYX);
            if (!(ZYX[0] >= id->idx->bounds[0] || ZYX[1] >= id->idx->bounds[1] || ZYX[2] >= id->idx->bounds[2]))
            {
              check_bit = 1, s = 0;    
              if (strcmp(var->type_name, FLOAT64) == 0)
              {
                dvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
                dvalue_2 = *(*((double**)var->hz_buffer[b]->buffer + i) + ((k * var->values_per_sample) + s));

                //printf("X[%lld]  %f %f\n", (long long)global_hz, dvalue_1, dvalue_2);
                check_bit = check_bit && (dvalue_1  == dvalue_2);
              }
              else if (strcmp(var->type_name, FLOAT32) == 0)
              {
                fvalue_1 = (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
                fvalue_2 = *(*((float**)var->hz_buffer[b]->buffer + i) + ((k * var->values_per_sample) + s));

                check_bit = check_bit && (fvalue_1  == fvalue_2);
              }
              else if (strcmp(var->type_name, FLOAT64_RGB) == 0)
              {
                for (s = 0; s < 3; s++)
                {
                  dvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);

                  memcpy(&dvalue_2, var->hz_buffer[b]->buffer[i] + ((k * 3) + s) * sizeof(double), sizeof(double));
                  //printf("Y[%d] %f -- %f\n", v, dvalue_1, dvalue_2);

                  check_bit = check_bit && (dvalue_1  == dvalue_2);
                }
              }
              if (strcmp(var->type_name, UINT64) == 0)
              {
                uvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
                uvalue_2 = *(*((unsigned long long**)var->hz_buffer[b]->buffer + i) + ((k * var->values_per_sample) + s));
                check_bit = check_bit && (uvalue_1  == uvalue_2);
              }
              else if (strcmp(var->type_name, UINT64_RGB) == 0)
              {
                for (s = 0; s < 3; s++)
                {
                  uvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_d->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
                  memcpy(&uvalue_2, var->hz_buffer[b]->buffer[i] + ((k * 3) + s) * sizeof(double), sizeof(double));

                  check_bit = check_bit && (uvalue_1  == uvalue_2);
                }
              }

              if (check_bit == 0)
              {
                //printf("%lld [HZ] %f %f (%lld :: %lld %lld %lld)\n", (long long)global_hz, fvalue_1, fvalue_2, (long long)global_hz, (long long)ZYX[0], (long long)ZYX[1], (long long)ZYX[2]);
                lost_element_count++;
              }
              else
              {
                //printf("HZ [%d] %f %f\n", rank, dvalue_1, dvalue_2);
                //printf("%lld [HZ] %f %f (%lld :: %lld %lld %lld)\n", (long long)global_hz, fvalue_1, fvalue_2, (long long)global_hz, (long long)ZYX[0], (long long)ZYX[1], (long long)ZYX[2]);
                element_count++;
              }
            }
          }
        }
      }
    }
  }
 
#if PIDX_HAVE_MPI
  int64_t global_volume = 0;
  if (id->idx_d->parallel_mode == 1)
    MPI_Allreduce(&element_count, &global_volume, 1, MPI_LONG_LONG, MPI_SUM, id->comm);
  else
    global_volume = element_count;
  
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
#endif
  return PIDX_success;
}
