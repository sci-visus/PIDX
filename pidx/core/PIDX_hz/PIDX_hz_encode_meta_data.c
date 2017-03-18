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



PIDX_return_code PIDX_hz_encode_meta_data_create(PIDX_hz_encode_id id)
{
  int **tpatch;
  int **allign_offset;
  int **allign_count;
  int j = 0, p = 0, d = 0, v = 0;

  int maxH = id->idx_d->maxh;

  tpatch = (int**) malloc(2 * sizeof (int*));
  memset(tpatch, 0, 2 * sizeof (int*));
  tpatch[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  tpatch[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[0], 0, PIDX_MAX_DIMENSIONS * sizeof (int));
  memset(tpatch[1], 0, PIDX_MAX_DIMENSIONS * sizeof (int));

  if(maxH <= 0)
  {
    fprintf(stderr, "[%s] [%d] maxH [%d] not set.\n", __FILE__, __LINE__, maxH);
    return PIDX_err_hz;
  }

  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    var->hz_buffer = malloc(sizeof(*(var->hz_buffer)) * var->patch_group_count);
    memset(var->hz_buffer, 0, sizeof(*(var->hz_buffer)) * var->patch_group_count);

    for (p = 0; p < var->patch_group_count; p++)
    {
      var->hz_buffer[p] = malloc(sizeof(*(var->hz_buffer[p])));
      memset(var->hz_buffer[p], 0, sizeof(*(var->hz_buffer[p])));

      HZ_buffer hz_buf = var->hz_buffer[p];
      if (var->chunk_patch_group[p]->type == 0)
      {
        unsigned long long buffer_size = 1;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          buffer_size = buffer_size * (var->chunk_patch_group[p]->patch[0]->size[d]/id->idx->chunk_size[d]);

        hz_buf->buffer_index = malloc(buffer_size * sizeof (unsigned long long));
        memset(hz_buf->buffer_index, 0, buffer_size * sizeof (unsigned long long));
      }

      hz_buf->type = var->chunk_patch_group[p]->type;

      hz_buf->start_hz_index = malloc(sizeof (unsigned long long) * maxH);
      hz_buf->end_hz_index = malloc(sizeof (unsigned long long) * maxH);
      memset(hz_buf->start_hz_index, 0, sizeof (unsigned long long) * maxH);
      memset(hz_buf->end_hz_index, 0, sizeof (unsigned long long) * maxH);

      hz_buf->compressed_buffer_size = malloc(sizeof (unsigned long long) * maxH);
      memset(hz_buf->compressed_buffer_size, 0, sizeof (unsigned long long) * maxH);

      hz_buf->nsamples_per_level = malloc(sizeof (int*) * maxH);
      if (hz_buf->nsamples_per_level == NULL)
      memset(hz_buf->nsamples_per_level, 0, sizeof (int*) * maxH);

      if (var->chunk_patch_group[p]->type == 0)
      {
        hz_buf->samples_per_level = malloc( maxH * sizeof (unsigned long long));
        memset(hz_buf->samples_per_level, 0, maxH * sizeof (unsigned long long));
      }

      allign_offset = malloc(sizeof (int*) * maxH);
      allign_count = malloc(sizeof (int*) * maxH);
      memset(allign_offset, 0, sizeof (int*) * maxH);
      memset(allign_count, 0, sizeof (int*) * maxH);

      for (j = 0; j < maxH; j++)
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

        Point3D startXYZ;
        startXYZ.x = allign_offset[j][0];
        startXYZ.y = allign_offset[j][1];
        startXYZ.z = allign_offset[j][2];
        hz_buf->start_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, startXYZ);

        Point3D endXYZ;
        endXYZ.x = allign_count[j][0];
        endXYZ.y = allign_count[j][1];
        endXYZ.z = allign_count[j][2];

        hz_buf->end_hz_index[j] = xyz_to_HZ(id->idx->bitPattern, maxH - 1, endXYZ);

        //if (rank == 0)
        //  printf("[%d] [%d] [%d %d %d - %d %d %d] : SE %d %d T %d %d %d\n", j, p, tpatch[0][0], tpatch[0][1], tpatch[0][2], tpatch[1][0], tpatch[1][1], tpatch[1][2], hz_buf->start_hz_index[j], hz_buf->end_hz_index[j], hz_buf->nsamples_per_level[j][0], hz_buf->nsamples_per_level[j][1], hz_buf->nsamples_per_level[j][2]);
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



PIDX_return_code PIDX_hz_encode_meta_data_destroy(PIDX_hz_encode_id id)
{
  int itr = 0, p = 0, v = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];
  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
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
