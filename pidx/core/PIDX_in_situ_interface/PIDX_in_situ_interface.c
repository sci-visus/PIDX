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

/**
 * \file PIDX_in_situ_interface.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_in_situ_interface.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_insitu_id PIDX_in_situ_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, idx_comm idx_c, idx_debug idx_dbg, int group_index, int var_start_index, int var_end_index)
{
  PIDX_insitu_id insitu_id;
  insitu_id = (PIDX_insitu_id)malloc(sizeof (*insitu_id));
  memset(insitu_id, 0, sizeof (*insitu_id));

  insitu_id->idx = idx_meta_data;
  insitu_id->idx_derived = idx_derived;
  insitu_id->idx_c = idx_c;
  insitu_id->idx_dbg = idx_dbg;

  insitu_id->group_index = group_index;
  insitu_id->first_index = var_start_index;
  insitu_id->last_index = var_end_index;

  return insitu_id;
}


PIDX_return_code PIDX_in_situ_finalize(PIDX_insitu_id insitu_id)
{
  free(insitu_id);
  insitu_id = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_in_situ_perform(PIDX_insitu_id id)
{
#if PIDX_HAVE_PMT
  int ret = 0;
  int i = 0;
  int color = 0;
  MPI_Comm insitu_comm;


  for (i = 0; i < id->idx->regridded_process_count[0] * id->idx->regridded_process_count[1] * id->idx->regridded_process_count[2]; i++)
  {
    if (id->idx_c->grank == id->idx->regridded_patch[i]->rank)
    {
      color = 1;
      break;
    }
  }

  //if (color == 1)
  //  fprintf(stderr, "[INSITU] [%d] --> %d (%d / %d) %d (%d / %d) %d (%d / %d) : %d\n", id->idx_c->grank, id->idx->regridded_process_count[0], id->idx->bounds[0], id->idx->reg_patch_size[0],  id->idx->regridded_process_count[1], id->idx->bounds[1], id->idx->reg_patch_size[1], id->idx->regridded_process_count[2], id->idx->bounds[2], id->idx->reg_patch_size[2], color);




//  if (id->idx->variable_grp[id->group_index]->variable[id->first_index]->patch_group_count != 0)
//    color = 1;

  ret = MPI_Comm_split(id->idx_c->global_comm, color, id->idx_c->grank, &(insitu_comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }



  if (color == 1)
  {
    //int sf = 2;
    Ndim_patch r_p = id->idx->variable_grp[id->group_index]->variable[id->first_index]->rst_patch_group[0]->reg_patch;

    /*
    int s_x = r_p->size[0];
    int s_y = r_p->size[1];
    int s_z = r_p->size[2];

    if (s_x % sf == 0)
      s_x = s_x / sf;
    else
      s_x = s_x / sf + 1;

    if (s_y % sf == 0)
      s_y = s_y / sf;
    else
      s_y = s_y / sf + 1;

    if (s_z % sf == 0)
      s_z = s_z / sf;
    else
      s_z = s_z / sf + 1;
    unsigned char* sub_sampled_buffer = malloc(sizeof(float) * s_x * s_y * s_z);
    */

#if 0
    int **tpatch;
    int **allign_offset;
    int **allign_count;
    int **nsamples_per_level;
    tpatch = (int**) malloc(2 * sizeof (int*));
    memset(tpatch, 0, 2 * sizeof (int*));
    tpatch[0] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
    tpatch[1] = (int*) malloc(PIDX_MAX_DIMENSIONS * sizeof (int));
    memset(tpatch[0], 0, PIDX_MAX_DIMENSIONS * sizeof (int));
    memset(tpatch[1], 0, PIDX_MAX_DIMENSIONS * sizeof (int));

    allign_offset = malloc(sizeof (int*) * id->idx_derived->maxh);
    allign_count = malloc(sizeof (int*) * id->idx_derived->maxh);
    memset(allign_offset, 0, sizeof (int*) * id->idx_derived->maxh);
    memset(allign_count, 0, sizeof (int*) * id->idx_derived->maxh);

    nsamples_per_level = malloc(sizeof (int*) * id->idx_derived->maxh);
    memset(nsamples_per_level, 0, sizeof (int*) * id->idx_derived->maxh);

    int j = 0, p = 0;
    for (j = 0; j < id->idx_derived->maxh; j++)
    {
      allign_offset[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(allign_offset[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

      allign_count[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(allign_count[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

      nsamples_per_level[j] = malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
      memset(nsamples_per_level[j], 0, sizeof (int) * PIDX_MAX_DIMENSIONS);
    }

    tpatch[0][0] = r_p->offset[0];
    tpatch[0][1] = r_p->offset[1];
    tpatch[0][2] = r_p->offset[2];

    tpatch[1][0] = r_p->offset[0] + r_p->size[0] - 1;
    tpatch[1][1] = r_p->offset[1] + r_p->size[1] - 1;
    tpatch[1][2] = r_p->offset[2] + r_p->size[2] - 1;

    Align((id->idx_derived->maxh - 1), 24, id->idx->bitPattern, tpatch, allign_offset, allign_count, nsamples_per_level);

#endif


    /*
    unsigned long long level_start_address = (unsigned long long)pow(2, 24);
    unsigned long long ZYX[PIDX_MAX_DIMENSIONS];
    Hz_to_xyz(id->idx->bitPattern, id->idx_derived->maxh - 1, level_start_address, ZYX);
    fprintf(stderr, "[%d %s] start %lld %lld %lld [%d %d %d]\n", id->idx->bitSequence, id->idx_derived->maxh, ZYX[0], ZYX[1], ZYX[2], r_p->size[0], r_p->size[1], r_p->size[2]);

    int i1 = 0, j1 = 0, k1 = 0;
    int count = 0;
    for (k1 = ZYX[2]; k1 < r_p->size[2]; k1 = k1 + sf)
    {
        for (j1 = ZYX[1]; j1 < r_p->size[1]; j1 = j1 + sf)
        {
            for (i1 = ZYX[0]; i1 < r_p->size[0]; i1 = i1 + sf)
            {
              // (r_p->size[0] * r_p->size[1] * k1) + (r_p->size[0] * j1) + i1
              ((float*)sub_sampled_buffer)[count] = ((float*)r_p->buffer)[(r_p->size[0] * r_p->size[1] * k1) + (r_p->size[0] * j1) + i1];
              assert (count < s_x * s_y * s_z);
              count++;
            }
        }
    }

    //

    r_p->offset[0] = r_p->offset[0] / sf;
    r_p->offset[1] = r_p->offset[1] / sf;
    r_p->offset[2] = r_p->offset[2] / sf;

    if (r_p->size[0] % sf == 0)
      r_p->size[0] = r_p->size[0] / sf;
    else
      r_p->size[0] = r_p->size[0] / sf + 1;

    if (r_p->size[1] % sf == 0)
      r_p->size[1] = r_p->size[1] / sf;
    else
      r_p->size[1] = r_p->size[1] / sf + 1;

    if (r_p->size[2] % sf == 0)
      r_p->size[2] = r_p->size[2] / sf;
    else
      r_p->size[2] = r_p->size[2] / sf + 1;
    */

    int low[3] = {(int)r_p->offset[0], (int)r_p->offset[1], (int)r_p->offset[2]};
    int high[3] = {(int)r_p->offset[0] + (int) r_p->size[0] - 1, (int)r_p->offset[1] + (int) r_p->size[1] - 1, (int)r_p->offset[2] + (int) r_p->size[2] - 1};
    int bound[3] = {(int)id->idx->bounds[0], (int)id->idx->bounds[1], (int)id->idx->bounds[2]};

    //bound[0] = bound[0] / sf;
    //bound[1] = bound[1] / sf;
    //bound[2] = bound[2] / sf;

    //if (id->idx_c->grank == 0)
    //{
    //  fprintf(stderr, "[xxxx %d] Low: %d %d %d High %d %d %d\n", id->idx_c->grank, (int)r_p->offset[0], (int)r_p->offset[1], (int)r_p->offset[2], (int) r_p->size[0], (int) r_p->size[1], (int) r_p->size[2]);
    //  fprintf(stderr, "PIDX %d %d %d\n", (int)id->idx->bounds[0], (int)id->idx->bounds[1], (int)id->idx->bounds[2]);
    //}

    //
    //float x1;
    //memcpy(&x1, r_p->buffer, sizeof(float));
    //fprintf(stderr, "x1 = %f\n", x1);

    //
    execute_pmt((char*)r_p->buffer,
                low,
                high,
                bound,
                (uint32_t*) id->idx->regridded_process_count,
                2, 2.9, "FINAL", 2,  insitu_comm);
    //

    //
    //double e_t = MPI_Wtime();
    //fprintf(stderr, "PIDX %d Time %f\n", id->idx_c->grank, e_t - s_t);
    //

    /*
    int local_bounds[6] = {(int)r_p->offset[0], (int)r_p->offset[0] + (int) r_p->size[0] - 1,
                           (int)r_p->offset[1], (int)r_p->offset[1] + (int) r_p->size[1] - 1,
                           (int)r_p->offset[2], (int)r_p->offset[2] + (int) r_p->size[2] - 1};

    int *image_bounds;
    unsigned char *image;
    unsigned char *zbuf;
    int rank;
    MPI_Comm_rank(insitu_comm, &rank);

    image = malloc(sizeof(*image) * 512 * 512 * 3);
    memset(image, 0, sizeof(*image) * 512 * 512 * 3);

    zbuf = malloc(sizeof(*zbuf) * 512 * 512);
    memset(zbuf, 0, sizeof(*zbuf) * 512 * 512);

    volume_render(local_bounds, (char*)r_p->buffer, 512, 512, bound, image, zbuf, rank);
    */
  }

  MPI_Comm_free(&insitu_comm);

#endif
  return PIDX_success;
}
