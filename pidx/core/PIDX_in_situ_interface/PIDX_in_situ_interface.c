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

  int ret = 0;
  int i = 0;
  int color = 0;
  MPI_Comm insitu_comm;


  for (i = 0; i < id->idx->number_processes[0] * id->idx->number_processes[1] * id->idx->number_processes[2]; i++)
  {
    if (id->idx_c->grank == id->idx->new_box_set[i]->rank)
    {
      color = 1;
      break;
    }
  }

  printf("[%d] --> %d (%d / %d) %d (%d / %d) %d (%d / %d) : %d\n", id->idx_c->grank, id->idx->number_processes[0], id->idx->bounds[0], id->idx->reg_patch_size[0],  id->idx->number_processes[1], id->idx->bounds[1], id->idx->reg_patch_size[1], id->idx->number_processes[2], id->idx->bounds[2], id->idx->reg_patch_size[2], color);




//  if (id->idx->variable_grp[id->group_index]->variable[id->first_index]->patch_group_count != 0)
//    color = 1;

  ret = MPI_Comm_split(id->idx_c->global_comm, color, id->idx_c->grank, &(insitu_comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }



  if (color == 1)
  {
    Ndim_patch r_p = id->idx->variable_grp[id->group_index]->variable[id->first_index]->rst_patch_group[0]->reg_patch;
    int low[3] = {(int)r_p->offset[0], (int)r_p->offset[1], (int)r_p->offset[2]};
    int high[3] = {(int)r_p->offset[0] + (int) r_p->size[0] - 1, (int)r_p->offset[1] + (int) r_p->size[1] - 1, (int)r_p->offset[2] + (int) r_p->size[2] - 1};
      //double s_t = MPI_Wtime();
    int bound[3] = {(int)id->idx->bounds[0], (int)id->idx->bounds[1], (int)id->idx->bounds[2]};
    printf("Color = %d\n", color);

    printf("xxxx Low: %d %d %d High %d %d %d\n", (int)r_p->offset[0], (int)r_p->offset[1], (int)r_p->offset[2], (int) r_p->size[0], (int) r_p->size[1], (int) r_p->size[2]);

    printf("PIDX %d %d %d\n", (int)id->idx->bounds[0], (int)id->idx->bounds[1], (int)id->idx->bounds[2]);

    float x1;
    memcpy(&x1, r_p->buffer, sizeof(float));
    printf("x1 = %f\n", x1);
    execute_pmt((char*)r_p->buffer,
                low,
                high,
                bound,
                (uint32_t*) id->idx->number_processes,
                2, 0, insitu_comm);

    //double e_t = MPI_Wtime();
    //printf("PIDX %d Time %f\n", id->idx_c->grank, e_t - s_t);
  }

  MPI_Comm_free(&insitu_comm);


  return PIDX_success;
}
