/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */

/**
 * \file PIDX_generic_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_generic_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_generic_rst_id PIDX_generic_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index)
{
  PIDX_generic_rst_id generic_rst_id;
  generic_rst_id = (PIDX_generic_rst_id)malloc(sizeof (*generic_rst_id));
  memset(generic_rst_id, 0, sizeof (*generic_rst_id));

  generic_rst_id->idx = idx_meta_data;
  generic_rst_id->idx_d = idx_d;
  generic_rst_id->idx_c = idx_c;
  generic_rst_id->idx_dbg = idx_dbg;

  generic_rst_id->group_index = 0;
  generic_rst_id->first_index = var_start_index;
  generic_rst_id->last_index = var_end_index;

  return (generic_rst_id);
}


#if 0
PIDX_return_code PIDX_generic_rst_perform(PIDX_generic_rst_id id)
{
  int i = 0, v = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->group_index];

  for (i = 0; i < id->idx->regridded_process_count[0] * id->idx->regridded_process_count[1] * id->idx->regridded_process_count[2]; i++)
  {
    if (id->idx_c->grank = id->idx->regridded_patch[i]->rank)
    {
      id->reg_patch_grp_count = 1;

        for (v = id->first_index; v <= id->last_index; v++)
        {
          PIDX_variable var = var_grp->variable[v];
          Ndim_patch_group patch_group = var->rst_patch_group[0];
          patch_group->data_source = 0;

          for(j = 0; j < id->reg_patch_grp[i]->count; j++)
          {
            patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

            if (patch_group->patch[j]->buffer == NULL)
              return PIDX_err_rst;
            memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
          }
        }
    }
    else
      id->reg_patch_grp_count = 0;

  }

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_generic_rst_finalize(PIDX_generic_rst_id generic_rst_id)
{  
  free(generic_rst_id);
  generic_rst_id = 0;
  
  return PIDX_success;
}
