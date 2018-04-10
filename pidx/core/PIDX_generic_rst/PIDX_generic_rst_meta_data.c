/*
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


#include "../../PIDX_inc.h"

#if 1

static int maximum_neighbor_count = 256;
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);
static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);


PIDX_return_code PIDX_generic_rst_meta_data_create(PIDX_generic_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  int j = 0;
  int v = 0;

  int r, d;
  unsigned long long i, patch_count;
  int reg_patch_count;


  var0->patch_group_count = 0;

  rst_id->rank_r_offset = malloc(sizeof (unsigned long long) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS);
  memset(rst_id->rank_r_offset, 0, (sizeof (unsigned long long) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS));

  rst_id->rank_r_count =  malloc(sizeof (unsigned long long) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS);
  memset(rst_id->rank_r_count, 0, (sizeof (unsigned long long) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS));

  MPI_Allgather(var0->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rst_id->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rst_id->idx_c->global_comm);

  MPI_Allgather(var0->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rst_id->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rst_id->idx_c->global_comm);

  /// STEP 1 : Compute the dimension of the regular patch
  memcpy(rst_id->reg_patch_size, rst_id->idx->reg_patch_size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  /// extents for the local process(rank)
  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_proc_patch->offset[d] = rst_id->rank_r_offset[PIDX_MAX_DIMENSIONS * rst_id->idx_c->grank + d];
    local_proc_patch->size[d] = rst_id->rank_r_count[PIDX_MAX_DIMENSIONS * rst_id->idx_c->grank + d];
  }

  unsigned long long adjusted_bounds[PIDX_MAX_DIMENSIONS];
  memcpy(adjusted_bounds, rst_id->idx->bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    adjusted_bounds[d] = rst_id->idx->box_bounds[d];
    if (rst_id->idx->box_bounds[d] % rst_id->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((rst_id->idx->box_bounds[d] / rst_id->idx->chunk_size[d]) + 1) * rst_id->idx->chunk_size[d];
  }

  rst_id->reg_patch_grp_count = 0;
  for (i = 0; i < rst_id->idx->regridded_process_count[0] * rst_id->idx->regridded_process_count[1] * rst_id->idx->regridded_process_count[2]; i++)
  {
    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    //Interior regular patches
    reg_patch->offset[0] = rst_id->idx->regridded_patch[i]->offset[0];
    reg_patch->offset[1] = rst_id->idx->regridded_patch[i]->offset[1];
    reg_patch->offset[2] = rst_id->idx->regridded_patch[i]->offset[2];
    reg_patch->size[0] = rst_id->idx->regridded_patch[i]->size[0];
    reg_patch->size[1] = rst_id->idx->regridded_patch[i]->size[1];
    reg_patch->size[2] = rst_id->idx->regridded_patch[i]->size[2];

    if (intersectNDChunk(reg_patch, local_proc_patch))
    {
      if (rst_id->idx_c->grank == rst_id->idx->regridded_patch[i]->rank)
      {
        free(reg_patch);
        continue;
      }

      rst_id->reg_patch_grp_count++;
    }

    free(reg_patch);
  }

  for (i = 0; i < rst_id->idx->regridded_process_count[0] * rst_id->idx->regridded_process_count[1] * rst_id->idx->regridded_process_count[2]; i++)
  {
    if (rst_id->idx_c->grank == rst_id->idx->regridded_patch[i]->rank)
      rst_id->reg_patch_grp_count++;
  }


  rst_id->reg_patch_grp = (Ndim_patch_group*)malloc(sizeof(*rst_id->reg_patch_grp) * rst_id->reg_patch_grp_count);
  memset(rst_id->reg_patch_grp, 0, sizeof(*rst_id->reg_patch_grp) * rst_id->reg_patch_grp_count);

  reg_patch_count = 0;
  /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with

  for (i = 0; i < rst_id->idx->regridded_process_count[0] * rst_id->idx->regridded_process_count[1] * rst_id->idx->regridded_process_count[2]; i++)
  {
    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    //Interior regular patches
    reg_patch->offset[0] = rst_id->idx->regridded_patch[i]->offset[0];
    reg_patch->offset[1] = rst_id->idx->regridded_patch[i]->offset[1];
    reg_patch->offset[2] = rst_id->idx->regridded_patch[i]->offset[2];
    reg_patch->size[0] = rst_id->idx->regridded_patch[i]->size[0];
    reg_patch->size[1] = rst_id->idx->regridded_patch[i]->size[1];
    reg_patch->size[2] = rst_id->idx->regridded_patch[i]->size[2];

    /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
    if (intersectNDChunk(reg_patch, local_proc_patch))
    {
      if (rst_id->idx_c->grank == rst_id->idx->regridded_patch[i]->rank)
      {
        free(reg_patch);
        continue;
      }

      rst_id->reg_patch_grp[reg_patch_count] = malloc(sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));
      memset(rst_id->reg_patch_grp[reg_patch_count], 0, sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));

      Ndim_patch_group patch_grp = rst_id->reg_patch_grp[reg_patch_count];

      patch_grp->source_patch_rank = (int*)malloc(sizeof(int) * maximum_neighbor_count);
      patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
      patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
      memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
      memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
      memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

      patch_count = 0;
      patch_grp->count = 0;
      patch_grp->type = rst_id->idx->regridded_patch[i]->edge;

      //Iterate through all processes
      for (r = 0; r < rst_id->idx_c->gnprocs; r++)
      {
        //Extent of process with rank r
        Ndim_patch rank_r_patch = malloc(sizeof (*rank_r_patch));
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));

        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          rank_r_patch->offset[d] = rst_id->rank_r_offset[PIDX_MAX_DIMENSIONS * r + d];
          rank_r_patch->size[d] = rst_id->rank_r_count[PIDX_MAX_DIMENSIONS * r + d];
        }

        //If process with rank r intersects with the regular patch, then calculate the offset, count and volume of the intersecting volume
        if (intersectNDChunk(reg_patch, rank_r_patch))
        {
          patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
          memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
            if (rank_r_patch->offset[d] <= reg_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) - reg_patch->offset[d] + 1;
            }
            else if (reg_patch->offset[d] <= rank_r_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - rank_r_patch->offset[d] + 1;
            }
            else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) && reg_patch->offset[d] >= rank_r_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
            }
            else if (( rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && rank_r_patch->offset[d] >= reg_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = rank_r_patch->size[d];
            }

            //offset and count of intersecting regular patch
            patch_grp->reg_patch->offset[d] = reg_patch->offset[d];
            patch_grp->reg_patch->size[d] = reg_patch->size[d];
          }

          patch_grp->source_patch_rank[patch_count] = r;
          patch_count++;

          if (patch_count >= maximum_neighbor_count)
          {
            maximum_neighbor_count = maximum_neighbor_count * 2;

            int *temp_buffer2 = realloc(patch_grp->source_patch_rank, maximum_neighbor_count * sizeof(int));
            if (temp_buffer2 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->source_patch_rank = temp_buffer2;

            Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
            if (temp_buffer3 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->patch = temp_buffer3;
          }

          patch_grp->count = patch_count;
        }
        free(rank_r_patch);
      }

      patch_grp->max_patch_rank = rst_id->idx->regridded_patch[i]->rank;
      if(rst_id->idx_c->grank == patch_grp->max_patch_rank)
      {
        var0->patch_group_count = var0->patch_group_count + 1;
        assert (var0->patch_group_count == 1);
      }

      reg_patch_count++;
    }
    free(reg_patch);
  }


  for (i = 0; i < rst_id->idx->regridded_process_count[0] * rst_id->idx->regridded_process_count[1] * rst_id->idx->regridded_process_count[2]; i++)
  {
    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    //Interior regular patches
    reg_patch->offset[0] = rst_id->idx->regridded_patch[i]->offset[0];// i;
    reg_patch->offset[1] = rst_id->idx->regridded_patch[i]->offset[1];// j;
    reg_patch->offset[2] = rst_id->idx->regridded_patch[i]->offset[2];// k;
    reg_patch->size[0] = rst_id->idx->regridded_patch[i]->size[0]; //rst_id->reg_patch_size[0];
    reg_patch->size[1] = rst_id->idx->regridded_patch[i]->size[1]; //rst_id->reg_patch_size[1];
    reg_patch->size[2] = rst_id->idx->regridded_patch[i]->size[2]; //rst_id->reg_patch_size[2];

    /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
    if (rst_id->idx_c->grank == rst_id->idx->regridded_patch[i]->rank)
    {
      rst_id->reg_patch_grp[reg_patch_count] = malloc(sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));
      memset(rst_id->reg_patch_grp[reg_patch_count], 0, sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));

      Ndim_patch_group patch_grp = rst_id->reg_patch_grp[reg_patch_count];

      patch_grp->source_patch_rank = (int*)malloc(sizeof(int) * maximum_neighbor_count);
      patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
      patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
      memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
      memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
      memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

      patch_count = 0;
      patch_grp->count = 0;
      patch_grp->type = rst_id->idx->regridded_patch[i]->edge;

      //Iterate through all processes
      for (r = 0; r < rst_id->idx_c->gnprocs; r++)
      {
        //Extent of process with rank r
        Ndim_patch rank_r_patch = malloc(sizeof (*rank_r_patch));
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));

        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          rank_r_patch->offset[d] = rst_id->rank_r_offset[PIDX_MAX_DIMENSIONS * r + d];
          rank_r_patch->size[d] = rst_id->rank_r_count[PIDX_MAX_DIMENSIONS * r + d];
        }

        //If process with rank r intersects with the regular patch, then calculate the offset, count and volume of the intersecting volume
        if (intersectNDChunk(reg_patch, rank_r_patch))
        {
          patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
          memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
            if (rank_r_patch->offset[d] <= reg_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) - reg_patch->offset[d] + 1;
            }
            else if (reg_patch->offset[d] <= rank_r_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - rank_r_patch->offset[d] + 1;
            }
            else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) && reg_patch->offset[d] >= rank_r_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
            }
            else if (( rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && rank_r_patch->offset[d] >= reg_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = rank_r_patch->size[d];
            }

            //offset and count of intersecting regular patch
            patch_grp->reg_patch->offset[d] = reg_patch->offset[d];
            patch_grp->reg_patch->size[d] = reg_patch->size[d];
          }

          patch_grp->source_patch_rank[patch_count] = r;
          patch_count++;

          if (patch_count >= maximum_neighbor_count)
          {
            maximum_neighbor_count = maximum_neighbor_count * 2;

            int *temp_buffer2 = realloc(patch_grp->source_patch_rank, maximum_neighbor_count * sizeof(int));
            if (temp_buffer2 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->source_patch_rank = temp_buffer2;

            Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
            if (temp_buffer3 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->patch = temp_buffer3;
          }

          patch_grp->count = patch_count;
        }
        free(rank_r_patch);
      }

      patch_grp->max_patch_rank = rst_id->idx->regridded_patch[i]->rank;
      if(rst_id->idx_c->grank == patch_grp->max_patch_rank)
      {
        var0->patch_group_count = var0->patch_group_count + 1;
        assert (var0->patch_group_count == 1);
      }

      reg_patch_count++;
    }
    free(reg_patch);
  }


  free(local_proc_patch);



  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    var->patch_group_count = var0->patch_group_count;

    var->rst_patch_group = malloc(var->patch_group_count * sizeof(*(var->rst_patch_group)));
    memset(var->rst_patch_group, 0, var->patch_group_count * sizeof(*(var->rst_patch_group)));
  }

  int cnt = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    cnt = 0;
    for (i = 0; i < rst_id->reg_patch_grp_count; i++)
    {
      if (rst_id->idx_c->grank == rst_id->reg_patch_grp[i]->max_patch_rank)
      {
        Ndim_patch_group patch_group = var->rst_patch_group;
        patch_group->count = rst_id->reg_patch_grp[i]->count;
        patch_group->type = rst_id->reg_patch_grp[i]->type;
        patch_group->patch = malloc(sizeof(*(patch_group->patch)) * rst_id->reg_patch_grp[i]->count);
        memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * rst_id->reg_patch_grp[i]->count);

        patch_group->reg_patch = malloc(sizeof(*(patch_group->reg_patch)));
        memset(patch_group->reg_patch, 0, sizeof(*(patch_group->reg_patch)));

        for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
        {
          patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
          memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));

          memcpy(patch_group->patch[j]->offset, rst_id->reg_patch_grp[i]->patch[j]->offset, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
          memcpy(patch_group->patch[j]->size, rst_id->reg_patch_grp[i]->patch[j]->size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
        }
        memcpy(patch_group->reg_patch->offset, rst_id->reg_patch_grp[i]->reg_patch->offset, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
        memcpy(patch_group->reg_patch->size, rst_id->reg_patch_grp[i]->reg_patch->size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

        cnt++;
        assert(cnt == 1);
      }
    }
    if (cnt != var->patch_group_count)
      return PIDX_err_rst;
  }



  free(rst_id->rank_r_offset);
  free(rst_id->rank_r_count);

  return PIDX_success;
}



PIDX_return_code PIDX_generic_rst_meta_data_write(PIDX_generic_rst_id rst_id)
{
  int *global_patch_offset;
  int *global_patch_size;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  int patch_count =var0->patch_group_count;

  int *local_patch_offset = malloc(sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_offset, 0, sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));

  int *local_patch_size = malloc(sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_size, 0, sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));

  int pcounter = 0;
  int d = 0;
  local_patch_offset[0] = (uint32_t)patch_count;
  local_patch_size[0] = (uint32_t)patch_count;
  if (patch_count != 0)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_patch_offset[PIDX_MAX_DIMENSIONS + d + 1] = (uint32_t)var0->rst_patch_group->reg_patch->offset[d];
      local_patch_size[PIDX_MAX_DIMENSIONS + d + 1] = (uint32_t)var0->rst_patch_group->reg_patch->size[d];
    }
    pcounter++;
  }

  int wc = (rst_id->idx_c->gnprocs * (PIDX_MAX_DIMENSIONS + 1) + 2);
  global_patch_offset = malloc(wc * sizeof(uint32_t));
  memset(global_patch_offset, 0,wc * sizeof(uint32_t));

  global_patch_size = malloc(wc * sizeof(uint32_t));
  memset(global_patch_size, 0, wc * sizeof(uint32_t));

  MPI_Allgather(local_patch_offset, PIDX_MAX_DIMENSIONS + 1, MPI_INT, global_patch_offset + 2, PIDX_MAX_DIMENSIONS + 1, MPI_INT, rst_id->idx_c->global_comm);
  MPI_Allgather(local_patch_size, PIDX_MAX_DIMENSIONS + 1, MPI_INT, global_patch_size + 2, PIDX_MAX_DIMENSIONS + 1, MPI_INT, rst_id->idx_c->global_comm);

  global_patch_size[0] = rst_id->idx_c->gnprocs;
  global_patch_offset[0] = rst_id->idx_c->gnprocs;
  global_patch_size[1] = 1;
  global_patch_offset[1] = 1;

  char *directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", directory_path);
  sprintf(size_path, "%s_SIZE", directory_path);
  free(directory_path);
  if (rst_id->idx_c->grank == 1 || rst_id->idx_c->gnprocs == 1)
  {
    if (rst_id->idx->flip_endian == 1)
    {
      int y = 0;
      int temp;
      int temp2;
      for (y = 0; y < wc; y++)
      {
        memcpy(&temp, global_patch_offset + y, sizeof(int));
        bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
        memcpy(global_patch_offset + y, &temp2, sizeof(int));
      }

      for (y = 0; y < wc; y++)
      {
        memcpy(&temp, global_patch_size + y, sizeof(int));
        bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
        memcpy(global_patch_size + y, &temp2, sizeof(int));
      }
    }

    int fp = open(offset_path, O_CREAT | O_WRONLY, 0664);
    ssize_t write_count = pwrite(fp, global_patch_offset, wc * sizeof(uint32_t), 0);
    if (write_count != wc * sizeof(uint32_t))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);

    fp = open(size_path, O_CREAT | O_WRONLY, 0664);
    write_count = pwrite(fp, global_patch_size, wc * sizeof(uint32_t), 0);
    if (write_count != wc * sizeof(uint32_t))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);
  }

  free(local_patch_offset);
  free(local_patch_size);

  free(global_patch_offset);
  global_patch_offset = 0;

  free(global_patch_size);
  global_patch_size = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_generic_rst_meta_data_destroy(PIDX_generic_rst_id rst_id)
{
  int j, v;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  if (var0->patch_group_count == 0)
    return PIDX_success;

  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    for(j = 0; j < var_grp->variable[v]->rst_patch_group->count; j++)
    {
      free(var->rst_patch_group->patch[j]);
      var->rst_patch_group->patch[j] = 0;
    }

    free(var->rst_patch_group->reg_patch);
    var->rst_patch_group->reg_patch = 0;

    free(var->rst_patch_group->patch);
    var->rst_patch_group->patch = 0;

    free(var->rst_patch_group);
    var->rst_patch_group = 0;
  }

  if (rst_id->idx->enable_rst == 1)
  {
    int i, j;
    for (i = 0; i < rst_id->reg_patch_grp_count; i++)
    {
      for (j = 0; j < rst_id->reg_patch_grp[i]->count ; j++ )
      {
        free(rst_id->reg_patch_grp[i]->patch[j]);
        rst_id->reg_patch_grp[i]->patch[j] = 0;
      }

      free(rst_id->reg_patch_grp[i]->source_patch_rank);
      rst_id->reg_patch_grp[i]->source_patch_rank = 0;

      free(rst_id->reg_patch_grp[i]->patch);
      rst_id->reg_patch_grp[i]->patch = 0;

      free(rst_id->reg_patch_grp[i]->reg_patch);
      rst_id->reg_patch_grp[i]->reg_patch = 0;

      free(rst_id->reg_patch_grp[i]);
      rst_id->reg_patch_grp[i] = 0;
    }

    free(rst_id->reg_patch_grp);
    rst_id->reg_patch_grp = 0;
  }

  return PIDX_success;
}


/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
  unsigned char *data = ((unsigned char *)val) + 3;
  unsigned char *out = (unsigned char *)outbuf;

  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out = *data;

  return;
}
#endif
