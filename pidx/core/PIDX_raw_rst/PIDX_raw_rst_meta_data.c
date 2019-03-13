/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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
 * \file PIDX_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_raw_rst.h
 *
 */

#include "../../PIDX_inc.h"


///
/// \brief intersectNDChunk
/// \param A
/// \param B
/// \return
///
static int intersectNDChunk(PIDX_patch A, PIDX_patch B);


///
/// \brief contains_patch
/// \param reg_patch
/// \param patches
/// \param count
/// \return
///
static int contains_patch(PIDX_patch reg_patch, PIDX_patch* patches, int count);


PIDX_return_code PIDX_raw_rst_meta_data_create(PIDX_raw_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];
  int p = 0, v = 0, j = 0;

  int r, d, c;
  uint64_t i, k, max_vol, patch_count, pc;
  uint64_t reg_patch_count = 0;
  int edge_case = 0;

  int start_var_index = rst_id->first_index;

  MPI_Allreduce(&rst_id->idx->variable[start_var_index]->sim_patch_count, &rst_id->sim_max_raw_io_restructured_super_patch_count, 1, MPI_INT, MPI_MAX, rst_id->idx_c->simulation_comm);

  rst_id->sim_raw_r_count = malloc(sizeof (uint64_t) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count);
  memset(rst_id->sim_raw_r_count, 0, (sizeof (uint64_t) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count));
  rst_id->sim_raw_r_offset = malloc(sizeof (uint64_t) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count);
  memset(rst_id->sim_raw_r_offset, 0, (sizeof (uint64_t) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count));

  for (pc=0; pc < rst_id->idx->variable[start_var_index]->sim_patch_count; pc++)
  {
    uint64_t* tempoff = rst_id->idx->variable[start_var_index]->sim_patch[pc]->offset;
    uint64_t* tempsize = rst_id->idx->variable[start_var_index]->sim_patch[pc]->size;

    uint64_t index = rst_id->idx_c->simulation_rank * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count) + pc*PIDX_MAX_DIMENSIONS;
    uint64_t* curr_patch_offset = &rst_id->sim_raw_r_offset[index];
    uint64_t* curr_patch_size = &rst_id->sim_raw_r_count[index];

    memcpy(curr_patch_offset, tempoff,sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(curr_patch_size, tempsize,sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  }

  uint64_t* count_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*count_buffer_copy));
  memset(count_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*count_buffer_copy));

  memcpy(count_buffer_copy, &rst_id->sim_raw_r_count[rst_id->idx_c->simulation_rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*count_buffer_copy));

  MPI_Allgather(count_buffer_copy, PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count*sizeof(uint64_t), MPI_BYTE, rst_id->sim_raw_r_count, PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(uint64_t), MPI_BYTE, rst_id->idx_c->simulation_comm);
  free(count_buffer_copy);

  uint64_t* offset_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*offset_buffer_copy));
  memset(offset_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*offset_buffer_copy));

  memcpy(offset_buffer_copy, &rst_id->sim_raw_r_offset[rst_id->idx_c->simulation_rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(*offset_buffer_copy));

  MPI_Allgather(offset_buffer_copy, PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(uint64_t), MPI_BYTE, rst_id->sim_raw_r_offset, PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count * sizeof(uint64_t), MPI_BYTE, rst_id->idx_c->simulation_comm);
  free(offset_buffer_copy);

  var0->raw_io_restructured_super_patch_count = 0;

  /// extents for the local process(rst_id->idx_c->simulation_rank)
  uint64_t adjusted_bounds[PIDX_MAX_DIMENSIONS];
  memcpy(adjusted_bounds, rst_id->idx->bounds, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  uint64_t max_found_reg_patches = 1;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    adjusted_bounds[d] = rst_id->idx->box_bounds[d];
    if (rst_id->idx->box_bounds[d] % rst_id->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((rst_id->idx->box_bounds[d] / rst_id->idx->chunk_size[d]) + 1) * rst_id->idx->chunk_size[d];
    max_found_reg_patches *= ceil((float)rst_id->idx->box_bounds[d]/(float)rst_id->reg_patch_size[d]);
  }

  //fprintf(stderr, "max_found_reg_patches = %d [%d %d %d]\n", max_found_reg_patches, rst_id->reg_patch_size[0], rst_id->reg_patch_size[1], rst_id->reg_patch_size[2]);
  rst_id->reg_raw_grp_count = 0;

  uint64_t pc0 = 0, d0 = 0;
  PIDX_patch* found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  int found_reg_patches_count = 0;
  for (pc0 = 0; pc0 < rst_id->idx->variable[start_var_index]->sim_patch_count; pc0++)
  {
    PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
    {
      local_proc_patch->offset[d0] = rst_id->idx->variable[start_var_index]->sim_patch[pc0]->offset[d0];
      local_proc_patch->size[d0] = rst_id->idx->variable[start_var_index]->sim_patch[pc0]->size[d0];
    }

    for (i = 0; i < adjusted_bounds[0]; i = i + rst_id->reg_patch_size[0])
      for (j = 0; j < adjusted_bounds[1]; j = j + rst_id->reg_patch_size[1])
        for (k = 0; k < adjusted_bounds[2]; k = k + rst_id->reg_patch_size[2])
        {
          PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
          memset(reg_patch, 0, sizeof (*reg_patch));

          //Interior regular patches
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = rst_id->reg_patch_size[0];
          reg_patch->size[1] = rst_id->reg_patch_size[1];
          reg_patch->size[2] = rst_id->reg_patch_size[2];

          //Edge regular patches
          if ((i + rst_id->reg_patch_size[0]) > adjusted_bounds[0])
            reg_patch->size[0] = adjusted_bounds[0] - i;
          if ((j + rst_id->reg_patch_size[1]) > adjusted_bounds[1])
            reg_patch->size[1] = adjusted_bounds[1] - j;
          if ((k + rst_id->reg_patch_size[2]) > adjusted_bounds[2])
            reg_patch->size[2] = adjusted_bounds[2] - k;

          if (intersectNDChunk(reg_patch, local_proc_patch))
          {
            if (!contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
            {
              found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
              memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));

              found_reg_patches_count++;
              rst_id->reg_raw_grp_count++;
            }
          }
          free(reg_patch);
        }
    free(local_proc_patch);
  }

  for (i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  found_reg_patches_count = 0;

  rst_id->reg_raw_grp = malloc(sizeof(*rst_id->reg_raw_grp) * rst_id->reg_raw_grp_count);
  memset(rst_id->reg_raw_grp, 0, sizeof(*rst_id->reg_raw_grp) * rst_id->reg_raw_grp_count);

  reg_patch_count = 0;

  /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with

  for (i = 0; i < adjusted_bounds[0]; i = i + rst_id->reg_patch_size[0])
    for (j = 0; j < adjusted_bounds[1]; j = j + rst_id->reg_patch_size[1])
      for (k = 0; k < adjusted_bounds[2]; k = k + rst_id->reg_patch_size[2])
      {
        PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));

        //Interior regular patches
        reg_patch->offset[0] = i;
        reg_patch->offset[1] = j;
        reg_patch->offset[2] = k;
        reg_patch->size[0] = rst_id->reg_patch_size[0];
        reg_patch->size[1] = rst_id->reg_patch_size[1];
        reg_patch->size[2] = rst_id->reg_patch_size[2];

        //Edge regular patches
        edge_case = 0;
        if ((i + rst_id->reg_patch_size[0]) > adjusted_bounds[0])
        {
          reg_patch->size[0] = adjusted_bounds[0] - i;
          edge_case = 1;
        }
        if ((j + rst_id->reg_patch_size[1]) > adjusted_bounds[1])
        {
          reg_patch->size[1] = adjusted_bounds[1] - j;
          edge_case = 1;
        }
        if ((k + rst_id->reg_patch_size[2]) > adjusted_bounds[2])
        {
          reg_patch->size[2] = adjusted_bounds[2] - k;
          edge_case = 1;
        }

        for (pc0 = 0; pc0 < rst_id->idx->variable[start_var_index]->sim_patch_count; pc0++)
        {
          PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
          memset(local_proc_patch, 0, sizeof (*local_proc_patch));
          for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
          {
            local_proc_patch->offset[d0] = rst_id->idx->variable[start_var_index]->sim_patch[pc0]->offset[d0];
            local_proc_patch->size[d0] = rst_id->idx->variable[start_var_index]->sim_patch[pc0]->size[d0];
          }

          /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
          if (intersectNDChunk(reg_patch, local_proc_patch) && !contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
          {
            found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
            memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));
            found_reg_patches_count++;

            rst_id->reg_raw_grp[reg_patch_count] = malloc(sizeof(*(rst_id->reg_raw_grp[reg_patch_count])));
            memset(rst_id->reg_raw_grp[reg_patch_count], 0, sizeof(*(rst_id->reg_raw_grp[reg_patch_count])));

            PIDX_super_patch patch_grp = rst_id->reg_raw_grp[reg_patch_count];

            patch_grp->source_patch = (PIDX_source_patch_index*)malloc(sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
            patch_grp->patch = malloc(sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
            patch_grp->restructured_patch = malloc(sizeof(*patch_grp->restructured_patch));
            memset(patch_grp->source_patch, 0, sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
            memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
            memset(patch_grp->restructured_patch, 0, sizeof(*patch_grp->restructured_patch));

            patch_count = 0;
            patch_grp->patch_count = 0;
            if (edge_case == 0)
              patch_grp->is_boundary_patch = 1;
            else
              patch_grp->is_boundary_patch = 2;

            //Iterate through all processes
            for (r = 0; r < rst_id->idx_c->simulation_nprocs; r++)
            {
              for (pc = 0; pc < rst_id->sim_max_raw_io_restructured_super_patch_count; pc++)
              {
                //Extent of process with rst_id->idx_c->simulation_rank r
                PIDX_patch curr_patch = malloc(sizeof (*curr_patch));
                memset(curr_patch, 0, sizeof (*curr_patch));

                for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                {
                  uint64_t index = r * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_raw_io_restructured_super_patch_count) + pc*PIDX_MAX_DIMENSIONS;
                  curr_patch->offset[d] = rst_id->sim_raw_r_offset[index+d];
                  curr_patch->size[d] = rst_id->sim_raw_r_count[index+d];
                }

                if (curr_patch->size[0] == 0)
                {
                  // not existing patch for current rank, skip
                  //fprintf(stderr, "%d: skipping not existing patch\n", rank);
                  free(curr_patch);
                  continue;
                }

                if (intersectNDChunk(reg_patch, curr_patch))
                {
                  patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
                  memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

                  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                  {
                    //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
                    if (curr_patch->offset[d] <= reg_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                    {
                      patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = (curr_patch->offset[d] + curr_patch->size[d] - 1) - reg_patch->offset[d] + 1;
                    }
                    else if (reg_patch->offset[d] <= curr_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                    {
                      patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - curr_patch->offset[d] + 1;
                    }
                    else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (curr_patch->offset[d] + curr_patch->size[d] - 1) && reg_patch->offset[d] >= curr_patch->offset[d])
                    {
                      patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
                    }
                    else if (( curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && curr_patch->offset[d] >= reg_patch->offset[d])
                    {
                      patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = curr_patch->size[d];
                    }

                    //offset and count of intersecting regular patch
                    patch_grp->restructured_patch->offset[d] = reg_patch->offset[d];
                    patch_grp->restructured_patch->size[d] = reg_patch->size[d];
                  }

                  patch_grp->source_patch[patch_count].rank = r;
                  patch_grp->source_patch[patch_count].index = pc;
                  patch_count++;

                  if (patch_count >= rst_id->maximum_neighbor_count)
                  {
                    rst_id->maximum_neighbor_count = rst_id->maximum_neighbor_count * 2;

                    PIDX_source_patch_index *temp_buffer2 = realloc(patch_grp->source_patch, rst_id->maximum_neighbor_count * sizeof(PIDX_source_patch_index));
                    if (temp_buffer2 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->source_patch = temp_buffer2;

                    PIDX_patch *temp_buffer3 = realloc(patch_grp->patch, rst_id->maximum_neighbor_count * sizeof(*patch_grp->patch));
                    if (temp_buffer3 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->patch = temp_buffer3;

                    if (rst_id->idx_c->simulation_rank == 0)
                      fprintf(stderr, "[ERROR] rst_id->maximum_neighbor_count needs to be increased\n");
                  }

                  patch_grp->patch_count = patch_count;
                }
                free(curr_patch);
              }
            }

            patch_grp->max_patch_rank = patch_grp->source_patch[0].rank;
            max_vol = 1;
            for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              max_vol = max_vol * patch_grp->patch[0]->size[d];

            uint64_t c_vol = 1;
            for (c = 1; c < patch_grp->patch_count; c++)
            {
              c_vol = 1;
              for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                c_vol = c_vol * patch_grp->patch[c]->size[d];
              if (c_vol > max_vol)
              {
                max_vol = c_vol;
                patch_grp->max_patch_rank = patch_grp->source_patch[c].rank;
              }
            }

            if (rst_id->idx_c->simulation_rank == patch_grp->max_patch_rank)
              var0->raw_io_restructured_super_patch_count = var0->raw_io_restructured_super_patch_count + 1;

            reg_patch_count++;
          }

          free(local_proc_patch);
        }
        free(reg_patch);
      }

  for (i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx->variable[v];
    var->raw_io_restructured_super_patch_count = var0->raw_io_restructured_super_patch_count;

    var->raw_io_restructured_super_patch_count = rst_id->idx->variable[rst_id->first_index]->raw_io_restructured_super_patch_count;

    var->raw_io_restructured_super_patch = malloc(var->raw_io_restructured_super_patch_count * sizeof(*(var->raw_io_restructured_super_patch)));
    memset(var->raw_io_restructured_super_patch, 0, var->raw_io_restructured_super_patch_count * sizeof(*(var->raw_io_restructured_super_patch)));
    for (p = 0; p < var->raw_io_restructured_super_patch_count; p++)
    {
      var->raw_io_restructured_super_patch[p] = malloc(sizeof(*(var->raw_io_restructured_super_patch[p])));
      memset(var->raw_io_restructured_super_patch[p], 0, sizeof(*(var->raw_io_restructured_super_patch[p])));
    }
  }

  j = 0;
  v = 0;
  p = 0;

  int cnt = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx->variable[v];
    cnt = 0;

    for (i = 0; i < rst_id->reg_raw_grp_count; i++)
    {
      if (rst_id->idx_c->simulation_rank == rst_id->reg_raw_grp[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->raw_io_restructured_super_patch[cnt]; // here use patch_group
        patch_group->patch_count = rst_id->reg_raw_grp[i]->patch_count;
        patch_group->is_boundary_patch = rst_id->reg_raw_grp[i]->is_boundary_patch;
        patch_group->patch = malloc(sizeof(*(patch_group->patch)) * rst_id->reg_raw_grp[i]->patch_count);
        memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * rst_id->reg_raw_grp[i]->patch_count);

        patch_group->restructured_patch = malloc(sizeof(*(patch_group->restructured_patch)));
        memset(patch_group->restructured_patch, 0, sizeof(*(patch_group->restructured_patch)));

        for (j = 0; j < rst_id->reg_raw_grp[i]->patch_count; j++)
        {
          patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
          memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));

          memcpy(patch_group->patch[j]->offset, rst_id->reg_raw_grp[i]->patch[j]->offset, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
          memcpy(patch_group->patch[j]->size, rst_id->reg_raw_grp[i]->patch[j]->size, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));
        }
        memcpy(patch_group->restructured_patch->offset, rst_id->reg_raw_grp[i]->restructured_patch->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
        memcpy(patch_group->restructured_patch->size, rst_id->reg_raw_grp[i]->restructured_patch->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
        cnt++;
      }
    }

    if (cnt != var->raw_io_restructured_super_patch_count) //TODO CHECK THIS
      return PIDX_err_rst;
  }

  free(rst_id->sim_raw_r_offset);
  free(rst_id->sim_raw_r_count);

  return PIDX_success;
}



PIDX_return_code PIDX_raw_rst_meta_data_write(PIDX_raw_rst_id rst_id)
{
  int *global_patch_offset;
  int *global_patch_size;
  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];
  int max_patch_count;
  int patch_count = var0->raw_io_restructured_super_patch_count;
  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, rst_id->idx_c->simulation_comm);

  int *local_patch_offset = malloc(sizeof(uint32_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_offset, 0, sizeof(uint32_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));

  int *local_patch_size = malloc(sizeof(uint32_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_size, 0, sizeof(uint32_t) * (max_patch_count * PIDX_MAX_DIMENSIONS + 1));

  int pcounter = 0;
  int i = 0, d = 0;
  local_patch_offset[0] = (uint32_t)patch_count;
  local_patch_size[0] = (uint32_t)patch_count;
  for (i = 0; i < patch_count; i++)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_patch_offset[i * PIDX_MAX_DIMENSIONS + d + 1] = (uint32_t)var0->raw_io_restructured_super_patch[i]->restructured_patch->offset[d];
      local_patch_size[i * PIDX_MAX_DIMENSIONS + d + 1] = (uint32_t)var0->raw_io_restructured_super_patch[i]->restructured_patch->size[d];
    }
    pcounter++;
  }

  global_patch_offset = malloc((rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));
  memset(global_patch_offset, 0,(rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));

  global_patch_size = malloc((rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));
  memset(global_patch_size, 0, (rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));

  MPI_Allgather(local_patch_offset, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_INT, global_patch_offset + 2, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_INT, rst_id->idx_c->simulation_comm);
  MPI_Allgather(local_patch_size, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_INT, global_patch_size + 2, PIDX_MAX_DIMENSIONS * max_patch_count + 1, MPI_INT, rst_id->idx_c->simulation_comm);

  global_patch_size[0] = rst_id->idx_c->simulation_nprocs;
  global_patch_offset[0] = rst_id->idx_c->simulation_nprocs;
  global_patch_size[1] = max_patch_count;
  global_patch_offset[1] = max_patch_count;

  char *directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", directory_path);
  sprintf(size_path, "%s_SIZE", directory_path);
  free(directory_path);
  if (rst_id->idx_c->simulation_rank == 1 || rst_id->idx_c->simulation_nprocs == 1)
  {
    int fp = open(offset_path, O_CREAT | O_WRONLY, 0664);
    uint64_t write_count = pwrite(fp, global_patch_offset, (rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t), 0);
    if (write_count != (rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);

    fp = open(size_path, O_CREAT | O_WRONLY, 0664);
    write_count = pwrite(fp, global_patch_size, (rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t), 0);
    if (write_count != (rst_id->idx_c->simulation_nprocs * (max_patch_count * PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t))
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


PIDX_return_code PIDX_raw_rst_meta_data_destroy(PIDX_raw_rst_id rst_id)
{
  for (int v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx->variable[v];
    for (int i = 0; i < rst_id->idx->variable[v]->raw_io_restructured_super_patch_count; i++)
    {
      for (int j = 0; j < rst_id->idx->variable[v]->raw_io_restructured_super_patch[i]->patch_count; j++)
        free(var->raw_io_restructured_super_patch[i]->patch[j]);

      free(var->raw_io_restructured_super_patch[i]->restructured_patch);
      free(var->raw_io_restructured_super_patch[i]->patch);
      free(var->raw_io_restructured_super_patch[i]);

    }
    free(var->raw_io_restructured_super_patch);
  }


  for (int i = 0; i < rst_id->reg_raw_grp_count; i++)
  {
    for (int j = 0; j < rst_id->reg_raw_grp[i]->patch_count ; j++ )
      free(rst_id->reg_raw_grp[i]->patch[j]);

    free(rst_id->reg_raw_grp[i]->source_patch);
    free(rst_id->reg_raw_grp[i]->patch);
    free(rst_id->reg_raw_grp[i]->restructured_patch);
    free(rst_id->reg_raw_grp[i]);
  }

  free(rst_id->reg_raw_grp);
  rst_id->reg_raw_grp = 0;

  return PIDX_success;
}

/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}


static int contains_patch(PIDX_patch reg_patch, PIDX_patch* patches, int count)
{
  int i=0;

  for (i=0; i<count; i++)
  {
    int d=0;
    int matches = 0;
    for (d=0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      if (reg_patch->offset[d] == patches[i]->offset[d] && reg_patch->size[d] == patches[i]->size[d])
        matches++;
    }

    if (matches == PIDX_MAX_DIMENSIONS)
      return 1;
  }

  return 0;
}
