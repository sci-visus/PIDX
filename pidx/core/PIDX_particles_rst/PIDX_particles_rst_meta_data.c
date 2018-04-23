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
 * \file PIDX_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"


static int intersectNDChunk(PIDX_patch A, PIDX_patch B);
static int pointInChunk(PIDX_patch p, const double *pos);
static int contains_patch(PIDX_patch reg_patch, PIDX_patch* patches, int count);
static void gather_all_patch_extents(PIDX_particles_rst_id rst_id);
static PIDX_return_code distribute_particle_info(PIDX_particles_rst_id rst_id);
static void calculate_number_of_intersecting_boxes(PIDX_particles_rst_id rst_id);
static PIDX_return_code populate_all_intersecting_restructured_super_patch_meta_data(PIDX_particles_rst_id rst_id);
static PIDX_return_code copy_reciever_patch_info(PIDX_particles_rst_id rst_id);
static void free_patch_extents(PIDX_particles_rst_id rst_id);


PIDX_return_code PIDX_particles_rst_meta_data_create(PIDX_particles_rst_id rst_id)
{
  // Gathers the extents of all patches
  // The outcome is stored in rst_id->sim_multi_patch_r_size and rst_id->sim_multi_patch_r_offset
  gather_all_patch_extents(rst_id);


  // Calculates the number of restructured box a process intersects with.
  // It can intersect either as a receiver or as a sender
  // The outcome is stored in rst_id->intersected_restructured_super_patch_count
  calculate_number_of_intersecting_boxes(rst_id);


  // Iterates through all the imposed restructured patches, and find all the patches that intersects with it
  // The outcome is stored in rst_id->intersected_restructured_super_patch
  populate_all_intersecting_restructured_super_patch_meta_data(rst_id);

  distribute_particle_info(rst_id);

  // If a processor is a reciever i.e. it holds a restructured patch, then copy all the relevant metadata from
  // rst_id->intersected_restructured_super_patch to var->restructured_super_patch
  copy_reciever_patch_info(rst_id);


  // Free the allocated patch extents
  free_patch_extents(rst_id);

  return PIDX_success;
}


static void gather_all_patch_extents(PIDX_particles_rst_id rst_id)
{
  uint64_t pc = 0;

  // rst_id->sim_max_patch_count will contain the maximum number of patch a process is holding
  MPI_Allreduce(&rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch_count, &rst_id->sim_max_patch_count, 1, MPI_INT, MPI_MAX, rst_id->idx_c->simulation_comm);

  // buffer to hold the offset and size information of all patches across all processes
  rst_id->sim_multi_patch_r_count = malloc(sizeof (double) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count);
  memset(rst_id->sim_multi_patch_r_count, 0, (sizeof (double) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count));
  rst_id->sim_multi_patch_r_offset = malloc(sizeof (double) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count);
  memset(rst_id->sim_multi_patch_r_offset, 0, (sizeof (double) * rst_id->idx_c->simulation_nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count));


  // every process populates rst_id->sim_multi_patch_r_count and rst_id->sim_multi_patch_r_offset with their local information
  for (pc=0; pc < rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch_count; pc++)
  {
    double* tempoff = rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc]->physical_offset;
    double* tempsize = rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc]->physical_size;

    uint64_t index = rst_id->idx_c->simulation_rank * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count) + pc*PIDX_MAX_DIMENSIONS;
    double* curr_patch_offset = &rst_id->sim_multi_patch_r_offset[index];
    double* curr_patch_size = &rst_id->sim_multi_patch_r_count[index];

    memcpy(curr_patch_offset, tempoff,sizeof(double) * PIDX_MAX_DIMENSIONS);
    memcpy(curr_patch_size, tempsize,sizeof(double) * PIDX_MAX_DIMENSIONS);
  }

  // making a local copy of the patch size since mpi does not allow the same buffer to be used as both a sender and a reciever
  double* count_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*count_buffer_copy));
  memset(count_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*count_buffer_copy));

  memcpy(count_buffer_copy, &rst_id->sim_multi_patch_r_count[rst_id->idx_c->simulation_rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*count_buffer_copy));

  // gather size of all patches across all processes
  MPI_Allgather(count_buffer_copy, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count*sizeof(double), MPI_BYTE, rst_id->sim_multi_patch_r_count, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count*sizeof(double), MPI_BYTE, rst_id->idx_c->simulation_comm);
  free(count_buffer_copy);

  // making a local copy of the patch size since mpi does not allow the same buffer to be used as both a sender and a reciever
  double* offset_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*offset_buffer_copy));
  memset(offset_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*offset_buffer_copy));

  memcpy(offset_buffer_copy, &rst_id->sim_multi_patch_r_offset[rst_id->idx_c->simulation_rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count * sizeof(*offset_buffer_copy));

  // gather size of all patches across all processes
  MPI_Allgather(offset_buffer_copy, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count*sizeof(double), MPI_BYTE, rst_id->sim_multi_patch_r_offset, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_count*sizeof(double), MPI_BYTE, rst_id->idx_c->simulation_comm);
  free(offset_buffer_copy);

  return;
}

static void calculate_number_of_intersecting_boxes(PIDX_particles_rst_id rst_id)
{
  int i = 0;
  uint64_t pc0 = 0, d0 = 0;
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  uint64_t *tpc = rst_id->restructured_grid->total_patch_count;
  int max_found_reg_patches = tpc[0] * tpc[1] * tpc[2];

  rst_id->intersected_restructured_super_patch_count = 0;

  PIDX_patch* found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  // If I am the reciever (holder) of a super patch
  for (i = 0; i < tpc[0] * tpc[1] * tpc[2]; i++)
  {
    if (rst_id->idx_c->simulation_rank == rst_id->restructured_grid->patch[i]->rank)
    {
      rst_id->intersected_restructured_super_patch_count++;
    }
  }

  // If I am a sender for a super patch
  int found_reg_patches_count = 0;
  for (pc0 = 0; pc0 < var0->sim_patch_count; pc0++)
  {
    PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
    {
      local_proc_patch->physical_offset[d0] = var0->sim_patch[pc0]->physical_offset[d0];
      local_proc_patch->physical_size[d0] = var0->sim_patch[pc0]->physical_size[d0];
    }

    for (i = 0; i < tpc[0] * tpc[1] * tpc[2]; i++)
    {
      PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
      memset(reg_patch, 0, sizeof (*reg_patch));

      reg_patch->physical_offset[0] = rst_id->restructured_grid->patch[i]->physical_offset[0];
      reg_patch->physical_offset[1] = rst_id->restructured_grid->patch[i]->physical_offset[1];
      reg_patch->physical_offset[2] = rst_id->restructured_grid->patch[i]->physical_offset[2];
      reg_patch->physical_size[0] = rst_id->restructured_grid->patch[i]->physical_size[0];
      reg_patch->physical_size[1] = rst_id->restructured_grid->patch[i]->physical_size[1];
      reg_patch->physical_size[2] = rst_id->restructured_grid->patch[i]->physical_size[2];

      if (intersectNDChunk(reg_patch, local_proc_patch))
      {
        // this is corresponding to the process being the reciever of a super patch which i already accounted for
        if (rst_id->idx_c->simulation_rank == rst_id->restructured_grid->patch[i]->rank)
        {
          free(reg_patch);
          continue;
        }

        //if (rst_id->idx_c->simulation_rank == 8)
        //  printf("should not be here R%d[%f %f %f and %f %f %f]  ||  8[%f %f %f and %f %f %f]\n", i, reg_patch->physical_offset[0], reg_patch->physical_offset[1], reg_patch->physical_offset[2],
        //          reg_patch->physical_size[0], reg_patch->physical_size[1], reg_patch->physical_size[2],
        //          local_proc_patch->physical_offset[0], local_proc_patch->physical_offset[1], local_proc_patch->physical_offset[2],
        //                            local_proc_patch->physical_size[0], local_proc_patch->physical_size[1], local_proc_patch->physical_size[2]);
        if (!contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
        {
          found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
          memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));

          found_reg_patches_count++;
          rst_id->intersected_restructured_super_patch_count++;
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

  // create buffer to hold meta data for the super patches a process intersects with
  rst_id->intersected_restructured_super_patch = (PIDX_super_patch*)malloc(sizeof(*rst_id->intersected_restructured_super_patch) * rst_id->intersected_restructured_super_patch_count);
  memset(rst_id->intersected_restructured_super_patch, 0, sizeof(*rst_id->intersected_restructured_super_patch) * rst_id->intersected_restructured_super_patch_count);

  return;
}



static PIDX_return_code populate_all_intersecting_restructured_super_patch_meta_data(PIDX_particles_rst_id rst_id)
{
  int i = 0, r = 0, d = 0;
  uint64_t patch_count = 0, pc = 0, pc0 = 0, d0 = 0;
  int found_reg_patches_count = 0, reg_patch_count = 0;

  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  // Patch group size is 1 for processes that holds a restructured patch or is 0 otherwise
  var0->restructured_super_patch_count = 0;
  uint64_t *tpc = rst_id->restructured_grid->total_patch_count;

  // Total number of restructured patch after the restructuring phase
  int max_found_reg_patches = tpc[0] * tpc[1] * tpc[2];

  PIDX_patch* found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  for (i = 0; i < tpc[0] * tpc[1] * tpc[2]; i++)
  {
    PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    Ndim_empty_patch ep = rst_id->restructured_grid->patch[i];

    // extents of the restrucuted super patches
    reg_patch->physical_offset[0] = ep->physical_offset[0];
    reg_patch->physical_offset[1] = ep->physical_offset[1];
    reg_patch->physical_offset[2] = ep->physical_offset[2];
    reg_patch->physical_size[0] = ep->physical_size[0];
    reg_patch->physical_size[1] = ep->physical_size[1];
    reg_patch->physical_size[2] = ep->physical_size[2];

    for (pc0 = 0; pc0 < rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch_count; pc0++)
    {
      PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
      memset(local_proc_patch, 0, sizeof (*local_proc_patch));
      for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
      {
        local_proc_patch->physical_offset[d0] = rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc0]->physical_offset[d0];
        local_proc_patch->physical_size[d0] = rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc0]->physical_size[d0];
      }

      // If local process intersects with regular patch, then find all other process that intersects with the regular patch.
      // (if a process is a holder of a super patch || if a process sends data to a super patch)
      if (rst_id->idx_c->simulation_rank == ep->rank || (rst_id->idx_c->simulation_rank != ep->rank && intersectNDChunk(reg_patch, local_proc_patch)))
      {
          if (rst_id->idx_c->simulation_rank == 8)
              printf("reg_patch_count %d [%f %f %f]\n", reg_patch_count, reg_patch->physical_offset[0], reg_patch->physical_offset[1], reg_patch->physical_offset[2]);
        if (!contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
        {
          found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
          memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));
          found_reg_patches_count++;

          rst_id->intersected_restructured_super_patch[reg_patch_count] = malloc(sizeof(*(rst_id->intersected_restructured_super_patch[reg_patch_count])));
          memset(rst_id->intersected_restructured_super_patch[reg_patch_count], 0, sizeof(*(rst_id->intersected_restructured_super_patch[reg_patch_count])));

          PIDX_super_patch patch_grp = rst_id->intersected_restructured_super_patch[reg_patch_count];



          patch_grp->source_patch = (PIDX_source_patch_index*)malloc(sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
          patch_grp->patch = malloc(sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
          patch_grp->restructured_patch = malloc(sizeof(*patch_grp->restructured_patch));
          memset(patch_grp->source_patch, 0, sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
          memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
          memset(patch_grp->restructured_patch, 0, sizeof(*patch_grp->restructured_patch));

          patch_count = 0;
          patch_grp->patch_count = 0;
          patch_grp->is_boundary_patch = ep->is_boundary_patch;

          // Iterate through all processes
          for (r = 0; r < rst_id->idx_c->simulation_nprocs; r++)
          {
            for (pc = 0; pc < rst_id->sim_max_patch_count; pc++)
            {
              // Extent of process with rst_id->idx_c->simulation_rank r
              PIDX_patch curr_patch = malloc(sizeof (*curr_patch));
              memset(curr_patch, 0, sizeof (*curr_patch));

              for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              {
                uint64_t index = r * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_count) + pc*PIDX_MAX_DIMENSIONS;
                curr_patch->physical_offset[d] = rst_id->sim_multi_patch_r_offset[index+d];
                curr_patch->physical_size[d] = rst_id->sim_multi_patch_r_count[index+d];
              }

              if (curr_patch->physical_size[0] == 0)
              {
                // not existing patch for current rank, skip
                free(curr_patch);
                continue;
              }

              if (intersectNDChunk(reg_patch, curr_patch))
              {
                patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
                memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

                if (r == rst_id->idx_c->simulation_rank)
                {
                  for (uint64_t i = 0; i < rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc]->particle_count; ++i)
                  {
                    PIDX_variable pos_var = rst_id->idx_metadata->variable[rst_id->first_index];
                    const uint64_t bytes_per_pos = pos_var->vps * pos_var->bpv/8;
                    if (pointInChunk(reg_patch, (double*)(rst_id->idx_metadata->variable[rst_id->first_index]->sim_patch[pc]->buffer + i * bytes_per_pos)))
                      patch_grp->patch[patch_count]->particle_count++;
                  }
                  //printf("My rank %d Particle intersection %d\n", r, patch_grp->patch[patch_count]->particle_count);
                }

                for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                {
                  //offset and count of intersecting regular patch
                  patch_grp->restructured_patch->physical_offset[d] = reg_patch->physical_offset[d];
                  patch_grp->restructured_patch->physical_size[d] = reg_patch->physical_size[d];
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

          patch_grp->max_patch_rank = ep->rank;
          if (rst_id->idx_c->simulation_rank == patch_grp->max_patch_rank)
          {
            var0->restructured_super_patch_count = var0->restructured_super_patch_count + 1;
            assert (var0->restructured_super_patch_count <= 1);
          }

          reg_patch_count++;
        }
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

  return PIDX_success;
}



static PIDX_return_code distribute_particle_info(PIDX_particles_rst_id rst_id)
{

  int req_count_estimate = 0;
  int req_count_actual = 0;
  MPI_Request *req;
  MPI_Status *status;

  for (int i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    for (int j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      req_count_estimate++;

  req = malloc(sizeof (*req) * req_count_estimate * 2);
  memset(req, 0, sizeof (*req) * req_count_estimate * 2);

  status = malloc(sizeof (*status) * req_count_estimate * 2);
  memset(status, 0, sizeof (*status) * req_count_estimate * 2);

  for (int i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    PIDX_super_patch patch_grp = rst_id->intersected_restructured_super_patch[i];

    if (rst_id->idx_c->simulation_rank == patch_grp->max_patch_rank)
    {
      for (int j = 0; j < patch_grp->patch_count; j++)
      {
        if (rst_id->idx_c->simulation_rank != patch_grp->source_patch[j].rank)
        {
          MPI_Irecv(&(patch_grp->patch[j]->particle_count), 1, MPI_INT, patch_grp->source_patch[j].rank, 123, rst_id->idx_c->simulation_comm, &req[req_count_actual]);
          req_count_actual++;
          //printf("My rank %d Receiving from %d\n", rst_id->idx_c->simulation_rank, patch_grp->source_patch[j].rank);
        }
      }
    }
    else
    {
      for (int j = 0; j < patch_grp->patch_count; j++)
      {
        if (rst_id->idx_c->simulation_rank == patch_grp->source_patch[j].rank)
        {
          MPI_Isend(&(patch_grp->patch[j]->particle_count), 1, MPI_INT, patch_grp->max_patch_rank, 123, rst_id->idx_c->simulation_comm, &req[req_count_actual]);
          req_count_actual++;
          //printf("My rank %d Sending to %d\n", rst_id->idx_c->simulation_rank, patch_grp->max_patch_rank);
        }
      }
    }
  }

  if (MPI_Waitall(req_count_actual, req, status) != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  free(req);
  free(status);


  for (int i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    if (rst_id->idx_c->simulation_rank == 8)
    {
    PIDX_super_patch patch_grp = rst_id->intersected_restructured_super_patch[i];
    printf("[%d] count %d\n", i, patch_grp->patch_count);
    for (int j = 0; j < patch_grp->patch_count; j++)
    {
        printf("[%d] PC %lld R %d MR %d\n", j, (unsigned long long)patch_grp->patch[j]->particle_count, patch_grp->source_patch[j].rank, patch_grp->max_patch_rank);
    }
    }
  }

  return PIDX_success;
}



static PIDX_return_code copy_reciever_patch_info(PIDX_particles_rst_id rst_id)
{
  int v = 0, j = 0;
  int cnt = 0;
  uint64_t i = 0;

  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    cnt = 0;

    PIDX_variable var = rst_id->idx_metadata->variable[v];
    var->restructured_super_patch_count = var0->restructured_super_patch_count;

    var->restructured_super_patch = malloc(var->restructured_super_patch_count * sizeof(*(var->restructured_super_patch)));
    memset(var->restructured_super_patch, 0, var->restructured_super_patch_count * sizeof(*(var->restructured_super_patch)));

    for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      PIDX_super_patch irsp = rst_id->intersected_restructured_super_patch[i];
      if (rst_id->idx_c->simulation_rank == irsp->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->restructured_super_patch;
        patch_group->patch_count = irsp->patch_count;
        patch_group->is_boundary_patch = irsp->is_boundary_patch;
        patch_group->patch = malloc(sizeof(*(patch_group->patch)) * irsp->patch_count);
        memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * irsp->patch_count);

        patch_group->restructured_patch = malloc(sizeof(*(patch_group->restructured_patch)));
        memset(patch_group->restructured_patch, 0, sizeof(*(patch_group->restructured_patch)));

        for (j = 0; j < irsp->patch_count; j++)
        {
          patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
          memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));

          patch_group->patch[j]->particle_count = irsp->patch[j]->particle_count;
        }
        memcpy(patch_group->restructured_patch->physical_offset, irsp->restructured_patch->physical_offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
        memcpy(patch_group->restructured_patch->physical_size, irsp->restructured_patch->physical_size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
        cnt++;
        assert(cnt == 1);
      }
    }

    if (cnt != var->restructured_super_patch_count)
      return PIDX_err_rst;
  }

  return PIDX_success;
}


static void free_patch_extents(PIDX_particles_rst_id rst_id)
{
  free(rst_id->sim_multi_patch_r_offset);
  free(rst_id->sim_multi_patch_r_count);

  return;
}


PIDX_return_code PIDX_particles_rst_meta_data_write(PIDX_particles_rst_id rst_id)
{
  int *global_patch_offset;
  int *global_patch_size;
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];
  int patch_count =var0->restructured_super_patch_count;

  int *local_patch_offset = malloc(sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_offset, 0, sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));

  int *local_patch_size = malloc(sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));
  memset(local_patch_size, 0, sizeof(uint32_t) * (PIDX_MAX_DIMENSIONS + 1));

  int d = 0;
  local_patch_offset[0] = (uint32_t)patch_count;
  local_patch_size[0] = (uint32_t)patch_count;

  if (patch_count != 0)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_patch_offset[d + 1] = (uint32_t)var0->restructured_super_patch->restructured_patch->offset[d];
      local_patch_size[d + 1] = (uint32_t)var0->restructured_super_patch->restructured_patch->size[d];
    }
  }

  global_patch_offset = malloc((rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));
  memset(global_patch_offset, 0,(rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));

  global_patch_size = malloc((rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));
  memset(global_patch_size, 0, (rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t));

  MPI_Allgather(local_patch_offset, PIDX_MAX_DIMENSIONS + 1, MPI_INT, global_patch_offset + 2, PIDX_MAX_DIMENSIONS + 1, MPI_INT, rst_id->idx_c->simulation_comm);
  MPI_Allgather(local_patch_size, PIDX_MAX_DIMENSIONS + 1, MPI_INT, global_patch_size + 2, PIDX_MAX_DIMENSIONS + 1, MPI_INT, rst_id->idx_c->simulation_comm);

  global_patch_size[0] = rst_id->idx_c->simulation_nprocs;
  global_patch_offset[0] = rst_id->idx_c->simulation_nprocs;
  global_patch_size[1] = 1;
  global_patch_offset[1] = 1;

  char *directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx_metadata->filename, strlen(rst_id->idx_metadata->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", directory_path);
  sprintf(size_path, "%s_SIZE", directory_path);
  free(directory_path);
  if (rst_id->idx_c->simulation_rank == 1 || rst_id->idx_c->simulation_nprocs == 1)
  {
    int fp = open(offset_path, O_CREAT | O_WRONLY, 0664);
    uint64_t write_count = pwrite(fp, global_patch_offset, (rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t), 0);
    if (write_count != (rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);

    fp = open(size_path, O_CREAT | O_WRONLY, 0664);
    write_count = pwrite(fp, global_patch_size, (rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t), 0);
    if (write_count != (rst_id->idx_c->simulation_nprocs * (PIDX_MAX_DIMENSIONS + 1) + 2) * sizeof(uint32_t))
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


PIDX_return_code PIDX_particles_rst_meta_data_destroy(PIDX_particles_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];
  if (var0->restructured_super_patch_count == 0)
      return PIDX_success;

  int i, j, v;
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];

    for (j = 0; j < rst_id->idx_metadata->variable[v]->restructured_super_patch->patch_count; j++)
    {
      free(var->restructured_super_patch->patch[j]);
      var->restructured_super_patch->patch[j] = 0;
    }

    free(var->restructured_super_patch->restructured_patch);
    var->restructured_super_patch->restructured_patch = 0;

    free(var->restructured_super_patch->patch);
    var->restructured_super_patch->patch = 0;

    free(var->restructured_super_patch);
    var->restructured_super_patch = 0;
  }

  for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    PIDX_super_patch irsp = rst_id->intersected_restructured_super_patch[i];
    for (j = 0; j < irsp->patch_count; j++ )
      free(irsp->patch[j]);

    free(irsp->source_patch);
    free(irsp->patch);
    free(irsp->restructured_patch);
    free(irsp);
  }

  free(rst_id->intersected_restructured_super_patch);
  rst_id->intersected_restructured_super_patch = 0;

  return PIDX_success;
}


// Function to check if patches A and B intersects
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int check_bit = 0;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
  {
    check_bit = check_bit
      || ((A->physical_offset[d] + A->physical_size[d]) <= B->physical_offset[d])
      || ((B->physical_offset[d] + B->physical_size[d]) <= A->physical_offset[d]);
  }

  //printf("CB %d %d\n", check_bit, !check_bit);
  return !check_bit;
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
      if (reg_patch->physical_offset[d] == patches[i]->physical_offset[d] && reg_patch->physical_size[d] == patches[i]->physical_size[d])
        matches++;
    }

    if (matches == PIDX_MAX_DIMENSIONS)
      return 1;
  }

  return 0;
}


static int pointInChunk(PIDX_patch p, const double *pos)
{
  int contains_point = 1;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    contains_point = contains_point
      && pos[d] >= p->physical_offset[d]
      && pos[d] <= p->physical_offset[d] + p->physical_size[d];

  return contains_point;
}
