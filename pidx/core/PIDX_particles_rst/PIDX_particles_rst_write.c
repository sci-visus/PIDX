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
 * \file PIDX_multi_patch_rst.c
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

#if 0
static int pointInChunk(PIDX_patch p, const double *pos);
#endif

PIDX_return_code PIDX_particles_rst_staged_write(PIDX_particles_rst_id rst_id)
{
  uint64_t req_count = 0;
  // creating ample requests and statuses
  for (uint64_t i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      req_count++;

#if 0
  unsigned char ***buffer;
  buffer = malloc(sizeof(*buffer) * rst_id->idx_metadata->variable_count);
#endif

  MPI_Request *req = malloc(sizeof (*req) * req_count * 2 * rst_id->idx_metadata->variable_count);
  MPI_Status *status = malloc(sizeof (*status) * req_count * 2 * rst_id->idx_metadata->variable_count);

  for (uint64_t i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    uint64_t req_counter = 0;

#if 0
    memset(buffer, 0, sizeof(*buffer) * rst_id->idx_metadata->variable_count);
#endif

    // If we're the receiver of this restructured patch, queue up recvs for everyone sending us data
    if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
    {
      for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      {
        for (uint64_t v = 0; v < rst_id->idx_metadata->variable_count; v++)
        {
          PIDX_variable var = rst_id->idx_metadata->variable[v];
          // TODO WILL: How is the length computed? How does the receiver of the restructured patch know
          // the number of particles they're getting?
          const int length = (int) (var->restructured_super_patch->patch[j]->particle_count) * ((var->vps * var->bpv) / CHAR_BIT);

          //printf("[%d] Receiving %d from %d\n", rst_id->idx_c->simulation_rank, length, rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank);
          const int ret = MPI_Irecv(var->restructured_super_patch->patch[j]->buffer, length, MPI_BYTE,
              rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank, 123,
              rst_id->idx_c->simulation_comm, &req[req_counter]);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
            return PIDX_err_mpi;
          }
          req_counter++;
        }
      }
    }

#if 0
    int buffer_count = 0;
    for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
    {
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
        buffer_count++;
    }

    for (int v = 0; v < rst_id->idx_metadata->variable_count; v++)
      buffer[v] = malloc(sizeof(*buffer[v]) * buffer_count);

    buffer_count = 0;
#endif


    // If we're a sender for any patches in this restructured patch, send the data to the owner
    for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
    {
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
      {
        const size_t particles_to_send = rst_id->intersected_restructured_super_patch[i]->patch[j]->particle_count;
#if 0
        PIDX_patch reg_patch = malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));
        for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          reg_patch->physical_offset[d] = rst_id->intersected_restructured_super_patch[i]->restructured_patch->physical_offset[d];
          reg_patch->physical_size[d] = rst_id->intersected_restructured_super_patch[i]->restructured_patch->physical_size[d];
        }

        PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];
        const uint64_t bytes_per_pos_v0 = ((var0->vps * var0->bpv) / CHAR_BIT);
        const int p_index = rst_id->intersected_restructured_super_patch[i]->source_patch[j].index;
        // Optional: Do a validation check that we're writing the right number of particles
        {
          int counter = 0;
          for (uint64_t p = 0; p < var0->sim_patch[p_index]->particle_count; ++p)
          {
            if (pointInChunk(reg_patch, (double*)(var0->sim_patch[p_index]->buffer + p * bytes_per_pos_v0)))
              counter++;
          }
          assert(counter == particles_to_send);
        }

        // Allocate buffers for each var to hold the data we're sending
        for (int v = 0; v < rst_id->idx_metadata->variable_count; v++)
        {
          PIDX_variable var = rst_id->idx_metadata->variable[v];
          const uint64_t bytes_per_var = ((var->vps * var->bpv) / CHAR_BIT);
          buffer[v][buffer_count] = malloc(particles_to_send * bytes_per_var);
        }


        // Populate each variable buffer we're going to send
        uint32_t pcount = 0;
        for (uint64_t p = 0; p < var0->sim_patch[p_index]->particle_count; ++p)
        {
          if (pointInChunk(reg_patch, (double*)(var0->sim_patch[p_index]->buffer + p * bytes_per_pos_v0)))
          {
            for (int v = 0; v < rst_id->idx_metadata->variable_count; v++)
            {
              PIDX_variable var = rst_id->idx_metadata->variable[v];
              const uint64_t bytes_per_var = ((var->vps * var->bpv) / CHAR_BIT);
              memcpy(buffer[v][buffer_count] + pcount * bytes_per_var,
                  var->sim_patch[p_index]->buffer + p * bytes_per_var, bytes_per_var);
            }
            pcount++;
          }
        }
        // Send the populated buffers over to the owner of this restructured patch
        assert(pcount == particles_to_send);
#endif

        for (int v = 0; v < rst_id->idx_metadata->variable_count; v++)
        {
          PIDX_variable var = rst_id->idx_metadata->variable[v];
          const uint64_t bytes_per_var = ((var->vps * var->bpv) / CHAR_BIT);

          //printf("[%d] Sending %d to %d\n", rst_id->idx_c->simulation_rank, (int) particles_to_send * bytes_per_var, rst_id->intersected_restructured_super_patch[i]->max_patch_rank);

          const int ret = MPI_Isend(rst_id->intersected_restructured_super_patch[i]->patch[j]->tbuffer[v], (int) particles_to_send * bytes_per_var, MPI_BYTE,
              rst_id->intersected_restructured_super_patch[i]->max_patch_rank, 123,
              rst_id->idx_c->simulation_comm, &req[req_counter]);
#if 0
          const int ret = MPI_Isend(buffer[v][buffer_count], (int) particles_to_send * bytes_per_var, MPI_BYTE,
              rst_id->intersected_restructured_super_patch[i]->max_patch_rank, 123,
              rst_id->idx_c->simulation_comm, &req[req_counter]);
#endif
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
            return PIDX_err_mpi;
          }
          req_counter++;
        }
#if 0
        free(reg_patch);
        buffer_count++;
#endif
      }
    }

    // Wait for the pending sends/recvs for this restructured patch to complete
    if (MPI_Waitall(req_counter, req, status) != MPI_SUCCESS)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

#if 0
    // Clean up the particle variable buffers
    for (int v = 0; v < rst_id->idx_metadata->variable_count; v++)
    {
      for (uint64_t j = 0; j < buffer_count; j++)
      {
        free(buffer[v][j]);
      }
      free(buffer[v]);
    }
#endif
  }

#if 0
  free(buffer);
#endif

  free(req);
  free(status);

  return PIDX_success;
}

#if 0
static int pointInChunk(PIDX_patch p, const double *pos)
{
  int contains_point = 1;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    contains_point = contains_point
        && pos[d] >= p->physical_offset[d]
        && pos[d] <= p->physical_offset[d] + p->physical_size[d];

  return contains_point;
}
#endif
