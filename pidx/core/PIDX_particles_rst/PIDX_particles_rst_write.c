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
static int pointInChunk(PIDX_patch p, const double *pos);

PIDX_return_code PIDX_particles_rst_staged_write(PIDX_particles_rst_id rst_id)
{
  uint64_t req_count = 0;
  uint64_t req_counter = 0;
  int ret = 0;

  MPI_Request *req;
  MPI_Status *status;

  //creating ample requests and statuses
  for (uint64_t i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      req_count++;

  int end_index = rst_id->idx_metadata->variable_count - 1;
  int start_index = 0;
  req_counter = 0;


  req = malloc(sizeof (*req) * req_count * 2 * (end_index - start_index + 1));
  if (!req)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  memset(req, 0, sizeof (*req) * req_count * 2 * (end_index - start_index + 1));

  status = malloc(sizeof (*status) * req_count * 2 * (end_index - start_index + 1));
  if (!status)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  memset(status, 0, sizeof (*status) * req_count * 2 * (end_index - start_index + 1));

  for (uint64_t i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
    {
      for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      {
        for (uint64_t v = start_index; v <= end_index; v++)
        {
          PIDX_variable var = rst_id->idx_metadata->variable[v];
          int length = (var->restructured_super_patch->patch[j]->particle_count) * var->vps * var->bpv/ CHAR_BIT;

          //printf("Receiving %d from %d\n", length, rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank);
          ret = MPI_Irecv(var->restructured_super_patch->patch[j]->buffer, length, MPI_BYTE, rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank, v, rst_id->idx_c->simulation_comm, &req[req_counter]);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
            return PIDX_err_mpi;
          }
          req_counter++;
        }
      }
    }

    for (uint64_t j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
    {
      if (rst_id->idx_c->simulation_rank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
      {
        PIDX_patch reg_patch = malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));
        for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          reg_patch->physical_offset[d] = rst_id->intersected_restructured_super_patch[i]->restructured_patch->physical_offset[d];
          reg_patch->physical_size[d] = rst_id->intersected_restructured_super_patch[i]->restructured_patch->physical_size[d];
        }

        int counter = 0;
        PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];
        const uint64_t bytes_per_pos_v0 = var0->vps * var0->bpv/ CHAR_BIT;
        int p_index = rst_id->intersected_restructured_super_patch[i]->source_patch[j].index;
        for (uint64_t p = 0; p < var0->sim_patch[p_index]->particle_count; ++p)
        {
          if (pointInChunk(reg_patch, (double*)(var0->sim_patch[p_index]->buffer + p * bytes_per_pos_v0)))
          {
            //patch_grp->patch[patch_count]->particle_count++;
            counter++;
          }
        }

        for (int v = start_index; v <= end_index; v++)
        {
          PIDX_variable var = rst_id->idx_metadata->variable[v];
          const uint64_t bytes_per_pos = var->vps * var->bpv/ CHAR_BIT;
          unsigned char *buffer = malloc(counter * bytes_per_pos);
          uint32_t pcount = 0;
          for (uint64_t p = 0; p < var0->sim_patch[p_index]->particle_count; ++p)
          {
            if (pointInChunk(reg_patch, (double*)(var0->sim_patch[p_index]->buffer + p * bytes_per_pos_v0)))
            {
              memcpy(buffer + pcount * bytes_per_pos, var->sim_patch[p_index]->buffer + p * bytes_per_pos, bytes_per_pos);
              double dx, dy, dz;
              memcpy(&dx, buffer + (pcount * 3 + 0) * sizeof(double), sizeof(double));
              memcpy(&dy, buffer + (pcount * 3 + 1) * sizeof(double), sizeof(double));
              memcpy(&dz, buffer + (pcount * 3 + 2) * sizeof(double), sizeof(double));
              //if (v == 0)
              //  printf("Buffer %f %f %f\n", dx, dy, dz);
              pcount++;
            }
          }
          assert(pcount == counter);

          //printf("[%d] Sending %d bytes to %d\n", rst_id->idx_c->simulation_rank, counter * bytes_per_pos, rst_id->intersected_restructured_super_patch[i]->max_patch_rank);
          if (MPI_Isend(buffer, counter * bytes_per_pos, MPI_BYTE, rst_id->intersected_restructured_super_patch[i]->max_patch_rank, v, rst_id->idx_c->simulation_comm, &req[req_counter]) != MPI_SUCCESS)
          {
            fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
            return PIDX_err_mpi;
          }
          free(buffer);
          req_counter++;
        }

        //if (rst_id->intersected_restructured_super_patch[i]->max_patch_rank == 0)
        //  printf("[%d] Counter %d index %d particle count %d -- %f %f %f - %f %f %f\n", rst_id->idx_c->simulation_rank, counter, p_index, var->sim_patch[p_index]->particle_count, reg_patch->physical_offset[0], reg_patch->physical_offset[1], reg_patch->physical_offset[2], reg_patch->physical_size[0], reg_patch->physical_size[1], reg_patch->physical_size[2]);

      }
    }
  }

  //printf("[%d] [%d - %d] Req count %d\n", rst_id->idx_c->simulation_rank, start_index, end_index, req_counter);
  if (MPI_Waitall(req_counter, req, status) != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }


  free(req);
  req = 0;
  free(status);
  status = 0;
  req_counter = 0;

  return PIDX_success;
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
