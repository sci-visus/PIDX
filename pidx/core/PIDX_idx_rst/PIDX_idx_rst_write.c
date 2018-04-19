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

PIDX_return_code PIDX_idx_rst_staged_write(PIDX_idx_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx_metadata->variable_grp[rst_id->group_index];

  unsigned long long k1 = 0, i1 = 0, j1 = 0;
  unsigned long long i, j, v, index, count1 = 0, req_count = 0;
  int *send_count, *send_offset;
  unsigned long long counter = 0, req_counter = 0, chunk_counter = 0;
  size_t send_c;
  off_t send_o;
  int ret = 0;

  MPI_Request *req;
  MPI_Status *status;
  MPI_Datatype *chunk_data_type;

  //creating ample requests and statuses
  for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      req_count++;

  int end_index = 0;
  int start_index = 0;
  for (start_index = rst_id->first_index; start_index < (rst_id->last_index + 1); start_index = start_index + 1)
  {
    send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
    end_index = ((start_index) >= (rst_id->last_index + 1)) ? (rst_id->last_index) : (start_index);

    req = malloc(sizeof (*req) * req_count * 2 * (end_index - start_index + 1));
    if (!req)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(req, 0, sizeof (*req) * req_count * 2 * (end_index - start_index + 1));

    status = malloc(sizeof (*status) * req_count * 2 * (end_index - start_index + 1));
    if (!status)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(status, 0, sizeof (*status) * req_count * 2 * (end_index - start_index + 1));

    chunk_data_type =  malloc(sizeof (*chunk_data_type) * req_count  * (end_index - start_index + 1));
    if (!chunk_data_type)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(chunk_data_type, 0, sizeof (*chunk_data_type) * req_count  * (end_index - start_index + 1));

    for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    {
      if (rst_id->idx_comm_metadata->simulation_rank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
      {
        for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          off_t *reg_patch_offset = rst_id->intersected_restructured_super_patch[i]->patch[j]->offset;
          size_t *reg_patch_count  = rst_id->intersected_restructured_super_patch[i]->patch[j]->size;

          if(rst_id->idx_comm_metadata->simulation_rank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
          {
            count1 = 0;
            int p_index = rst_id->intersected_restructured_super_patch[i]->source_patch[j].index;
            off_t *sim_patch_offset = var_grp->variable[start_index]->sim_patch[p_index]->offset;
            size_t *sim_patch_count = var_grp->variable[start_index]->sim_patch[p_index]->size;

            for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
              for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                {
                  index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                          (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                          (i1 - sim_patch_offset[0]);

                  for(v = start_index; v <= end_index; v++)
                  {
                    PIDX_variable var = var_grp->variable[v];
                    send_o = index * var->vps;
                    send_c = reg_patch_count[0] * var->vps;

                    if (rst_id->idx_debug_metadata->state_dump != PIDX_NO_IO_AND_META_DATA_DUMP)
                      memcpy(var->restructured_super_patch->patch[j]->buffer + (count1 * send_c * var->bpv/8), var->sim_patch[p_index]->buffer + send_o * var->bpv/8, send_c * var->bpv/8);

                    if (rst_id->idx_debug_metadata->state_dump == PIDX_META_DATA_DUMP_ONLY || rst_id->idx_debug_metadata->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
                    {
                      fprintf(rst_id->idx_debug_metadata->local_dump_fp, "[M] [%lld] Dest offset %lld Dest size %lld Source offset %lld Source size %lld\n", v, (unsigned long long)(count1 * send_c * var->bpv/8), (unsigned long long)(send_c * var->bpv/8), (unsigned long long)(send_o * var->bpv/8), (unsigned long long)(send_c * var->bpv/8));
                      fflush(rst_id->idx_debug_metadata->local_dump_fp);
                    }
                  }
                  count1++;
                }
          }
          else
          {
            for(v = start_index; v <= end_index; v++)
            {
              PIDX_variable var = var_grp->variable[v];
              int length = (reg_patch_count[0] * reg_patch_count[1] * reg_patch_count[2]) * var->vps * var->bpv/8;
              if (rst_id->idx_debug_metadata->state_dump != PIDX_NO_IO_AND_META_DATA_DUMP)
              {
                ret = MPI_Irecv(var->restructured_super_patch->patch[j]->buffer, length, MPI_BYTE, rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank, 123, rst_id->idx_comm_metadata->simulation_comm, &req[req_counter]);
                if (ret != MPI_SUCCESS)
                {
                  fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                  return PIDX_err_mpi;
                }
              }

              if (rst_id->idx_debug_metadata->state_dump == PIDX_META_DATA_DUMP_ONLY || rst_id->idx_debug_metadata->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
              {
                fprintf(rst_id->idx_debug_metadata->mpi_dump_fp, "[N REC] [%lld] Dest offset 0 Dest size %d My rank %d Source rank %d\n", v, length, rst_id->idx_comm_metadata->simulation_rank,  rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank);
                fflush(rst_id->idx_debug_metadata->mpi_dump_fp);
              }

              req_counter++;
            }
          }
        }
        counter++;
        assert(counter == 1);
      }
      else
      {
        for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
        {
          if(rst_id->idx_comm_metadata->simulation_rank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
          {
            for(v = start_index; v <= end_index; v++)
            {
              PIDX_variable var = var_grp->variable[v];

              size_t *reg_patch_count = rst_id->intersected_restructured_super_patch[i]->patch[j]->size;
              off_t *reg_patch_offset = rst_id->intersected_restructured_super_patch[i]->patch[j]->offset;

              send_offset = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
              if (!send_offset)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
              memset(send_offset, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

              send_count = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
              if (!send_count)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
              memset(send_count, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

              count1 = 0;
              int p_index =  rst_id->intersected_restructured_super_patch[i]->source_patch[j].index;
              size_t *sim_patch_count  = var_grp->variable[start_index]->sim_patch[p_index]->size;
              off_t *sim_patch_offset = var_grp->variable[start_index]->sim_patch[p_index]->offset;

              int total_send_count = 0;
              for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                  for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                  {
                    index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                            (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                            (i1 - sim_patch_offset[0]);
                    send_offset[count1] = index * var->vps * var->bpv/8;
                    send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;
                    total_send_count = total_send_count + send_count[count1];

                    count1++;
                  }

              MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type[chunk_counter]);
              MPI_Type_commit(&chunk_data_type[chunk_counter]);

              if (rst_id->idx_debug_metadata->state_dump != PIDX_NO_IO_AND_META_DATA_DUMP)
              {
                ret = MPI_Isend(var->sim_patch[p_index]->buffer, 1, chunk_data_type[chunk_counter], rst_id->intersected_restructured_super_patch[i]->max_patch_rank, 123, rst_id->idx_comm_metadata->simulation_comm, &req[req_counter]);
                if (ret != MPI_SUCCESS)
                {
                  fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                  return PIDX_err_mpi;
                }
              }

              if (rst_id->idx_debug_metadata->state_dump == PIDX_META_DATA_DUMP_ONLY || rst_id->idx_debug_metadata->state_dump == PIDX_NO_IO_AND_META_DATA_DUMP)
              {
                fprintf(rst_id->idx_debug_metadata->mpi_dump_fp, "[N SND] [%lld] Source offset 0 Source size %d My rank %d Dest rank %d\n", v, total_send_count, rst_id->idx_comm_metadata->simulation_rank,  rst_id->intersected_restructured_super_patch[i]->max_patch_rank);
                fflush(rst_id->idx_debug_metadata->mpi_dump_fp);
              }

              req_counter++;
              chunk_counter++;

              free(send_offset);
              free(send_count);
            }
          }
        }
      }
    }

    if (rst_id->idx_debug_metadata->state_dump != PIDX_NO_IO_AND_META_DATA_DUMP)
    {
      ret = MPI_Waitall(req_counter, req, status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
        return (-1);
      }
    }

    for (i = 0; i < chunk_counter; i++)
      MPI_Type_free(&chunk_data_type[i]);
    free(chunk_data_type);
    chunk_data_type = 0;

    free(req);
    req = 0;
    free(status);
    status = 0;
    req_counter = 0;
  }

  return PIDX_success;
}
