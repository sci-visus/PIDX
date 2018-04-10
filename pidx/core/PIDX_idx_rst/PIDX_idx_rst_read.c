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


PIDX_return_code PIDX_idx_rst_read(PIDX_idx_rst_id rst_id)
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

  //fprintf(stderr, "rst_id->intersected_restructured_super_patch_count = %d\n", rst_id->intersected_restructured_super_patch_count);
  for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
    for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      req_count++;

  //creating ample requests and statuses
  req = (MPI_Request*) malloc(sizeof (*req) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));
  if (!req)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(req, 0, sizeof (*req) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));

  status = (MPI_Status*) malloc(sizeof (*status) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));
  if (!status)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(status, 0, sizeof (*status) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));

  chunk_data_type = malloc(sizeof (*chunk_data_type) * req_count * (rst_id->last_index - rst_id->first_index + 1));
  if (!chunk_data_type)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(chunk_data_type, 0, sizeof (*chunk_data_type) * req_count * (rst_id->last_index - rst_id->first_index + 1));

  for (i = 0; i < rst_id->intersected_restructured_super_patch_count; i++)
  {
    if (rst_id->idx_comm_metadata->grank == rst_id->intersected_restructured_super_patch[i]->max_patch_rank)
    {
      for(j = 0; j < rst_id->intersected_restructured_super_patch[i]->patch_count; j++)
      {
        off_t *reg_patch_offset = rst_id->intersected_restructured_super_patch[i]->patch[j]->offset;
        size_t *reg_patch_count  = rst_id->intersected_restructured_super_patch[i]->patch[j]->size;

        if(rst_id->idx_comm_metadata->grank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
        {
          count1 = 0;
          int p_index = rst_id->intersected_restructured_super_patch[i]->source_patch[j].index;

          off_t *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->offset;
          size_t *sim_patch_count = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->size;

          for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
            for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
              for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
              {
                index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                    (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                    (i1 - sim_patch_offset[0]);

                for(v = rst_id->first_index; v <= rst_id->last_index; v++)
                {
                  PIDX_variable var = var_grp->variable[v];
                  send_o = index * var->vps;
                  send_c = reg_patch_count[0] * var->vps;

                  memcpy(var->sim_patch[p_index]->buffer + send_o * var->bpv/8, var->restructured_super_patch->patch[j]->buffer + (count1 * send_c * var->bpv/8), send_c * var->bpv/8);

                  //double a1;
                  //memcpy(&a1, var->restructured_super_patch->patch[j]->buffer + (count1 * send_c * var->bpv/8), sizeof(double));
                  //fprintf(stderr, "a1 = %f\n", a1);
                }
                count1++;
              }
        }
        else
        {
          for(v = rst_id->first_index; v <= rst_id->last_index; v++)
          {
            PIDX_variable var = var_grp->variable[v];

            int length = (reg_patch_count[0] * reg_patch_count[1] * reg_patch_count[2]) * var->vps * var->bpv/8;

            ret = MPI_Isend(var->restructured_super_patch->patch[j]->buffer, length, MPI_BYTE, rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank, 123, rst_id->idx_comm_metadata->global_comm, &req[req_counter]);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
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
        if(rst_id->idx_comm_metadata->grank == rst_id->intersected_restructured_super_patch[i]->source_patch[j].rank)
        {
          for(v = rst_id->first_index; v <= rst_id->last_index; v++)
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

            size_t *sim_patch_count  = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->size;
            off_t *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->offset;

            for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
              for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                {

                  index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                      (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                      (i1 - sim_patch_offset[0]);
                  send_offset[count1] = index * var->vps * var->bpv/8;
                  send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;

                  count1++;
                }

            MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type[chunk_counter]);
            MPI_Type_commit(&chunk_data_type[chunk_counter]);

            ret = MPI_Irecv(var->sim_patch[p_index]->buffer, 1, chunk_data_type[chunk_counter], rst_id->intersected_restructured_super_patch[i]->max_patch_rank, 123, rst_id->idx_comm_metadata->global_comm, &req[req_counter]);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
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

  ret = MPI_Waitall(req_counter, req, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_mpi;
  }

  for (i = 0; i < chunk_counter; i++)
    MPI_Type_free(&chunk_data_type[i]);
  free(chunk_data_type);
  chunk_data_type = 0;

  free(req);
  req = 0;
  free(status);
  status = 0;

  return PIDX_success;
}
