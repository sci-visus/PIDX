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


#include "../../PIDX_inc.h"

PIDX_return_code PIDX_rst_write(PIDX_rst_id rst_id)
{
  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(rst_id->comm,  &rank);
#endif

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];

  if (rst_id->idx->enable_rst != 1)
  {
    int v = 0, j = 0, p = 0;
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p];
        for(j = 0; j < patch_group->count; j++)
          memcpy(patch_group->patch[j]->buffer, var->sim_patch[p]->buffer, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps));
      }
    }
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  unsigned long long k1 = 0, i1 = 0, j1 = 0;
  unsigned long long i, j, v, index, count1 = 0, req_count = 0;
  int *send_count, *send_offset;
  unsigned long long send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
  int ret = 0;

  MPI_Request *req;
  MPI_Status *status;
  MPI_Datatype *chunk_data_type;

  for (i = 0; i < rst_id->reg_patch_grp_count; i++)
    for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
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

  chunk_data_type =  malloc(sizeof (*chunk_data_type) * req_count  * (rst_id->last_index - rst_id->first_index + 1));
  if (!chunk_data_type)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(chunk_data_type, 0, sizeof (*chunk_data_type) * req_count  * (rst_id->last_index - rst_id->first_index + 1));

  for (i = 0; i < rst_id->reg_patch_grp_count; i++)
  {
    if (rank == rst_id->reg_patch_grp[i]->max_patch_rank)
    {
      for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
      {
        unsigned long long *reg_patch_offset = rst_id->reg_patch_grp[i]->patch[j]->offset;
        unsigned long long *reg_patch_count  = rst_id->reg_patch_grp[i]->patch[j]->size;

        if(rank == rst_id->reg_patch_grp[i]->source_patch_rank[j])
        {
          count1 = 0;
              for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                  for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                  {
                    unsigned long long *sim_patch_offset = var0->sim_patch[0]->offset;
                    unsigned long long *sim_patch_count = var0->sim_patch[0]->size;

                    index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                            (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                            (i1 - sim_patch_offset[0]);


                    for(v = rst_id->first_index; v <= rst_id->last_index; v++)
                    {
                      PIDX_variable var = var_grp->variable[v];
                      send_o = index * var->vps;
                      send_c = reg_patch_count[0] * var->vps;
#if !SIMULATE_IO
                      memcpy(var->rst_patch_group[counter]->patch[j]->buffer + (count1 * send_c * var->bpv/8), var->sim_patch[0]->buffer + send_o * var->bpv/8, send_c * var->bpv/8);
#endif
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

#if !SIMULATE_IO
            ret = MPI_Irecv(var->rst_patch_group[counter]->patch[j]->buffer, length, MPI_BYTE, rst_id->reg_patch_grp[i]->source_patch_rank[j], 123, rst_id->comm, &req[req_counter]);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
#endif
            req_counter++;
          }
        }
      }
      counter++;
    }
    else
    {
      for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
      {
        if(rank == rst_id->reg_patch_grp[i]->source_patch_rank[j])
        {
          for(v = rst_id->first_index; v <= rst_id->last_index; v++)
          {
            PIDX_variable var = var_grp->variable[v];

            unsigned long long *reg_patch_count = rst_id->reg_patch_grp[i]->patch[j]->size;
            unsigned long long *reg_patch_offset = rst_id->reg_patch_grp[i]->patch[j]->offset;


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
                for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                  for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                    for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                    {
                      unsigned long long *sim_patch_count  = var0->sim_patch[0]->size;
                      unsigned long long *sim_patch_offset = var0->sim_patch[0]->offset;

                      index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                              (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                              (i1 - sim_patch_offset[0]);
                      send_offset[count1] = index * var->vps * var->bpv/8;
                      send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;

                      count1++;
                    }


            //MPI_Datatype chunk_data_type;
            MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type[chunk_counter]);
            MPI_Type_commit(&chunk_data_type[chunk_counter]);

#if !SIMULATE_IO
            ret = MPI_Isend(var->sim_patch[0]->buffer, 1, chunk_data_type[chunk_counter], rst_id->reg_patch_grp[i]->max_patch_rank, 123, rst_id->comm, &req[req_counter]);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
#endif

            req_counter++;
            chunk_counter++;

            free(send_offset);
            free(send_count);
          }
        }
      }
    }
  }

#if !SIMULATE_IO
  ret = MPI_Waitall(req_counter, req, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
#endif

  for (i = 0; i < chunk_counter; i++)
    MPI_Type_free(&chunk_data_type[i]);
  free(chunk_data_type);
  chunk_data_type = 0;

  free(req);
  req = 0;
  free(status);
  status = 0;

  return PIDX_success;
#else
  if (rst_id->idx->enable_rst == 1)
    return PIDX_err_rst;
  else
    return PIDX_success;
#endif
}



PIDX_return_code PIDX_rst_staged_write(PIDX_rst_id rst_id)
{
  int rank = 0;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(rst_id->comm,  &rank);
#endif

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  if (rst_id->idx->enable_rst != 1)
  {
    int v = 0, j = 0, p = 0;
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p];
        for(j = 0; j < patch_group->count; j++)
          memcpy(patch_group->patch[j]->buffer, var->sim_patch[p]->buffer, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps));
      }
    }
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  unsigned long long k1 = 0, i1 = 0, j1 = 0;
  unsigned long long i, j, v, index, count1 = 0, req_count = 0;
  int *send_count, *send_offset;
  unsigned long long send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
  int ret = 0;
  int pipe_length = 0;

  MPI_Request *req;
  MPI_Status *status;
  MPI_Datatype *chunk_data_type;


  //creating ample requests and statuses
  for (i = 0; i < rst_id->reg_patch_grp_count; i++)
    for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
      req_count++;

  int end_index = 0;
  int start_index = 0;
  for (start_index = rst_id->first_index; start_index < (rst_id->last_index + 1); start_index = start_index + pipe_length + 1)
  {
    send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
    end_index = ((start_index + pipe_length) >= (rst_id->last_index + 1)) ? (rst_id->last_index) : (start_index + pipe_length);

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


    for (i = 0; i < rst_id->reg_patch_grp_count; i++)
    {
      if (rank == rst_id->reg_patch_grp[i]->max_patch_rank)
      {
        for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
        {
          unsigned long long *reg_patch_offset = rst_id->reg_patch_grp[i]->patch[j]->offset;
          unsigned long long *reg_patch_count  = rst_id->reg_patch_grp[i]->patch[j]->size;

          if(rank == rst_id->reg_patch_grp[i]->source_patch_rank[j])
          {
            count1 = 0;
                for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                  for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                    for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                    {
                      unsigned long long *sim_patch_offset = var_grp->variable[start_index]->sim_patch[0]->offset;
                      unsigned long long *sim_patch_count = var_grp->variable[start_index]->sim_patch[0]->size;

                      index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                              (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                              (i1 - sim_patch_offset[0]);



                      for(v = start_index; v <= end_index; v++)
                      {
                        PIDX_variable var = var_grp->variable[v];
                        send_o = index * var->vps;
                        send_c = reg_patch_count[0] * var->vps;
#if !SIMULATE_IO
                        memcpy(var->rst_patch_group[counter]->patch[j]->buffer + (count1 * send_c * var->bpv/8), var->sim_patch[0]->buffer + send_o * var->bpv/8, send_c * var->bpv/8);
#endif
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

#if !SIMULATE_IO
              ret = MPI_Irecv(var->rst_patch_group[counter]->patch[j]->buffer, length, MPI_BYTE, rst_id->reg_patch_grp[i]->source_patch_rank[j], 123, rst_id->comm, &req[req_counter]);
              if (ret != MPI_SUCCESS)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
#endif
              req_counter++;
            }
          }
        }
        counter++;
      }
      else
      {
        for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
        {
          if(rank == rst_id->reg_patch_grp[i]->source_patch_rank[j])
          {
            for(v = start_index; v <= end_index; v++)
            {
              PIDX_variable var = var_grp->variable[v];

              unsigned long long *reg_patch_count = rst_id->reg_patch_grp[i]->patch[j]->size;
              unsigned long long *reg_patch_offset = rst_id->reg_patch_grp[i]->patch[j]->offset;

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

                  for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                    for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                      for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                      {
                        unsigned long long *sim_patch_count  = var_grp->variable[start_index]->sim_patch[0]->size;
                        unsigned long long *sim_patch_offset = var_grp->variable[start_index]->sim_patch[0]->offset;

                        index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                                (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                                (i1 - sim_patch_offset[0]);
                        send_offset[count1] = index * var->vps * var->bpv/8;
                        send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;

                        count1++;
                      }


              MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type[chunk_counter]);
              MPI_Type_commit(&chunk_data_type[chunk_counter]);

#if !SIMULATE_IO
              ret = MPI_Isend(var->sim_patch[0]->buffer, 1, chunk_data_type[chunk_counter], rst_id->reg_patch_grp[i]->max_patch_rank, 123, rst_id->comm, &req[req_counter]);
              if (ret != MPI_SUCCESS)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
#endif
              req_counter++;
              chunk_counter++;

              free(send_offset);
              free(send_count);

            }
          }
        }
      }
    }


#if !SIMULATE_IO
    ret = MPI_Waitall(req_counter, req, status);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
#endif

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
#else
  if (rst_id->idx->enable_rst == 1)
    return PIDX_err_rst;
  else
    return PIDX_success;
#endif
}
