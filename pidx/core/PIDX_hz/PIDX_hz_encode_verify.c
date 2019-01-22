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


#include "../../PIDX_inc.h"


PIDX_return_code HELPER_Hz_encode(PIDX_hz_encode_id id)
{
  int i = 0, k = 0, v = 0;
  uint64_t global_hz, element_count = 0, lost_element_count = 0;
  uint64_t ZYX[PIDX_MAX_DIMENSIONS];
  int check_bit = 1, s = 0;
  double dvalue_1 = 0, dvalue_2 = 0;
  float fvalue_1 = 0, fvalue_2 = 0;
  uint64_t uvalue_1 = 0, uvalue_2 = 0;

  for (v = id->first_index; v <= id->last_index; v++)
  {
    PIDX_variable var = id->idx->variable[v];

    for (i = 0; i < id->idx->maxh; i++)
    {
      if (var->hz_buffer->nsamples_per_level[i][0] * var->hz_buffer->nsamples_per_level[i][1] * var->hz_buffer->nsamples_per_level[i][2] != 0)
      {
        for (k = 0; k <= (var->hz_buffer->end_hz_index[i] - var->hz_buffer->start_hz_index[i]) * 1; k++)
        {
          global_hz = var->hz_buffer->start_hz_index[i] + k;
          Hz_to_xyz(id->idx->bitPattern, id->idx->maxh - 1, global_hz, ZYX);
          if ((ZYX[0] < id->idx->box_bounds[0] && ZYX[1] < id->idx->box_bounds[1] && ZYX[2] < id->idx->box_bounds[2]))
          {
            check_bit = 1, s = 0;
            if (strcmp(var->type_name, PIDX_DType.FLOAT64) == 0)
            {
              dvalue_1 = 100 + v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0];
              dvalue_2 = *(*((double**)var->hz_buffer->buffer + i) + ((k * var->vps) + s));

              //fprintf(stderr, "X[%lld]  %f %f\n", (long long)global_hz, dvalue_1, dvalue_2);
              check_bit = check_bit && (dvalue_1  == dvalue_2);
            }
            else if (strcmp(var->type_name, PIDX_DType.FLOAT32) == 0)
            {
              fvalue_1 = 100 + v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0];
              fvalue_2 = *(*((float**)var->hz_buffer->buffer + i) + ((k * var->vps) + s));

              //if (file->idx_c->partition_rank == 0)
              //fprintf(stderr, "%f %f\n", fvalue_1, fvalue_2);
              check_bit = check_bit && (fvalue_1 == fvalue_2);
            }
            else if (strcmp(var->type_name, PIDX_DType.FLOAT64_RGB) == 0)
            {
              for (s = 0; s < 3; s++)
              {
                dvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0];
                memcpy(&dvalue_2, var->hz_buffer->buffer[i] + ((k * 3) + s) * sizeof(double), sizeof(double));
                //fprintf(stderr, "Y[%d] %f -- %f\n", v, dvalue_1, dvalue_2);

                check_bit = check_bit && (dvalue_1  == dvalue_2);
              }
            }
            if (strcmp(var->type_name, PIDX_DType.UINT64) == 0)
            {
              uvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_c->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
              uvalue_2 = *(*((uint64_t**)var->hz_buffer->buffer + i) + ((k * var->vps) + s));
              check_bit = check_bit && (uvalue_1  == uvalue_2);
            }
            else if (strcmp(var->type_name, PIDX_DType.UINT64_RGB) == 0)
            {
              for (s = 0; s < 3; s++)
              {
                uvalue_1 = v + s + (id->idx->bounds[0] * id->idx->bounds[1]*(ZYX[2]))+(id->idx->bounds[0]*(ZYX[1])) + ZYX[0] + (id->idx_c->color * id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]);
                memcpy(&uvalue_2, var->hz_buffer->buffer[i] + ((k * 3) + s) * sizeof(double), sizeof(double));

                check_bit = check_bit && (uvalue_1  == uvalue_2);
              }
            }

            if (check_bit == 1 && fvalue_1 != 0 && fvalue_2 != 0)
              element_count++;

            else
              lost_element_count++;
          }
        }
      }
    }
  }

  uint64_t global_volume = 0;
  MPI_Allreduce(&element_count, &global_volume, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, id->idx_c->partition_comm);

  //fprintf(stderr, "[HZ] Volume [%lld] and Volume [%lld]\n", global_volume, (uint64_t)(id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1)));

  if (global_volume != (uint64_t) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1))
  {
    if (id->idx_c->partition_rank == 0)
      fprintf(stderr, "[HZ Debug FAILED!!!!] [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", id->idx_c->color, (long long) global_volume, (long long) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1));

    if (id->idx_c->partition_rank == 0)
      fprintf(stderr, "[HZ]  file->idx_c->partition_rank %d Color %d [LOST ELEMENT COUNT %lld] [FOUND ELEMENT COUNT %lld] [TOTAL ELEMNTS %lld] \n", id->idx_c->partition_rank,  id->idx_c->color, (long long) lost_element_count, (long long) element_count, (long long) (id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2]) * (id->last_index - id->first_index + 1));

    return PIDX_err_hz;
  }
  else
  {
    if (id->idx_c->partition_rank == 0)
      fprintf(stderr, "[HZ Debug PASSED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld] [Lost Elemet COunt %lld]\n", id->idx_c->color, (long long) global_volume, (long long) id->idx->bounds[0] * id->idx->bounds[1] * id->idx->bounds[2] * (id->last_index - id->first_index + 1), (long long) lost_element_count);
  }

  return PIDX_success;
}
