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
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_rst.h
 *
 */

#include "../../PIDX_inc.h"

PIDX_return_code HELPER_idx_rst(PIDX_idx_rst_id rst_id)
{
  int i, j, k, v = 0, s = 0, n, bytes_for_datatype;
  uint64_t element_count = 0;
  uint64_t lost_element_count = 0;

  float fvalue_1, fvalue_2;
  double dvalue_1, dvalue_2;
  uint64_t uvalue_1, uvalue_2;
  int ivalue_1, ivalue_2;
  int vol = 0;
  uint64_t global_volume;

  uint64_t *bounds = rst_id->idx_metadata->bounds;
  PIDX_variable var0 = rst_id->idx_metadata->variable[rst_id->first_index];

  if (var0->restructured_super_patch_count == 0)
    goto skip_verify;

  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx_metadata->variable[v];
    bytes_for_datatype = var->bpv / 8;

    for (n = 0; n < var->restructured_super_patch->patch_count; n++)
    {
      uint64_t *count_ptr = var->restructured_super_patch->patch[n]->size;
      uint64_t *offset_ptr = var->restructured_super_patch->patch[n]->offset;
      vol = vol + (count_ptr[0] * count_ptr[1] * count_ptr[2]);

      for (k = 0; k < count_ptr[2]; k++)
        for (j = 0; j < count_ptr[1]; j++)
          for (i = 0; i < count_ptr[0]; i++)
          {
            uint64_t index = (count_ptr[0] * count_ptr[1] * k) + (count_ptr[0] * j) + i;
            int check_bit = 1;
            for (s = 0; s < var->vps; s++)
            {
              if (strcmp(var->type_name, PIDX_DType.FLOAT32) == 0)
              {
                fvalue_1 = 100 + v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i;
                memcpy(&fvalue_2, var->restructured_super_patch->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);
                check_bit = check_bit && (fvalue_1 == fvalue_2);
              }
              else if (strcmp(var->type_name, PIDX_DType.FLOAT64) == 0)
              {
                dvalue_1 = 100 + v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_c->color * bounds[0] * bounds[1] * bounds[2]);
                memcpy(&dvalue_2, var->restructured_super_patch->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                check_bit = check_bit && (dvalue_1 == dvalue_2);
              }
              else if (strcmp(var->type_name, PIDX_DType.FLOAT64_RGB) == 0)
              {
                for (s = 0; s < 3; s++)
                {
                  dvalue_1 = 100 + v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_c->color * bounds[0] * bounds[1] * bounds[2]);

                  memcpy(&dvalue_2, var->restructured_super_patch->patch[n]->buffer + ((index * 3) + s) * sizeof(double), sizeof(double));
                  check_bit = check_bit && (dvalue_1  == dvalue_2);
                }
              }
              else if (strcmp(var->type_name, PIDX_DType.UINT64) == 0)
              {
                uvalue_1 = v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_c->color * bounds[0] * bounds[1] * bounds[2]);

                memcpy(&uvalue_2, var->restructured_super_patch->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                check_bit = check_bit && (uvalue_1 == uvalue_2);
              }
              else if (strcmp(var->type_name, PIDX_DType.INT32) == 0)
              {
                ivalue_1 = 100 + v + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_c->color * bounds[0] * bounds[1] * bounds[2]);

                memcpy(&ivalue_2, var->restructured_super_patch->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                check_bit = check_bit && (ivalue_1 == ivalue_2);
              }
            }

            if (check_bit == 0)
            {
              lost_element_count++;
            }
            else
            {
              element_count++;
            }
          }
    }

  }

skip_verify:
  MPI_Allreduce(&element_count, &global_volume, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, rst_id->idx_c->simulation_comm);

  if (global_volume != (uint64_t) bounds[0] * bounds[1] * bounds[2] * (rst_id->last_index - rst_id->first_index + 1))
  {
    if (rst_id->idx_c->simulation_rank == 0)
      fprintf(stderr, "[RST Debug FAILED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", rst_id->idx_c->color, (long long) global_volume, (long long) bounds[0] * bounds[1] * bounds[2]  * (rst_id->last_index - rst_id->first_index + 1));

    if (rst_id->idx_c->simulation_rank == 0)
      fprintf(stderr, "[RST]  Rank %d Color %d [LOST ELEMENT COUNT %lld] [FOUND ELEMENT COUNT %lld] [TOTAL ELEMNTS %lld] [LV %d]\n", rst_id->idx_c->simulation_rank, rst_id->idx_c->color, (long long) lost_element_count, (long long) element_count, (long long) (bounds[0] * bounds[1] * bounds[2]) * (rst_id->last_index - rst_id->first_index + 1), vol);

    return PIDX_err_rst;
  }
  else
    if (rst_id->idx_c->simulation_rank == 0)
      fprintf(stderr, "[RST Debug PASSED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", rst_id->idx_c->color, (long long) global_volume, (long long) bounds[0] * bounds[1] * bounds[2]  * (rst_id->last_index - rst_id->first_index + 1));


  return PIDX_success;
}
