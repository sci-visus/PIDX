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
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_rst.h
 *
 */

#include "../../PIDX_inc.h"

PIDX_return_code HELPER_rst(PIDX_rst_id rst_id)
{
#if !SIMULATE_IO
  int i, j, k, v = 0, s = 0, m, n, bytes_for_datatype;
  unsigned long long element_count = 0;
  unsigned long long lost_element_count = 0;

  float fvalue_1, fvalue_2;
  double dvalue_1, dvalue_2;
  unsigned long long uvalue_1, uvalue_2;
  int vol = 0;

  unsigned long long *bounds = rst_id->idx->bounds;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    bytes_for_datatype = var->bpv / 8;

    for (m = 0; m < var->patch_group_count; m++)
    {
      for(n = 0; n < var->rst_patch_group[m]->count; n++)
      {
        unsigned long long *count_ptr = var->rst_patch_group[m]->patch[n]->size;
        unsigned long long *offset_ptr = var->rst_patch_group[m]->patch[n]->offset;
        vol = vol + (count_ptr[0] * count_ptr[1] * count_ptr[2]);

        for (k = 0; k < count_ptr[2]; k++)
          for (j = 0; j < count_ptr[1]; j++)
            for (i = 0; i < count_ptr[0]; i++)
            {
              unsigned long long index = (count_ptr[0] * count_ptr[1] * k) + (count_ptr[0] * j) + i;
              int check_bit = 1;
              for (s = 0; s < var->vps; s++)
              {
                if (strcmp(var->type_name, FLOAT32) == 0)
                {
                  fvalue_1 = 100 + v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i;
                  memcpy(&fvalue_2, var->rst_patch_group[m]->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);
                  //fprintf(stderr, "VAL: %f %f\n", fvalue_1, fvalue_2);
                  check_bit = check_bit && (fvalue_1 == fvalue_2);
                }
                else if (strcmp(var->type_name, FLOAT64) == 0)
                {
                  dvalue_1 = 100 + v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_d->color * bounds[0] * bounds[1] * bounds[2]);
                  memcpy(&dvalue_2, var->rst_patch_group[m]->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                  check_bit = check_bit && (dvalue_1 == dvalue_2);
                }
                else if (strcmp(var->type_name, FLOAT64_RGB) == 0)
                {
                  for (s = 0; s < 3; s++)
                  {
                    dvalue_1 = v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_d->color * bounds[0] * bounds[1] * bounds[2]);

                    memcpy(&dvalue_2, var->rst_patch_group[m]->patch[n]->buffer + ((index * 3) + s) * sizeof(double), sizeof(double));
                    check_bit = check_bit && (dvalue_1  == dvalue_2);
                  }
                }
                else if (strcmp(var->type_name, UINT64) == 0)
                {
                  uvalue_1 = v + s + (bounds[0] * bounds[1] * (offset_ptr[2] + k)) + (bounds[0] * (offset_ptr[1] + j)) + offset_ptr[0] + i + ( rst_id->idx_d->color * bounds[0] * bounds[1] * bounds[2]);

                  memcpy(&uvalue_2, var->rst_patch_group[m]->patch[n]->buffer + ((index * var->vps) + s) * bytes_for_datatype, bytes_for_datatype);

                  check_bit = check_bit && (uvalue_1 == uvalue_2);
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
  }

  unsigned long long global_volume;
#if PIDX_HAVE_MPI
  if (rst_id->idx_d->parallel_mode == 1)
    MPI_Allreduce(&element_count, &global_volume, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, rst_id->idx_c->global_comm);
  else
    global_volume = element_count;
#else
  global_volume = element_count;
#endif

  if (global_volume != (unsigned long long) bounds[0] * bounds[1] * bounds[2] * (rst_id->last_index - rst_id->first_index + 1))
  {
    if (rst_id->idx_c->grank == 0)
      fprintf(stderr, "[RST Debug FAILED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", rst_id->idx_d->color, (long long) global_volume, (long long) bounds[0] * bounds[1] * bounds[2]  * (rst_id->last_index - rst_id->first_index + 1));

    if (rst_id->idx_c->grank == 0)
      fprintf(stderr, "[RST]  Rank %d Color %d [LOST ELEMENT COUNT %lld] [FOUND ELEMENT COUNT %lld] [TOTAL ELEMNTS %lld] [LV %d]\n", rst_id->idx_c->grank,  rst_id->idx_d->color, (long long) lost_element_count, (long long) element_count, (long long) (bounds[0] * bounds[1] * bounds[2]) * (rst_id->last_index - rst_id->first_index + 1), vol);

    return PIDX_err_rst;
  }
  else
    if (rst_id->idx_c->grank == 0)
      fprintf(stderr, "[RST Debug PASSED!!!!]  [Color %d] [Recorded Volume %lld] [Actual Volume %lld]\n", rst_id->idx_d->color, (long long) global_volume, (long long) bounds[0] * bounds[1] * bounds[2]  * (rst_id->last_index - rst_id->first_index + 1));
#endif
  return PIDX_success;
}
