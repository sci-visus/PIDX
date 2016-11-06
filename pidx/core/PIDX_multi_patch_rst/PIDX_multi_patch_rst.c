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

static int getPowerOftwo(int x);


PIDX_multi_patch_rst_id PIDX_multi_patch_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, int var_start_index, int var_end_index)
{
  //Creating the restructuring ID
  PIDX_multi_patch_rst_id multi_patch_rst_id;
  multi_patch_rst_id = (PIDX_multi_patch_rst_id)malloc(sizeof (*multi_patch_rst_id));
  memset(multi_patch_rst_id, 0, sizeof (*multi_patch_rst_id));

  multi_patch_rst_id->idx = idx_meta_data;
  multi_patch_rst_id->idx_derived = idx_derived;

  multi_patch_rst_id->group_index = 0;
  multi_patch_rst_id->first_index = var_start_index;
  multi_patch_rst_id->last_index = var_end_index;

  return (multi_patch_rst_id);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_multi_patch_rst_set_communicator(PIDX_multi_patch_rst_id multi_patch_rst_id, MPI_Comm comm)
{
  if (multi_patch_rst_id == NULL)
  return PIDX_err_id;

  multi_patch_rst_id->comm = comm;

  return PIDX_success;
}
#endif



PIDX_return_code PIDX_multi_patch_rst_set_reg_patch_size(PIDX_multi_patch_rst_id rst_id, PIDX_point reg_box)
{
  memcpy(rst_id->reg_patch_size, reg_box, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  return PIDX_success;
}


PIDX_return_code PIDX_multi_patch_rst_set_reg_patch_size_from_bit_string(PIDX_multi_patch_rst_id rst_id)
{
  int ncores;
  MPI_Comm_size(rst_id->comm, &ncores);
  int bits = log2(getPowerOf2(ncores));
  int counter = 1;
  unsigned long long power_two_bound[PIDX_MAX_DIMENSIONS];
  power_two_bound[0] = rst_id->idx_derived->partition_count[0] * rst_id->idx_derived->partition_size[0];
  power_two_bound[1] = rst_id->idx_derived->partition_count[1] * rst_id->idx_derived->partition_size[1];
  power_two_bound[2] = rst_id->idx_derived->partition_count[2] * rst_id->idx_derived->partition_size[2];

  memcpy(rst_id->idx->reg_patch_size, power_two_bound, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  while (bits != 0)
  {
    if (rst_id->idx->bitSequence[counter] == '0')
      rst_id->idx->reg_patch_size[0] = rst_id->idx->reg_patch_size[0] / 2;

    else if (rst_id->idx->bitSequence[counter] == '1')
      rst_id->idx->reg_patch_size[1] = rst_id->idx->reg_patch_size[1] / 2;

    else if (rst_id->idx->bitSequence[counter] == '2')
      rst_id->idx->reg_patch_size[2] = rst_id->idx->reg_patch_size[2] / 2;

    counter++;
    bits--;
  }

  memcpy(rst_id->reg_patch_size, rst_id->idx->reg_patch_size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  return PIDX_success;
}



PIDX_return_code PIDX_multi_patch_rst_auto_set_reg_patch_size(PIDX_multi_patch_rst_id multi_patch_rst_id, int factor)
{
  PIDX_variable_group var_grp = multi_patch_rst_id->idx->variable_grp[multi_patch_rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[multi_patch_rst_id->first_index];

  multi_patch_rst_id->reg_patch_size[0] = factor * getPowerOftwo(var0->sim_patch[0]->size[0]);
  multi_patch_rst_id->reg_patch_size[1] = factor * getPowerOftwo(var0->sim_patch[0]->size[1]);
  multi_patch_rst_id->reg_patch_size[2] = factor * getPowerOftwo(var0->sim_patch[0]->size[2]);

  memcpy(multi_patch_rst_id->idx->reg_patch_size, multi_patch_rst_id->reg_patch_size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);


  return PIDX_success;
}



PIDX_return_code PIDX_multi_patch_rst_finalize(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  free(multi_patch_rst_id);
  multi_patch_rst_id = 0;

  return PIDX_success;
}



/// Function to find the power of 2 of an integer value (example 5->8)
static int getPowerOftwo(int x)
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}
