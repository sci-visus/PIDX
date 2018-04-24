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
static PIDX_return_code populate_restructured_grid(PIDX_io file);
static void guess_restructured_box_size(PIDX_io file, int svi);
static void adjust_restructured_box_size(PIDX_io file);
static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file);


PIDX_return_code set_rst_box_size_for_write(PIDX_io file, int svi)
{
  PIDX_time time = file->time;
  time->set_reg_box_start = PIDX_get_time();

  // Guess the initial restrucuting box size
  // Compute the largest box length in each dimension (requires MPI reduce)
  // The restructuring box size is the closest power of two of the largest box
  guess_restructured_box_size(file, svi);

  // The box size needs to be adjusted so that you have total
  // number of box always less than or equal to total number of processes
  adjust_restructured_box_size(file);

  // Assign rank to each of the restructured super patch
  int ret = populate_restructured_grid(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}



PIDX_return_code set_rst_box_size_for_read(PIDX_io file, int svi)
{
  PIDX_time time = file->time;
  time->set_reg_box_start = PIDX_get_time();

  // The bitstring is used to populate the restructuring box size
  // we use log_2(ncores) bits of the bitstring
  // This is the initial guess, we can use more bits to make the restructuring
  // box bigger to make sure that a process holding only one patch at a time.
  set_reg_patch_size_from_bit_string(file);

  populate_restructured_grid(file);

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}


static void guess_restructured_box_size(PIDX_io file, int svi)
{
  int patch_size_x = 0, patch_size_y = 0, patch_size_z = 0;
  int max_patch_size_x = 0, max_patch_size_y = 0, max_patch_size_z = 0;

  // check if the procee at least have on patch
  // TODO for multi patch runs, pick the patch with maximum length
  if (file->idx->variable[svi]->sim_patch_count != 0)
  {
    patch_size_x = file->idx->variable[svi]->sim_patch[0]->size[0];
    patch_size_y = file->idx->variable[svi]->sim_patch[0]->size[1];
    patch_size_z = file->idx->variable[svi]->sim_patch[0]->size[2];
  }

  // compute the maximum box length in each dimension
  MPI_Allreduce(&patch_size_x, &max_patch_size_x, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_y, &max_patch_size_y, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_z, &max_patch_size_z, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);

  // closest power two (since HZ is a multiresolution LOD scheme that works effeciently only with power of two domain datasets)
  file->restructured_grid->patch_size[0] = getPowerOf2(max_patch_size_x);
  file->restructured_grid->patch_size[1] = getPowerOf2(max_patch_size_y);
  file->restructured_grid->patch_size[2] = getPowerOf2(max_patch_size_z);

  return;
}


static void adjust_restructured_box_size(PIDX_io file)
{
  // The logic is to increase (double it) the restructuring box size (super patch size) untill the
  // number of super patches is less than equal to total number of processes.
  // This constrain ensures that not more than one patch per process is assigned crucial for partitioning.

    int box_size_factor_x = 1;
  int box_size_factor_y = 1;
  int box_size_factor_z = 1;
  int counter = 0;

  uint64_t *ps = file->restructured_grid->patch_size;

  recaliberate:

  ps[0] = ps[0] * box_size_factor_x;
  ps[1] = ps[1] * box_size_factor_y;
  ps[2] = ps[2] * box_size_factor_z;

  // total number of super patches
  uint64_t *tpc = file->restructured_grid->total_patch_count;
  tpc[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
  tpc[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
  tpc[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

  if (tpc[0] * tpc[1] * tpc[2] > file->idx_c->simulation_nprocs)
  {
    if (counter % 3 == 0)
      box_size_factor_x = box_size_factor_x * 2;
    else if (counter % 3 == 1)
      box_size_factor_y = box_size_factor_y * 2;
    else if (counter % 3 == 2)
      box_size_factor_z = box_size_factor_z * 2;

    counter++;
    goto recaliberate;
  }

  return;
}



static PIDX_return_code populate_restructured_grid(PIDX_io file)
{
  // TODO: cache computation
  // This function assigns the logical extents (offset and size) of each of the super patch in the restructured grid
  // and also assigns ranks to them. Rank can be assigned in many ways to the super patches (currently using row-major order).

  uint64_t *rgp = file->restructured_grid->total_patch_count;
  uint64_t total_patch_count = rgp[0] * rgp[1] * rgp[2];

  file->restructured_grid->patch = malloc(total_patch_count * sizeof(*file->restructured_grid->patch));
  memset(file->restructured_grid->patch, 0, (total_patch_count * sizeof(*file->restructured_grid->patch)));

  for (uint64_t i = 0; i < total_patch_count; i++)
  {
    file->restructured_grid->patch[i] = malloc(sizeof(*(file->restructured_grid->patch[i])));
    memset(file->restructured_grid->patch[i], 0, sizeof(*(file->restructured_grid->patch[i])));

    file->restructured_grid->patch[i]->rank = -1;
  }

  uint32_t rank_count = 0;
  uint64_t index = 0;
  uint64_t *ps = file->restructured_grid->patch_size;
  for (uint64_t k = 0; k < file->idx->box_bounds[2]; k = k + ps[2])
    for (uint64_t j = 0; j < file->idx->box_bounds[1]; j = j + ps[1])
      for (uint64_t i = 0; i < file->idx->box_bounds[0]; i = i + ps[0])
      {
        // access super patches is row order
        index = ((k / ps[2]) * rgp[0] * rgp[1]) + ((j / ps[1]) * rgp[0]) + (i / ps[0]);
        assert(index < total_patch_count);

        // assign offset and size of the super patches
        Ndim_empty_patch *patch = file->restructured_grid->patch;
        patch[index]->offset[0] = i;
        patch[index]->offset[1] = j;
        patch[index]->offset[2] = k;

        patch[index]->size[0] = ps[0];
        patch[index]->size[1] = ps[1];
        patch[index]->size[2] = ps[2];

        // default to non boundary (1 is for non-boundary and 2 is boundary)
        patch[index]->is_boundary_patch = 1;

        //Edge regular patches
        if ((i + ps[0]) > file->idx->box_bounds[0])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->size[0] = file->idx->box_bounds[0] - i;
        }

        if ((j + ps[1]) > file->idx->box_bounds[1])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->size[1] = file->idx->box_bounds[1] - j;
        }

        if ((k + ps[2]) > file->idx->box_bounds[2])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->size[2] = file->idx->box_bounds[2] - k;
        }

        // adjust the patch sizes if zfp compression is applied
        if (patch[index]->size[0] % file->idx->chunk_size[0] != 0)
          patch[index]->size[0] = (ceil((float)patch[index]->size[0] / file->idx->chunk_size[0])) * file->idx->chunk_size[0];
        if (patch[index]->size[1] % file->idx->chunk_size[1] != 0)
          patch[index]->size[1] = (ceil((float)patch[index]->size[1] / file->idx->chunk_size[1])) * file->idx->chunk_size[1];
        if (patch[index]->size[2] % file->idx->chunk_size[2] != 0)
          patch[index]->size[2] = (ceil((float)patch[index]->size[2] / file->idx->chunk_size[2])) * file->idx->chunk_size[2];

        // assign rank to super patch
        patch[index]->rank = rank_count * (file->idx_c->simulation_nprocs / (total_patch_count));
        rank_count++;
      }

#if 0
  if (file->idx_c->simulation_rank == 0)
  {
    for (k = 0; k < rgp[2]; k++)
      for (j = 0; j < rgp[1]; j++)
        for (i = 0; i < rgp[0]; i++)
        {
            index = (k * rgp[0] * rgp[1]) + (j * rgp[0]) + i;
            fprintf(stderr, "Index %d\n", index);
            Ndim_empty_patch *patch = file->restructured_grid->patch;
            fprintf(stderr, "patch %d %d %d - %d %d %d R %d\n", patch[index]->offset[0], patch[index]->offset[1], patch[index]->offset[2], patch[index]->size[0], patch[index]->size[1], patch[index]->size[2], patch[index]->rank);
        }
  }
#endif
  assert(rank_count <= total_patch_count);

  return PIDX_success;
}


static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file)
{
  // If we do not follow the bitstring to compute the restructuring box size, then we will
  // run into interleaved read accesses during file io, while using a fewer processes than what
  // was used to write the idx file.
  // With this scheme we ensure that aggregators access large contiguous chunks of data.

  // initially use log_2(cores) bits to estimate the size of the restructuring box
  uint32_t bits = (int)log2(getPowerOf2(file->idx_c->simulation_nprocs));
  uint32_t counter = 1;

  uint64_t power_two_bound[PIDX_MAX_DIMENSIONS];
  power_two_bound[0] = getPowerOf2(file->idx->box_bounds[0]);
  power_two_bound[1] = getPowerOf2(file->idx->box_bounds[1]);
  power_two_bound[2] = getPowerOf2(file->idx->box_bounds[2]);

  increase_box_size:
  // initialize the restructuring box size to the bounds of the data
  // and then iteratively half the resolution in each dimension following the bit sequence
  // of the "bits" number of bits of the bitsequence
  memcpy(file->restructured_grid->patch_size, power_two_bound, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);

  uint32_t  temp_bits = bits;

  counter = 1;
  uint64_t *ps = file->restructured_grid->patch_size;
  while (temp_bits != 0)
  {
    if (file->idx->bitSequence[counter] == '0')
      ps[0] = ps[0] / 2;

    else if (file->idx->bitSequence[counter] == '1')
      ps[1] = ps[1] / 2;

    else if (file->idx->bitSequence[counter] == '2')
      ps[2] = ps[2] / 2;

    counter++;
    temp_bits--;
  }

  file->restructured_grid->total_patch_count[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
  file->restructured_grid->total_patch_count[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
  file->restructured_grid->total_patch_count[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

  // making sure that the total number of restructured super patches are not more than the total number of processes
  // if total number of patches are more than to total number of patches, use one more bit of the bit string, make the
  // restructuring box size larger, creating fewer super patches.
  if (file->restructured_grid->total_patch_count[0] * file->restructured_grid->total_patch_count[1] * file->restructured_grid->total_patch_count[2] > file->idx_c->simulation_nprocs)
  {
    bits = bits - 1;
    goto increase_box_size;
  }

  return PIDX_success;
}
