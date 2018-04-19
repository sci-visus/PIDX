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

#define ACTUAL_IO 1
#include "../../PIDX_inc.h"

static int cvi = 0;
static int lgi = 0;

static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file);
static PIDX_return_code populate_restructured_grid(PIDX_io file, int gi, int svi);

static void guess_restructured_box_size(PIDX_io file, int gi, int svi);
static void adjust_restructured_box_size(PIDX_io file);

static PIDX_return_code free_idx_rst_box(PIDX_io file);

// Initialiazation and creation of buffers for restructuring phase
PIDX_return_code idx_restructure_setup(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;
  cvi = svi;
  lgi = gi;

#if 0
  time->set_reg_box_start = MPI_Wtime();
  if (mode == PIDX_WRITE)
  {
    if (set_rst_box_size_for_write(file, gi, svi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  else
  {
    ret = set_rst_box_size_for_read(file, gi, svi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  time->set_reg_box_end = MPI_Wtime();
#endif

  // Initialize the restructuring phase
  time->rst_init_start[lgi][cvi] = PIDX_get_time();
  file->idx_rst_id = PIDX_idx_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
  time->rst_init_end[lgi][cvi] = PIDX_get_time();


  // Populates the relevant meta-data
  time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_idx_rst_meta_data_create(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
  time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();


  // Creating the buffers required for restructurig
  time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_idx_rst_buf_create(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}


  // Aggregating the aligned small buffers after restructuring into one single buffer
  ret = PIDX_idx_rst_aggregate_buf_create(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
  time->rst_buffer_end[lgi][cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code idx_restructure(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Perform data restructuring
      time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_staged_write(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[lgi][cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        ret = HELPER_idx_rst(file->idx_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }

      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregate(file->idx_rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_destroy(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
  }

  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregate(file->idx_rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();

      if (file->idx_dbg->debug_rst == 1)
      {
        ret = HELPER_idx_rst(file->idx_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }

      // Perform data restructuring
      time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_read(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[lgi][cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_destroy(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code idx_restructure_io(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
    {
      // Write out restructured data
      time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregated_write(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
    }
  }
  else if (mode == PIDX_READ)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Read restructured data
      time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_idx_rst_buf_aggregated_read(file->idx_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code idx_restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  // Destroy buffers allocated during restructuring phase
  time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_idx_rst_aggregate_buf_destroy(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_idx_rst_meta_data_destroy(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  // Deleting the restructuring ID
  PIDX_idx_rst_finalize(file->idx_rst_id);
  time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();

  free_idx_rst_box(file);

  return PIDX_success;
}


PIDX_return_code idx_restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;

  file->idx_rst_id = PIDX_idx_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);

  ret = PIDX_idx_rst_forced_raw_read(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_idx_rst_finalize(file->idx_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  return PIDX_success;
}


PIDX_return_code idx_restructure_rst_comm_create(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  MPI_Comm_rank(file->idx_c->simulation_comm, &(file->idx_c->simulation_rank));

  MPI_Comm_split(file->idx_c->simulation_comm, var0->restructured_super_patch_count, file->idx_c->simulation_rank, &(file->idx_c->rst_comm));

  MPI_Comm_rank(file->idx_c->partition_comm, &(file->idx_c->partition_rank));
  MPI_Comm_size(file->idx_c->partition_comm, &(file->idx_c->partition_nprocs));

  MPI_Comm_rank(file->idx_c->rst_comm, &(file->idx_c->rrank));
  MPI_Comm_size(file->idx_c->rst_comm, &(file->idx_c->rnprocs));

  return PIDX_success;
}

PIDX_return_code idx_restructure_copy_rst_comm_to_partition_comm(PIDX_io file, int gi, int svi)
{
  file->idx_c->partition_comm = file->idx_c->rst_comm;

  MPI_Comm_rank(file->idx_c->partition_comm, &(file->idx_c->partition_rank));
  MPI_Comm_size(file->idx_c->partition_comm, &(file->idx_c->partition_nprocs));

  return PIDX_success;
}


PIDX_return_code free_restructured_communicators(PIDX_io file, int gi)
{
  MPI_Comm_free(&(file->idx_c->rst_comm));

  return PIDX_success;
}



PIDX_return_code set_rst_box_size_for_write(PIDX_io file, int gi, int svi)
{
  PIDX_time time = file->idx_d->time;
  time->set_reg_box_start = PIDX_get_time();

  guess_restructured_box_size(file, gi, svi);

  adjust_restructured_box_size(file);

  int ret = populate_restructured_grid(file, gi, svi);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}


static void guess_restructured_box_size(PIDX_io file, int gi, int svi)
{
  int patch_size_x = 0, patch_size_y = 0, patch_size_z = 0;
  int max_patch_size_x = 0, max_patch_size_y = 0, max_patch_size_z = 0;

  if (file->idx->variable_grp[gi]->variable[svi]->sim_patch_count != 0)
  {
    patch_size_x = file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[0];
    patch_size_y = file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[1];
    patch_size_z = file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[2];
  }

  MPI_Allreduce(&patch_size_x, &max_patch_size_x, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_y, &max_patch_size_y, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_z, &max_patch_size_z, 1, MPI_INT, MPI_MAX, file->idx_c->simulation_comm);

  file->idx_d->restructured_grid->patch_size[0] = getPowerOf2(max_patch_size_x);
  file->idx_d->restructured_grid->patch_size[1] = getPowerOf2(max_patch_size_y);
  file->idx_d->restructured_grid->patch_size[2] = getPowerOf2(max_patch_size_z);

  //fprintf(stderr, "Guess %d %d %d\n", file->idx_d->restructured_grid->patch_size[0], file->idx_d->restructured_grid->patch_size[1], file->idx_d->restructured_grid->patch_size[2]);
  return;
}

static void adjust_restructured_box_size(PIDX_io file)
{
  int box_size_factor_x = 1;
  int box_size_factor_y = 1;
  int box_size_factor_z = 1;
  int counter = 0;

  size_t *ps = file->idx_d->restructured_grid->patch_size;

  recaliberate:

  ps[0] = ps[0] * box_size_factor_x;
  ps[1] = ps[1] * box_size_factor_y;
  ps[2] = ps[2] * box_size_factor_z;

  int *tpc = file->idx_d->restructured_grid->total_patch_count;
  tpc[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
  tpc[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
  tpc[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

  //fprintf(stderr, "BB %d %d %d -> tpc %d %d %d\n", file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], tpc[0], tpc[1], tpc[2]);

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



static PIDX_return_code populate_restructured_grid(PIDX_io file, int gi, int svi)
{
  // TODO: cache computation
  int *rgp = file->idx_d->restructured_grid->total_patch_count;
  int total_patch_count = rgp[0] * rgp[1] * rgp[2];

  file->idx_d->restructured_grid->patch = malloc(total_patch_count * sizeof(*file->idx_d->restructured_grid->patch));
  memset(file->idx_d->restructured_grid->patch, 0, (total_patch_count * sizeof(*file->idx_d->restructured_grid->patch)));

  int i = 0, j = 0, k = 0;
  for (i = 0; i < total_patch_count; i++)
  {
    file->idx_d->restructured_grid->patch[i] = malloc(sizeof(*(file->idx_d->restructured_grid->patch[i])));
    memset(file->idx_d->restructured_grid->patch[i], 0, sizeof(*(file->idx_d->restructured_grid->patch[i])));

    file->idx_d->restructured_grid->patch[i]->rank = -1;
  }

  int rank_count = 0;
  int index = 0;
  size_t *ps = file->idx_d->restructured_grid->patch_size;
  for (k = 0; k < file->idx->box_bounds[2]; k = k + ps[2])
    for (j = 0; j < file->idx->box_bounds[1]; j = j + ps[1])
      for (i = 0; i < file->idx->box_bounds[0]; i = i + ps[0])
      {
        //Interior regular patches
        index = ((k / ps[2]) * rgp[0] * rgp[1]) + ((j / ps[1]) * rgp[0]) + (i / ps[0]);

        if (index >= total_patch_count)
        {
          fprintf(stderr, "[%d %d %d] -- %d [%d %d %d] [%d %d %d]\n", rgp[2], rgp[1], rgp[0], index, i, j, k, (int)(k / ps[2]), (int)(j / ps[1]), (int)(i / ps[0]));
        }

        assert(index < total_patch_count);

        Ndim_empty_patch *patch = file->idx_d->restructured_grid->patch;
        patch[index]->offset[0] = i;
        patch[index]->offset[1] = j;
        patch[index]->offset[2] = k;

        patch[index]->size[0] = ps[0];
        patch[index]->size[1] = ps[1];
        patch[index]->size[2] = ps[2];

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

        if (patch[index]->size[0] % file->idx->chunk_size[0] != 0)
          patch[index]->size[0] = (ceil((float)patch[index]->size[0] / file->idx->chunk_size[0])) * file->idx->chunk_size[0];
        if (patch[index]->size[1] % file->idx->chunk_size[1] != 0)
          patch[index]->size[1] = (ceil((float)patch[index]->size[1] / file->idx->chunk_size[1])) * file->idx->chunk_size[1];
        if (patch[index]->size[2] % file->idx->chunk_size[2] != 0)
          patch[index]->size[2] = (ceil((float)patch[index]->size[2] / file->idx->chunk_size[2])) * file->idx->chunk_size[2];

        patch[index]->rank = rank_count * (file->idx_c->simulation_nprocs / (total_patch_count));
        rank_count++;
      }
#if 0
  int nx = 0, ny = 0, nz = 0;
  int int_x = file->idx_c->gnproc_x / rgp[0];
  int int_y = file->idx_c->gnproc_y / rgp[1];
  int int_z = file->idx_c->gnproc_z / rgp[2];

  if (file->idx_c->gnproc_x != -1 && file->idx_c->gnproc_y != -1 && file->idx_c->gnproc_z != -1)
  {
    for (nz = 0; nz < file->idx_c->gnproc_z; nz = nz + int_z)
      for (ny = 0; ny < file->idx_c->gnproc_y; ny = ny + int_y)
        for (nx = 0; nx < file->idx_c->gnproc_x; nx = nx + int_x)
        {
          index = ((nz / int_z) * rgp[0] * rgp[1]) + ((ny / int_y) * rgp[0]) + (nx / int_x);

          patch[index]->rank = (nz * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (ny * file->idx_c->gnproc_x) + nx;

          if (file->idx_c->simulation_rank == patch[index]->rank)
          {
            file->idx_c->simulation_rank_x = nx;
            file->idx_c->simulation_rank_y = ny;
            file->idx_c->simulation_rank_z = nz;
            //fprintf(stderr, "patch[index]->rank %d [%d %d %d]\n", patch[index]->rank, nx, ny, nz);
          }
        }
  }
#endif

#if 0
  if (file->idx_c->simulation_rank == 0)
  {
    for (k = 0; k < rgp[2]; k++)
      for (j = 0; j < rgp[1]; j++)
        for (i = 0; i < rgp[0]; i++)
        {
            index = (k * rgp[0] * rgp[1]) + (j * rgp[0]) + i;
            fprintf(stderr, "Index %d\n", index);
            Ndim_empty_patch *patch = file->idx_d->restructured_grid->patch;
            fprintf(stderr, "patch %d %d %d - %d %d %d R %d\n", patch[index]->offset[0], patch[index]->offset[1], patch[index]->offset[2], patch[index]->size[0], patch[index]->size[1], patch[index]->size[2], patch[index]->rank);
        }
  }
#endif
  assert(rank_count <= total_patch_count);

  return PIDX_success;

}



PIDX_return_code set_rst_box_size_for_read(PIDX_io file, int gi, int svi)
{
  PIDX_time time = file->idx_d->time;
  time->set_reg_box_start = PIDX_get_time();

  set_reg_patch_size_from_bit_string(file);
  populate_restructured_grid(file, gi, svi);

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}


static PIDX_return_code free_idx_rst_box(PIDX_io file)
{
  int *rgp = file->idx_d->restructured_grid->total_patch_count;
  int total_patch_count = rgp[0] * rgp[1] * rgp[2];

  int i = 0;
  for (i = 0; i < total_patch_count; i++)
    free(file->idx_d->restructured_grid->patch[i]);

  free(file->idx_d->restructured_grid->patch);

  return PIDX_success;
}



static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file)
{
  int core = (int)log2(getPowerOf2(file->idx_c->simulation_nprocs));
  int bits;
  int counter = 1;
  size_t power_two_bound[PIDX_MAX_DIMENSIONS];
  power_two_bound[0] = getPowerOf2(file->idx->box_bounds[0]);//file->idx_d->partition_count[0] * file->idx_d->partition_size[0];
  power_two_bound[1] = getPowerOf2(file->idx->box_bounds[1]);//file->idx_d->partition_count[1] * file->idx_d->partition_size[1];
  power_two_bound[2] = getPowerOf2(file->idx->box_bounds[2]);//file->idx_d->partition_count[2] * file->idx_d->partition_size[2];

  //fprintf(stderr, "PTB %d %d %d BB %d %d %d\n", power_two_bound[0], power_two_bound[1], power_two_bound[2], file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2]);

  increase_box_size:
  memcpy(file->idx_d->restructured_grid->patch_size, power_two_bound, sizeof(size_t) * PIDX_MAX_DIMENSIONS);

  bits = core;

#if 1
  counter = 1;
  size_t *ps = file->idx_d->restructured_grid->patch_size;
  int np[3] = {1,1,1};
  while (bits != 0)
  //while (np[0] * np[1] * np[2] < file->idx_c->simulation_nprocs)
  {
    if (file->idx->bitSequence[counter] == '0')
      ps[0] = ps[0] / 2;

    else if (file->idx->bitSequence[counter] == '1')
      ps[1] = ps[1] / 2;

    else if (file->idx->bitSequence[counter] == '2')
      ps[2] = ps[2] / 2;

    //np[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
    //np[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
    //np[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

    //if (file->idx_c->simulation_rank == 0)
    //  fprintf(stderr, "[%c] : %d %d %d -> %d (np: %d) BB %d %d %d\n", file->idx->bitSequence[counter], ps[0], ps[1], ps[2], bits, file->idx_c->simulation_nprocs, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2]);

    counter++;
    bits--;
  }
#endif

  np[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
  np[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
  np[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

  //if (file->idx_c->simulation_rank == 0)
  //  fprintf(stderr, "np %d (%d / %d) %d %d\n", np[0], file->idx->box_bounds[0], ps[0], np[1], np[2]);

  if (np[0] * np[1] * np[2] > file->idx_c->simulation_nprocs)
  {
    core = core - 1;
    goto increase_box_size;
  }

  file->idx_d->restructured_grid->total_patch_count[0] = ceil((float)file->idx->box_bounds[0] / ps[0]);
  file->idx_d->restructured_grid->total_patch_count[1] = ceil((float)file->idx->box_bounds[1] / ps[1]);
  file->idx_d->restructured_grid->total_patch_count[2] = ceil((float)file->idx->box_bounds[2] / ps[2]);

  //fprintf(stderr, "TPC: %d %d %d\n", file->idx_d->restructured_grid->total_patch_count[0], file->idx_d->restructured_grid->total_patch_count[1], file->idx_d->restructured_grid->total_patch_count[2]);

  return PIDX_success;
}
