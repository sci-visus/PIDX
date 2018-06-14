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
static PIDX_return_code group_meta_data_init(PIDX_io file, int svi, int evi);
static void log_status(char* log_message, int step, int line_number, MPI_Comm comm);

PIDX_return_code PIDX_particle_rst_write(PIDX_io file, int svi, int evi)
{
  log_status("[Particle restructured io Step 0]: Entering\n", 0, __LINE__, file->idx_c->simulation_comm);

  if (particles_set_rst_box_size_for_raw_write(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  log_status("[Particle restructured io Step 1]: Setting restructuring box\n", 1, __LINE__, file->idx_c->simulation_comm);

#if 0
  if (file->idx_c->simulation_rank == 0)
  {
    printf("RST grid size %f %f %f\n", file->restructured_grid->physical_patch_size[0], file->restructured_grid->physical_patch_size[1], file->restructured_grid->physical_patch_size[2]);
    printf("Total patches %d %d %d\n", file->restructured_grid->total_patch_count[0], file->restructured_grid->total_patch_count[1], file->restructured_grid->total_patch_count[2]);

    for (int i = 0; i < file->restructured_grid->total_patch_count[0] * file->restructured_grid->total_patch_count[1] * file->restructured_grid->total_patch_count[2]; i++)
    {
      printf("[%d] Offset %f %f %f Count %f %f %f Rank %d\n", i, file->restructured_grid->patch[i]->physical_offset[0],
                                                         file->restructured_grid->patch[i]->physical_offset[1],
                                                         file->restructured_grid->patch[i]->physical_offset[2],
                                                         file->restructured_grid->patch[i]->physical_size[0],
                                                         file->restructured_grid->patch[i]->physical_size[1],
                                                         file->restructured_grid->patch[i]->physical_size[2],
                                                         file->restructured_grid->patch[i]->rank);
    }
  }
#endif


  if (group_meta_data_init(file, svi, evi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  log_status("[Particle restructured io Step 2]: Creating the folder structure\n", 2, __LINE__, file->idx_c->simulation_comm);

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (int si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    int ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    if (particles_restructure_setup(file, si, ei) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    log_status("[Particle restructured io Step 2]: Particle restructuring setup\n", 3, __LINE__, file->idx_c->simulation_comm);

    // Step 2: Perform data restructuring
    if (particles_restructure(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    log_status("[Particle restructured io Step 2]: Particle restructuring\n", 4, __LINE__, file->idx_c->simulation_comm);

    // Step 3: Write out restructured data
    if (particles_restructure_io(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    log_status("[Particle restructured io Step 2]: Particle restructuring IO\n", 5, __LINE__, file->idx_c->simulation_comm);

    // Step 4: Cleanup all buffers and ids
    if (particles_restructure_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    log_status("[Particle restructured io Step 2]: Particle restructuring cleanup\n", 6, __LINE__, file->idx_c->simulation_comm);
  }

  return PIDX_success;
}



static PIDX_return_code group_meta_data_init(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;

  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  int ret = init_raw_headers_layout(file, svi, evi, file->idx->filename);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();

  return PIDX_success;
}



static void log_status(char* log_message, int step, int line_number, MPI_Comm comm)
{
  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Barrier(comm);

  if (rank == 0)
    fprintf(stderr, "[nprocs %d] R%d S%d [%d] Log message: %s", size, rank, step, line_number, log_message);

  return;
}
