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

#define ACTUAL_IO 1
#include "../../PIDX_inc.h"

static int cvi = 0;
static void guess_restructured_box_size(PIDX_io file, int svi);
static void adjust_restructured_box_size(PIDX_io file);
static PIDX_return_code populate_restructured_grid(PIDX_io file);
static PIDX_return_code free_particles_rst_box(PIDX_io file);
static void log_status(char* log_message, int step, int line_number, MPI_Comm comm);


// Initialiazation and creation of buffers for restructuring phase
PIDX_return_code particles_restructure_setup(PIDX_io file, int svi, int evi)
{
  PIDX_time time = file->time;
  cvi = svi;

  // Initialize the restructuring phase
  time->rst_init_start[cvi] = PIDX_get_time();
  file->particles_rst_id = PIDX_particles_rst_init(file->idx, file->idx_c, file->idx_dbg, file->restructured_grid, svi, evi);
  log_status("[Restructuring Step 0]: Init phase\n", 0, __LINE__, file->idx_c->simulation_comm);
  time->rst_init_end[cvi] = PIDX_get_time();


  // Populates the relevant meta-data
  time->rst_meta_data_create_start[cvi] = PIDX_get_time();
  if (PIDX_particles_rst_meta_data_create(file->particles_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  log_status("[Restructuring Step 1]: Metadata create phase\n", 1, __LINE__, file->idx_c->simulation_comm);
  time->rst_meta_data_create_end[cvi] = PIDX_get_time();


  // Creating the buffers required for restructurig
  time->rst_buffer_start[cvi] = PIDX_get_time();
  if (PIDX_particles_rst_buf_create(file->particles_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  log_status("[Restructuring Step 2]: Buffer create phase\n", 2, __LINE__, file->idx_c->simulation_comm);

  // Aggregating the aligned small buffers after restructuring into one single buffer
  if (PIDX_particles_rst_aggregate_buf_create(file->particles_rst_id) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  log_status("[Restructuring Step 3]: Aggregation buffer create phase\n", 3, __LINE__, file->idx_c->simulation_comm);
  time->rst_buffer_end[cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code particles_restructure(PIDX_io file)
{
  PIDX_time time = file->time;

  if (file->idx_dbg->debug_do_rst == 1)
  {
    // Perform data restructuring
    time->rst_write_read_start[cvi] = PIDX_get_time();
    if (PIDX_particles_rst_staged_write(file->particles_rst_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    log_status("[Restructuring Step 4]: Actual restructuring phase\n", 2, __LINE__, file->idx_c->simulation_comm);
    time->rst_write_read_end[cvi] = PIDX_get_time();

    // Aggregating in memory restructured buffers into one large buffer
    time->rst_buff_agg_start[cvi] = PIDX_get_time();
    if (PIDX_particles_rst_buf_aggregate(file->particles_rst_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    log_status("[Restructuring Step 5]: Restructuring buffer aggregation phase\n", 2, __LINE__, file->idx_c->simulation_comm);
    time->rst_buff_agg_end[cvi] = PIDX_get_time();

    // Saving the metadata info needed for reading back the data.
    // Especially when number of cores is different from number of cores
    // used to create the dataset
    time->rst_meta_data_io_start[cvi] = PIDX_get_time();
    if (PIDX_particles_rst_meta_data_write(file->particles_rst_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    log_status("[Restructuring Step 6]: Restructuring meta data create phase\n", 2, __LINE__, file->idx_c->simulation_comm);
    time->rst_meta_data_io_end[cvi] = PIDX_get_time();

    // Destroying the restructure buffers (as they are now combined into one large buffer)
    time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
    if (PIDX_particles_rst_buf_destroy(file->particles_rst_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    log_status("[Restructuring Step 7]: Restructuring buffer destroy phase\n", 2, __LINE__, file->idx_c->simulation_comm);
    time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
  }

  return PIDX_success;
}



PIDX_return_code particles_restructure_io(PIDX_io file)
{
  PIDX_time time = file->time;

  if (file->idx_dbg->debug_do_rst == 1)
  {
    // Write out restructured data
    time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
    if (PIDX_particles_rst_buf_aggregated_write(file->particles_rst_id) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    log_status("[Restructuring Step 8]: Restructuring IO phase\n", 2, __LINE__, file->idx_c->simulation_comm);
    time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
  }

  return PIDX_success;
}



PIDX_return_code particles_restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->time;

  // Destroy buffers allocated during restructuring phase
  time->rst_cleanup_start[cvi] = PIDX_get_time();
  ret = PIDX_particles_rst_aggregate_buf_destroy(file->particles_rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  log_status("[Restructuring Step 9]: Restructuring aggregation buffer destroy\n", 2, __LINE__, file->idx_c->simulation_comm);

  ret = PIDX_particles_rst_meta_data_destroy(file->particles_rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }
  log_status("[Restructuring Step 10]: Restructuring meta data destroy\n", 2, __LINE__, file->idx_c->simulation_comm);

  // Deleting the restructuring ID
  PIDX_particles_rst_finalize(file->particles_rst_id);
  free_particles_rst_box(file);
  log_status("[Restructuring Step 11]: Restructuring finalize\n", 2, __LINE__, file->idx_c->simulation_comm);
  time->rst_cleanup_end[cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code particles_set_rst_box_size_for_raw_write(PIDX_io file, int svi)
{
  PIDX_time time = file->time;
  time->set_reg_box_start = PIDX_get_time();

  guess_restructured_box_size(file, svi);

  adjust_restructured_box_size(file);

  int ret = populate_restructured_grid(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}



static void guess_restructured_box_size(PIDX_io file, int svi)
{
  double max_patch_size_x = file->idx->variable[svi]->sim_patch[0]->physical_size[0];
  double max_patch_size_y = file->idx->variable[svi]->sim_patch[0]->physical_size[1];
  double max_patch_size_z = file->idx->variable[svi]->sim_patch[0]->physical_size[2];

#ifdef PARTICLE_OPTIMIZED
  double patch_size_x = 0, patch_size_y = 0, patch_size_z = 0;
  if (file->idx->variable[svi]->sim_patch_count != 0)
  {
    patch_size_x = file->idx->variable[svi]->sim_patch[0]->physical_size[0];
    patch_size_y = file->idx->variable[svi]->sim_patch[0]->physical_size[1];
    patch_size_z = file->idx->variable[svi]->sim_patch[0]->physical_size[2];
  }

  MPI_Allreduce(&patch_size_x, &max_patch_size_x, 1, MPI_DOUBLE, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_y, &max_patch_size_y, 1, MPI_DOUBLE, MPI_MAX, file->idx_c->simulation_comm);
  MPI_Allreduce(&patch_size_z, &max_patch_size_z, 1, MPI_DOUBLE, MPI_MAX, file->idx_c->simulation_comm);
#endif



  // This scaling factor will set how many process are aggregated to one rank
  file->restructured_grid->physical_patch_size[0] = max_patch_size_x * 2;
  file->restructured_grid->physical_patch_size[1] = max_patch_size_y * 2;
  file->restructured_grid->physical_patch_size[2] = max_patch_size_z * 2;

  return;
}



static void adjust_restructured_box_size(PIDX_io file)
{
  int box_size_factor_x = 1;
  int box_size_factor_y = 1;
  int box_size_factor_z = 1;
  int counter = 0;

  double *ps = file->restructured_grid->physical_patch_size;

  recaliberate:

  ps[0] = ps[0] * box_size_factor_x;
  ps[1] = ps[1] * box_size_factor_y;
  ps[2] = ps[2] * box_size_factor_z;

  uint64_t *tpc = file->restructured_grid->total_patch_count;
  tpc[0] = ceil((float)file->idx->physical_box_bounds[0] / ps[0]);
  tpc[1] = ceil((float)file->idx->physical_box_bounds[1] / ps[1]);
  tpc[2] = ceil((float)file->idx->physical_box_bounds[2] / ps[2]);

  //fprintf(stderr, "BB %f %f %f -> tpc %d %d %d\n", file->idx->physical_box_bounds[0], file->idx->physical_box_bounds[1], file->idx->physical_box_bounds[2], tpc[0], tpc[1], tpc[2]);

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

  int rank_count = 0;
  int index = 0;
  double *ps = file->restructured_grid->physical_patch_size;
  for (double k = 0; k < file->idx->physical_box_bounds[2]; k = k + ps[2])
    for (double j = 0; j < file->idx->physical_box_bounds[1]; j = j + ps[1])
      for (double i = 0; i < file->idx->physical_box_bounds[0]; i = i + ps[0])
      {
        //Interior regular patches
        index = ((k / ps[2]) * rgp[0] * rgp[1]) + ((j / ps[1]) * rgp[0]) + (i / ps[0]);
        assert(index < total_patch_count);

        Ndim_empty_patch *patch = file->restructured_grid->patch;
        patch[index]->physical_offset[0] = i;
        patch[index]->physical_offset[1] = j;
        patch[index]->physical_offset[2] = k;

        patch[index]->physical_size[0] = ps[0];
        patch[index]->physical_size[1] = ps[1];
        patch[index]->physical_size[2] = ps[2];

        patch[index]->is_boundary_patch = 1;


        //Edge regular patches
        if ((i + ps[0]) > file->idx->physical_box_bounds[0])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->physical_size[0] = file->idx->physical_box_bounds[0] - i;
        }

        if ((j + ps[1]) > file->idx->physical_box_bounds[1])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->physical_size[1] = file->idx->physical_box_bounds[1] - j;
        }

        if ((k + ps[2]) > file->idx->physical_box_bounds[2])
        {
          patch[index]->is_boundary_patch = 2;
          patch[index]->physical_size[2] = file->idx->physical_box_bounds[2] - k;
        }

        patch[index]->rank = rank_count * (file->idx_c->simulation_nprocs / (total_patch_count));
        rank_count++;
      }

  assert(rank_count <= total_patch_count);

  return PIDX_success;

}



static PIDX_return_code free_particles_rst_box(PIDX_io file)
{
  uint64_t *rgp = file->restructured_grid->total_patch_count;
  uint64_t total_patch_count = rgp[0] * rgp[1] * rgp[2];

  for (uint64_t i = 0; i < total_patch_count; i++)
    free(file->restructured_grid->patch[i]);

  free(file->restructured_grid->patch);

  return PIDX_success;
}



static void log_status(char* log_message, int step, int line_number, MPI_Comm comm)
{
  /*
  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Barrier(comm);

  if (rank == 0)
    fprintf(stderr, "[nprocs %d] R%d S%d [%d] Log message: %s", size, rank, step, line_number, log_message);
  */
  return;
}
