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
static int lgi = 0;


// Initialiazation and creation of buffers for restructuring phase
PIDX_return_code particles_restructure_setup(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;
  cvi = svi;
  lgi = gi;

  // Initialize the restructuring phase
  time->rst_init_start[lgi][cvi] = PIDX_get_time();
  file->particles_rst_id = PIDX_particles_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
  time->rst_init_end[lgi][cvi] = PIDX_get_time();


  // Populates the relevant meta-data
  time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_particles_rst_meta_data_create(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
  time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();


  // Saving the metadata info needed for reading back the data.
  // Especially when number of cores is different from number of cores
  // used to create the dataset
  time->rst_meta_data_io_start[lgi][cvi] = PIDX_get_time();
  //
  //if (file->idx->cached_ts == file->idx->current_time_step)
  {
    if (mode == PIDX_WRITE)
    {
      ret = PIDX_particles_rst_meta_data_write(file->particles_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    }
  }
  //
  time->rst_meta_data_io_end[lgi][cvi] = PIDX_get_time();


  // Creating the buffers required for restructurig
  time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_particles_rst_buf_create(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}


  // Aggregating the aligned small buffers after restructuring into one single buffer
  ret = PIDX_particles_rst_aggregate_buf_create(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
  time->rst_buffer_end[lgi][cvi] = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code particles_restructure(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1)
    {
      // Perform data restructuring
      time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_staged_write(file->particles_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[lgi][cvi] = PIDX_get_time();


      // Aggregating in memory restructured buffers into one large buffer
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_buf_aggregate(file->particles_rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_buf_destroy(file->particles_rst_id);
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
      ret = PIDX_particles_rst_buf_aggregate(file->particles_rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();

      // Perform data restructuring
      time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_read(file->particles_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_write_read_end[lgi][cvi] = PIDX_get_time();

      // Destroying the restructure buffers (as they are now combined into one large buffer)
      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_buf_destroy(file->particles_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code particles_restructure_io(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (mode == PIDX_WRITE)
  {
    if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
    {
      // Write out restructured data
      time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
      ret = PIDX_particles_rst_buf_aggregated_write(file->particles_rst_id);
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
      ret = PIDX_particles_rst_buf_aggregated_read(file->particles_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
    }
  }

  return PIDX_success;
}



PIDX_return_code particles_restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  // Destroy buffers allocated during restructuring phase
  time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
  ret = PIDX_particles_rst_aggregate_buf_destroy(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_particles_rst_meta_data_destroy(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  // Deleting the restructuring ID
  PIDX_particles_rst_finalize(file->particles_rst_id);
  time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();

  return PIDX_success;
}


PIDX_return_code particles_restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;

  file->particles_rst_id = PIDX_particles_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);

  ret = PIDX_particles_rst_forced_raw_read(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_particles_rst_finalize(file->particles_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  return PIDX_success;
}




PIDX_return_code particles_set_rst_box_size_for_raw_write(PIDX_io file, int gi, int svi)
{
  PIDX_time time = file->idx_d->time;
  time->set_reg_box_start = PIDX_get_time();

  if (file->idx_d->restructured_grid->patch_size[0] == -1 || file->idx_d->restructured_grid->patch_size[1] == -1 || file->idx_d->restructured_grid->patch_size[2] == -1)
  {
    fprintf(stderr,"Restructuring box needs to be set File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  time->set_reg_box_end = MPI_Wtime();

  return PIDX_success;
}
