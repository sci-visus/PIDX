#define ACTUAL_IO 1
#include "../PIDX_inc.h"

 // 0 for patch per process and 1 for multi patch per process
static int rst_case_type = 2;
static int cvi = 0;
static int lgi = 0;


// Initialiazation and creation of buffers for restructuring phase
// 1. Find if this is a patch per process problem or  multi patch per process problem
// 2. If the restructuring box size is provided then use that else find restructuring box size
PIDX_return_code restructure_setup(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  int patch_count = var0->sim_patch_count;
  int max_patch_count = 0;
  PIDX_time time = file->idx_d->time;
  cvi = svi;
  lgi = gi;

  // if using analysis then use generic restructuring
  if (file->idx->reg_box_set == PIDX_UNIFORMLY_DISTRIBUTED_BOX)
    rst_case_type = 2;
  else if (file->idx->reg_box_set == PIDX_WAVELET_BOX)
    rst_case_type = 3;
  else
  {
    MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);
    if (max_patch_count > 1)
      rst_case_type = 1;
  }

  if (rst_case_type == 0)
  {
    time->rst_init_start[lgi][cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
    time->rst_init_end[lgi][cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
    // Creating the metadata to perform retructuring
    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();

#if ACTUAL_IO
    if (file->idx->io_type == PIDX_RAW_IO)
    {
      if (file->idx->cached_ts == file->idx->current_time_step)
      {
        if (mode == PIDX_WRITE)
        {
          time->rst_meta_data_io_start[lgi][cvi] = PIDX_get_time();
          // Saving the metadata info needed for reading back the data.
          // Especially when number of cores is different from number of cores
          // used to create the dataset
          ret = PIDX_rst_meta_data_write(file->rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
          time->rst_meta_data_io_end[lgi][cvi] = PIDX_get_time();
        }
      }
    }
#endif


    time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_rst_aggregate_buf_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[lgi][cvi] = PIDX_get_time();
  }

  // Case for more than one patch per process
  else if (rst_case_type == 1)
  {
    time->rst_init_start[lgi][cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->multi_patch_rst_id = PIDX_multi_patch_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
    time->rst_init_end[lgi][cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
    ret = PIDX_multi_patch_rst_meta_data_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();


    if (file->idx->cached_ts == file->idx->current_time_step)
    {
      if (mode == PIDX_WRITE)
      {
        time->rst_meta_data_io_start[lgi][cvi] = PIDX_get_time();
        // Saving the metadata info needed for reading back the data.
        // Especially when number of cores is different from number of cores
        // used to create the dataset
        ret = PIDX_multi_patch_rst_meta_data_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_meta_data_io_end[lgi][cvi] = PIDX_get_time();
      }
    }


    time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_multi_patch_rst_buf_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_multi_patch_rst_aggregate_buf_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[lgi][cvi] = PIDX_get_time();
  }

  // Boxes are distributed uniformly
  else if (rst_case_type == 2)
  {
    time->rst_init_start[lgi][cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->generic_rst_id = PIDX_generic_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
    time->rst_init_end[lgi][cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
    // Creating the metadata to perform retructuring
    ret = PIDX_generic_rst_meta_data_create(file->generic_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();


    if (file->idx->io_type == PIDX_RAW_IO)
    {
      if (file->idx->cached_ts == file->idx->current_time_step)
      {
        if (mode == PIDX_WRITE)
        {
          time->rst_meta_data_io_start[lgi][cvi] = PIDX_get_time();
          // Saving the metadata info needed for reading back the data.
          // Especially when number of cores is different from number of cores
          // used to create the dataset
          ret = PIDX_generic_rst_meta_data_write(file->generic_rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
          time->rst_meta_data_io_end[lgi][cvi] = PIDX_get_time();
        }
      }
    }


    time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_generic_rst_buf_create(file->generic_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_generic_rst_aggregate_buf_create(file->generic_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[lgi][cvi] = PIDX_get_time();
  }

  // Boxes are distributed uniformly for wavelet computation
  else if (rst_case_type == 3)
  {
    time->rst_init_start[lgi][cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->wavelet_rst_id = PIDX_wavelet_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
    time->rst_init_end[lgi][cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
    // Creating the metadata to perform retructuring
    ret = PIDX_wavelet_rst_meta_data_create(file->wavelet_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();

#if 1
    time->rst_buffer_start[lgi][cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_wavelet_rst_buf_create(file->wavelet_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_wavelet_rst_aggregate_buf_create(file->wavelet_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[lgi][cvi] = PIDX_get_time();
#endif
  }

  return PIDX_success;
}



PIDX_return_code restructure(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;
  if (rst_case_type == 0)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_rst_staged_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        //
        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_rst(file->rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        //
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();
      }

#if ACTUAL_IO
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      // Aggregating in memory restructured buffers into one large buffer
      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();
#endif

      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
    else if (mode == PIDX_READ)
    {
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      // Perform data restructuring
      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();


      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_rst(file->rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        ret = PIDX_rst_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();
      }


      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
  }

  // Case for more than one patch per process
  else if (rst_case_type == 1)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_c->grank == 0)
        fprintf(stderr, "Multi patch case\n");

      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_multi_patch_rst_staged_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();


        time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
        // Aggregating in memory restructured buffers into one large buffer
        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_WRITE);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();


        time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
        // Destroying the restructure buffers (as they are now combined into one large buffer)
        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
        // Aggregating in memory restructured buffers into one large buffer
        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_READ);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();


        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_multi_patch_rst_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();


        time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
        // Destroying the restructure buffers (as they are now combined into one large buffer)
        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
      }
    }
  }

  else if (rst_case_type == 2)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_generic_rst_staged_write(file->generic_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_generic_rst(file->generic_rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();
      }

#if ACTUAL_IO
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      // Aggregating in memory restructured buffers into one large buffer
      ret = PIDX_generic_rst_buf_aggregate(file->generic_rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();
#endif

      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_generic_rst_buf_destroy(file->generic_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
    else if (mode == PIDX_READ)
    {
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      // Perform data restructuring
      ret = PIDX_generic_rst_buf_aggregate(file->generic_rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();


      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_generic_rst(file->generic_rst_id);
          if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        ret = PIDX_generic_rst_read(file->generic_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();
      }


      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_generic_rst_buf_destroy(file->generic_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
    }
  }

  else if (rst_case_type == 3)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[lgi][cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_wavelet_rst_staged_write(file->wavelet_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        //if (file->idx_dbg->debug_rst == 1)
        //{
        // ret = HELPER_rst(file->rst_id);
          //if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        //}
        time->rst_write_read_end[lgi][cvi] = PIDX_get_time();
      }
#if 1
#if ACTUAL_IO
      time->rst_buff_agg_start[lgi][cvi] = PIDX_get_time();
      // Aggregating in memory restructured buffers into one large buffer
      ret = PIDX_wavelet_rst_buf_aggregate(file->wavelet_rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[lgi][cvi] = PIDX_get_time();
#endif

      time->rst_buff_agg_free_start[lgi][cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_wavelet_rst_buf_destroy(file->wavelet_rst_id);
      if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[lgi][cvi] = PIDX_get_time();
#endif
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure_io(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (rst_case_type == 0)
  {
#if ACTUAL_IO
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_rst_buf_aggregated_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Read restructured data
        ret = PIDX_rst_buf_aggregated_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
#endif
  }

  // Case for more than one patch per process
  else if (rst_case_type == 1)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_multi_patch_rst_buf_aggregated_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Read restructured data
        ret = PIDX_multi_patch_rst_buf_aggregated_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
  }

  else if (rst_case_type == 2)
  {
#if ACTUAL_IO
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_generic_rst_buf_aggregated_write(file->generic_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Read restructured data
        ret = PIDX_generic_rst_buf_aggregated_read(file->generic_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
#endif
  }

  else if (rst_case_type == 3)
  {
#if ACTUAL_IO
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1 && file->idx_dbg->simulate_rst_io != PIDX_NO_IO_AND_META_DATA_DUMP)
      {
        time->rst_buff_agg_io_start[lgi][cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_wavelet_rst_buf_aggregated_write(file->wavelet_rst_id);
        if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[lgi][cvi] = PIDX_get_time();
      }
    }
#endif
  }

  return PIDX_success;
}



PIDX_return_code restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (rst_case_type == 0)
  {
    time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_rst_finalize(file->rst_id);
    time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();
  }

  // Case for more than one patch per process
  else if (rst_case_type == 1)
  {
    time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_multi_patch_rst_aggregate_buf_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_multi_patch_rst_meta_data_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_multi_patch_rst_finalize(file->multi_patch_rst_id);
    time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();
  }

  else if (rst_case_type == 2)
  {
    time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_generic_rst_aggregate_buf_destroy(file->generic_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_generic_rst_meta_data_destroy(file->generic_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_generic_rst_finalize(file->generic_rst_id);
    time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();
  }

  else if (rst_case_type == 3)
  {
    time->rst_cleanup_start[lgi][cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_wavelet_rst_aggregate_buf_destroy(file->wavelet_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_wavelet_rst_meta_data_destroy(file->wavelet_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_wavelet_rst_finalize(file->wavelet_rst_id);
    time->rst_cleanup_end[lgi][cvi] = PIDX_get_time();
  }

  return PIDX_success;
}


PIDX_return_code restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;

  file->rst_id = PIDX_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);

  ret = PIDX_rst_forced_raw_read(file->rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_rst_finalize(file->rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  return PIDX_success;
}
