#include "../PIDX_inc.h"

 // 0 for patch per process and 1 for multi patch per process
static int rst_case_type = 0;
static int cvi = 0;


// Initialiazation and creation of buffers for restructuring phase
// 1. Find if this is a patch per process problem or  multi patch per process problem
// 2. If the restructuring box size is provided then use that else find restructuring box size
PIDX_return_code restructure_setup(PIDX_io file, int gi, int svi, int evi)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  int patch_count = var0->sim_patch_count;
  int max_patch_count = 0;
  PIDX_time time = file->idx_d->time;
  cvi = svi;

  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->comm);
  if (max_patch_count > 1)
    rst_case_type = 1;

  if (rst_case_type == 0)
  {
    time->rst_init_start[cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, svi, evi);

    // attach communicator to the restructuring phase
    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_init_end[cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[cvi] = PIDX_get_time();
    // Creating the metadata to perform retructuring
    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[cvi] = PIDX_get_time();


    time->rst_meta_data_io_start[cvi] = PIDX_get_time();
    // Saving the metadata info needed for reading back the data.
    // Especially when number of cores is different from number of cores
    // used to create the dataset
    ret = PIDX_rst_meta_data_write(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_io_end[cvi] = PIDX_get_time();


    time->rst_buffer_start[cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_rst_aggregate_buf_create(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[cvi] = PIDX_get_time();
  }

  // Case for more than one patch per process
  else
  {
    time->rst_init_start[cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->multi_patch_rst_id = PIDX_multi_patch_rst_init(file->idx, file->idx_d, svi, evi);

    // attach communicator to the restructuring phase
    ret = PIDX_multi_patch_rst_set_communicator(file->multi_patch_rst_id, file->comm);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_init_end[cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[cvi] = PIDX_get_time();
    ret = PIDX_multi_patch_rst_meta_data_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[cvi] = PIDX_get_time();


    time->rst_meta_data_io_start[cvi] = PIDX_get_time();
    // Saving the metadata info needed for reading back the data.
    // Especially when number of cores is different from number of cores
    // used to create the dataset
    ret = PIDX_multi_patch_rst_meta_data_write(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_io_end[cvi] = PIDX_get_time();


    time->rst_buffer_start[cvi] = PIDX_get_time();
    // Creating the buffers required for restructurig
    ret = PIDX_multi_patch_rst_buf_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Aggregating the aligned small buffers after restructuring into one single buffer
    ret = PIDX_multi_patch_rst_aggregate_buf_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_buffer_end[cvi] = PIDX_get_time();
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
        time->rst_write_read_start[cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_rst_staged_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_rst(file->rst_id);
          if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        time->rst_write_read_end[cvi] = PIDX_get_time();
      }


      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      // Aggregating in memory restructured buffers into one large buffer
      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_WRITE);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[cvi] = PIDX_get_time();


      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
    else if (mode == PIDX_READ)
    {
      time->rst_buff_agg_start[cvi] = PIDX_get_time();
      // Perform data restructuring
      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_READ);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_end[cvi] = PIDX_get_time();


      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[cvi] = PIDX_get_time();
        if (file->idx_dbg->debug_rst == 1)
        {
          ret = HELPER_rst(file->rst_id);
          if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
        ret = PIDX_rst_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[cvi] = PIDX_get_time();
      }


      time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
      // Destroying the restructure buffers (as they are now combined into one large buffer)
      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
    }
  }

  // Case for more than one patch per process
  else
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_write_read_start[cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_multi_patch_rst_staged_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[cvi] = PIDX_get_time();


        time->rst_buff_agg_start[cvi] = PIDX_get_time();
        // Aggregating in memory restructured buffers into one large buffer
        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_WRITE);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_end[cvi] = PIDX_get_time();


        time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
        // Destroying the restructure buffers (as they are now combined into one large buffer)
        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_start[cvi] = PIDX_get_time();
        // Aggregating in memory restructured buffers into one large buffer
        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_READ);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_end[cvi] = PIDX_get_time();


        time->rst_write_read_start[cvi] = PIDX_get_time();
        // Perform data restructuring
        ret = PIDX_multi_patch_rst_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_write_read_end[cvi] = PIDX_get_time();


        time->rst_buff_agg_free_start[cvi] = PIDX_get_time();
        // Destroying the restructure buffers (as they are now combined into one large buffer)
        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_free_end[cvi] = PIDX_get_time();
      }
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
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_rst_buf_aggregated_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
        // Read restructured data
        ret = PIDX_rst_buf_aggregated_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
      }
    }
  }

  // Case for more than one patch per process
  else
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
        // Write out restructured data
        ret = PIDX_multi_patch_rst_buf_aggregated_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
      }
    }
    else if (mode == PIDX_READ)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        time->rst_buff_agg_io_start[cvi] = PIDX_get_time();
        // Read restructured data
        ret = PIDX_multi_patch_rst_buf_aggregated_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        time->rst_buff_agg_io_end[cvi] = PIDX_get_time();
      }
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  if (rst_case_type == 0)
  {
    time->rst_cleanup_start[cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_rst_finalize(file->rst_id);
    time->rst_cleanup_end[cvi] = PIDX_get_time();
  }

  // Case for more than one patch per process
  else
  {
    time->rst_cleanup_start[cvi] = PIDX_get_time();
    // Destroy buffers allocated during restructuring phase
    ret = PIDX_multi_patch_rst_aggregate_buf_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_multi_patch_rst_meta_data_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    // Deleting the restructuring ID
    PIDX_multi_patch_rst_finalize(file->multi_patch_rst_id);
    time->rst_cleanup_end[cvi] = PIDX_get_time();
  }

  return PIDX_success;
}


PIDX_return_code restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;

  file->rst_id = PIDX_rst_init(file->idx, file->idx_d, svi, evi);

  ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_rst_forced_raw_read(file->rst_id);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_rst_finalize(file->rst_id);
  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  return PIDX_success;
}
