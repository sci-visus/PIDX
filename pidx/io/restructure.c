#define ACTUAL_IO 1
#include "../PIDX_inc.h"

 // 0 for patch per process and 1 for multi patch per process
static int rst_case_type = 2;
static int cvi = 0;
static int lgi = 0;

static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file);
static PIDX_return_code calculate_rank_mapping(PIDX_io file, int gi, int svi);


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


  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);
  if (max_patch_count > 1)
    rst_case_type = 1;


  // Case for more than one patch per process
  if (rst_case_type == 1)
  {
    time->rst_init_start[lgi][cvi] = PIDX_get_time();
    // Initialize the restructuring phase
    file->multi_patch_rst_id = PIDX_multi_patch_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);
    time->rst_init_end[lgi][cvi] = PIDX_get_time();


    time->rst_meta_data_create_start[lgi][cvi] = PIDX_get_time();
    ret = PIDX_multi_patch_rst_meta_data_create(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    time->rst_meta_data_create_end[lgi][cvi] = PIDX_get_time();


    //if (file->idx->cached_ts == file->idx->current_time_step)
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
      //if (file->idx->cached_ts == file->idx->current_time_step)
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


  return PIDX_success;
}



PIDX_return_code restructure(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  // Case for more than one patch per process
  if (rst_case_type == 1)
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

        //if (file->idx_dbg->debug_rst == 1)
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

  return PIDX_success;
}



PIDX_return_code restructure_io(PIDX_io file, int mode)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  // Case for more than one patch per process
  if (rst_case_type == 1)
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


  return PIDX_success;
}



PIDX_return_code restructure_cleanup(PIDX_io file)
{
  int ret = 0;
  PIDX_time time = file->idx_d->time;

  // Case for more than one patch per process
  if (rst_case_type == 1)
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


  return PIDX_success;
}


PIDX_return_code restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;

  file->generic_rst_id = PIDX_generic_rst_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, svi, evi);

  ret = PIDX_generic_rst_forced_raw_read(file->generic_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  ret = PIDX_generic_rst_finalize(file->generic_rst_id);
  if (ret != PIDX_success) {fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

  return PIDX_success;
}


PIDX_return_code create_restructured_communicators(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  int temp_rank = 0;
  MPI_Comm_rank(file->idx_c->global_comm, &(temp_rank));

  MPI_Comm_split(file->idx_c->global_comm, var0->patch_group_count, file->idx_c->grank, &(file->idx_c->rst_comm));
  file->idx_c->global_comm = file->idx_c->rst_comm;
  file->idx_c->local_comm = file->idx_c->rst_comm;

  MPI_Comm_rank(file->idx_c->global_comm, &(file->idx_c->grank));
  MPI_Comm_rank(file->idx_c->local_comm, &(file->idx_c->lrank));
  MPI_Comm_size(file->idx_c->global_comm, &(file->idx_c->gnprocs));
  MPI_Comm_size(file->idx_c->local_comm, &(file->idx_c->lnprocs));

  //if (var0->patch_group_count != 0)
  //  printf("[%d] -> %d %d\n", file->idx_c->gnprocs, temp_rank, file->idx_c->grank);

  return PIDX_success;
}



PIDX_return_code set_rst_box_size_for_write(PIDX_io file, int gi, int svi)
{
  if (file->idx->reg_patch_size[0] == -1 || file->idx->reg_patch_size[1] == -1 || file->idx->reg_patch_size[2] == -1)
  {
    file->idx->reg_patch_size[0] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[0]);
    file->idx->reg_patch_size[1] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[1]);
    file->idx->reg_patch_size[2] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[2]);
  }

  calculate_rank_mapping(file, gi, svi);

  return PIDX_success;
}



PIDX_return_code set_rst_box_size_for_read(PIDX_io file, int gi, int svi)
{
  set_reg_patch_size_from_bit_string(file);
  calculate_rank_mapping(file, gi, svi);

  return PIDX_success;
}


static PIDX_return_code calculate_rank_mapping(PIDX_io file, int gi, int svi)
{
  int box_size_factor_x = 1;
  int box_size_factor_y = 1;
  int box_size_factor_z = 1;
  int counter = 0;

  int *rgp = file->idx->regridded_process_count;

  recaliberate:

  file->idx->reg_patch_size[0] = file->idx->reg_patch_size[0] * box_size_factor_x;
  file->idx->reg_patch_size[1] = file->idx->reg_patch_size[1] * box_size_factor_y;
  file->idx->reg_patch_size[2] = file->idx->reg_patch_size[2] * box_size_factor_z;

  rgp[0] = ceil((float)file->idx->box_bounds[0] / file->idx->reg_patch_size[0]);
  rgp[1] = ceil((float)file->idx->box_bounds[1] / file->idx->reg_patch_size[1]);
  rgp[2] = ceil((float)file->idx->box_bounds[2] / file->idx->reg_patch_size[2]);

  if (rgp[0] * rgp[1] * rgp[2] > file->idx_c->gnprocs)
  {
    if (counter % 3 == 0)
      box_size_factor_x = box_size_factor_x * 2;
    else if (counter % 3 == 1)
      box_size_factor_y = box_size_factor_y * 2;
    else if (counter % 3 == 2)
      box_size_factor_z = box_size_factor_z * 2;

    goto recaliberate;
  }

  file->idx->regridded_patch = malloc(rgp[0] * rgp[1] * rgp[2] * sizeof(*file->idx->regridded_patch));
  memset(file->idx->regridded_patch, 0, (rgp[0] * rgp[1] * rgp[2] * sizeof(*file->idx->regridded_patch)));

  int i = 0, j = 0, k = 0;
  for (i = 0; i < rgp[0] * rgp[1] * rgp[2]; i++)
  {
    file->idx->regridded_patch[i] = malloc(sizeof(*(file->idx->regridded_patch[i])));
    memset(file->idx->regridded_patch[i], 0, sizeof(*(file->idx->regridded_patch[i])));

    file->idx->regridded_patch[i]->rank = -1;
  }

  int rank_count = 0;
  int index = 0;
  for (k = 0; k < file->idx->box_bounds[2]; k = k + file->idx->reg_patch_size[2])
    for (j = 0; j < file->idx->box_bounds[1]; j = j + file->idx->reg_patch_size[1])
      for (i = 0; i < file->idx->box_bounds[0]; i = i + file->idx->reg_patch_size[0])
      {
        //Interior regular patches
        index = ((k / file->idx->reg_patch_size[2]) * rgp[0] * rgp[1]) + ((j / file->idx->reg_patch_size[1]) * rgp[0]) + (i / file->idx->reg_patch_size[0]);

        if (index >= rgp[0] * rgp[1] * rgp[2])
        {
          fprintf(stderr, "[%d %d %d] -- %d [%d %d %d] [%d %d %d]\n", rgp[2], rgp[1], rgp[0], index, i, j, k, (int)(k / file->idx->reg_patch_size[2]), (int)(j / file->idx->reg_patch_size[1]), (int)(i / file->idx->reg_patch_size[0]));
        }

        assert(index < rgp[0] * rgp[1] * rgp[2]);

        file->idx->regridded_patch[index]->offset[0] = i;
        file->idx->regridded_patch[index]->offset[1] = j;
        file->idx->regridded_patch[index]->offset[2] = k;

        file->idx->regridded_patch[index]->size[0] = file->idx->reg_patch_size[0];
        file->idx->regridded_patch[index]->size[1] = file->idx->reg_patch_size[1];
        file->idx->regridded_patch[index]->size[2] = file->idx->reg_patch_size[2];

        file->idx->regridded_patch[index]->edge = 1;


        //Edge regular patches
        if ((i + file->idx->reg_patch_size[0]) > file->idx->box_bounds[0])
        {
          file->idx->regridded_patch[index]->edge = 2;
          file->idx->regridded_patch[index]->size[0] = file->idx->box_bounds[0] - i;
        }

        if ((j + file->idx->reg_patch_size[1]) > file->idx->box_bounds[1])
        {
          file->idx->regridded_patch[index]->edge = 2;
          file->idx->regridded_patch[index]->size[1] = file->idx->box_bounds[1] - j;
        }

        if ((k + file->idx->reg_patch_size[2]) > file->idx->box_bounds[2])
        {
          file->idx->regridded_patch[index]->edge = 2;
          file->idx->regridded_patch[index]->size[2] = file->idx->box_bounds[2] - k;
        }

        file->idx->regridded_patch[index]->rank = rank_count * (file->idx_c->gnprocs / (rgp[0] * rgp[1] * rgp[2]));
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

          file->idx->regridded_patch[index]->rank = (nz * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (ny * file->idx_c->gnproc_x) + nx;

          if (file->idx_c->grank == file->idx->regridded_patch[index]->rank)
          {
            file->idx_c->grank_x = nx;
            file->idx_c->grank_y = ny;
            file->idx_c->grank_z = nz;
            //fprintf(stderr, "file->idx->regridded_patch[index]->rank %d [%d %d %d]\n", file->idx->regridded_patch[index]->rank, nx, ny, nz);
          }
        }
  }
#endif

  //assert(rank_count == file->idx_c->gnprocs);
  assert(rank_count == rgp[0] * rgp[1] * rgp[2]);

  return PIDX_success;

}


static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file)
{
  int core = log2(getPowerOf2(file->idx_c->gnprocs));
  int bits;
  int counter = 1;
  unsigned long long power_two_bound[PIDX_MAX_DIMENSIONS];
  power_two_bound[0] = file->idx_d->partition_count[0] * file->idx_d->partition_size[0];
  power_two_bound[1] = file->idx_d->partition_count[1] * file->idx_d->partition_size[1];
  power_two_bound[2] = file->idx_d->partition_count[2] * file->idx_d->partition_size[2];

  increase_box_size:

  bits = core;
  counter = 1;

  memcpy(file->idx->reg_patch_size, power_two_bound, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  while (bits != 0)
  {
    if (file->idx->bitSequence[counter] == '0')
      file->idx->reg_patch_size[0] = file->idx->reg_patch_size[0] / 2;

    else if (file->idx->bitSequence[counter] == '1')
      file->idx->reg_patch_size[1] = file->idx->reg_patch_size[1] / 2;

    else if (file->idx->bitSequence[counter] == '2')
      file->idx->reg_patch_size[2] = file->idx->reg_patch_size[2] / 2;

    counter++;
    bits--;
  }

  int np[3];
  np[0] = ceil((float)file->idx->box_bounds[0] / file->idx->reg_patch_size[0]);
  np[1] = ceil((float)file->idx->box_bounds[1] / file->idx->reg_patch_size[1]);
  np[2] = ceil((float)file->idx->box_bounds[2] / file->idx->reg_patch_size[2]);

  if (np[0] * np[1] * np[2] > file->idx_c->gnprocs)
  {
    core = core - 1;
    goto increase_box_size;
  }

  return PIDX_success;
}
