#include "../PIDX_inc.h"

static int hz_from_non_shared = 0;
static int hz_to_non_shared = 0;
//static int agg_l_nshared = 0;
static PIDX_return_code find_agg_level(PIDX_io file, int gi);
static PIDX_return_code select_io_mode(PIDX_io file, int gi);
static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);
static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi);

// IDX Write Steps
/********************************************************
*  Step 0: Setup Group and IDX related meta-data        *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform data Restructuring                   *
*  Step 3: Setup HZ encoding Phase                      *
*  Step 4: Perform HZ encoding                          *
*  Step 5: Setup aggregation buffers                    *
*  Step 6: Perform data aggregation                     *
*  Step 7: Perform actual file IO                       *
*  Step 8: cleanup for Steps 1, 3, 5                    *
*                                                       *
*  Step 9: Cleanup the group and IDX related meta-data  *
*********************************************************/

PIDX_return_code PIDX_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    ret = restructure_setup(file, gi, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 2: Perform data restructuring
    ret = restructure(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Setup HZ buffers
    ret = hz_encode_setup(file, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4: Perform HZ encoding
    ret = hz_encode(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (hz_io(file, gi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 5: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, 0, 0, var_grp->agg_l_nshared, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 6: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, 0, 0, var_grp->agg_l_nshared, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 7: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[li] = PIDX_get_time();
      create_async_buffers(file, gi, 0, 0, var_grp->agg_l_nshared);

      ret = data_io(file, gi, li, 0, 0, var_grp->agg_l_nshared, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_async_buffers(file, gi, 0, 0, var_grp->agg_l_nshared);
      finalize_aggregation(file, gi, li, 0, 0, var_grp->agg_l_nshared);
      time->io_end[li] = PIDX_get_time();
    }

    // Step 8: Cleanup all buffers and ids
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = restructure_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 9
  ret = group_meta_data_finalize(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



// IDX Read Steps
/********************************************************
*  Step 0: Setup Group and IDX related meta-data        *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Setup HZ encoding Phase                      *
*  Step 3: Setup aggregation buffers                    *
*  Step 4: Perform actual file IO                       *
*  Step 5: Perform data aggregation                     *
*  Step 6: Perform HZ encoding                          *
*  Step 7: Perform data Restructuring                   *
*  Step 8: cleanup for Steps 1, 3, 5                    *
*                                                       *
*  Step 9: Cleanup the group and IDX related meta-data  *
*********************************************************/

PIDX_return_code PIDX_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_count + 1))
  {
    ei = ((si + file->idx->variable_count) >= (evi)) ? (evi - 1) : (si + file->idx->variable_count);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    ret = restructure_setup(file, gi, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 2: Setup HZ buffers
    ret = hz_encode_setup(file, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (hz_io(file, gi, PIDX_READ) != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, 0, 0, var_grp->agg_l_nshared, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 4: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[li] = PIDX_get_time();
      ret = data_io(file, gi, li, 0, 0, var_grp->agg_l_nshared, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[li] = PIDX_get_time();
    }

    // Step 5: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, 0, 0, var_grp->agg_l_nshared, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      finalize_aggregation(file, gi, li, 0, 0, var_grp->agg_l_nshared);
    }

    // Step 6: Perform HZ encoding
    ret = hz_encode(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 7: Perform data restructuring
    ret = restructure(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 8: Cleanup all buffers and ids
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = restructure_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 9
  ret = group_meta_data_finalize(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


static PIDX_return_code find_agg_level(PIDX_io file, int gi)
{
  int i = 0;
  int no_of_aggregators = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  if (file->idx->enable_agg == 0)
    var_grp->agg_l_nshared = var_grp->nshared_start_layout_index;
  else
  {
    for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->nshared_block_layout_by_level[i - var_grp->nshared_start_layout_index]->efc;
      if (no_of_aggregators <= file->idx_c->gnprocs)
        var_grp->agg_l_nshared = i + 1;
    }
  }

  return PIDX_success;
}


static PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  var_grp->f0_start_layout_index = 0;
  var_grp->f0_end_layout_index = 0;
  var_grp->shared_start_layout_index = 0;
  var_grp->shared_end_layout_index = 0;

  hz_from_non_shared = 0;
  hz_to_non_shared =  idx->maxh;

  return PIDX_success;
}


static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->init_start = PIDX_get_time();
  ret = idx_init(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->init_end = PIDX_get_time();


  time->set_reg_box_start = PIDX_get_time();
  ret = set_rst_box_size(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->set_reg_box_end = MPI_Wtime();


  time->bit_string_start = PIDX_get_time();
  // calculates maxh and bitstring
  ret = populate_bit_string(file, mode);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // selects layout levels based on maxh
  select_io_mode(file, gi);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layout, given this is pure IDX only non-share block layout is populated
  ret = populate_block_layouts(file, gi, svi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared, PIDX_IDX_IO);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Creates the agg and io ids
  ret = create_agg_io_buffer(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, file->idx_d->perm_layout_count);
    file->one_time_initializations = 1;
  }
  time->layout_end = PIDX_get_time();


  if (mode == PIDX_WRITE)
  {
    time->header_io_start = PIDX_get_time();
    // Creates the file heirarchy and writes the header info for all binary files
    ret = write_headers(file, gi, svi, evi, PIDX_IDX_IO);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->header_io_end = PIDX_get_time();
  }

  return PIDX_success;
}


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->group_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = delete_block_layout(file, gi, 0, 0, 0, 0, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  free(var_grp->rank_buffer);
  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
