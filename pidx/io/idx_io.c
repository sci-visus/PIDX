#include "../PIDX_inc.h"

static int hz_from_non_shared = 0, hz_from_shared = 0, hz_from_file_zero = 0;
static int hz_to_non_shared = 0, hz_to_shared = 0, hz_to_file_zero = 0;
static PIDX_return_code find_agg_level(PIDX_io file, int gi);
static PIDX_return_code select_io_mode(PIDX_io file, int gi);
static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);
static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);

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

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 1: Setup restructuring buffers
  ret = restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2: Perform data restructuring
  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    //ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);

    if ((si + file->idx->variable_pipe_length) >= evi)
      file->idx->variable_pipe_length = evi - (si + 1);

    ei = si + file->idx->variable_pipe_length;

    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    //fprintf(stderr, "[%d %d] ---- %d %d pipe %d\n", svi, evi, si, ei, file->idx->variable_pipe_length);

    // Step 3: Setup HZ buffers
    ret = hz_encode_setup(file, gi, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4: Perform HZ encoding
    ret = hz_encode(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (hz_io(file, gi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


#if 1
    // Step 5: Setup aggregation buffers
    //for (li = si; li <= ei; li = li + 1)
    //{
      ret = data_aggregate(file, gi, si, li, ei, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    //}

    // Step 6: Performs data aggregation
    //for (li = si; li <= ei; li = li + 1)
    //{
      ret = data_aggregate(file, gi, si, li, ei, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    //}


    // Step 7: Performs actual file io
    //for (li = si; li <= ei; li = li + 1)
    //{
      time->io_start[gi][li] = PIDX_get_time();
      create_async_buffers(file, gi);

      ret = data_io(file, gi, si, li, ei, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_async_buffers(file, gi);
      finalize_aggregation(file, gi, li, si);
      time->io_end[gi][li] = PIDX_get_time();
    //}
#else

    int aei = 0;
    int agg_pipe_length = 2;
    for (li = si; li <= ei; li = li + (agg_pipe_length + 1))
    {
      aei = ((li + agg_pipe_length) >= (ei)) ? (ei) : (li + agg_pipe_length);
      ret = data_aggregate(file, gi, si, li, aei, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 6: Performs data aggregation
    for (li = si; li <= ei; li = li + (agg_pipe_length + 1))
    {
      aei = ((li + agg_pipe_length) >= (ei)) ? (ei) : (li + agg_pipe_length);
      ret = data_aggregate(file, gi, si, li, ei, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }


    // Step 7: Performs actual file io
    for (li = si; li <= ei; li = li + (agg_pipe_length + 1))
    {
      aei = ((li + agg_pipe_length) >= (ei)) ? (ei) : (li + agg_pipe_length);
      time->io_start[gi][li] = PIDX_get_time();
      create_async_buffers(file, gi);

      ret = data_io(file, gi, si, li, ei, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_async_buffers(file, gi);
      finalize_aggregation(file, gi, li, si);
      time->io_end[gi][li] = PIDX_get_time();
    }

#endif

    // Step 8: Cleanup all buffers and ids
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  ret = restructure_cleanup(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 9
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
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

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_count + 1))
  {
    ei = ((si + file->idx->variable_count) >= (evi)) ? (evi - 1) : (si + file->idx->variable_count);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    ret = restructure_setup(file, gi, si, ei, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 2: Setup HZ buffers
    ret = hz_encode_setup(file, gi, si, ei);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    //PIDX_variable_group var_grp = file->idx->variable_grp[gi];
    //printf("X %d %d %d Y %d %d %d\n", var_grp->nshared_start_layout_index, var_grp->agg_l_nshared, var_grp->nshared_end_layout_index, var_grp->shared_start_layout_index, var_grp->agg_l_shared, var_grp->shared_end_layout_index);


    if (hz_io(file, gi, PIDX_READ) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, si, li, li, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 4: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[gi][li] = PIDX_get_time();
      ret = data_io(file, gi, si, li, ei, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[gi][li] = PIDX_get_time();
    }

    // Step 5: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, si, li, li, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      finalize_aggregation(file, gi, li, si);
    }

    // Step 6: Perform HZ encoding
    ret = hz_encode(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 7: Perform data restructuring
    ret = restructure(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 8: Cleanup all buffers and ids
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = restructure_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 9
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}

#if 0
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

  if (file->idx_c->grank == 0)
    fprintf(stderr, "%d - %d %d %d\n", var_grp->agg_l_nshared, var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index, var_grp->nshared_layout_count);

  return PIDX_success;
}


static PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  var_grp->agg_l_f0 = 0;
  var_grp->f0_start_layout_index = 0;
  var_grp->f0_end_layout_index = 0;
  var_grp->f0_layout_count = 0;

  var_grp->agg_l_shared = 0;
  var_grp->shared_start_layout_index = 0;
  var_grp->shared_end_layout_index = 0;
  var_grp->shared_layout_count = 0;

  hz_from_non_shared = 0;
  hz_to_non_shared =  idx->maxh;

  return PIDX_success;
}
#else

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
      if (no_of_aggregators <= file->idx_c->lnprocs)
        var_grp->agg_l_nshared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    var_grp->agg_l_shared = var_grp->shared_start_layout_index;
  else
  {
    for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->shared_block_layout_by_level[i - var_grp->shared_start_layout_index]->efc;
      if (no_of_aggregators <= file->idx_c->lnprocs)
        var_grp->agg_l_shared = i + 1;
    }
  }

  if (file->idx->enable_agg == 0)
    var_grp->agg_l_f0 = var_grp->f0_start_layout_index;
  else
  {
    for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
    {
      no_of_aggregators = var_grp->f0_block_layout_by_level[i - var_grp->f0_start_layout_index]->efc;
      if (no_of_aggregators <= file->idx_c->lnprocs)
        var_grp->agg_l_f0 = i + 1;
    }
  }

  return PIDX_success;
}


static PIDX_return_code select_io_mode(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  hz_from_file_zero = 0;
  hz_to_file_zero =  0;

  hz_from_shared = 0;
  hz_to_shared = idx->total_partiton_level;

  hz_from_non_shared = idx->total_partiton_level;
  hz_to_non_shared =  idx->maxh;

  if (hz_from_file_zero == hz_to_file_zero)
  {
    var_grp->f0_start_layout_index = 0;
    var_grp->f0_end_layout_index = 0;
  }

  if (hz_from_shared == hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
  }

  if (hz_from_non_shared == hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
  }

  return PIDX_success;
}

#endif

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->init_start = PIDX_get_time();
  ret = idx_init(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->init_end = PIDX_get_time();

  time->set_reg_box_start = PIDX_get_time();
  if (mode == PIDX_WRITE)
  {
    ret = set_rst_box_size_for_write(file, gi, svi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  else if (mode == PIDX_READ)
  {
    ret = set_rst_box_size_for_read(file, gi, svi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  time->set_reg_box_end = MPI_Wtime();


  time->bit_string_start = PIDX_get_time();
  // calculates maxh and bitstring
  ret = populate_global_bit_string(file, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // selects layout levels based on maxh
  select_io_mode(file, gi);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layout, given this is pure IDX only non-share block layout is populated
  ret = populate_block_layouts(file, gi, svi, 0, 0, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_IDX_IO);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Creates the agg and io ids
  ret = create_agg_io_buffer(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->layout_end = PIDX_get_time();


  if (mode == PIDX_WRITE)
  {
    time->header_io_start = PIDX_get_time();
    // Creates the file heirarchy and writes the header info for all binary files
    ret = write_headers(file, gi, svi, evi, PIDX_IDX_IO);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->header_io_end = PIDX_get_time();
  }

  return PIDX_success;
}


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->group_cleanup_start = PIDX_get_time();

  ret = destroy_agg_io_buffer(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = delete_block_layout(file, gi, 0, 0, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = idx_finalize(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
