#include "../../../PIDX_inc.h"

static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code populate_block_layout_and_buffers(PIDX_io file, int gi, int svi, int evi, int mode);

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

  if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 1: Setup restructuring buffers
  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2: Perform data restructuring
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if 1
  if (idx_restructure_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    if (populate_block_layout_and_buffers(file, gi, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    file->idx->variable_pipe_length = file->idx->variable_count;
    for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      if ((si + file->idx->variable_pipe_length) >= evi)
        file->idx->variable_pipe_length = evi - (si + 1);

      ei = si + file->idx->variable_pipe_length;

      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

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


      // Step 5: Setup aggregation buffers
      //for (li = si; li <= ei; li = li + 1)
      //{
      //printf("[PL %d] [%d %d %d]\n", file->idx->variable_pipe_length, si, li, ei);
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

        ret = data_io(file, gi, si, li, ei, PIDX_WRITE);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        finalize_aggregation(file, gi, li, si);
        time->io_end[gi][li] = PIDX_get_time();
      //}
#if 0
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

    // Step 9
    ret = group_meta_data_finalize(file, gi, svi, evi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  free_restructured_communicators(file, gi);

  ret = idx_restructure_cleanup(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif
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

  // Step 1: Setup restructuring buffers
  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    if (populate_block_layout_and_buffers(file, gi, svi, evi, PIDX_READ) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    file->idx->variable_pipe_length = file->idx->variable_count;
    for (si = svi; si < evi; si = si + (file->idx->variable_count + 1))
    {
      ei = ((si + file->idx->variable_count) >= (evi)) ? (evi - 1) : (si + file->idx->variable_count);
      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      // Step 2: Setup HZ buffers
      ret = hz_encode_setup(file, gi, si, ei);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }


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

      // Step 8: Cleanup all buffers and ids
      ret = hz_encode_cleanup(file);
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
  }

  free_restructured_communicators(file, gi);

  // Step 7: Perform data restructuring
  ret = idx_restructure(file, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = idx_restructure_cleanup(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


static PIDX_return_code populate_block_layout_and_buffers(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

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
  ret = populate_rst_block_layouts(file, gi, svi, file->hz_from_shared, file->hz_to_non_shared);
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



  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  ret = write_headers(file, gi, svi, evi, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();


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

  ret = delete_rst_block_layout(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}