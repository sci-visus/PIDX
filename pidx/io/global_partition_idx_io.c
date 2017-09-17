#include "../PIDX_inc.h"

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);


/// Global Partitioned IDX Write Steps
/*********************************************************
*  Step 0:  group and IDX related meta data              *
*                                                        *
*  Step 1:  Restrucure setup                             *
*  Step 2:  Restrucure                                   *
*  Step 3:  Partition                                    *
*                                                        *
*  Step 4:  Post partition group meta data               *
*                                                        *
*  Step 5:  Setup HZ encoding Phase                      *
*  Step 6:  Perform HZ encoding                          *
*  Step 7:  Setup aggregation buffers                    *
*  Step 8:  Perform data aggregation                     *
*  Step 9:  Perform actual file IO                       *
*  Step 10: cleanup for Steps 6                          *
*                                                        *
*  Step 11: Cleanup the group and IDX related meta-data  *
*                                                        *
*  Step 12: Partition cleanup                            *
*  Step 13: Restructuring cleanup                        *
**********************************************************/

PIDX_return_code PIDX_global_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 0:  group and IDX related meta data
  if (group_meta_data_init(file, gi, svi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 1:  Restrucure setup
  if (restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2:  Restrucure
  if (restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 3:  Partition
  if (partition(file, gi, svi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->idx->variable_grp[gi]->variable[svi]->restructured_super_patch_count == 0)
    goto cleanup;

  // Step 4:  Post partition group meta data
  ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 5:  Setup HZ encoding Phase
    if (hz_encode_setup(file, gi, si, ei) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 6: Perform HZ encoding
    if (hz_encode(file, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (hz_io(file, gi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Setup 7: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, si, li, li, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Setup 8: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, si, li, li, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
#if 1
    // Setup 9: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[gi][li] = PIDX_get_time();
      create_async_buffers(file, gi);

      ret = data_io(file, gi, si, li, li, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_async_buffers(file, gi);
      finalize_aggregation(file, gi, li, si);
      time->io_end[gi][li] = PIDX_get_time();
    }

    // Step 10: Cleanup for step 6
    if (hz_encode_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#endif
  }
#if 1
  // Step 11: Cleanup the group and IDX related meta-data
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif
  // Step 12: Partition cleanup
  cleanup:
  time->partition_cleanup_start = MPI_Wtime();
  if (destroy_local_comm(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_cleanup_end = MPI_Wtime();

  // Step 13: Restructuring cleanup
  if (restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


/// Global Partitioned IDX Read Steps
/*********************************************************
*  Step 0:  group and IDX related meta data              *
*                                                        *
*  Step 1:  Restrucure setup                             *
*  Step 2:  Partition                                    *
*                                                        *
*  Step 3:  Post partition group meta data               *
*                                                        *
*  Step 4:  Setup HZ encoding Phase                      *
*  Step 5:  Setup aggregation buffers                    *
*  Step 6:  Perform actual file IO                       *
*  Step 7:  Perform data aggregation                     *
*  Step 8:  Perform HZ encoding                          *
*  Step 9:  cleanup for Steps 6                          *
*                                                        *
*  Step 10: Cleanup the group and IDX related meta-data  *
*                                                        *
*  Step 11: Partition cleanup                            *
*  Step 12: Restrucure                                   *
*  Step 13: Restructuring cleanup                        *
**********************************************************/

PIDX_return_code PIDX_global_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 0:  group and IDX related meta data
  ret = group_meta_data_init(file, gi, svi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 1:  Restrucure setup
  if (restructure_setup(file, gi, svi, evi - 1, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2:  Partition
  ret = partition(file, gi, svi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->idx->variable_grp[gi]->variable[svi]->restructured_super_patch_count == 0)
    goto cleanup;

  // Step 3:  Post partition group meta data
  ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 4:  Setup HZ encoding Phase
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

    // Setup 5: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, si, li, li, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Setup 6: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, si, li, li, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[gi][li] = PIDX_get_time();
    }

    //
    // Setup 7: Performs data aggregation
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
    //

    // Step 8: Perform HZ encoding
    ret = hz_encode(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 9: Cleanup for step 6
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 10: Cleanup the group and IDX related meta-data
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 11: Partition cleanup
  cleanup:
  time->partition_cleanup_start = MPI_Wtime();
  if (destroy_local_comm(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_cleanup_end = MPI_Wtime();

  // Step 12:  Restrucure
  if (restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 13: Restructuring cleanup
  if (restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}





static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->partition_start = MPI_Wtime();
  // Calculates the number of partititons
  if (mode == PIDX_WRITE)
  {
    ret = find_partition_count(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // assign same color to processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Splits the global communicator into local communicators
  ret = create_local_comm(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_end = MPI_Wtime();

  return PIDX_success;
}



static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->init_start = MPI_Wtime();
  time->init_end = MPI_Wtime();

  time->set_reg_box_start = MPI_Wtime();
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

  return PIDX_success;
}


static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
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

  // selects levels based on mode and maxh
  select_io_mode(file, gi);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layout, given this is pure IDX only non-share block layout is populated
  ret = populate_block_layouts(file, gi, svi, file->hz_from_shared, file->hz_to_shared, file->hz_from_non_shared, file->hz_to_non_shared);
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
    ret = write_headers(file, gi, svi, evi, PIDX_GLOBAL_PARTITION_IDX_IO);
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

  ret = delete_block_layout(file, gi, file->hz_from_shared, file->hz_to_shared, file->hz_from_non_shared, file->hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
