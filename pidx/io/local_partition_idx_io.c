#include "../PIDX_inc.h"

static int hz_from_file_zero = 0, hz_from_shared = 0, hz_from_non_shared = 0;
static int hz_to_file_zero = 0, hz_to_shared = 0, hz_to_non_shared = 0;
//static int agg_l_nshared = 0, agg_l_shared = 0, agg_l_f0 = 0;
static PIDX_return_code find_agg_level(PIDX_io file, int gi);
static PIDX_return_code select_io_mode(PIDX_io file, int gi);
static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi);
static PIDX_return_code re_adjust_offsets(PIDX_io file, int gi, int svi);
static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);


/// local Partitioned IDX Write Steps
/*********************************************************
*  Step 0:  group and IDX related meta data              *
*                                                        *
*  Step 1:  Restrucure setup                             *
*  Step 2:  Restrucure                                   *
*  Step 3:  Partition                                    *
*                                                        *
*  Step 4:  Adjust offset                                *
*                                                        *
*  Step 5:  Post partition group meta data               *
*                                                        *
*  Step 6:  Setup HZ encoding Phase                      *
*  Step 7:  Perform HZ encoding                          *
*  Step 8:  Setup aggregation buffers                    *
*  Step 9:  Perform data aggregation                     *
*  Step 10: Perform actual file IO                       *
*  Step 11: cleanup for Steps 6                          *
*                                                        *
*  Step 12: Cleanup the group and IDX related meta-data  *
*                                                        *
*  Step 13: Partition cleanup                            *
*  Step 14: Restructuring cleanup                        *
**********************************************************/

PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
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

  // Step 4: Adjust per process offsets and global bounds as per the partition
  if (adjust_offsets(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->idx->variable_grp[gi]->variable[svi]->patch_group_count == 0)
    goto cleanup;

  // Step 5:  Post partition group meta data
  ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#if 1
  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 6:  Setup HZ encoding Phase
    if (hz_encode_setup(file, gi, si, ei) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 7: Perform HZ encoding
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

    // Setup 8: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, si, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Setup 9: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, si, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Setup 10: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[gi][li] = PIDX_get_time();
      create_async_buffers(file, gi);

      ret = data_io(file, gi, li, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      wait_and_destroy_async_buffers(file, gi);
      finalize_aggregation(file, gi, li);
      time->io_end[gi][li] = PIDX_get_time();
    }

    // Step 11: Cleanup for step 6
    if (hz_encode_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 12: Cleanup the group and IDX related meta-data
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif
  // Step 13: Partition cleanup
  cleanup:
  time->partition_cleanup_start = MPI_Wtime();
  if (destroy_local_comm(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_cleanup_end = MPI_Wtime();

  // Step 14: Restructuring cleanup
  if (restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


/// local Partitioned IDX Read Steps
/*********************************************************
*  Step 0:  group and IDX related meta data              *
*                                                        *
*  Step 1:  Restrucure setup                             *
*  Step 2:  Partition                                    *
*                                                        *
*  Step 3:  Adjust offsets                               *
*                                                        *
*  Step 4:  Post partition group meta data               *
*                                                        *
*  Step 5:  Setup HZ encoding Phase                      *
*  Step 6:  Setup aggregation buffers                    *
*  Step 7:  Perform actual file IO                       *
*  Step 8:  Perform data aggregation                     *
*  Step 9:  Perform HZ encoding                          *
*  Step 10:  cleanup for Steps 6                         *
*                                                        *
*  Step 11: Cleanup the group and IDX related meta-data  *
*                                                        *
*  Step 12: Partition cleanup                            *
*  Step 13: Restrucure                                   *
*  Step 14: Restructuring cleanup                        *
**********************************************************/

PIDX_return_code PIDX_local_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
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
  ret = partition(file, gi, svi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 3: Adjust per process offsets and global bounds as per the partition
  if (adjust_offsets(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (file->idx->variable_grp[gi]->variable[svi]->patch_group_count == 0)
    goto cleanup;

  // Step 4:  Post partition group meta data
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

    // Step 5:  Setup HZ encoding Phase
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

    // Setup 6: Setup aggregation buffers
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, si, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Setup 7: Performs actual file io
    for (li = si; li <= ei; li = li + 1)
    {
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, li, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[gi][li] = PIDX_get_time();
    }

    //
    // Setup 8: Performs data aggregation
    for (li = si; li <= ei; li = li + 1)
    {
      ret = data_aggregate(file, gi, li, si, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      finalize_aggregation(file, gi, li);
    }
    //

    // Step 9: Perform HZ encoding
    ret = hz_encode(file, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 10: Cleanup for step 6
    ret = hz_encode_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Step 11: Cleanup the group and IDX related meta-data
  ret = group_meta_data_finalize(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 12: Partition cleanup
  cleanup:
  time->partition_cleanup_start = MPI_Wtime();
  if (destroy_local_comm(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_cleanup_end = MPI_Wtime();

  // Step 3: Adjust per process offsets and global bounds as per the partition
  if (re_adjust_offsets(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 13:  Restrucure
  if (restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 14: Restructuring cleanup
  if (restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
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
  hz_to_shared =  idx->total_partiton_level;

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

  // Splits the local communicator into local communicators
  ret = create_local_comm(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_end = MPI_Wtime();

  return PIDX_success;
}



static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi)
{
  int i = 0, p = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[svi];

  if (var->patch_group_count > 1)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (p = 0 ; p < var->patch_group_count; p++)
  {
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      var->rst_patch_group[0]->reg_patch->offset[i] = var->rst_patch_group[0]->reg_patch->offset[i] - file->idx_d->partition_offset[i];
    }
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx_d->partition_offset[i] + file->idx_d->partition_size[i] <= file->idx->box_bounds[i])
      file->idx->box_bounds[i] = file->idx_d->partition_size[i];
    else
      file->idx->box_bounds[i] = file->idx->box_bounds[i] - file->idx_d->partition_offset[i];

    if (getPowerOf2(file->idx->box_bounds[i]) < file->idx->reg_patch_size[i])
      file->idx->box_bounds[i] = file->idx->reg_patch_size[i];//(file->idx->box_bounds[i] / file->idx->reg_patch_size[i] + 1) * file->idx->reg_patch_size[i];
  }

  return PIDX_success;
}


static PIDX_return_code re_adjust_offsets(PIDX_io file, int gi, int svi)
{
  int i = 0, p = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[svi];

  if (var->patch_group_count > 1)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (p = 0 ; p < var->patch_group_count; p++)
  {
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      var->rst_patch_group[0]->reg_patch->offset[i] = var->rst_patch_group[0]->reg_patch->offset[i] + file->idx_d->partition_offset[i];
    }
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    //if (file->idx_c->grank == 0)
    //  fprintf(stderr, "%d: %d\n", i, file->idx->bounds[i]);

    file->idx->box_bounds[i] = file->idx->bounds[i];
  }

  // TODO
  /*
  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx_d->partition_offset[i] + file->idx_d->partition_size[i] <= file->idx->bounds[i])
      file->idx->bounds[i] = file->idx_d->partition_size[i];
    else
      file->idx->bounds[i] = file->idx->bounds[i] - file->idx_d->partition_offset[i];

    //file->idx->bounds[i] = 256;
  }
  */

  return PIDX_success;
}



static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->init_start = MPI_Wtime();
  ret = idx_init(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->init_end = MPI_Wtime();

  time->set_reg_box_start = MPI_Wtime();
  ret = set_rst_box_size(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
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
  ret = populate_local_bit_string(file, mode);
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
  ret = populate_block_layouts(file, gi, svi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared, PIDX_LOCAL_PARTITION_IDX_IO);
  if (ret != PIDX_success)
  {
     fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
     return PIDX_err_file;
  }
#if 1
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
    ret = write_headers(file, gi, svi, evi, PIDX_LOCAL_PARTITION_IDX_IO);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->header_io_end = PIDX_get_time();
  }
#endif
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

  ret = delete_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero, hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
