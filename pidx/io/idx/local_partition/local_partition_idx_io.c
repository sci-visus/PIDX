#include "../../../PIDX_inc.h"

static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);


/// local Partitioned IDX Write Steps
/*********************************************************
*  Step 0:  group and IDX related meta data        *
*                            *
*  Step 1:  Restrucure setup               *
*  Step 2:  Restrucure                   *
*  Step 3:  Partition                  *
*                            *
*  Step 4:  Adjust offset                *
*                            *
*  Step 5:  Post partition group meta data         *
*                            *
*  Step 6:  Setup HZ encoding Phase            *
*  Step 7:  Perform HZ encoding              *
*  Step 8:  Setup aggregation buffers          *
*  Step 9:  Perform data aggregation           *
*  Step 10: Perform actual file IO             *
*  Step 11: cleanup for Steps 6              *
*                            *
*  Step 12: Cleanup the group and IDX related meta-data  *
*                            *
*  Step 13: Partition cleanup              *
*  Step 14: Restructuring cleanup            *
**********************************************************/

PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 1:  Restrucure setup
  set_rst_box_size_for_write(file, gi, svi);

  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2:  Restrucure
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_rst_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    // Step 3:  Partition
    if (partition(file, gi, svi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = populate_bit_string(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, gi, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 5:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if 1
    for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
    {
      ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
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
      ret = data_aggregate(file, gi, si, ei, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 9: Performs data aggregation
      ret = data_aggregate(file, gi, si, ei, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 10: Performs actual file io
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, si, ei, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      finalize_aggregation(file, gi, si);
      time->io_end[gi][li] = PIDX_get_time();


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

    // Step 13: Partition cleanup
    time->partition_cleanup_start = MPI_Wtime();
    if (destroy_local_comm(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->partition_cleanup_end = MPI_Wtime();
#endif
  }

  //free_restructured_communicators(file, gi);

  // Step 14: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
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

  // Step 1:  Restrucure setup

  set_rst_box_size_for_write(file, gi, svi);

  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_rst_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    // Step 2:  Partition
    ret = partition(file, gi, svi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, gi, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
    {
      ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
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
      ret = data_aggregate(file, gi, si, ei, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }


      // Setup 7: Performs actual file io
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, si, ei, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[gi][li] = PIDX_get_time();

#if 1
      //
      // Setup 8: Performs data aggregation
      ret = data_aggregate(file, gi, si, ei, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      finalize_aggregation(file, gi, si);

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
      #endif
    }

    // Step 11: Cleanup the group and IDX related meta-data
    ret = group_meta_data_finalize(file, gi, svi, evi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 12: Partition cleanup
    time->partition_cleanup_start = MPI_Wtime();
    if (destroy_local_comm(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->partition_cleanup_end = MPI_Wtime();

    int i = 0;
    PIDX_variable var = var_grp->variable[svi];
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] + file->idx_d->partition_offset[i];

  }

#if 1
  // Step 13:  Restrucure
  if (idx_restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 14: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_local_partition_mapped_idx_read(PIDX_io file, int gi, int svi, int evi)
{

}




static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->partition_start = MPI_Wtime();
  // Calculates the number of partititons
  ret = find_partition_count(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
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



static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi, int evi)
{
  int i = 0, v = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  fprintf(stderr, "SE %d %d\n", svi, evi);
  for (v = svi; v < evi; v++)
  {
  PIDX_variable var = var_grp->variable[v];

  if (var->restructured_super_patch_count != 1)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] - file->idx_d->partition_offset[i];
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx_d->partition_offset[i] + file->idx_d->partition_size[i] <= file->idx->box_bounds[i])
      file->idx->box_bounds[i] = file->idx_d->partition_size[i];
    else
      file->idx->box_bounds[i] = file->idx->box_bounds[i] - file->idx_d->partition_offset[i];
  }

  memcpy(file->idx->bounds, file->idx->box_bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  //fprintf(stderr, "%d - %d %d %d -- PO %d %d %d\n", file->idx_c->grank, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], file->idx_d->partition_offset[0], file->idx_d->partition_offset[1], file->idx_d->partition_offset[2]);

  return PIDX_success;
}



static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  //fprintf(stderr, "%d %d %d\n", file->idx_d->restructured_grid->patch_size[0], file->idx_d->restructured_grid->patch_size[1], file->idx_d->restructured_grid->patch_size[2]);

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
  ret = populate_rst_block_layouts(file, gi, svi, file->hz_from_shared, file->hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi, svi, evi);
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
