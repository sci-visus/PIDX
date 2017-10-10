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

PIDX_return_code PIDX_serial_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int bytes_for_datatype;
  int i = 0, j = 0, k = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  char file_name[PATH_MAX];
  PIDX_time time = file->idx_d->time;
  MPI_File fp = 0;
  MPI_Status status;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

    if (populate_block_layout_and_buffers(file, gi, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    file->idx->variable_pipe_length = file->idx->variable_count;
    for (si = svi; si < evi; si++)
    {
      bytes_for_datatype = ((var_grp->variable[si]->bpv / 8) * var_grp->variable[si]->vps);

      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      unsigned long long index = 0;
      unsigned long long hz;
      unsigned long long xyz[PIDX_MAX_DIMENSIONS];
      unsigned char* block_buffer = malloc(file->idx_d->samples_per_block * bytes_for_datatype);

      for (i = 0; i < file->idx_d->max_file_count; i++)
      {

        if (generate_file_name(file->idx->blocks_per_file, file->idx->filename_template_partition, i, file_name, PATH_MAX) == 1)
        {
          fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }

        if (MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
          return PIDX_err_io;
        }

        //
        for (j = 0; j < file->idx->blocks_per_file; j++)
        {
          if (file->idx_d->block_bitmap[i][j] != 0)
          {
            for (k = 0; k < file->idx_d->samples_per_block; k++)
            {
              hz = (i * file->idx->blocks_per_file * file->idx_d->samples_per_block) + (j * file->idx_d->samples_per_block) + k;
              Hz_to_xyz(file->idx->bitPattern, file->idx_d->maxh, hz, xyz);

              index = (var_grp->variable[si]->sim_patch[0]->size[0] * var_grp->variable[si]->sim_patch[0]->size[1] * xyz[2])
                      + (var_grp->variable[si]->sim_patch[0]->size[0] * xyz[1])
                      + xyz[0];

              if (xyz[0] >= file->idx->bounds[0] || xyz[1] >= file->idx->bounds[1] || xyz[2] >= file->idx->bounds[2])
                  continue;

              //printf("[%d %d %d] HZ %lld -> %lld %lld %lld [%d] bdt %d\n", i, j, k, hz, xyz[0], xyz[1], xyz[2], index, bytes_for_datatype);
              memcpy(block_buffer + (k * bytes_for_datatype),
                   var_grp->variable[si]->sim_patch[0]->buffer + (index * bytes_for_datatype),
                   bytes_for_datatype);

            }

            if (MPI_File_write_at(fp, file->idx_d->block_offset_bitmap[i][j], block_buffer, file->idx_d->samples_per_block * bytes_for_datatype, MPI_BYTE, &status) != MPI_SUCCESS)
            {
              fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
              return PIDX_err_io;
            }
          }
        }
        //
        MPI_File_close(&fp);

      }
      free(block_buffer);
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

PIDX_return_code PIDX_serial_idx_read(PIDX_io file, int gi, int svi, int evi)
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
  ret = populate_sim_block_layouts(file, gi, svi, file->hz_from_shared, file->hz_to_non_shared);
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
  ret = delete_sim_block_layout(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
