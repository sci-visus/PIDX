#include "../../../PIDX_inc.h"

static PIDX_return_code read_block(PIDX_io file, int gi, int vi, int p, int block_number, unsigned long long *patch_offset, unsigned long long *patch_size, unsigned char* patch_buffer);
static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code populate_block_layout_and_buffers(PIDX_io file, int gi, int svi, int evi, int mode);
static PIDX_return_code parse_local_partition_idx_file(PIDX_io file, int partition_index);

PIDX_return_code PIDX_serial_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int bytes_for_datatype;
  int i = 0, j = 0, k = 0;
  int si = 0;
  PIDX_return_code ret;
  char file_name[PATH_MAX];
  MPI_File fp = 0;
  MPI_Status status;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  set_rst_box_size_for_write(file, gi, svi);

  if (populate_block_layout_and_buffers(file, gi, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

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

      for (j = 0; j < file->idx->blocks_per_file; j++)
      {
        if (file->idx_d->block_bitmap[i][j] != 0)
        {
          printf("File number %d Block number %d\n", i, j);
          for (k = 0; k < file->idx_d->samples_per_block; k++)
          {
            hz = (i * file->idx->blocks_per_file * file->idx_d->samples_per_block) + (j * file->idx_d->samples_per_block) + k;
            Hz_to_xyz(file->idx->bitPattern, file->idx_d->maxh, hz, xyz);

            index = (var_grp->variable[si]->sim_patch[0]->size[0] * var_grp->variable[si]->sim_patch[0]->size[1] * xyz[2])
                + (var_grp->variable[si]->sim_patch[0]->size[0] * xyz[1])
                + xyz[0];

            if (xyz[0] >= file->idx->bounds[0] || xyz[1] >= file->idx->bounds[1] || xyz[2] >= file->idx->bounds[2])
              continue;

            memcpy(block_buffer + (k * bytes_for_datatype),
                 var_grp->variable[si]->sim_patch[0]->buffer + (index * bytes_for_datatype),
                bytes_for_datatype);
          }

          if (MPI_File_write_at(fp, file->idx_d->block_offset_bitmap[si][i][j], block_buffer, file->idx_d->samples_per_block * bytes_for_datatype, MPI_BYTE, &status) != MPI_SUCCESS)
          {
            fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
        }
      }
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



PIDX_return_code PIDX_parallel_local_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  int i = 0, j = 0, p = 0;
  int si = 0;
  int ret = 0;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  ret = populate_global_bit_string(file, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  for (si = svi; si < evi; si++)
  {
    PIDX_variable var = var_grp->variable[si];

    for (p = 0; p < var_grp->variable[si]->sim_patch_count; p++)
    {
      int par = 0;
      for (par = 0; par < file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]; par++)
      {
        if (parse_local_partition_idx_file(file, par) == -1) continue;
        char dirname[1024], basename[1024];
        VisusSplitFilename(file->idx->filename_template_partition, dirname, basename);
        sprintf(file->idx->filename_template_partition, "%s/time%09d/%s", dirname, file->idx->current_time_step, basename );

        //printf("Partition %d ----> offset %d %d %d size %d %d %d -> %s\n", par, file->idx_d->partition_offset[0], file->idx_d->partition_offset[1], file->idx_d->partition_offset[2], file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2], file->idx->filename_template_partition);

        int d = 0, check_bit = 0;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          check_bit = check_bit || (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) < file->idx_d->partition_offset[d] || (file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1) < var->sim_patch[p]->offset[d];

        if (!check_bit)
        {
          unsigned long long intersected_box_offset[PIDX_MAX_DIMENSIONS];
          unsigned long long intersected_box_size[PIDX_MAX_DIMENSIONS];

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            if (var->sim_patch[p]->offset[d] <= file->idx_d->partition_offset[d] && (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) <= (file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1))
            {
              intersected_box_offset[d] = file->idx_d->partition_offset[d];
              intersected_box_size[d] = (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) - file->idx_d->partition_offset[d] + 1;
            }
            else if (file->idx_d->partition_offset[d] <= var->sim_patch[p]->offset[d] && (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) >= (file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1))
            {
              intersected_box_offset[d] = var->sim_patch[p]->offset[d];
              intersected_box_size[d] = (file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1) - var->sim_patch[p]->offset[d] + 1;
            }
            else if (( file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1) <= (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) && file->idx_d->partition_offset[d] >= var->sim_patch[p]->offset[d])
            {
              intersected_box_offset[d] = file->idx_d->partition_offset[d];
              intersected_box_size[d] = file->idx_d->partition_size[d];
            }
            else if (( var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) <= (file->idx_d->partition_offset[d] + file->idx_d->partition_size[d] - 1) && var->sim_patch[p]->offset[d] >= file->idx_d->partition_offset[d])
            {
              intersected_box_offset[d] = var->sim_patch[p]->offset[d];
              intersected_box_size[d] = var->sim_patch[p]->size[d];
            }

            intersected_box_offset[d] = intersected_box_offset[d] - file->idx_d->partition_offset[d];
          }

          for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
          {
            bounding_box[0][i] = intersected_box_offset[i];
            bounding_box[1][i] = intersected_box_offset[i] + intersected_box_size[i];
          }
          int bytes_for_datatype = ((var_grp->variable[si]->bpv / 8) * var_grp->variable[si]->vps);
          unsigned char* intersected_box_buffer = malloc(intersected_box_size[0] * intersected_box_size[1] * intersected_box_size[2] * bytes_for_datatype);

          PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
          memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
          ret = PIDX_blocks_initialize_layout(per_patch_local_block_layout, 0, file->idx_d->maxh, file->idx_d->maxh, file->idx->bits_per_block);
          if (ret != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          ret = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
          if (ret != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          //PIDX_blocks_print_layout(per_patch_local_block_layout, file->idx->bits_per_block);
          int ctr = 1;
          int block_number = 0;
          read_block(file, gi, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer);
          for (i = file->idx->bits_per_block + 1 ; i < per_patch_local_block_layout->resolution_to ; i++)
          {
            for (j = 0 ; j < ctr ; j++)
            {
              if(per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
              {
                block_number = per_patch_local_block_layout->hz_block_number_array[i][j];
                read_block(file, gi, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer);
              }
            }
            ctr = ctr * 2;
          }

          PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx_d->maxh, per_patch_local_block_layout);
          free(per_patch_local_block_layout);

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            intersected_box_offset[d] = intersected_box_offset[d] + file->idx_d->partition_offset[d];

          int k1, j1, i1, index = 0, recv_o = 0, send_o = 0, send_c = 0;
          for (k1 = intersected_box_offset[2]; k1 < intersected_box_offset[2] + intersected_box_size[2]; k1++)
          {
            for (j1 = intersected_box_offset[1]; j1 < intersected_box_offset[1] + intersected_box_size[1]; j1++)
            {
              for (i1 = intersected_box_offset[0]; i1 < intersected_box_offset[0] + intersected_box_size[0]; i1 = i1 + intersected_box_size[0])
              {
                index = ((intersected_box_size[0])* (intersected_box_size[1]) * (k1 - intersected_box_offset[2])) + ((intersected_box_size[0]) * (j1 - intersected_box_offset[1])) + (i1 - intersected_box_offset[0]);
                send_o = index * var->vps * (var->bpv/8);
                send_c = (intersected_box_size[0]);
                recv_o = (var_grp->variable[si]->sim_patch[p]->size[0] * var_grp->variable[si]->sim_patch[p]->size[1] * (k1 - var_grp->variable[si]->sim_patch[p]->offset[2])) + (var_grp->variable[si]->sim_patch[p]->size[0] * (j1 - var_grp->variable[si]->sim_patch[p]->offset[1])) + (i1 - var_grp->variable[si]->sim_patch[p]->offset[0]);

                memcpy(var_grp->variable[si]->sim_patch[p]->buffer + (recv_o * var->vps * (var->bpv/8)), intersected_box_buffer + send_o, send_c * var->vps * (var->bpv/8));
              }
            }
          }
          free(intersected_box_buffer);
        }
      }
    }
  }

  return PIDX_success;
}



static PIDX_return_code read_block(PIDX_io file, int gi, int vi, int p, int block_number, unsigned long long* patch_offset, unsigned long long* patch_size, unsigned char* patch_buffer)
{
  int k = 0, ret = 0;
  MPI_File fp = 0;
  MPI_Status status;
  unsigned long long index = 0;
  unsigned long long hz;
  unsigned long long xyz[PIDX_MAX_DIMENSIONS];
  char file_name[PATH_MAX];

  int file_number = block_number / file->idx->blocks_per_file;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  int bytes_for_datatype = ((var_grp->variable[vi]->bpv / 8) * var_grp->variable[vi]->vps);

  unsigned char* block_buffer = malloc(file->idx_d->samples_per_block * bytes_for_datatype);

  //generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename_partition, file->idx->current_time_step, file->idx->filename_template_partition);

  if (generate_file_name(file->idx->blocks_per_file, file->idx->filename_template_partition, file_number, file_name, PATH_MAX) == 1)
  {
    fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  if (MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_open() block number %d file number %d filename %s failed.\n", __FILE__, __LINE__, block_number, file_number, file_name);
    return PIDX_err_io;
  }

  uint32_t *headers;
  int total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  headers = malloc(total_header_size);
  memset(headers, 0, total_header_size);

  ret = MPI_File_read_at(fp, 0, headers, total_header_size , MPI_BYTE, &status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Data offset = [%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }
  int read_count = 0;
  MPI_Get_count(&status, MPI_BYTE, &read_count);
  if (read_count != total_header_size)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %dd\n", __FILE__, __LINE__, read_count, total_header_size);
    return PIDX_err_io;
  }

  off_t data_offset = htonl(headers[12 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);
  size_t data_size = htonl(headers[14 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);

  //printf("File number %d Block number %d\n", file_number, block_number);
  assert (data_size != 0);

  ret = MPI_File_read_at(fp, data_offset, block_buffer, data_size, MPI_BYTE, &status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }

  for (k = 0; k < file->idx_d->samples_per_block; k++)
  {
    hz = (block_number * file->idx_d->samples_per_block) + k;
    Hz_to_xyz(file->idx->bitPattern, file->idx_d->maxh, hz, xyz);

    if ( ((xyz[0] < patch_offset[0] || xyz[0] >= patch_offset[0] + patch_size[0]) ||
        (xyz[1] < patch_offset[1] || xyz[1] >= patch_offset[1] + patch_size[1]) ||
        (xyz[2] < patch_offset[2] || xyz[2] >= patch_offset[2] + patch_size[2]) ) )
      continue;

    index = (patch_size[0] * patch_size[1] * (xyz[2] - patch_offset[2]))
        + (patch_size[0] * (xyz[1] - patch_offset[1]))
        + (xyz[0] - patch_offset[0]);

    memcpy(patch_buffer + (index * bytes_for_datatype),
         block_buffer + (k * bytes_for_datatype),
         bytes_for_datatype);

  }

  MPI_File_close(&fp);

  free(headers);
  free(block_buffer);

  return PIDX_success;
}


static PIDX_return_code parse_local_partition_idx_file(PIDX_io file, int partition_index)
{

    char file_name_skeleton[1024];
    strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
    file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';
    sprintf(file->idx->filename_partition, "%s_%d.idx", file_name_skeleton, partition_index);

    char *pch;
    int count = 0;
    char line [ 512 ];

    FILE *fp = fopen(file->idx->filename_partition, "r");
    if (fp == NULL)  return -1;


    while (fgets(line, sizeof (line), fp) != NULL)
    {
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "(partition index)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        file->idx_d->color = atoi(line);
      }

      if (strcmp(line, "(partition size)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          file->idx_d->partition_size[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
      }

      if (strcmp(line, "(bits)") == 0)
      {
        int i = 0;
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        strcpy(file->idx->bitSequence, line);
        file->idx_d->maxh = strlen(file->idx->bitSequence);
        for (i = 0; i <= file->idx_d->maxh; i++)
          file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);
        //printf("BS %s MH %d\n", file->idx->bitSequence, file->idx_d->maxh);
      }

      if (strcmp(line, "(filename_template)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        memset(file->idx->filename_template_partition, 0, 1024 * sizeof(char));
        strcpy(file->idx->filename_template_partition, line);
      }

      if (strcmp(line, "(partition offset)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          file->idx_d->partition_offset[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
      }
    }
    fclose(fp);

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
