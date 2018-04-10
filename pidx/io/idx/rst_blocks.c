#include "../../PIDX_inc.h"

static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int group_index, int start_index, int hz_level_from, int hz_level_to);

static PIDX_return_code populate_idx_layout(PIDX_io file, int gi, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level);

static PIDX_return_code destroy_block_layout(PIDX_io file, int gi);


PIDX_return_code populate_rst_block_layouts(PIDX_io file, int gi, int svi, int hz_from_shared, int hz_to_non_shared)
{
  int i = 0, v = 0, ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  file->idx_d->block_bitmap = malloc(file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  memset(file->idx_d->block_bitmap, 0, file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    file->idx_d->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
    memset(file->idx_d->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
  }

  file->idx_d->block_offset_bitmap = malloc(file->idx->variable_count * sizeof (*file->idx_d->block_offset_bitmap));
  memset(file->idx_d->block_offset_bitmap, 0, file->idx->variable_count * sizeof (*file->idx_d->block_offset_bitmap));
  for (v = 0; v < file->idx->variable_count; v++)
  {
    file->idx_d->block_offset_bitmap[v] = malloc(file->idx_d->max_file_count * sizeof (*(file->idx_d->block_offset_bitmap[v])));
    memset(file->idx_d->block_offset_bitmap[v], 0, file->idx_d->max_file_count * sizeof (*(file->idx_d->block_offset_bitmap[v])));

    for (i = 0; i < file->idx_d->max_file_count; i++)
    {
      file->idx_d->block_offset_bitmap[v][i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_offset_bitmap[v][i]));
      memset(file->idx_d->block_offset_bitmap[v][i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_offset_bitmap[v][i]));
    }
  }

  int total_layout_count = var_grp->shared_layout_count + var_grp->nshared_layout_count;
  var_grp->block_layout = malloc(sizeof (*var_grp->block_layout));
  memset(var_grp->block_layout, 0, sizeof (*var_grp->block_layout));

  var_grp->block_layout_by_level = malloc(sizeof (*var_grp->block_layout_by_level) * total_layout_count);
  memset(var_grp->block_layout_by_level, 0, sizeof (*var_grp->block_layout_by_level) * total_layout_count);

  for (i = 0; i < total_layout_count ; i++)
  {
    var_grp->block_layout_by_level[i] = malloc(sizeof(*(var_grp->block_layout_by_level[i])));
    memset(var_grp->block_layout_by_level[i], 0, sizeof(*(var_grp->block_layout_by_level[i])));
  }

  ret = populate_idx_block_layout(file,
                  var_grp->block_layout, var_grp->block_layout_by_level,
                  var_grp->shared_start_layout_index, var_grp->nshared_end_layout_index,
                  var_grp->shared_layout_count + var_grp->nshared_layout_count,
                  gi,
                  svi,
                  hz_from_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


static PIDX_return_code populate_idx_layout(PIDX_io file, int gi, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level)
{
  int i, j, ctr = 1;
  PIDX_return_code ret_code;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  int lvi = start_var_index;


  PIDX_block_layout all_patch_local_block_layout = malloc(sizeof (*all_patch_local_block_layout));
  memset(all_patch_local_block_layout, 0, sizeof (*all_patch_local_block_layout));
  ret_code = PIDX_blocks_initialize_layout(all_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[lvi];

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    bounding_box[0][i] = var->restructured_super_patch->restructured_patch->offset[i];
    bounding_box[1][i] = var->restructured_super_patch->restructured_patch->size[i] + var->restructured_super_patch->restructured_patch->offset[i];

    bounding_box[0][i] = (bounding_box[0][i] / file->idx->chunk_size[i]);

    if (bounding_box[1][i] % file->idx->chunk_size[i] == 0)
      bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]);
    else
      bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]) + 1;
  }

  PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
  memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
  ret_code = PIDX_blocks_initialize_layout(per_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  //fprintf(stderr, "MH %d BS %s F %d T %d [%d %d %d]\n", file->idx_d->maxh, file->idx->bitSequence, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2]);

  ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (all_patch_local_block_layout->resolution_from <= file->idx->bits_per_block)
  {
    for (i = all_patch_local_block_layout->resolution_from ; i <=   file->idx->bits_per_block ; i++)
    {
      if (per_patch_local_block_layout->hz_block_number_array[i][0] == 0)
      {
        all_patch_local_block_layout->hz_block_number_array[i][0] = per_patch_local_block_layout->hz_block_number_array[i][0];
        break;
      }
    }

    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
    {
      for (j = 0 ; j < ctr ; j++)
      {
        if(per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
          all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
    {
      if (i >= all_patch_local_block_layout->resolution_from)
      {
        for (j = 0 ; j < ctr ; j++)
        {
          if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
            all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
        }
      }
      ctr = ctr * 2;
    }
  }

  PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx_d->maxh, per_patch_local_block_layout);
  free(per_patch_local_block_layout);
  per_patch_local_block_layout = 0;

  //fprintf(stderr, "RF %d BPB %d\n", block_layout->resolution_from, file->idx->bits_per_block);
  if (block_layout->resolution_from <= file->idx->bits_per_block)
  {
    int level_count = 1;
    for (i = block_layout->resolution_from; i <= file->idx->bits_per_block; i++)
    {
      //if (i >= file->idx_d->maxh)
      //  continue;

      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
    }

    for (i = file->idx->bits_per_block + 1; i < (block_layout->resolution_to); i++)
    {
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
      level_count = level_count * 2;
    }
  }
  else
  {
    int level_count = 1;
    for (i =   file->idx->bits_per_block + 1; i < (block_layout->resolution_to); i++)
    {
      if (i >= block_layout->resolution_from)
      {
        MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
      }
      level_count = level_count * 2;
    }
  }

  //if (rank == 4)
  //  PIDX_blocks_print_layout(all_patch_local_block_layout);

  PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx_d->maxh, all_patch_local_block_layout);
  free(all_patch_local_block_layout);
  all_patch_local_block_layout = 0;

  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->bcpf = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->bcpf, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->lbi = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->lbi, 0, sizeof(int) * (file->idx_d->max_file_count));

  int file_number = 0;
  if (block_layout->resolution_from <=   file->idx->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->bcpf[file_number]++;
        block_layout->lbi[file_number] = 0;
        break;
      }
    }

    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->bcpf[file_number]++;
          block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
        }
      }
      ctr = ctr * 2;
    }

    //if (file->idx_c->lnprocsrank == 4)
    //  PIDX_blocks_print_layout(block_layout);
  }
  else
  {
    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->bcpf[file_number]++;
            block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
          }
        }
      }
      ctr = ctr * 2;
    }
  }
  block_layout->efc = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->efc++;

  block_layout->existing_file_index = (int*) malloc(block_layout->efc * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->efc * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;

      count++;
    }
  }

  return PIDX_success;
}



static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout block_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int gi, int si, int hz_level_from, int hz_level_to)
{
  if (hz_level_from == 0 && hz_level_to == 0)
    return PIDX_success;

  PIDX_return_code ret_code;

  int i = 0, j = 0, ctr, file_number = 0;

  int lower_hz_level = 0, higher_hz_level = 0;
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  lower_hz_level = hz_level_from;
  higher_hz_level = hz_level_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (start_layout_index == 0)
  {
    lower_level_low_layout = 0;
    higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

    if (higher_level_low_layout >= higher_hz_level)
      higher_level_low_layout = higher_hz_level;

    ret_code = PIDX_blocks_initialize_layout(layout_by_level[0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    //fprintf(stderr, "LL HL %d %d\n", lower_level_low_layout, higher_level_low_layout);

    ret_code = populate_idx_layout(file, gi, si, layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
    {
      if (j >= file->idx_d->maxh)
        continue;

      memcpy(block_layout->hz_block_number_array[j], layout_by_level[0]->hz_block_number_array[j], sizeof(int));
    }

    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
    {
      memcpy(block_layout->hz_block_number_array[j], layout_by_level[0]->hz_block_number_array[j], sizeof(int) * ctr);
      ctr = ctr * 2;
    }

    for (i = 1; i < layout_count; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = populate_idx_layout(file, gi, si, layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], layout_by_level[i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
      ctr = ctr * 2;

    ctr = (int)pow(2, start_layout_index - 1) * file->idx->blocks_per_file;
    for (i = start_layout_index; i < end_layout_index; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(layout_by_level[i - start_layout_index], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = populate_idx_layout(file, gi, si, layout_by_level[i - start_layout_index], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], layout_by_level[i - start_layout_index]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }

  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->bcpf = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->bcpf, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->lbi = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->lbi, 0, sizeof(int) * (file->idx_d->max_file_count));

  if (block_layout->resolution_from <=   file->idx->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        file->idx_d->block_bitmap[file_number][block_layout->hz_block_number_array[i][0] % file->idx->blocks_per_file] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->bcpf[file_number]++;
        break;
      }
    }

    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          file->idx_d->block_bitmap[file_number][block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->bcpf[file_number]++;
          block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
        }
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            file->idx_d->block_bitmap[file_number][block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->bcpf[file_number]++;
            block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
          }
        }
      }
      ctr = ctr * 2;
    }
  }

  block_layout->efc = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->efc++;

  block_layout->existing_file_index = (int*) malloc(block_layout->efc * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->efc * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  //if (file->idx_c->grank == 0)
  //  PIDX_blocks_print_layout(block_layout, file->idx->bits_per_block);

  return PIDX_success;
}


PIDX_return_code delete_rst_block_layout(PIDX_io file, int gi)
{
  int i;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  PIDX_free_layout(var_grp->block_layout);
  PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx_d->maxh, var_grp->block_layout);

  for (i = 0; i < var_grp->shared_layout_count + var_grp->nshared_layout_count ; i++)
  {
    PIDX_free_layout(var_grp->block_layout_by_level[i]);
    PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx_d->maxh, var_grp->block_layout_by_level[i]);
  }
  destroy_block_layout(file, gi);

  int v = 0;
  for (v = 0; v < file->idx->variable_count; v++)
  {
    for (i = 0; i < file->idx_d->max_file_count; i++)
      free(file->idx_d->block_offset_bitmap[v][i]);
    free(file->idx_d->block_offset_bitmap[v]);
  }
  free(file->idx_d->block_offset_bitmap);

  for (i = 0; i < file->idx_d->max_file_count; i++)
    free(file->idx_d->block_bitmap[i]);
  free(file->idx_d->block_bitmap);

  return PIDX_success;
}

PIDX_return_code destroy_block_layout(PIDX_io file, int gi)
{
  int i = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (i = 0; i < var_grp->shared_layout_count + var_grp->nshared_layout_count; i++)
    free(var_grp->block_layout_by_level[i]);

  free(var_grp->block_layout);
  free(var_grp->block_layout_by_level);

  return PIDX_success;
}
