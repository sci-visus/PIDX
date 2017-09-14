#include "../PIDX_inc.h"

static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int group_index, int start_index, int hz_level_from, int hz_level_to);

static PIDX_return_code populate_idx_layout(PIDX_io file, int gi, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level);

static PIDX_return_code destroy_block_layout(PIDX_io file, int gi);

static PIDX_return_code create_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to);
static PIDX_return_code create_non_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to);

PIDX_return_code populate_global_bit_string(PIDX_io file, int mode)
{
  int i = 0;
  unsigned long long cb[PIDX_MAX_DIMENSIONS];
  unsigned long long* cs = file->idx->chunk_size;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx->bounds[i] % file->idx->chunk_size[i] == 0)
      cb[i] = (int) file->idx->bounds[i] / file->idx->chunk_size[i];
    else
      cb[i] = (int) (file->idx->bounds[i] / file->idx->chunk_size[i]) + 1;
  }

  if (mode == PIDX_WRITE)
  {
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];
    char partition_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx_d->restructured_grid->patch_size[0] / cs[0];
    rpp.y = (int) file->idx_d->restructured_grid->patch_size[1] / cs[1];
    rpp.z = (int) file->idx_d->restructured_grid->patch_size[2] / cs[2];
    guess_bit_string_ZYX(reg_patch_bs, rpp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[1] %s : %d %d %d\n", reg_patch_bs, rpp.x, rpp.y, rpp.z);
#endif

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) file->idx_d->partition_size[0] / file->idx_d->restructured_grid->patch_size[0];
    prcp.y = (int) file->idx_d->partition_size[1] / file->idx_d->restructured_grid->patch_size[1];
    prcp.z = (int) file->idx_d->partition_size[2] / file->idx_d->restructured_grid->patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[2] %s : %d %d %d\n", process_bs, prcp.x, prcp.y, prcp.z);
#endif

    // Last part of the bitstring
    Point3D pcp;
    pcp.x = (int) file->idx_d->partition_count[0];
    pcp.y = (int) file->idx_d->partition_count[1];
    pcp.z = (int) file->idx_d->partition_count[2];
    guess_bit_string(partition_bs, pcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[3] %s : %d %d %d\n", partition_bs, pcp.x, pcp.y, pcp.z);
#endif

    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    strcpy(file->idx->bitSequence, partition_bs);
    strcat(file->idx->bitSequence, temp_bs + 1);
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

#if DETAIL_OUTPUT
  if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
    fprintf(stderr, "Bitstring %s maxh %d\n", file->idx->bitSequence, file->idx_d->maxh);
#endif

  unsigned long long total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  unsigned long long max_sample_per_file = (unsigned long long) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

  int partion_level = (int) log2(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  return PIDX_success;
}

PIDX_return_code populate_local_bit_string(PIDX_io file, int mode)
{
  int i = 0;
  unsigned long long cb[PIDX_MAX_DIMENSIONS];

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx->box_bounds[i] % file->idx->chunk_size[i] == 0)
      cb[i] = (int) file->idx->box_bounds[i] / file->idx->chunk_size[i];
    else
      cb[i] = (int) (file->idx->box_bounds[i] / file->idx->chunk_size[i]) + 1;
  }

  if (mode == PIDX_WRITE)
  {
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx_d->restructured_grid->patch_size[0];
    rpp.y = (int) file->idx_d->restructured_grid->patch_size[1];
    rpp.z = (int) file->idx_d->restructured_grid->patch_size[2];
    guess_bit_string_ZYX(reg_patch_bs, rpp);
    //if (file->idx_c->lrank == 0)
    //  fprintf(stderr, "[1X %d] %s : %d %d %d\n", file->idx_d->color, reg_patch_bs, rpp.x, rpp.y, rpp.z);

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) getPowerOf2(file->idx->box_bounds[0]) / file->idx_d->restructured_grid->patch_size[0];
    prcp.y = (int) getPowerOf2(file->idx->box_bounds[1]) / file->idx_d->restructured_grid->patch_size[1];
    prcp.z = (int) getPowerOf2(file->idx->box_bounds[2]) / file->idx_d->restructured_grid->patch_size[2];
    //prcp.x = (int) file->idx_d->partition_size[0] / file->idx->reg_patch_size[0];
    //prcp.y = (int) file->idx_d->partition_size[1] / file->idx->reg_patch_size[1];
    //prcp.z = (int) file->idx_d->partition_size[2] / file->idx->reg_patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);
    //if (file->idx_c->lrank == 0)
    //  fprintf(stderr, "[2Y %d] %s : %d %d %d\n", file->idx_d->color, process_bs, prcp.x, prcp.y, prcp.z);


    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    //strcpy(file->idx->bitSequence, partition_bs);
    //strcat(file->idx->bitSequence, temp_bs + 1);
    strcpy(file->idx->bitSequence, temp_bs);
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

  //if (file->idx_c->lrank == 0)
  //  fprintf(stderr, "%d Bitstring %s maxh %d\n", file->idx_d->color, file->idx->bitSequence, file->idx_d->maxh);

  unsigned long long total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  unsigned long long max_sample_per_file = (unsigned long long) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

  //if (file->idx_c->lrank == 0)
  //fprintf(stderr, "[%d] MFC %d : %d %d %d (%d %d %d)\n", file->idx_d->color, file->idx_d->max_file_count, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2]);
  file->idx_d->block_bitmap = malloc(file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  memset(file->idx_d->block_bitmap, 0, file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    file->idx_d->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
    memset(file->idx_d->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
  }

  int partion_level = (int) log2(/*file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]*/1);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  if (cb[0] == 0 && cb[1] == 0 && cb[2] == 0)
  {
    file->idx_d->maxh = 0;
    file->idx_d->max_file_count = 0;
  }

  return PIDX_success;
}


PIDX_return_code populate_block_layouts(PIDX_io file, int gi, int svi, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared)
{
  int i = 0, ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  file->idx_d->block_bitmap = malloc(file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  memset(file->idx_d->block_bitmap, 0, file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    file->idx_d->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
    memset(file->idx_d->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
  }

  if (hz_from_shared == hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
    var_grp->shared_layout_count = 0;
  }
  else
    create_shared_block_layout(file, gi, hz_from_shared, hz_to_shared);


  if (hz_from_non_shared == hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
    var_grp->nshared_layout_count = 0;
  }
  else
    create_non_shared_block_layout(file, gi, hz_from_non_shared, hz_to_non_shared);


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

  int lvi = start_var_index;//file->local_variable_index;


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

  if (block_layout->resolution_from <= file->idx->bits_per_block)
  {
    int level_count = 1;
    for (i = block_layout->resolution_from; i <= file->idx->bits_per_block; i++)
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);

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



PIDX_return_code create_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to)
{
  int lower_hz_level = hz_level_from;
  int higher_hz_level = hz_level_to;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  var_grp->shared_start_layout_index = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->shared_start_layout_index <= 0)
    var_grp->shared_start_layout_index = 0;

  var_grp->shared_end_layout_index = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->shared_end_layout_index <= 0)
    var_grp->shared_end_layout_index = 1;

  var_grp->shared_layout_count = var_grp->shared_end_layout_index - var_grp->shared_start_layout_index;

  return PIDX_success;
}


PIDX_return_code create_non_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to)
{
  int lower_hz_level = hz_level_from;
  int higher_hz_level = hz_level_to;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  var_grp->nshared_start_layout_index = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->nshared_start_layout_index <= 0)
    var_grp->nshared_start_layout_index = 0;

  var_grp->nshared_end_layout_index = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->nshared_end_layout_index <= 0)
    var_grp->nshared_end_layout_index = 1;

  var_grp->nshared_layout_count = var_grp->nshared_end_layout_index - var_grp->nshared_start_layout_index;

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

    ret_code = populate_idx_layout(file, gi, si, layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
      memcpy(block_layout->hz_block_number_array[j], layout_by_level[0]->hz_block_number_array[j], sizeof(int));

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

  //if (file->idx_c->grank == 32)
  //  PIDX_blocks_print_layout(block_layout);

  return PIDX_success;
}


PIDX_return_code delete_block_layout(PIDX_io file, int gi, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared)
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
