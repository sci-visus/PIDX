#include "../PIDX_inc.h"

static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int group_index, int start_index, int hz_level_from, int hz_level_to, int io_type);

static PIDX_return_code populate_idx_layout(PIDX_io file, int gi, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level, int io_type);

static PIDX_return_code create_file_zero_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to);
static PIDX_return_code create_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to);
static PIDX_return_code create_non_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to);

static PIDX_return_code destroy_file_zero_block_layout(PIDX_io file, int gi);
static PIDX_return_code destroy_shared_block_layout(PIDX_io file, int gi);
static PIDX_return_code destroy_non_shared_block_layout(PIDX_io file, int gi);

static PIDX_return_code PIDX_file_initialize_time_step(PIDX_io file, char* filename, char* filename_template, int current_time_step);

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
#if 1
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];
    char partition_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx->reg_patch_size[0] / cs[0];
    rpp.y = (int) file->idx->reg_patch_size[1] / cs[1];
    rpp.z = (int) file->idx->reg_patch_size[2] / cs[2];
    guess_bit_string_ZYX(reg_patch_bs, rpp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      printf("[1] %s : %d %d %d\n", reg_patch_bs, rpp.x, rpp.y, rpp.z);
#endif

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) file->idx_d->partition_size[0] / file->idx->reg_patch_size[0];
    prcp.y = (int) file->idx_d->partition_size[1] / file->idx->reg_patch_size[1];
    prcp.z = (int) file->idx_d->partition_size[2] / file->idx->reg_patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    if (file->idx->bitsequence_type == 0)
      guess_bit_string_Z(process_bs, prcp);
    else if (file->idx->bitsequence_type == 1)
      guess_bit_string_Y(process_bs, prcp);
    else if (file->idx->bitsequence_type == 2)
      guess_bit_string_X(process_bs, prcp);
    else if (file->idx->bitsequence_type == 3)
      guess_bit_string_XZY(process_bs, prcp);
    else if (file->idx->bitsequence_type == 4)
      guess_bit_string_YXZ(process_bs, prcp);
    else if (file->idx->bitsequence_type == 5)
      guess_bit_string_ZYX(process_bs, prcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      printf("[2] %s : %d %d %d\n", process_bs, prcp.x, prcp.y, prcp.z);
#endif

    // Last part of the bitstring
    Point3D pcp;
    pcp.x = (int) file->idx_d->partition_count[0];
    pcp.y = (int) file->idx_d->partition_count[1];
    pcp.z = (int) file->idx_d->partition_count[2];
    guess_bit_string(partition_bs, pcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      printf("[3] %s : %d %d %d\n", partition_bs, pcp.x, pcp.y, pcp.z);
#endif

    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    strcpy(file->idx->bitSequence, partition_bs);
    strcat(file->idx->bitSequence, temp_bs + 1);

#else
    Point3D pcp;
    pcp.x = (int) file->idx->bounds[0];
    pcp.y = (int) file->idx->bounds[1];
    pcp.z = (int) file->idx->bounds[2];
    GuessBitmaskPattern(file->idx->bitSequence, pcp);
    //guess_bit_string(file->idx->bitSequence, pcp);
    //if (file->idx_c->grank == 0)
    //  printf("[BS] %s : %d %d %d\n", file->idx->bitSequence, pcp.x, pcp.y, pcp.z);
#endif
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

#if DETAIL_OUTPUT
  if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
    printf("Bitstring %s maxh %d\n", file->idx->bitSequence, file->idx_d->maxh);
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

  file->idx_d->block_bitmap = malloc(file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  memset(file->idx_d->block_bitmap, 0, file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    file->idx_d->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
    memset(file->idx_d->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
  }

  file->idx_d->shared_block_level = (int)log2(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]) + file->idx->bits_per_block + 1;
  if (file->idx_d->shared_block_level >= file->idx_d->maxh)
    file->idx_d->shared_block_level = file->idx_d->maxh;

  int partion_level = (int) log2(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;



  file->idx->random_agg_list = malloc(sizeof(*file->idx->random_agg_list) * file->idx_d->max_file_count * file->idx->variable_count);
  memset(file->idx->random_agg_list, 0, sizeof(*file->idx->random_agg_list) * file->idx_d->max_file_count * file->idx->variable_count);

  file->idx->random_agg_counter = 0;

  if (file->idx_c->grank == 0)
  {
    time_t t;
    srand((unsigned) time(&t));


    int M = file->idx_d->max_file_count * file->idx->variable_count;
    int N = file->idx_c->gnprocs - 1;

    //
    unsigned char *is_used;
    is_used = malloc(sizeof(*is_used) * N);
    memset(is_used, 0, sizeof(*is_used) * N);

    int in, im;
    im = 0;

    for (in = N - M; in < N && im < M; ++in)
    {
      int r = rand() % (in + 1);
      if (is_used[r])
      {
        //printf("RANDOM %d %d ", r, in);
        r = in;
      }

      assert(!is_used[r]);
      file->idx->random_agg_list[im++] = r;
        is_used[r] = 1;
    }

    assert(im == M);
    //

    /*
    for (i = 0; i < file->idx_d->max_file_count * file->idx->variable_count; i++)
      file->idx->random_agg_list[i] = rand();
    for (i = 0; i < file->idx_d->max_file_count * file->idx->variable_count; i++)
      file->idx->random_agg_list[i] = file->idx->random_agg_list[i] % file->idx_c->gnprocs;
    */

    printf("\nAggs: ");
    for (i = 0; i < file->idx_d->max_file_count * file->idx->variable_count; i++)
      printf ("%d ", file->idx->random_agg_list[i]);
    printf("\n");
  }
  MPI_Bcast(file->idx->random_agg_list, (file->idx_d->max_file_count * file->idx->variable_count), MPI_INT, 0, file->idx_c->global_comm);


  if (cb[0] == 0 && cb[1] == 0 && cb[2] == 0)
  {
    file->idx_d->maxh = 0;
    file->idx_d->max_file_count = 0;
  }

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
#if 1
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx->reg_patch_size[0];
    rpp.y = (int) file->idx->reg_patch_size[1];
    rpp.z = (int) file->idx->reg_patch_size[2];
    guess_bit_string_ZYX(reg_patch_bs, rpp);
    //if (file->idx_c->lrank == 0)
    //  printf("[1X %d] %s : %d %d %d\n", file->idx_d->color, reg_patch_bs, rpp.x, rpp.y, rpp.z);

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) getPowerOf2(file->idx->box_bounds[0]) / file->idx->reg_patch_size[0];
    prcp.y = (int) getPowerOf2(file->idx->box_bounds[1]) / file->idx->reg_patch_size[1];
    prcp.z = (int) getPowerOf2(file->idx->box_bounds[2]) / file->idx->reg_patch_size[2];
    //prcp.x = (int) file->idx_d->partition_size[0] / file->idx->reg_patch_size[0];
    //prcp.y = (int) file->idx_d->partition_size[1] / file->idx->reg_patch_size[1];
    //prcp.z = (int) file->idx_d->partition_size[2] / file->idx->reg_patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);
    //if (file->idx_c->lrank == 0)
    //  printf("[2Y %d] %s : %d %d %d\n", file->idx_d->color, process_bs, prcp.x, prcp.y, prcp.z);


    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    //strcpy(file->idx->bitSequence, partition_bs);
    //strcat(file->idx->bitSequence, temp_bs + 1);
    strcpy(file->idx->bitSequence, temp_bs);

#else
    Point3D pcp;
    pcp.x = (int) file->idx->bounds[0];
    pcp.y = (int) file->idx->bounds[1];
    pcp.z = (int) file->idx->bounds[2];
    GuessBitmaskPattern(file->idx->bitSequence, pcp);
    //guess_bit_string(file->idx->bitSequence, pcp);
    //if (file->idx_c->grank == 0)
    //  printf("[BS] %s : %d %d %d\n", file->idx->bitSequence, pcp.x, pcp.y, pcp.z);
#endif
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

  //if (file->idx_c->lrank == 0)
  //  printf("%d Bitstring %s maxh %d\n", file->idx_d->color, file->idx->bitSequence, file->idx_d->maxh);

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
  //printf("[%d] MFC %d : %d %d %d (%d %d %d)\n", file->idx_d->color, file->idx_d->max_file_count, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2]);
  file->idx_d->block_bitmap = malloc(file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  memset(file->idx_d->block_bitmap, 0, file->idx_d->max_file_count * sizeof (*file->idx_d->block_bitmap));
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    file->idx_d->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
    memset(file->idx_d->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_d->block_bitmap[i]));
  }

  file->idx_d->shared_block_level = (int)log2(/*file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]*/1) + file->idx->bits_per_block + 1;
  if (file->idx_d->shared_block_level >= file->idx_d->maxh)
    file->idx_d->shared_block_level = file->idx_d->maxh;

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


PIDX_return_code populate_block_layouts(PIDX_io file, int gi, int svi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared, int io_type)
{
  int ret;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->filename_template, file->idx->current_time_step);
  PIDX_file_initialize_time_step(file, file->idx->filename_global, file->idx->filename_template_global, file->idx->current_time_step);
  PIDX_file_initialize_time_step(file, file->idx->filename_partition, file->idx->filename_template_partition, file->idx->current_time_step);
  PIDX_file_initialize_time_step(file, file->idx->filename_file_zero, file->idx->filename_template_file_zero, file->idx->current_time_step);

  if (hz_from_file_zero == hz_to_file_zero)
  {
    var_grp->f0_start_layout_index = 0;
    var_grp->f0_end_layout_index = 0;
    var_grp->f0_layout_count = 0;
  }
  else
  {
    create_file_zero_block_layout(file, gi, hz_from_file_zero, hz_to_file_zero);
    ret = populate_idx_block_layout(file,
                               var_grp->f0_block_layout,
                               var_grp->f0_block_layout_by_level,
                               var_grp->f0_start_layout_index,
                               var_grp->f0_end_layout_index,
                               var_grp->f0_layout_count,
                               gi,
                               svi,
                               hz_from_file_zero, hz_to_file_zero,
                               io_type);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
#if 1
  if (hz_from_shared == hz_to_shared)
  {
    var_grp->shared_start_layout_index = 0;
    var_grp->shared_end_layout_index = 0;
    var_grp->shared_layout_count = 0;
  }
  else
  {
    create_shared_block_layout(file, gi, hz_from_shared, hz_to_shared);
    ret = populate_idx_block_layout(file,
                               var_grp->shared_block_layout, var_grp->shared_block_layout_by_level,
                               var_grp->shared_start_layout_index, var_grp->shared_end_layout_index,
                               var_grp->shared_layout_count,
                               gi,
                               svi,
                               hz_from_shared, hz_to_shared,
                               io_type);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
#endif
#if 1

  //if (file->idx_c->grank == 32)
  //  printf("HZ %d %d -- %d %d\n", hz_from_shared, hz_to_shared, hz_from_non_shared, hz_to_non_shared);

  if (hz_from_non_shared == hz_to_non_shared)
  {
    var_grp->nshared_start_layout_index = 0;
    var_grp->nshared_end_layout_index = 0;
    var_grp->nshared_layout_count = 0;
  }
  else
  {
    create_non_shared_block_layout(file, gi, hz_from_non_shared, hz_to_non_shared);
    ret = populate_idx_block_layout(file,
                               var_grp->nshared_block_layout, var_grp->nshared_block_layout_by_level,
                               var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index,
                               var_grp->nshared_layout_count,
                               gi,
                               svi,
                               hz_from_non_shared, hz_to_non_shared,
                               io_type);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    } 
  }
#endif

  return PIDX_success;
}


static PIDX_return_code populate_idx_layout(PIDX_io file, int gi, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level, int io_type)
{
  int i, j;
  int p = 0, ctr = 1;
  PIDX_return_code ret_code;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };
#if 1
  int lvi = start_var_index;//file->local_variable_index;
#if 1
  if (file->idx_d->parallel_mode == 1 /*&& file->idx->compression_type == PIDX_NO_COMPRESSION*/)
  {
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

    if (io_type == PIDX_IDX_IO)
    {
      for (p = 0 ; p < var->sim_patch_count ; p++)
      //for (p = 0 ; p < var->patch_group_count ; p++)
      {
        for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
        {
          bounding_box[0][i] = var->sim_patch[p]->offset[i];
          bounding_box[1][i] = var->sim_patch[p]->size[i] + var->sim_patch[p]->offset[i];

          //bounding_box[0][i] = var->rst_patch_group[p]->reg_patch->offset[i];
          //bounding_box[1][i] = var->rst_patch_group[p]->reg_patch->size[i] + var->rst_patch_group[p]->reg_patch->offset[i];


          bounding_box[0][i] = (bounding_box[0][i] / file->idx->chunk_size[i]);

          if (bounding_box[1][i] % file->idx->chunk_size[i] == 0)
            bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]);
          else
            bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]) + 1;
        }
        //printf("BB ---- [%d %d]: %d %d %d - %d %d %d\n", rank, var->sim_patch_count, bounding_box[0][0], bounding_box[0][1], bounding_box[0][2], bounding_box[1][0], bounding_box[1][1], bounding_box[1][2]);

        PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
        memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
        ret_code = PIDX_blocks_initialize_layout(per_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
        if (ret_code != PIDX_success)
        {
          fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bitPattern, per_patch_local_block_layout, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
        if (ret_code != PIDX_success)
        {
          fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        if (all_patch_local_block_layout->resolution_from <= all_patch_local_block_layout->bits_per_block)
        {
          for (i = all_patch_local_block_layout->resolution_from ; i <= all_patch_local_block_layout->bits_per_block ; i++)
          {
            if (per_patch_local_block_layout->hz_block_number_array[i][0] == 0)
            {
              all_patch_local_block_layout->hz_block_number_array[i][0] = per_patch_local_block_layout->hz_block_number_array[i][0];
              break;
            }
          }

          ctr = 1;
          for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
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
          for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
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

        PIDX_blocks_free_layout(per_patch_local_block_layout);
        free(per_patch_local_block_layout);
        per_patch_local_block_layout = 0;
      }
    }
    else
    {
      assert(var->patch_group_count <= 1);
      //for (p = 0 ; p < var->sim_patch_count ; p++)
      for (p = 0 ; p < var->patch_group_count ; p++)
      {
        for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
        {
          //bounding_box[0][i] = var->sim_patch[p]->offset[i];
          //bounding_box[1][i] = var->sim_patch[p]->size[i] + var->sim_patch[p]->offset[i];

          bounding_box[0][i] = var->rst_patch_group[0]->reg_patch->offset[i];
          bounding_box[1][i] = var->rst_patch_group[0]->reg_patch->size[i] + var->rst_patch_group[0]->reg_patch->offset[i];


          bounding_box[0][i] = (bounding_box[0][i] / file->idx->chunk_size[i]);

          if (bounding_box[1][i] % file->idx->chunk_size[i] == 0)
            bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]);
          else
            bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]) + 1;
        }
        //if (file->idx_d->color == 1)
        //  printf("A: [%d %d]: %d %d %d - %d %d %d\n", file->idx_d->color, var->sim_patch_count, bounding_box[0][0], bounding_box[0][1], bounding_box[0][2], bounding_box[1][0], bounding_box[1][1], bounding_box[1][2]);

        PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
        memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
        ret_code = PIDX_blocks_initialize_layout(per_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
        if (ret_code != PIDX_success)
        {
          fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bitPattern, per_patch_local_block_layout, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
        if (ret_code != PIDX_success)
        {
          fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        if (all_patch_local_block_layout->resolution_from <= all_patch_local_block_layout->bits_per_block)
        {
          for (i = all_patch_local_block_layout->resolution_from ; i <= all_patch_local_block_layout->bits_per_block ; i++)
          {
            if (per_patch_local_block_layout->hz_block_number_array[i][0] == 0)
            {
              all_patch_local_block_layout->hz_block_number_array[i][0] = per_patch_local_block_layout->hz_block_number_array[i][0];
              break;
            }
          }

          ctr = 1;
          for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
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
          for (i = all_patch_local_block_layout->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
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

        PIDX_blocks_free_layout(per_patch_local_block_layout);
        free(per_patch_local_block_layout);
        per_patch_local_block_layout = 0;
      }
    }
#if 1
    if (block_layout->resolution_from <= block_layout->bits_per_block)
    {
      int level_count = 1;
      for (i = block_layout->resolution_from; i <= block_layout->bits_per_block; i++)
      {
#if PIDX_HAVE_MPI
        if (file->idx_d->parallel_mode == 1)
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
        else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
        memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
      }

      for (i = block_layout->bits_per_block + 1; i < (block_layout->resolution_to); i++)
      {
#if PIDX_HAVE_MPI
        if (file->idx_d->parallel_mode == 1)
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
        else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
        memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
        level_count = level_count * 2;
      }
    }
    else
    {
#if 1
      int level_count = 1;
      for (i = block_layout->bits_per_block + 1; i < (block_layout->resolution_to); i++)
      {
        if (i >= block_layout->resolution_from)
        {
#if PIDX_HAVE_MPI
          //if (file->idx_c->lrank == 0)
          //{
          //  int k = 0;
          //  for (k = 0; k < level_count; k++)
          //    printf("i = %d %d\n", i, block_layout->hz_block_number_array[i][k], all_patch_local_block_layout->hz_block_number_array[i][k]);
          //}
          if (file->idx_d->parallel_mode == 1)
            MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->local_comm);
          else
            memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
        }
        level_count = level_count * 2;
      }
#endif
    }

    /*
    if (rank == 4)
    {
      printf("[XXXX] Final Block Bitmap\n");
      PIDX_blocks_print_layout(all_patch_local_block_layout);
    }
    */

    PIDX_blocks_free_layout(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
#endif
  }
  else
#endif
  {
    unsigned long long cb[PIDX_MAX_DIMENSIONS];
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      if (file->idx->bounds[i] % file->idx->chunk_size[i] == 0)
        cb[i] = (int) file->idx->bounds[i] / file->idx->chunk_size[i];
      else
        cb[i] = (int) (file->idx->bounds[i] / file->idx->chunk_size[i]) + 1;
    }

    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      bounding_box[0][i] = 0;
      bounding_box[1][i] = cb[i];
    }

    if (file->idx_c->grank == 0)
      printf("BB %d %d %d\n", cb[0], cb[1], cb[2]);

    ret_code = PIDX_blocks_create_layout (bounding_box, file->idx_d->maxh, file->idx->bitPattern, block_layout, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
      return PIDX_err_file;
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

  int file_number = 0;
  //printf("[AAAA] RES from TO to : %d %d\n", block_layout->resolution_from, block_layout->resolution_to);
  if (block_layout->resolution_from <= block_layout->bits_per_block)
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
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
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

    /*
    if (file->idx_c->lnprocsrank == 4)
    {
      printf("[XXXX] Final Block Bitmap\n");
      PIDX_blocks_print_layout(block_layout);
    }
    */
  }
  else
  {
#if 1
    ctr = 1;
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
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
#endif
  }
#if 1
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
#endif
#endif
  return PIDX_success;
}


PIDX_return_code create_file_zero_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to)
{
  int i_1 = 0, i = 0;
  int lower_hz_level = hz_level_from;
  int higher_hz_level = hz_level_to;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  var_grp->f0_start_layout_index = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->f0_start_layout_index <= 0)
    var_grp->f0_start_layout_index = 0;

  var_grp->f0_end_layout_index = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (var_grp->f0_end_layout_index <= 0)
    var_grp->f0_end_layout_index = 1;

  var_grp->f0_layout_count = var_grp->f0_end_layout_index - var_grp->f0_start_layout_index;

  var_grp->f0_block_layout = malloc(sizeof (*var_grp->f0_block_layout));
  memset(var_grp->f0_block_layout, 0, sizeof (*var_grp->f0_block_layout));

  var_grp->f0_block_layout_by_level = malloc(sizeof (*var_grp->f0_block_layout_by_level) * var_grp->f0_layout_count);
  memset(var_grp->f0_block_layout_by_level, 0, sizeof (*var_grp->f0_block_layout_by_level) * var_grp->f0_layout_count);

  for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
  {
    i_1 = i - var_grp->f0_start_layout_index;
    var_grp->f0_block_layout_by_level[i_1] = malloc(sizeof(*(var_grp->f0_block_layout_by_level[i_1])));
    memset(var_grp->f0_block_layout_by_level[i_1], 0, sizeof(*(var_grp->f0_block_layout_by_level[i_1])));
  }

  return PIDX_success;
}


PIDX_return_code create_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to)
{
  int i_1 = 0, i = 0;
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

  var_grp->shared_block_layout = malloc(sizeof (*var_grp->shared_block_layout));
  memset(var_grp->shared_block_layout, 0, sizeof (*var_grp->shared_block_layout));

  var_grp->shared_block_layout_by_level = malloc(sizeof (*var_grp->shared_block_layout_by_level) * var_grp->shared_layout_count);
  memset(var_grp->shared_block_layout_by_level, 0, sizeof (*var_grp->shared_block_layout_by_level) * var_grp->shared_layout_count);

  for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
  {
    i_1 = i - var_grp->shared_start_layout_index;
    var_grp->shared_block_layout_by_level[i_1] = malloc(sizeof(*(var_grp->shared_block_layout_by_level[i_1])));
    memset(var_grp->shared_block_layout_by_level[i_1], 0, sizeof(*(var_grp->shared_block_layout_by_level[i_1])));
  }

  return PIDX_success;
}


PIDX_return_code create_non_shared_block_layout(PIDX_io file, int gi, int hz_level_from, int hz_level_to)
{
  int i_1 = 0, i = 0;
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

  var_grp->nshared_block_layout = malloc(sizeof (*var_grp->nshared_block_layout));
  memset(var_grp->nshared_block_layout, 0, sizeof (*var_grp->nshared_block_layout));

  var_grp->nshared_block_layout_by_level = malloc(sizeof (*var_grp->nshared_block_layout_by_level) * var_grp->nshared_layout_count);
  memset(var_grp->nshared_block_layout_by_level, 0, sizeof (*var_grp->nshared_block_layout_by_level) * var_grp->nshared_layout_count);

  for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
  {
    i_1 = i - var_grp->nshared_start_layout_index;
    var_grp->nshared_block_layout_by_level[i_1] = malloc(sizeof(*(var_grp->nshared_block_layout_by_level[i_1])));
    memset(var_grp->nshared_block_layout_by_level[i_1], 0, sizeof(*(var_grp->nshared_block_layout_by_level[i_1])));
  }

  return PIDX_success;
}



static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int gi, int si, int hz_level_from, int hz_level_to, int io_type)
{
  PIDX_return_code ret_code;

  int i = 0, j = 0, ctr;
  int file_number = 0;
  int lower_hz_level = 0, higher_hz_level = 0;
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  PIDX_block_layout block_layout = global_layout;

  lower_hz_level = hz_level_from;
  higher_hz_level = hz_level_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (hz_level_from == 0 && hz_level_to == 0)
    return PIDX_success;

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

    ret_code = populate_idx_layout(file, gi, si, layout_by_level[0], lower_level_low_layout, higher_level_low_layout, io_type);
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

      ret_code = populate_idx_layout(file, gi, si, layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout, io_type);
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

      ret_code = populate_idx_layout(file, gi, si, layout_by_level[i - start_layout_index], lower_level_higher_layout, higher_level_higher_layout, io_type);
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

  if (block_layout->resolution_from <= block_layout->bits_per_block)
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
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
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
    for (i = block_layout->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          //if (rank == 0)
          //  printf("ctr = %d %d\n", ctr, block_layout->hz_block_number_array[i][j]);
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            //if (file->idx_c->lrank == 0)
            //  printf("%d file number %d %d \n", file->idx_c->grank, block_layout->hz_block_number_array[i][j], file_number);
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
      //if (rank == 0)
      //  printf("[%d %d] BPF %d = %d %d\n", block_layout->resolution_from, block_layout->resolution_to, i, block_layout->bcpf[i], count);
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  //if (file->idx_c->grank == 32)
  //{
  //  PIDX_blocks_print_layout(block_layout);
  //}

  return PIDX_success;
}


PIDX_return_code delete_block_layout(PIDX_io file, int gi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared)
{
  int i, i_1;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  if (hz_from_file_zero != hz_to_file_zero)
  {
    PIDX_free_layout(var_grp->f0_block_layout);
    PIDX_blocks_free_layout(var_grp->f0_block_layout);

    for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
    {
      i_1 = i - var_grp->f0_start_layout_index;
      PIDX_free_layout(var_grp->f0_block_layout_by_level[i_1]);
      PIDX_blocks_free_layout(var_grp->f0_block_layout_by_level[i_1]);
    }
    destroy_file_zero_block_layout(file, gi);
  }

  if (hz_from_shared != hz_to_shared)
  {
    PIDX_free_layout(var_grp->shared_block_layout);
    PIDX_blocks_free_layout(var_grp->shared_block_layout);

    for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
    {
      i_1 = i - var_grp->shared_start_layout_index;
      PIDX_free_layout(var_grp->shared_block_layout_by_level[i_1]);
      PIDX_blocks_free_layout(var_grp->shared_block_layout_by_level[i_1]);
    }
    destroy_shared_block_layout(file, gi);
  }

  if (hz_from_non_shared != hz_to_non_shared)
  {
    PIDX_free_layout(var_grp->nshared_block_layout);
    PIDX_blocks_free_layout(var_grp->nshared_block_layout);

    for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
    {
      i_1 = i - var_grp->nshared_start_layout_index;
      PIDX_free_layout(var_grp->nshared_block_layout_by_level[i_1]);
      PIDX_blocks_free_layout(var_grp->nshared_block_layout_by_level[i_1]);
    }
    destroy_non_shared_block_layout(file, gi);
  }

  return PIDX_success;
}

PIDX_return_code destroy_file_zero_block_layout(PIDX_io file, int gi)
{
  int i_1 = 0, i = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (i = var_grp->f0_start_layout_index; i < var_grp->f0_end_layout_index ; i++)
  {
    i_1 = i - var_grp->f0_start_layout_index;
    free(var_grp->f0_block_layout_by_level[i_1]);
  }

  free(var_grp->f0_block_layout);
  free(var_grp->f0_block_layout_by_level);

  return PIDX_success;
}


PIDX_return_code destroy_shared_block_layout(PIDX_io file, int gi)
{
  int i_1 = 0, i = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (i = var_grp->shared_start_layout_index; i < var_grp->shared_end_layout_index ; i++)
  {
    i_1 = i - var_grp->shared_start_layout_index;
    free(var_grp->shared_block_layout_by_level[i_1]);
  }

  free(var_grp->shared_block_layout);
  free(var_grp->shared_block_layout_by_level);

  return PIDX_success;
}


PIDX_return_code destroy_non_shared_block_layout(PIDX_io file, int gi)
{
  int i_1 = 0, i = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (i = var_grp->nshared_start_layout_index; i < var_grp->nshared_end_layout_index ; i++)
  {
    i_1 = i - var_grp->nshared_start_layout_index;
    free(var_grp->nshared_block_layout_by_level[i_1]);
  }

  free(var_grp->nshared_block_layout);
  free(var_grp->nshared_block_layout_by_level);

  return PIDX_success;
}


/// TODO: get rid of this function
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_io file, char* filename, char* filename_template, int current_time_step)
{
  int N;
  char dirname[1024], basename[1024];
  int nbits_blocknumber;
  char *directory_path;
  char *data_set_path;

  data_set_path = malloc(sizeof(*data_set_path) * 1024);
  memset(data_set_path, 0, sizeof(*data_set_path) * 1024);

  directory_path = malloc(sizeof(*directory_path) * 1024);
  memset(directory_path, 0, sizeof(*directory_path) * 1024);

  strncpy(directory_path, filename, strlen(filename) - 4);
  sprintf(data_set_path, "%s/time%09d.idx", directory_path, current_time_step);
  free(directory_path);

  nbits_blocknumber = (file->idx_d->maxh - file->idx->bits_per_block - 1);
  VisusSplitFilename(data_set_path, dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--)
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

#if 0
  //if i put . as the first character, if I move files VisusOpen can do path remapping
  sprintf(filename_template, "./%s", basename);
#endif
  //pidx does not do path remapping
  strcpy(filename_template, data_set_path);
  for (N = strlen(filename_template) - 1; N >= 0; N--)
  {
    int ch = filename_template[N];
    filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0)
    strcat(filename_template, "/%01x.bin");

  else
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4)
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8)
      strcat(filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12)
      strcat(filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16)
      strcat(filename_template, "/%04x.bin"); //no directories, 65536  files
    else
    {
      while (nbits_blocknumber > 16)
      {
        strcat(filename_template, "/%02x"); //256 subdirectories
        nbits_blocknumber -= 8;
      }
      strcat(filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }

  free(data_set_path);
  return PIDX_success;
}
