#include "../PIDX_io.h"

//static int regular_bounds[PIDX_MAX_DIMENSIONS] = {256, 256, 128, 1, 1};
static PIDX_return_code populate_idx_layout(PIDX_global_idx_io file, int start_var_index, int end_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level);
//static PIDX_return_code delete_idx_dataset(PIDX_global_idx_io file, int start_var_index, int end_var_index);
//static PIDX_return_code delete_idx_dataset(PIDX_global_idx_io file, int start_index, int end_index, int start_hz_level, int end_hz_level);
//static PIDX_return_code populate_idx_dataset(PIDX_global_idx_io file, int start_var_index, int end_var_index, int start_layout_index, int end_layout_index);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_global_idx_io file, char* file_name, int current_time_step);
static PIDX_return_code partition_setup(PIDX_global_idx_io file, int start_var_index, int end_var_index, int pipe_length);
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);
static void set_default_patch_size(PIDX_global_idx_io file, int64_t* process_bounds, int nprocs);
static int getPowerOftwo(int x);

struct PIDX_global_idx_io_descriptor
{

#if PIDX_HAVE_MPI
  MPI_Comm global_comm;                               ///< MPI sub-communicator (including all processes per IDX file)
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
#endif

  PIDX_header_io_id header_io_id;              ///< IDX metadata id
  PIDX_rst_id rst_id;                          ///< Restructuring phase id
  PIDX_chunk_id chunk_id;              ///< Block restructuring id (prepration for compression)
  PIDX_comp_id comp_id;          ///< Compression (lossy and lossless) id
  PIDX_hz_encode_id hz_id;                     ///< HZ encoding phase id
  PIDX_agg_id agg_id;                          ///< Aggregation phase id
  PIDX_agg_id** tagg_id;                          ///< Aggregation phase id
  PIDX_file_io_id io_id;                            ///< IO phase id
  PIDX_file_io_id** tio_id;                            ///< IO phase id

  //int flush_used;
  //int write_on_close;                          ///< HPC Writes
  int one_time_initializations;                ///<


  int small_agg_comm;

  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;


  //PIDX_time time;
};

static int getPowerOftwo(int x)
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}



static PIDX_return_code populate_idx_file_structure(PIDX_global_idx_io file, int partition_level)
{
  PointND bounds_point;
  int d = 0, i = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    if (file->idx->bounds[d] % file->idx->chunk_size[d] == 0)
      file->idx->chunked_bounds[d] = (int) file->idx->bounds[d] / file->idx->chunk_size[d];
    else
      file->idx->chunked_bounds[d] = (int) (file->idx->bounds[d] / file->idx->chunk_size[d]) + 1;
  }
  //int rank;
  //MPI_Comm_rank(file->comm, &rank);

  int64_t* cb = file->idx->chunked_bounds;
  bounds_point.x = (int) cb[0];
  bounds_point.y = (int) cb[1];
  bounds_point.z = (int) cb[2];
  bounds_point.u = (int) cb[3];
  bounds_point.v = (int) cb[4];
  /*
  GuessBitmaskPattern(file->idx->bitSequence, bounds_point);
  printf("Initial %s %d\n", file->idx->bitSequence, strlen(file->idx->bitSequence));
  */

  if (partition_level == 0)
  {
    Point3D idx_g_point;
    idx_g_point.x = (int) file->idx_d->idx_count[0];
    idx_g_point.y = (int) file->idx_d->idx_count[1];
    idx_g_point.z = (int) file->idx_d->idx_count[2];
    //GuessBitmaskPattern(file->idx->idx_cg_bitSequence, idx_g_point);
    guess_bit_string(file->idx->idx_cg_bitSequence, idx_g_point);
    //if (rank == 0)
    //  printf("PARTITION %s\n", file->idx->idx_cg_bitSequence);

#if 0
    PointND idx_l_point;
    idx_l_point.x = (int) file->idx_d->idx_size[0];
    idx_l_point.y = (int) file->idx_d->idx_size[1];
    idx_l_point.z = (int) file->idx_d->idx_size[2];
    idx_l_point.u = (int) file->idx_d->idx_size[3];
    idx_l_point.v = (int) file->idx_d->idx_size[4];
    GuessBitmaskPattern(file->idx->idx_cl_bitSequence, idx_l_point);
    //printf("Local %s\n", file->idx->idx_cl_bitSequence);

    //merge_bitmask(file->idx->idx_cg_bitSequence, file->idx->idx_cl_bitSequence, file->idx->bitSequence);
    strcpy(file->idx->bitSequence, file->idx->idx_cg_bitSequence);
    strcat(file->idx->bitSequence, file->idx->idx_cl_bitSequence + 1);
    //printf("Final %s %d\n", file->idx->bitSequence, strlen(file->idx->bitSequence));
#else
    Point3D idx_l1_point;
    idx_l1_point.x = (int) file->idx_d->idx_size[0] / file->idx->reg_patch_size[0];
    idx_l1_point.y = (int) file->idx_d->idx_size[1] / file->idx->reg_patch_size[1];
    idx_l1_point.z = (int) file->idx_d->idx_size[2] / file->idx->reg_patch_size[2];
    //GuessBitmaskPattern(file->idx->idx_cl1_bitSequence, idx_l1_point);
    guess_bit_string2(file->idx->idx_cl1_bitSequence, idx_l1_point);
    //if (rank == 0)
    //  printf("PC %d %d %d: %s\n", idx_l1_point.x, idx_l1_point.y, idx_l1_point.z, file->idx->idx_cl1_bitSequence);

    Point3D idx_l2_point;
    idx_l2_point.x = (int) file->idx->reg_patch_size[0];
    idx_l2_point.y = (int) file->idx->reg_patch_size[1];
    idx_l2_point.z = (int) file->idx->reg_patch_size[2];
    //GuessBitmaskPattern(file->idx->idx_cl2_bitSequence, idx_l2_point);
    guess_bit_string(file->idx->idx_cl2_bitSequence, idx_l2_point);

    //if (rank == 0)
    //  printf("CC %d %d %d: %s\n", idx_l2_point.x, idx_l2_point.y, idx_l2_point.z, file->idx->idx_cl2_bitSequence);

    strcpy(file->idx->idx_cl_bitSequence, file->idx->idx_cl1_bitSequence);
    strcat(file->idx->idx_cl_bitSequence, file->idx->idx_cl2_bitSequence + 1);
    //printf("FB %s\n", file->idx->idx_cl_bitSequence);

    int tempH = strlen(file->idx->idx_cl_bitSequence);
    for (i = 0; i <= tempH; i++)
      file->idx->idx_cl_bitPattern[i] = RegExBitmaskBit(file->idx->idx_cl_bitSequence, i);

    strcpy(file->idx->bitSequence, file->idx->idx_cg_bitSequence);
    strcat(file->idx->bitSequence, file->idx->idx_cl_bitSequence + 1);
    //printf("CB %s\n", file->idx->bitSequence);
#endif
  }
  else
  {
    //file->idx->bitSequence = "";
    memset(file->idx->bitSequence, 0, 512);
    memset(file->idx->bitPattern, 0, 512);
    GuessBitmaskPattern(file->idx->bitSequence, bounds_point);
    //printf("Initial %s %d\n", file->idx->bitSequence, strlen(file->idx->bitSequence));
  }


  file->idx_d->maxh = strlen(file->idx->bitSequence);

  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

  int64_t total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]) * getPowerOf2(cb[3]) * getPowerOf2(cb[4]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int64_t max_sample_per_file = (uint64_t) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

  if (cb[0] == 0 && cb[1] == 0 && cb[2] == 0)
  {
    file->idx_d->maxh = 0;
    file->idx_d->max_file_count = 0;
  }

  return PIDX_success;
}


static PIDX_return_code populate_idx_layout(PIDX_global_idx_io file, int start_var_index, int end_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level)
{
  int i, j;
  int p = 0, ctr = 1;
  PIDX_return_code ret_code;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  int lvi = start_var_index;//file->local_variable_index;

  if (file->idx_d->parallel_mode == 1 && file->idx->compression_type == PIDX_NO_COMPRESSION)
  {
    PIDX_block_layout all_patch_local_block_layout = malloc(sizeof (*all_patch_local_block_layout));
    memset(all_patch_local_block_layout, 0, sizeof (*all_patch_local_block_layout));
    ret_code = PIDX_blocks_initialize_layout(all_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (p = 0 ; p < file->idx->variable[lvi]->sim_patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
        bounding_box[0][i] = file->idx->variable[lvi]->sim_patch[p]->offset[i];
        bounding_box[1][i] = file->idx->variable[lvi]->sim_patch[p]->size[i] + file->idx->variable[lvi]->sim_patch[p]->offset[i];

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

    if (block_layout->resolution_from <= block_layout->bits_per_block)
    {
      int level_count = 1;
      for (i = block_layout->resolution_from; i <= block_layout->bits_per_block; i++)
      {
#if PIDX_HAVE_MPI
        if (file->idx_d->parallel_mode == 1)
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
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
          MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
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
      int level_count = 1;
      for (i = block_layout->bits_per_block + 1; i < (block_layout->resolution_to); i++)
      {
        if (i >= block_layout->resolution_from)
        {
#if PIDX_HAVE_MPI
          if (file->idx_d->parallel_mode == 1)
            MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->comm);
          else
            memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#else
          memcpy(block_layout->hz_block_number_array[i], all_patch_local_block_layout->hz_block_number_array[i], level_count * sizeof(int));
#endif
        }
        level_count = level_count * 2;
      }
    }

    PIDX_blocks_free_layout(all_patch_local_block_layout);
    free(all_patch_local_block_layout);
    all_patch_local_block_layout = 0;
  }
  else
  {
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    {
      bounding_box[0][i] = 0;
      bounding_box[1][i] = file->idx->chunked_bounds[i];
    }

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

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));


  int file_number = 0;
  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
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
          block_layout->block_count_per_file[file_number]++;
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
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }

  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

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


#if 0
static PIDX_return_code delete_idx_dataset(PIDX_global_idx_io file, int start_index, int end_index, int start_hz_level, int end_hz_level)
{
  int lvi = start_index;//file->local_variable_index;
  PIDX_variable var = file->idx->variable[lvi];

  PIDX_blocks_free_layout(file->idx->variable[lvi]->global_block_layout);
  PIDX_free_layout(file->idx->variable[lvi]->global_block_layout);

  free(file->idx->variable[lvi]->global_block_layout);
  file->idx->variable[lvi]->global_block_layout = 0;

  int i = 0;

  for (i = file->idx_d->start_layout_index; i < file->idx_d->end_layout_index ; i++)
  {
    PIDX_blocks_free_layout(var->block_layout_by_level[i - file->idx_d->start_layout_index]);
    PIDX_free_layout(var->block_layout_by_level[i - file->idx_d->start_layout_index]);

    free(var->block_layout_by_level[i - file->idx_d->start_layout_index]);
    var->block_layout_by_level[i - file->idx_d->start_layout_index] = 0;
  }


  free(var->block_layout_by_level);
  var->block_layout_by_level = 0;

  return PIDX_success;
}
#endif


static PIDX_return_code delete_idx_dataset_shared(PIDX_global_idx_io file, int start_index, int end_index, int start_hz_level, int end_hz_level)
{
  int lvi = start_index;//file->local_variable_index;
  PIDX_variable var = file->idx->variable[lvi];

  if (file->idx_d->start_layout_index_shared == file->idx_d->end_layout_index_shared)
    return PIDX_success;


  PIDX_blocks_free_layout(file->idx->variable[lvi]->global_block_layout_files[0]);
  PIDX_free_layout(file->idx->variable[lvi]->global_block_layout_files[0]);

  free(file->idx->variable[lvi]->global_block_layout_files[0]);
  file->idx->variable[lvi]->global_block_layout_files[0] = 0;


  int i = 0;
  for (i = file->idx_d->start_layout_index_shared; i < file->idx_d->end_layout_index_shared ; i++)
  {
    PIDX_blocks_free_layout(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared]);
    PIDX_free_layout(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared]);

    free(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared]);
    var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared] = 0;
  }

  free(var->block_layout_by_level_files[0]);
  var->block_layout_by_level_files[0] = 0;

  return PIDX_success;
}


static PIDX_return_code delete_idx_dataset_non_shared(PIDX_global_idx_io file, int start_index, int end_index, int start_hz_level, int end_hz_level)
{
  int lvi = start_index;//file->local_variable_index;
  PIDX_variable var = file->idx->variable[lvi];

  if (file->idx_d->start_layout_index_non_shared == file->idx_d->end_layout_index_non_shared)
    return PIDX_success;

  PIDX_blocks_free_layout(file->idx->variable[lvi]->global_block_layout_files[1]);
  PIDX_free_layout(file->idx->variable[lvi]->global_block_layout_files[1]);

  free(file->idx->variable[lvi]->global_block_layout_files[1]);
  file->idx->variable[lvi]->global_block_layout_files[1] = 0;

  int i = 0;

  for (i = file->idx_d->start_layout_index_non_shared; i < file->idx_d->end_layout_index_non_shared ; i++)
  {
    PIDX_blocks_free_layout(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared]);
    PIDX_free_layout(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared]);

    free(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared]);
    var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared] = 0;
  }

  free(var->block_layout_by_level_files[1]);
  var->block_layout_by_level_files[1] = 0;

  return PIDX_success;
}




/// TODO: get rid of this function
PIDX_return_code PIDX_file_initialize_time_step(PIDX_global_idx_io file, char* filename, int current_time_step)
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
  sprintf(pidx->filename_template, "./%s", basename);
#endif
  //pidx does not do path remapping
  strcpy(file->idx->filename_template, data_set_path);
  for (N = strlen(file->idx->filename_template) - 1; N >= 0; N--)
  {
    int ch = file->idx->filename_template[N];
    file->idx->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0)
    strcat(file->idx->filename_template, "/%01x.bin");

  else
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4)
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8)
      strcat(file->idx->filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12)
      strcat(file->idx->filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16)
      strcat(file->idx->filename_template, "/%04x.bin"); //no directories, 65536  files
    else
    {
      while (nbits_blocknumber > 16)
      {
    strcat(file->idx->filename_template, "/%02x"); //256 subdirectories
    nbits_blocknumber -= 8;
      }
      strcat(file->idx->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }

  free(data_set_path);
  return PIDX_success;
}


static PIDX_return_code populate_idx_dataset_shared(PIDX_global_idx_io file, int start_index, int end_index, int hz_level_from, int hz_level_to)
{
  int i = 0, j = 0, ctr;
  int file_number = 0;

  int rank = 0;
  int nprocs = 1;

  if (hz_level_from == hz_level_to)
  {
    file->idx_d->start_layout_index_shared = 0;
    file->idx_d->end_layout_index_shared = 0;
    //file->idx_d->max_file_count = 0;

    file->idx_d->layout_count_shared = 0;
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm, &nprocs);
  }
#endif

  PIDX_return_code ret_code;

  int lvi = start_index;//file->local_variable_index;
  int lower_hz_level = 0, higher_hz_level = 0;

  PIDX_variable var = file->idx->variable[lvi];
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  var->global_block_layout_files[0] = malloc(sizeof (*var->global_block_layout_files[0]));
  memset(var->global_block_layout_files[0], 0, sizeof (*var->global_block_layout_files[0]));
  PIDX_block_layout block_layout = var->global_block_layout_files[0];


  lower_hz_level = hz_level_from;//0;//file->idx_d->reduced_res_from;
  higher_hz_level = hz_level_to;// file->idx_d->maxh;// - file->idx_d->reduced_res_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (lower_hz_level == 0 && higher_hz_level == 0)
    return PIDX_success;


#if 1
  file->idx_d->start_layout_index_shared = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->start_layout_index_shared <= 0)
    file->idx_d->start_layout_index_shared = 0;
#else
  file->idx_d->end_layout_index_shared = 0;
#endif

#if 1
  file->idx_d->end_layout_index_shared = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->end_layout_index_shared <= 0)
    file->idx_d->end_layout_index_shared = 1;
#else
  file->idx_d->end_layout_index_shared = 1;
#endif


  file->idx_d->layout_count_shared = file->idx_d->end_layout_index_shared - file->idx_d->start_layout_index_shared;

  //PIDX_block_layout layout_by_level = file->idx->variable[lvi]->block_layout_by_level_files[0];
  var->block_layout_by_level_files[0] = malloc(sizeof(*(var->block_layout_by_level_files[0])) * file->idx_d->layout_count_shared);
  memset(var->block_layout_by_level_files[0], 0, sizeof(*(var->block_layout_by_level_files[0])) * file->idx_d->layout_count_shared);
  //for (i = 0; i < file->idx_d->layout_count ; i++)
  for (i = file->idx_d->start_layout_index_shared; i < file->idx_d->end_layout_index_shared ; i++)
  {
    var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared] = malloc(sizeof(*(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared])));
    memset(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared], 0, sizeof(*(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared])));
  }

  if (file->idx_d->start_layout_index_shared == 0)
  {
    lower_level_low_layout = 0;
    higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

    if (higher_level_low_layout >= higher_hz_level)
      higher_level_low_layout = higher_hz_level;

    ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[0][0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[0][0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
      memcpy(block_layout->hz_block_number_array[j], var->block_layout_by_level_files[0][0]->hz_block_number_array[j], sizeof(int));

    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
    {
      memcpy(block_layout->hz_block_number_array[j], var->block_layout_by_level_files[0][0]->hz_block_number_array[j], sizeof(int) * ctr);
      ctr = ctr * 2;
    }

    //printf("file->idx_d->layout_count = %d\n", file->idx_d->layout_count);
    for (i = 1; i < file->idx_d->layout_count_shared; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[0][i], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[0][i], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], var->block_layout_by_level_files[0][i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
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


    ctr = (int)pow(2, file->idx_d->start_layout_index_shared - 1) * file->idx->blocks_per_file;
    for (i = file->idx_d->start_layout_index_shared; i < file->idx_d->end_layout_index_shared; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], var->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }


  /*
  if (rank == 0)// && nprocs == 2)
  {
    printf("[B] Final Block Bitmap\n");
    PIDX_blocks_print_layout(block_layout);
  }
  */


  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
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
          block_layout->block_count_per_file[file_number]++;
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
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }


  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      //if (rank == 0)
      //  printf("[%d %d] BPF %d = %d\n", block_layout->resolution_from, block_layout->resolution_to, i, block_layout->block_count_per_file[i]);
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  return PIDX_success;
}


static PIDX_return_code populate_idx_dataset_non_shared(PIDX_global_idx_io file, int start_index, int end_index, int hz_level_from, int hz_level_to)
{
  int i = 0, j = 0, ctr;
  int file_number = 0;

  int rank = 0;
  int nprocs = 1;

  if (hz_level_from == hz_level_to)
  {
    file->idx_d->start_layout_index_non_shared = 0;
    file->idx_d->end_layout_index_non_shared = 0;
    //file->idx_d->max_file_count = 0;

    file->idx_d->layout_count_non_shared = 0;
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm, &nprocs);
  }
#endif

  PIDX_return_code ret_code;

  int lvi = start_index;//file->local_variable_index;
  int lower_hz_level = 0, higher_hz_level = 0;

  PIDX_variable var = file->idx->variable[lvi];
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  var->global_block_layout_files[1] = malloc(sizeof (*var->global_block_layout_files[1]));
  memset(var->global_block_layout_files[1], 0, sizeof (*var->global_block_layout_files[1]));
  PIDX_block_layout block_layout = var->global_block_layout_files[1];

  lower_hz_level = hz_level_from;//0;//file->idx_d->reduced_res_from;
  higher_hz_level = hz_level_to;// file->idx_d->maxh;// - file->idx_d->reduced_res_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (lower_hz_level == 0 && higher_hz_level == 0)
    return PIDX_success;


#if 1
  file->idx_d->start_layout_index_non_shared = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->start_layout_index_non_shared <= 0)
    file->idx_d->start_layout_index_non_shared = 0;
#else
  file->idx_d->end_layout_index_non_shared = 0;
#endif

#if 1
  file->idx_d->end_layout_index_non_shared = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->end_layout_index_non_shared <= 0)
    file->idx_d->end_layout_index_non_shared = 1;
#else
  file->idx_d->end_layout_index_non_shared = 1;
#endif


  file->idx_d->layout_count_non_shared = file->idx_d->end_layout_index_non_shared - file->idx_d->start_layout_index_non_shared;

  var->block_layout_by_level_files[1] = malloc(sizeof(*(var->block_layout_by_level_files[1])) * file->idx_d->layout_count_non_shared);
  memset(var->block_layout_by_level_files[1], 0, sizeof(*(var->block_layout_by_level_files[1])) * file->idx_d->layout_count_non_shared);

  for (i = file->idx_d->start_layout_index_non_shared; i < file->idx_d->end_layout_index_non_shared ; i++)
  {
    var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared] = malloc(sizeof(*(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared])));
    memset(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared], 0, sizeof(*(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared])));
  }

  if (file->idx_d->start_layout_index_non_shared == 0)
  {
    lower_level_low_layout = 0;
    higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

    if (higher_level_low_layout >= higher_hz_level)
      higher_level_low_layout = higher_hz_level;

    ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[1][0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[1][0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
      memcpy(block_layout->hz_block_number_array[j], var->block_layout_by_level_files[1][0]->hz_block_number_array[j], sizeof(int));

    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
    {
      memcpy(block_layout->hz_block_number_array[j], var->block_layout_by_level_files[1][0]->hz_block_number_array[j], sizeof(int) * ctr);
      ctr = ctr * 2;
    }

    //printf("file->idx_d->layout_count = %d\n", file->idx_d->layout_count);
    for (i = 1; i < file->idx_d->layout_count_non_shared; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[1][i], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[1][i], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], var->block_layout_by_level_files[1][i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
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


    ctr = (int)pow(2, file->idx_d->start_layout_index_non_shared - 1) * file->idx->blocks_per_file;
    for (i = file->idx_d->start_layout_index_non_shared; i < file->idx_d->end_layout_index_non_shared; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], var->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }


  /*
  if (rank == 0)// && nprocs == 2)
  {
    printf("[B] Final Block Bitmap\n");
    PIDX_blocks_print_layout(block_layout);
  }
  */


  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
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
          block_layout->block_count_per_file[file_number]++;
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
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }


  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      //if (rank == 0)
      //  printf("[%d %d] BPF %d = %d\n", block_layout->resolution_from, block_layout->resolution_to, i, block_layout->block_count_per_file[i]);
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  return PIDX_success;
}


#if 0
static PIDX_return_code populate_idx_dataset(PIDX_global_idx_io file, int start_index, int end_index, int hz_level_from, int hz_level_to)
{
  int i = 0, j = 0, ctr;
  int file_number = 0;

  int rank = 0;
  int nprocs = 1;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(file->comm, &rank);
    MPI_Comm_size(file->comm, &nprocs);
  }
#endif

  PIDX_return_code ret_code;

  if (hz_level_from == hz_level_to)
  {
    file->idx_d->start_layout_index = 0;
    file->idx_d->end_layout_index = 0;

    file->idx_d->layout_count = 0;
    return PIDX_success;
  }

  int lvi = start_index;//file->local_variable_index;
  int lower_hz_level = 0, higher_hz_level = 0;

  PIDX_variable var = file->idx->variable[lvi];
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  file->idx->variable[lvi]->global_block_layout = malloc(sizeof (*file->idx->variable[lvi]->global_block_layout));
  memset(file->idx->variable[lvi]->global_block_layout, 0, sizeof (*file->idx->variable[lvi]->global_block_layout));
  PIDX_block_layout block_layout = file->idx->variable[lvi]->global_block_layout;

  lower_hz_level = hz_level_from;//0;//file->idx_d->reduced_res_from;
  higher_hz_level = hz_level_to;// file->idx_d->maxh;// - file->idx_d->reduced_res_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx_d->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (lower_hz_level == 0 && higher_hz_level == 0)
    return PIDX_success;


#if 1
  file->idx_d->start_layout_index = (lower_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->start_layout_index <= 0)
    file->idx_d->start_layout_index = 0;
#else
  file->idx_d->end_layout_index = 0;
#endif

#if 1
  file->idx_d->end_layout_index = (higher_hz_level - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->end_layout_index <= 0)
    file->idx_d->end_layout_index = 1;
#else
  file->idx_d->end_layout_index = 1;
#endif


  file->idx_d->layout_count = file->idx_d->end_layout_index - file->idx_d->start_layout_index;

  var->block_layout_by_level = malloc(sizeof(*(var->block_layout_by_level)) * file->idx_d->layout_count);
  memset(var->block_layout_by_level, 0, sizeof(*(var->block_layout_by_level)) * file->idx_d->layout_count);
  //for (i = 0; i < file->idx_d->layout_count ; i++)
  for (i = file->idx_d->start_layout_index; i < file->idx_d->end_layout_index ; i++)
  {
    var->block_layout_by_level[i - file->idx_d->start_layout_index] = malloc(sizeof(*(var->block_layout_by_level[i - file->idx_d->start_layout_index])));
    memset(var->block_layout_by_level[i - file->idx_d->start_layout_index], 0, sizeof(*(var->block_layout_by_level[i - file->idx_d->start_layout_index])));
  }

  if (file->idx_d->start_layout_index == 0)
  {
    lower_level_low_layout = 0;
    higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

    if (higher_level_low_layout >= higher_hz_level)
      higher_level_low_layout = higher_hz_level;

    ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout, file->idx_d->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret_code = populate_idx_layout(file, start_index, end_index, file->idx->variable[lvi]->block_layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
      memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int));

    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
    {
      memcpy(block_layout->hz_block_number_array[j], file->idx->variable[lvi]->block_layout_by_level[0]->hz_block_number_array[j], sizeof(int) * ctr);
      ctr = ctr * 2;
    }

    //printf("file->idx_d->layout_count = %d\n", file->idx_d->layout_count);
    for (i = 1; i < file->idx_d->layout_count; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, file->idx->variable[lvi]->block_layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], file->idx->variable[lvi]->block_layout_by_level[i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
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


    ctr = (int)pow(2, file->idx_d->start_layout_index - 1) * file->idx->blocks_per_file;
    for (i = file->idx_d->start_layout_index; i < file->idx_d->end_layout_index; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(file->idx->variable[lvi]->block_layout_by_level[i - file->idx_d->start_layout_index], lower_level_higher_layout, higher_level_higher_layout, file->idx_d->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      ret_code = populate_idx_layout(file, start_index, end_index, file->idx->variable[lvi]->block_layout_by_level[i - file->idx_d->start_layout_index], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], file->idx->variable[lvi]->block_layout_by_level[i - file->idx_d->start_layout_index]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }


  /*
  if (rank == 0)// && nprocs == 2)
  {
    printf("[B] Final Block Bitmap\n");
    PIDX_blocks_print_layout(block_layout);
  }
  */


  block_layout->file_bitmap = malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx_d->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx_d->max_file_count));

  block_layout->block_count_per_file = malloc(sizeof(int) * (file->idx_d->max_file_count));
  memset(block_layout->block_count_per_file, 0, sizeof(int) * (file->idx_d->max_file_count));

  if (block_layout->resolution_from <= block_layout->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->block_count_per_file[file_number]++;
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
          block_layout->block_count_per_file[file_number]++;
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
            block_layout->file_index[file_number] = 1;
            block_layout->block_count_per_file[file_number]++;
          }
        }
      }
      ctr = ctr * 2;
    }
  }


  block_layout->existing_file_count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->existing_file_count++;

  block_layout->existing_file_index = (int*) malloc(block_layout->existing_file_count * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->existing_file_count * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx_d->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx_d->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx_d->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      //if (rank == 0)
      //  printf("[%d %d] BPF %d = %d\n", block_layout->resolution_from, block_layout->resolution_to, i, block_layout->block_count_per_file[i]);
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  return PIDX_success;
}
#endif



PIDX_global_idx_io PIDX_global_idx_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_global_idx_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_d = idx_derived_ptr;
  idx_io_id->idx_dbg = idx_dbg;

  return (idx_io_id);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_global_idx_io_set_communicator(PIDX_global_idx_io id, MPI_Comm comm)
{
  if (id == NULL)
    return PIDX_err_id;

  id->global_comm = comm;

  return PIDX_success;
}
#endif


static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}


static PIDX_return_code partition_destroy(PIDX_global_idx_io file)
{
  int ret = 0;

  /* Destroy buffers allocated during restructuring phase */
  ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  ret = PIDX_rst_meta_data_destroy(file->rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  /* Deleting the restructuring ID */
  PIDX_rst_finalize(file->rst_id);


  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;

  return PIDX_success;
}


static PIDX_return_code partition_setup(PIDX_global_idx_io file, int start_var_index, int end_var_index, int pipe_len)
{
  int ret = 0, nprocs = 1;
  int start_index = 0, end_index = 0;

  //PIDX_time time = file->idx_d->time;

#if PIDX_HAVE_MPI
  file->comm = file->global_comm;

  if (file->idx_d->parallel_mode == 1)
    MPI_Comm_size(file->comm,  &nprocs);
#endif

  file->idx_d->rank_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_offset, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

  file->idx_d->rank_r_count =  malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS);
  memset(file->idx_d->rank_r_count, 0, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS));

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);

    MPI_Allgather(file->idx->variable[start_var_index]->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->idx_d->rank_r_count, PIDX_MAX_DIMENSIONS, MPI_LONG_LONG, file->comm);
  }
  else
  {
    memcpy(file->idx_d->rank_r_offset, file->idx->variable[start_var_index]->sim_patch[0]->offset, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    memcpy(file->idx_d->rank_r_count, file->idx->variable[start_var_index]->sim_patch[0]->size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  }
#endif

  set_default_patch_size(file, file->idx_d->rank_r_count, nprocs);

  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (pipe_len + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

#if PIDX_HAVE_MPI
    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_rst;
#endif

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Creating the buffers required for restructurig */
    ret = PIDX_rst_buf_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Perform data restructuring */
    if (file->idx_dbg->debug_do_rst == 1)
    {
      ret = PIDX_rst_write(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_WRITE);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      ret = PIDX_rst_buf_destroy(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;
    }
  }

  return PIDX_success;
}


static PIDX_return_code partition(PIDX_global_idx_io file, int start_var_index, int end_var_index)
{
  int *colors;
  int index_i = 0, index_j = 0, index_k = 0;
  int i = 0, j = 0, k = 0, d = 0;
  //PIDX_return_code ret;

  int rank = 0;
  MPI_Comm_rank(file->global_comm, &rank);

  colors = malloc(sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  memset(colors, 0, sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  file->idx_d->color = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]) + 1;//MPI_UNDEFINED;

  for (k = 0; k < file->idx_d->idx_count[2]; k++)
    for (j = 0; j < file->idx_d->idx_count[1]; j++)
      for (i = 0; i < file->idx_d->idx_count[0]; i++)
      {
        colors[(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i] = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i;
      }

  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));

  PIDX_variable var = file->idx->variable[start_var_index];
  if (var->patch_group_count == 1)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = file->idx->variable[start_var_index]->rst_patch_group[0]->reg_patch->offset[d];
      local_proc_patch->size[d] = file->idx->variable[start_var_index]->rst_patch_group[0]->reg_patch->size[d];
    }
    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    PointND bounds_point;
    int maxH = 0;
    bounds_point.x = (int) file->idx_d->idx_count[0];
    bounds_point.y = (int) file->idx_d->idx_count[1];
    bounds_point.z = (int) file->idx_d->idx_count[2];
    bounds_point.u = (int) 1;
    bounds_point.v = (int) 1;
    char bitSequence[512];
    char bitPattern[512];
    GuessBitmaskPattern(bitSequence, bounds_point);
    maxH = strlen(bitSequence);

    //printf("IDX C %d %d %d mh %d\n", file->idx_d->idx_count[0], file->idx_d->idx_count[1], file->idx_d->idx_count[2], maxH);
    for (i = 0; i <= maxH; i++)
      bitPattern[i] = RegExBitmaskBit(bitSequence, i);

    //printf("maxH = %d\n", maxH);
    int z_order = 0;
    int number_levels = maxH - 1;

    if ((file->idx_d->idx_size[0]) > file->idx->bounds[0])
      file->idx_d->idx_size[0] = file->idx->bounds[0];
    if ((file->idx_d->idx_size[1]) > file->idx->bounds[1])
      file->idx_d->idx_size[1] = file->idx->bounds[1];
    if ((file->idx_d->idx_size[2]) > file->idx->bounds[2])
      file->idx_d->idx_size[2] = file->idx->bounds[2];

    for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx_d->idx_size[0], index_i++)
    {
      for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx_d->idx_size[1], index_j++)
      {
        for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx_d->idx_size[2], index_k++)
        {
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = file->idx_d->idx_size[0];
          reg_patch->size[1] = file->idx_d->idx_size[1];
          reg_patch->size[2] = file->idx_d->idx_size[2];

         if (intersectNDChunk(reg_patch, local_proc_patch))
          {
            PointND xyzuv_Index;
            xyzuv_Index.x = index_i;
            xyzuv_Index.y = index_j;
            xyzuv_Index.z = index_k;
            xyzuv_Index.u = 0;
            xyzuv_Index.v = 0;

            z_order = 0;
            PointND zero;
            zero.x = 0;
            zero.y = 0;
            zero.z = 0;
            zero.u = 0;
            zero.v = 0;
            memset(&zero, 0, sizeof (PointND));

            int cnt = 0;
            for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (PointND)); cnt++, number_levels--)
            {
              int bit = bitPattern[number_levels];
              z_order |= ((int64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            file->idx_d->color = colors[z_order];
            //printf("[%d] ---> %d\n", rank, file->idx_d->color);

            /*
            distance_x = index_i * file->idx_d->idx_size[0];
            distance_y = index_j * file->idx_d->idx_size[1];
            distance_z = index_k * file->idx_d->idx_size[2];

            var->rst_patch_group[0]->reg_patch->offset[0] = var->rst_patch_group[0]->reg_patch->offset[0] - distance_x;
            var->rst_patch_group[0]->reg_patch->offset[1] = var->rst_patch_group[0]->reg_patch->offset[1] - distance_y;
            var->rst_patch_group[0]->reg_patch->offset[2] = var->rst_patch_group[0]->reg_patch->offset[2] - distance_z;

            var->sim_patch_count = 1;
            var->sim_patch[0]->offset[0] = var->rst_patch_group[0]->reg_patch->offset[0];
            var->sim_patch[0]->offset[1] = var->rst_patch_group[0]->reg_patch->offset[1];
            var->sim_patch[0]->offset[2] = var->rst_patch_group[0]->reg_patch->offset[2];
            var->sim_patch[0]->offset[3] = 0;
            var->sim_patch[0]->offset[4] = 0;

            var->sim_patch[0]->size[0] = var->rst_patch_group[0]->reg_patch->size[0];
            var->sim_patch[0]->size[1] = var->rst_patch_group[0]->reg_patch->size[1];
            var->sim_patch[0]->size[2] = var->rst_patch_group[0]->reg_patch->size[2];
            var->sim_patch[0]->size[3] = 1;
            var->sim_patch[0]->size[4] = 1;

            file->idx->bounds[0] = reg_patch->size[0];
            file->idx->bounds[1] = reg_patch->size[1];
            file->idx->bounds[2] = reg_patch->size[2];
            */

            //printf("[UO %d] : %d %d %d\n", rank, var->rst_patch_group[0]->reg_patch->offset[0], var->rst_patch_group[0]->reg_patch->offset[1], var->rst_patch_group[0]->reg_patch->offset[2]);
            //printf("[B %d] : %d %d %d\n", rank, file->idx->bounds[0], file->idx->bounds[1], file->idx->bounds[2]);

            break;
          }
        }
      }
    }
    free(reg_patch);
  }
  else
  {
    file->idx->bounds[0] = 0;//reg_patch->size[0];
    file->idx->bounds[1] = 0;//reg_patch->size[1];
    file->idx->bounds[2] = 0;//reg_patch->size[2];

    var->sim_patch_count = 0;
    var->sim_patch[0]->offset[0] = 0;
    var->sim_patch[0]->offset[1] = 0;
    var->sim_patch[0]->offset[2] = 0;
    var->sim_patch[0]->offset[3] = 0;
    var->sim_patch[0]->offset[4] = 0;

    var->sim_patch[0]->size[0] = 0;//-1;
    var->sim_patch[0]->size[1] = 0;//-1;
    var->sim_patch[0]->size[2] = 0;//-1;
    var->sim_patch[0]->size[3] = 0;//-1;
    var->sim_patch[0]->size[4] = 0;//-1;
  }

  MPI_Comm_split(file->global_comm, file->idx_d->color, rank, &(file->comm));
  free(colors);

  //
  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';

  if (file->idx_d->idx_count[0] != 1 || file->idx_d->idx_count[1] != 1 || file->idx_d->idx_count[2] != 1)
    sprintf(file->idx->filename, "%s_%d.idx", file_name_skeleton, file->idx_d->color);
  //

  //printf("File name = %s\n", file->idx->filename);
  free(local_proc_patch);
  local_proc_patch = 0;

  return PIDX_success;
}


static void set_default_patch_size(PIDX_global_idx_io file, int64_t* process_bounds, int nprocs)
{
  int i = 0, j = 0;
  int64_t average_count = 0;
  int check_bit = 0;
  int64_t max_dim_length[PIDX_MAX_DIMENSIONS] = {0, 0, 0, 0, 0};
  int equal_partiton = 0;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * 0 + i];
    for (j = 0; j < nprocs; j++)
    {
      if (max_dim_length[i] <= process_bounds[PIDX_MAX_DIMENSIONS * j + i])
        max_dim_length[i] = process_bounds[PIDX_MAX_DIMENSIONS * j + i];
    }
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    average_count = average_count + max_dim_length[i];
  }
  average_count = average_count / PIDX_MAX_DIMENSIONS;
  average_count = getPowerOftwo(average_count);

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    check_bit = check_bit || ((double) file->idx->bounds[i] / average_count > (double) file->idx->bounds[i] / max_dim_length[i]);

  while (check_bit)
  {
    average_count = average_count * 2;
    check_bit = 0;
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      check_bit = check_bit || ((double) file->idx->bounds[i] / average_count > (double) file->idx->bounds[i] / max_dim_length[i]);
  }
  int64_t reg_patch_size[PIDX_MAX_DIMENSIONS];
  //reg_patch_size =  average_count;
  if (equal_partiton == 1)
  {
    reg_patch_size[0] = average_count / 1;
    reg_patch_size[1] = average_count / 1;
    reg_patch_size[2] = average_count / 1;
    reg_patch_size[3] = 1;
    reg_patch_size[4] = 1;
  }
  else
  {
    reg_patch_size[0] = getPowerOftwo(max_dim_length[0]) * 1;
    reg_patch_size[1] = getPowerOftwo(max_dim_length[1]) * 1;
    reg_patch_size[2] = getPowerOftwo(max_dim_length[2]) * 1;
    //printf("Box size = %d %d %d\n", reg_patch_size[0], reg_patch_size[1], reg_patch_size[2]);
    //reg_patch_size[0] = getPowerOftwo(process_bounds[0]) * 1;
    //reg_patch_size[1] = getPowerOftwo(process_bounds[1]) * 1;
    //reg_patch_size[2] = getPowerOftwo(process_bounds[2]) * 1;
    reg_patch_size[3] = 1;//getPowerOftwo(process_bounds[3]) * 1;
    reg_patch_size[4] = 1;//getPowerOftwo(process_bounds[4]) * 1;
  }

  memcpy(file->idx->reg_patch_size, reg_patch_size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
  //reg_patch_size = reg_patch_size * 4;
}



static PIDX_return_code write_headers_shared(PIDX_global_idx_io file, int start_var_index, int end_var_index, int layout_type)
{
  int ret = 0;
#if !SIMULATE_IO
  if (file->idx_dbg->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif
    ret = PIDX_header_io_file_create(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[0]);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (file->idx->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[0],  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (file->idx->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[0], 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code write_headers_non_shared(PIDX_global_idx_io file, int start_var_index, int end_var_index, int layout_type)
{
  int ret = 0;
#if !SIMULATE_IO
  if (file->idx_dbg->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif
    ret = PIDX_header_io_file_create(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[1]);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (file->idx->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[1],  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (file->idx->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout_files[1], 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
#endif

  return PIDX_success;
}

#if 0
static PIDX_return_code write_headers(PIDX_global_idx_io file, int start_var_index, int end_var_index, int layout_type)
{
  int ret = 0;
#if !SIMULATE_IO
  if (file->idx_dbg->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, start_var_index, end_var_index);
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      ret = PIDX_header_io_set_communicator(file->header_io_id, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
#endif
    ret = PIDX_header_io_file_create(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (file->idx->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout,  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (file->idx->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code PIDX_global_write_io(PIDX_global_idx_io file, int start_var_index, int end_var_index)
{
  int i, j;
  int rank = 0, nprocs = 1;
  int ret = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;
  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    time->startup_start[start_index] = PIDX_get_time();
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    int agg_io_level = 0, no_of_aggregators = 0;
    if (file->idx_d->agg_type != 0)
    {
      for (i = file->idx_d->start_layout_index; i < file->idx_d->end_layout_index ; i++)
      {
        no_of_aggregators = file->idx->variable[start_var_index]->block_layout_by_level[i - file->idx_d->start_layout_index]->existing_file_count;
        if (no_of_aggregators <= nprocs)
          agg_io_level = i;
      }
      agg_io_level = agg_io_level + 1;
    }
    else
    {
      no_of_aggregators = file->idx->variable[start_var_index]->global_block_layout->existing_file_count;
      if (no_of_aggregators <= nprocs)
       agg_io_level = file->idx_d->end_layout_index;
     else
        agg_io_level = file->idx_d->start_layout_index;//0;
    }

    if (file->idx->enable_agg == 0)
      agg_io_level = file->idx_d->start_layout_index;//0;

    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
    /* Create the aggregation ID */
    /* Create the I/O ID */
    file->tagg_id = malloc(sizeof(*(file->tagg_id)) * file->idx->variable_count);
    memset(file->tagg_id, 0, sizeof(*(file->tagg_id)) * file->idx->variable_count);

    file->tio_id = malloc(sizeof(*(file->tio_id)) * file->idx->variable_count);
    memset(file->tio_id, 0, sizeof(*(file->tio_id)) * file->idx->variable_count);

    int agg_var_pipe = 0;
    int agg_end_index = 0;
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      file->tagg_id[i] = malloc(sizeof(*(file->tagg_id[i])) * file->idx_d->layout_count);
      memset(file->tagg_id[i], 0, sizeof(*(file->tagg_id[i])) * file->idx_d->layout_count);

      file->tio_id[i] = malloc(sizeof(*(file->tio_id[i])) * file->idx_d->layout_count);
      memset(file->tio_id[i], 0, sizeof(*(file->tio_id[i])) * file->idx_d->layout_count);
    }

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      agg_end_index = ((i + agg_var_pipe) >= (end_index + 1)) ? (end_index) : (i + agg_var_pipe);
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        file->tagg_id[i][j - file->idx_d->start_layout_index] = PIDX_agg_init(file->idx, file->idx_d, start_var_index, i, agg_end_index);

      for(j = file->idx_d->start_layout_index ; j < file->idx_d->end_layout_index; j++)
        file->tio_id[i][j - file->idx_d->start_layout_index] = PIDX_file_io_init(file->idx, file->idx_d, start_var_index, i, agg_end_index);
    }

    file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
    memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);

    int j_1 = 0;
    for (i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      file->idx_d->agg_buffer[i] = malloc(sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      memset(file->idx_d->agg_buffer[i], 0, sizeof(*(file->idx_d->agg_buffer[i])) * agg_io_level);
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      {
        j_1 = j - file->idx_d->start_layout_index;
        file->idx_d->agg_buffer[i][j_1] = malloc(sizeof(*(file->idx_d->agg_buffer[i][j_1])) );
        memset(file->idx_d->agg_buffer[i][j_1], 0, sizeof(*(file->idx_d->agg_buffer[i][j_1])) );

        file->idx_d->agg_buffer[i][j_1]->file_number = -1;
        file->idx_d->agg_buffer[i][j_1]->var_number = -1;
        file->idx_d->agg_buffer[i][j_1]->sample_number = -1;

        file->idx_d->agg_buffer[i][j_1]->no_of_aggregators = 0;
        file->idx_d->agg_buffer[i][j_1]->aggregator_interval = 0;
        file->idx_d->agg_buffer[i][j_1]->aggregation_factor = 1;
      }
      //file->idx_d->aggregator_multiplier = 1;
      for (j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        file->idx_d->agg_buffer[i][j_1]->aggregation_factor = file->idx_d->aggregator_multiplier;//(int)pow(2, (agg_io_level - j));
    }
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      /* Attaching the communicator to the aggregation phase */
      /* Attaching the communicator to the I/O phase */
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        {
          j_1 = j - file->idx_d->start_layout_index;
          ret = PIDX_agg_set_global_communicator(file->tagg_id[i][j_1], file->global_comm);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_agg;
          }

          ret = PIDX_agg_set_communicator(file->tagg_id[i][j_1], file->comm);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_agg;
          }
        }

        for(j = file->idx_d->start_layout_index ; j < file->idx_d->end_layout_index; j++)
        {
          j_1 = j - file->idx_d->start_layout_index;
          ret = PIDX_file_io_set_communicator(file->tio_id[i][j_1], file->comm);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
        }
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/

    time->startup_end[start_index] = PIDX_get_time();

#if PIDX_DEBUG_OUTPUT
    l_init = 1;
    MPI_Allreduce(&l_init, &g_init, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_init == nprocs)
      printf("All Modules Initialized (Variable index [%d %d] Variable Pipe length %d)\n", start_index, end_index, file->idx_d->var_pipe_length);
#endif

    /*----------------------------------------------Agg [staart]-----------------------------------------------*/

    /* Creating the buffers required for Aggregation */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      {
        j_1 = j - file->idx_d->start_layout_index;
        PIDX_agg_id temp_id = file->tagg_id[i][j_1];
        Agg_buffer temp_agg = file->idx_d->agg_buffer[i][j_1];
        PIDX_block_layout temp_layout = file->idx->variable[start_var_index]->block_layout_by_level[j_1];
        PIDX_block_layout temp_global_layout = file->idx->variable[start_var_index]->global_block_layout;

        time->agg_meta_start[i][j] = PIDX_get_time();
        /* Creating the buffers required for Aggregation */
        ret = PIDX_agg_meta_data_create(temp_id, temp_agg, temp_global_layout, temp_layout);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_meta_end[i][j] = PIDX_get_time();

        time->agg_buf_start[i][j] = PIDX_get_time();
        ret = PIDX_agg_buf_create(temp_id, temp_agg, temp_layout, temp_global_layout, i, j);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_buf_end[i][j] = PIDX_get_time();

#if PIDX_DEBUG_OUTPUT
    l_agg_buf = 1;
    MPI_Allreduce(&l_agg_buf, &g_agg_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg_buf == nprocs)
      printf("[A] Aggregation Buffer Created\n");
#endif

        if (file->idx_dbg->debug_do_agg == 1)
        {
          time->agg_start[i][j] = PIDX_get_time();
          ret = PIDX_agg(temp_id, temp_agg, j, temp_layout, PIDX_WRITE, i, j_1);

          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_rst;
          }
          time->agg_end[i][j] = PIDX_get_time();
        }

#if PIDX_DEBUG_OUTPUT
    l_agg = 1;
    MPI_Allreduce(&l_agg, &g_agg, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg == nprocs)
      printf("[A] Aggregation Completed\n");
#endif

        if (file->idx_dbg->debug_do_io == 1)
        {
          time->io_start[i][j] = PIDX_get_time();
          ret = PIDX_aggregated_io(file->tio_id[i][j_1], temp_agg, temp_layout, PIDX_WRITE);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
          time->io_end[i][j] = PIDX_get_time();
        }

#if PIDX_DEBUG_OUTPUT
    l_io = 1;
    MPI_Allreduce(&l_io, &g_io, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_io == nprocs)
      printf("[I] Agg I/O completed\n");
#endif

        ret = PIDX_agg_buf_destroy(temp_id, temp_agg);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_agg;
        }

        free(temp_agg);
        temp_agg = 0;

        ret = PIDX_agg_meta_data_destroy(temp_id, temp_layout, temp_global_layout);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }

        PIDX_agg_finalize(temp_id);
      }

      free(file->idx_d->agg_buffer[i]);
      file->idx_d->agg_buffer[i] = 0;
    }
    free(file->idx_d->agg_buffer);\
    file->idx_d->agg_buffer = 0;

    if (file->idx_dbg->debug_do_io == 1)
    {
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = agg_io_level ; j < file->idx_d->end_layout_index; j++)
        {
          time->io_per_process_start[i][j] = PIDX_get_time();
          ret = PIDX_file_io_per_process(file->tio_id[i][j - agg_io_level], file->idx->variable[start_var_index]->block_layout_by_level[j - agg_io_level], PIDX_WRITE);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
          time->io_per_process_end[i][j] = PIDX_get_time();
        }
      }
    }

    /*--------------------------------------------Agg [end]--------------------------------------------------*/


    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[start_index] = PIDX_get_time();

    /* Deleting the I/O ID */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      for(j = file->idx_d->start_layout_index ; j < file->idx_d->end_layout_index; j++)
        PIDX_file_io_finalize(file->tio_id[i][j - file->idx_d->start_layout_index]);

    for(i = 0 ; i < file->idx->variable_count ; i++)
    {
      free(file->tagg_id[i]);
      file->tagg_id[i] = 0;

      free(file->tio_id[i]);
      file->tio_id[i] = 0;
    }
    free(file->tagg_id);
    file->tagg_id = 0;

    free(file->tio_id);
    file->tio_id = 0;

    time->finalize_end[start_index] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/
  }

  return PIDX_success;
}
#endif



static PIDX_return_code PIDX_global_io(PIDX_global_idx_io file, int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level)
{
  int j;
  int rank = 0, nprocs = 1;
  int ret = 0;
  int j_1 = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;

  /*------------------------------------Create ALL the IDs [start]---------------------------------------*/

  for(j = layout_start ; j < layout_end; j++)
  {
    j_1 = j - layout_start;
    file->tio_id[var_index][j] = PIDX_file_io_init(file->idx, file->idx_d, init_index, var_index, var_index);

    ret = PIDX_file_io_set_communicator(file->tio_id[var_index][j], file->comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  for(j = layout_start ; j < agg_io_level; j++)
  {
    j_1 = j - layout_start;
    Agg_buffer temp_agg = file->idx_d->agg_buffer[var_index][j];
    PIDX_block_layout temp_layout = file->idx->variable[init_index]->block_layout_by_level_files[index][j_1];

    if (file->idx_dbg->debug_do_io == 1)
    {
      time->io_start[var_index][j] = PIDX_get_time();
      ret = PIDX_aggregated_io(file->tio_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      time->io_end[var_index][j] = PIDX_get_time();
    }

    ret = PIDX_agg_buf_destroy(file->tagg_id[var_index][j], temp_agg);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(temp_agg);
    temp_agg = 0;
  }
  //free(file->idx_d->agg_buffer[var_index]);
  //file->idx_d->agg_buffer[var_index] = 0;

  if (file->idx_dbg->debug_do_io == 1)
  {
    for(j = agg_io_level ; j < layout_end; j++)
    {
      time->io_per_process_start[var_index][j] = PIDX_get_time();
      ret = PIDX_file_io_per_process(file->tio_id[var_index][j - agg_io_level], file->idx->variable[init_index]->block_layout_by_level_files[index][j - agg_io_level], PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      time->io_per_process_end[var_index][j] = PIDX_get_time();
    }
  }

  for(j = layout_start ; j < layout_end; j++)
    PIDX_file_io_finalize(file->tio_id[var_index][j]);

  for(j = layout_start ; j < agg_io_level; j++)
    PIDX_agg_finalize(file->tagg_id[var_index][j]);

  /*
  free(file->tio_id[var_index]);
  file->tio_id[var_index] = 0;

  free(file->tagg_id[var_index]);
  file->tagg_id[var_index] = 0;
  */

  return PIDX_success;
}




static PIDX_return_code PIDX_global_async_io(PIDX_global_idx_io file, int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level)
{
  int j;
  int rank = 0, nprocs = 1;
  int ret = 0;
  int j_1 = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;

  /*------------------------------------Create ALL the IDs [start]---------------------------------------*/

  for(j = layout_start ; j < layout_end; j++)
  {
    j_1 = j - layout_start;
    file->tio_id[var_index][j] = PIDX_file_io_init(file->idx, file->idx_d, init_index, var_index, var_index);

    ret = PIDX_file_io_set_communicator(file->tio_id[var_index][j], file->comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
  }

  //MPI_Request *req = *request;
  for(j = layout_start ; j < agg_io_level; j++)
  {
    j_1 = j - layout_start;
    Agg_buffer temp_agg = file->idx_d->agg_buffer[var_index][j];
    PIDX_block_layout temp_layout = file->idx->variable[init_index]->block_layout_by_level_files[index][j_1];

    if (file->idx_dbg->debug_do_io == 1)
    {
      time->io_start[var_index][j] = PIDX_get_time();
      ret = PIDX_async_aggregated_io(file->tio_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE, /*&(req[j_1])*/&(file->idx_d->request[j_1]), &(file->idx_d->fp[j_1]));
      //ret = PIDX_aggregated_io(file->tio_id[var_index][j], temp_agg, temp_layout, PIDX_WRITE);
      //if (ret != PIDX_success)
      //{
      //  fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      //  return PIDX_err_io;
      //}
      time->io_end[var_index][j] = PIDX_get_time();
    }
    //free(temp_agg);
    //temp_agg = 0;
  }
  //free(file->idx_d->agg_buffer[var_index]);
  //file->idx_d->agg_buffer[var_index] = 0;

  if (file->idx_dbg->debug_do_io == 1)
  {
    for(j = agg_io_level ; j < layout_end; j++)
    {
      time->io_per_process_start[var_index][j] = PIDX_get_time();
      ret = PIDX_file_io_per_process(file->tio_id[var_index][j - agg_io_level], file->idx->variable[init_index]->block_layout_by_level_files[index][j - agg_io_level], PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      time->io_per_process_end[var_index][j] = PIDX_get_time();
    }
  }


  /*
  for(j = layout_start ; j < layout_end; j++)
    PIDX_file_io_finalize(file->tio_id[var_index][j]);

  for(j = layout_start ; j < agg_io_level; j++)
    PIDX_agg_finalize(file->tagg_id[var_index][j]);

  free(file->tio_id[var_index]);
  file->tio_id[var_index] = 0;

  free(file->tagg_id[var_index]);
  file->tagg_id[var_index] = 0;
  */

  return PIDX_success;
}



static PIDX_return_code PIDX_global_aggregate(PIDX_global_idx_io file, int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level)
{
  int j;
  int rank = 0, nprocs = 1;
  int ret = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  PIDX_time time = file->idx_d->time;

  /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
  /* Create the aggregation ID */
  /* Create the I/O ID */
  int j_1 = 0;

  for(j = layout_start ; j < agg_io_level; j++)
  {
    j_1 = j - layout_start;

    file->tagg_id[var_index][j] = PIDX_agg_init(file->idx, file->idx_d, init_index, var_index, var_index);

    file->idx_d->agg_buffer[var_index][j] = malloc(sizeof(*(file->idx_d->agg_buffer[var_index][j])) );
    memset(file->idx_d->agg_buffer[var_index][j], 0, sizeof(*(file->idx_d->agg_buffer[var_index][j])) );

    file->idx_d->agg_buffer[var_index][j]->file_number = -1;
    file->idx_d->agg_buffer[var_index][j]->var_number = -1;
    file->idx_d->agg_buffer[var_index][j]->sample_number = -1;

    file->idx_d->agg_buffer[var_index][j]->no_of_aggregators = 0;
    file->idx_d->agg_buffer[var_index][j]->aggregator_interval = 0;
    file->idx_d->agg_buffer[var_index][j]->aggregation_factor = 1;
    file->idx_d->agg_buffer[var_index][j]->aggregation_factor = file->idx_d->aggregator_multiplier;//(int)pow(2, (agg_io_level - j));

    ret = PIDX_agg_set_global_communicator(file->tagg_id[var_index][j], file->global_comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    ret = PIDX_agg_set_communicator(file->tagg_id[var_index][j], file->comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    PIDX_agg_id temp_id = file->tagg_id[var_index][j];
    Agg_buffer temp_agg = file->idx_d->agg_buffer[var_index][j];
    PIDX_block_layout temp_layout = file->idx->variable[init_index]->block_layout_by_level_files[index][j_1];
    PIDX_block_layout temp_global_layout = file->idx->variable[init_index]->global_block_layout_files[index];

    time->agg_meta_start[var_index][j] = PIDX_get_time();
    /* Creating the buffers required for Aggregation */
    ret = PIDX_agg_meta_data_create(temp_id, temp_agg, temp_global_layout, temp_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    time->agg_meta_end[var_index][j] = PIDX_get_time();

    time->agg_buf_start[var_index][j] = PIDX_get_time();
    ret = PIDX_agg_buf_create(temp_id, temp_agg, temp_layout, temp_global_layout, var_index, j);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    time->agg_buf_end[var_index][j] = PIDX_get_time();

    if (file->idx_dbg->debug_do_agg == 1)
    {
      time->agg_start[var_index][j] = PIDX_get_time();
      ret = PIDX_agg(temp_id, temp_agg, j, temp_layout, PIDX_WRITE, var_index, j_1);

      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_end[var_index][j] = PIDX_get_time();
    }

    ret = PIDX_agg_meta_data_destroy(temp_id, temp_layout, temp_global_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
  }

  return PIDX_success;
}



static PIDX_return_code one_time_initialize(PIDX_global_idx_io file)
{
  int ret = 0;
  int total_header_size;
  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, file->idx_d->perm_layout_count);

    ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;

    file->one_time_initializations = 1;
  }

  return PIDX_success;
}


PIDX_return_code destroy_hz_buffers(PIDX_global_idx_io file, int start_var_index, int end_var_index)
{
  PIDX_return_code ret;

  PIDX_time time = file->idx_d->time;
  int start_index = 0;//, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    time->startup_start[start_index] = PIDX_get_time();
    //end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[start_index] = PIDX_get_time();

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    /* Deleting the compression ID */
    PIDX_compression_finalize(file->comp_id);

    /* Deleting the chunking ID */
    PIDX_chunk_finalize(file->chunk_id);

    time->finalize_end[start_index] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);
  }

  return PIDX_success;
}



PIDX_return_code create_hz_buffers(PIDX_global_idx_io file, int start_var_index, int end_var_index)
{
  int rank = 0, nprocs = 1;
  int ret = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  //PIDX_time time = file->idx_d->time;
  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    //time->startup_start[start_index] = PIDX_get_time();
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);


    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
    /* Create the chunking ID */
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the compression ID */
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the HZ encoding ID */
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      /* Attaching the communicator to the chunking phase */
      ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }

      /* Attaching the communicator to the compression phase */
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }

      /* Attaching the communicator to the HZ encodig phase phase */
      PIDX_hz_encode_set_global_communicator(file->hz_id, file->comm);
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/

    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_hz_encode_set_resolution(file->hz_id, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    //time->startup_end[start_index] = PIDX_get_time();

#if PIDX_DEBUG_OUTPUT
    l_init = 1;
    MPI_Allreduce(&l_init, &g_init, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_init == nprocs)
      printf("All Modules Initialized (Variable index [%d %d] Variable Pipe length %d)\n", start_index, end_index, file->idx_d->var_pipe_length);
#endif


    /*----------------------------------------Chunking [start]------------------------------------------------*/
    //time->chunk_start[start_index] = PIDX_get_time();

    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

#if PIDX_DEBUG_OUTPUT
    l_chunk_buf = 1;
    MPI_Allreduce(&l_chunk_buf, &g_chunk_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk_buf == nprocs)
      printf("[C] Chunking Buffer Created\n");
#endif

    /* Perform Chunking */
    if (file->idx_dbg->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_chunk = 1;
    MPI_Allreduce(&l_chunk, &g_chunk, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_chunk == nprocs)
      printf("[C] Chunking Completed\n");
#endif

    //time->chunk_end[start_index] = PIDX_get_time();
    /*-----------------------------------------Chunking [end]------------------------------------------------*/


    /*----------------------------------------Compression [start]--------------------------------------------*/
    //time->compression_start[start_index] = PIDX_get_time();

    /* Perform Compression */
    if (file->idx_dbg->debug_do_compress == 1)
    {
#if !SIMULATE_IO
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }
#endif
    }

#if PIDX_DEBUG_OUTPUT
    l_cmp = 1;
    MPI_Allreduce(&l_cmp, &g_cmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_cmp == nprocs)
      printf("[CMP] Compression Completed\n");
#endif

    //time->compression_end[start_index] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/



    /*---------------------------------------------HZ [start]------------------------------------------------*/
    //time->hz_start[start_index] = PIDX_get_time();

    /* Creating the buffers required for HZ encoding */
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

#if PIDX_DEBUG_OUTPUT
    l_hz_buf = 1;
    MPI_Allreduce(&l_hz_buf, &g_hz_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz_buf == nprocs)
      printf("[H] HZ Buffer Created\n");
#endif

    /* Perform HZ encoding */
    if (file->idx_dbg->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_hz = 1;
    MPI_Allreduce(&l_hz, &g_hz, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_hz == nprocs)
      printf("[H] HZ Encoding Finished\n");
#endif

    /* Verify the HZ encoding */
    if(file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    /* Destroy buffers allocated during chunking phase */
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

    //time->hz_end[start_index] = PIDX_get_time();
    /*----------------------------------------------HZ [end]-------------------------------------------------*/
  }

  return PIDX_success;
}

PIDX_return_code init_agg_io(PIDX_global_idx_io file, int start_index)
{
  PIDX_variable var = file->idx->variable[start_index];
  var->global_block_layout_files = malloc(2 * sizeof (*var->global_block_layout_files));
  memset(var->global_block_layout_files, 0, 2 * sizeof (*var->global_block_layout_files));

  var->block_layout_by_level_files = malloc(2 * sizeof (*var->block_layout_by_level_files));
  memset(var->block_layout_by_level_files, 0, 2 * sizeof (*var->block_layout_by_level_files));

  file->tagg_id = malloc(sizeof(*(file->tagg_id)) * file->idx->variable_count);
  memset(file->tagg_id, 0, sizeof(*(file->tagg_id)) * file->idx->variable_count);

  file->tio_id = malloc(sizeof(*(file->tio_id)) * file->idx->variable_count);
  memset(file->tio_id, 0, sizeof(*(file->tio_id)) * file->idx->variable_count);

  file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
  memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);

  int v = 0;
  for (v = 0; v < file->idx->variable_count; v++)
  {
    file->idx_d->agg_buffer[v] = malloc(sizeof(*(file->idx_d->agg_buffer[v])) * file->idx_d->perm_layout_count);
    memset(file->idx_d->agg_buffer[v], 0, sizeof(*(file->idx_d->agg_buffer[v])) * file->idx_d->perm_layout_count);

    file->tagg_id[v] = malloc(sizeof(*(file->tagg_id[v])) * file->idx_d->perm_layout_count);
    memset(file->tagg_id[v], 0, sizeof(*(file->tagg_id[v])) * file->idx_d->perm_layout_count);

    file->tio_id[v] = malloc(sizeof(*(file->tio_id[v])) * file->idx_d->perm_layout_count);
    memset(file->tio_id[v], 0, sizeof(*(file->tio_id[v])) * file->idx_d->perm_layout_count);
  }

  return PIDX_success;
}

PIDX_return_code finalize_agg_io(PIDX_global_idx_io file, int start_index)
{
  PIDX_variable var = file->idx->variable[start_index];
  free(var->global_block_layout_files);
  free(var->block_layout_by_level_files);

  int v = 0;
  for (v = 0; v < file->idx->variable_count; v++)
  {
    free (file->idx_d->agg_buffer[v]);
    file->idx_d->agg_buffer[v] = 0;
    free (file->tagg_id[v]);
    free(file->tio_id[v]);
  }

  free(file->idx_d->agg_buffer);
  free (file->tagg_id);
  free(file->tio_id);

  return PIDX_success;
}


PIDX_return_code PIDX_global_idx_write(PIDX_global_idx_io file, int start_var_index, int end_var_index)
{
  file->idx_d->var_pipe_length = file->idx->variable_count - 1;
  if (file->idx_d->var_pipe_length == 0)
    file->idx_d->var_pipe_length = 1;

#if PIDX_DEBUG_OUTPUT
  unsigned long long l_populate = 0, g_populate = 0;
  unsigned long long l_filec = 0, g_filec = 0;
  unsigned long long l_init = 0, g_init = 0;
  unsigned long long l_rst_buf = 0, g_rst_buf = 0;
  unsigned long long l_chunk_buf = 0, g_chunk_buf = 0;
  unsigned long long l_rst = 0, g_rst = 0;
  unsigned long long l_chunk = 0, g_chunk = 0;
  unsigned long long l_cmp = 0, g_cmp = 0;
  unsigned long long l_hz_buf = 0, g_hz_buf = 0;
  unsigned long long l_hz = 0, g_hz = 0;
  unsigned long long l_agg_buf = 0, g_agg_buf = 0;
  unsigned long long l_agg = 0, g_agg = 0;
  unsigned long long l_io = 0, g_io = 0;
  unsigned long long l_pidx = 0, g_pidx = 0;
#endif


  PIDX_return_code ret;
  //static int header_io = 0;
  int nprocs = 1, rank = 0;

  PIDX_time time = file->idx_d->time;
  time->populate_idx_start_time = PIDX_get_time();
#if 1

#if PIDX_HAVE_MPI
  ret = partition_setup(file, start_var_index, end_var_index, file->idx_d->var_pipe_length);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif

  int d = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
   {
    file->idx_d->idx_count[d] = file->idx->bounds[d] / file->idx_d->idx_size[d];
    if (file->idx->bounds[d] % file->idx_d->idx_size[d] != 0)
      file->idx_d->idx_count[d]++;

    file->idx_d->idx_count[d] = pow(2, (int)ceil(log2(file->idx_d->idx_count[d])));
  }
  //printf("IDX Count = %d %d %d\n", file->idx_d->idx_count[0], file->idx_d->idx_count[1], file->idx_d->idx_count[2]);

  ret = partition(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int grank = 0, gnprocs = 1;
#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);

    MPI_Comm_rank(file->global_comm, &grank);
    MPI_Comm_size(file->global_comm, &gnprocs);
  }
#endif

  file->idx_d->rank_buffer = malloc(gnprocs * sizeof(*file->idx_d->rank_buffer));
  memset(file->idx_d->rank_buffer, 0, gnprocs * sizeof(*file->idx_d->rank_buffer));
  MPI_Allgather(&rank, 1, MPI_INT, file->idx_d->rank_buffer, 1, MPI_INT, file->global_comm);

  ret = populate_idx_file_structure(file, 0);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int partion_level = (int) log2(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  int total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  //int remainder_level = file->idx_d->maxh - total_partiton_level;

  file->idx_d->perm_layout_count = (file->idx_d->maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->perm_layout_count <= 0)
    file->idx_d->perm_layout_count = 1;

  one_time_initialize(file);

  if (total_partiton_level >= file->idx_d->maxh)
    total_partiton_level = file->idx_d->maxh;

  int file_zero_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1;
  if (file_zero_level >= file->idx_d->maxh)
    file_zero_level = file->idx_d->maxh;

  ret = init_agg_io(file, start_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  //printf("total_partiton_level = %d maxh = %d\n", total_partiton_level, file->idx_d->maxh);

#if 1
  //printf("[0 %d %d]\n", total_partiton_level, file->idx_d->maxh);
  int hz_from = 0;//file_zero_level;//total_partiton_level;
  int hz_to =  total_partiton_level;//file->idx_d->maxh;//file_zero_level;// + 1;//total_partiton_level;//file->idx_d->maxh;
  if (hz_from == hz_to)
  {
    file->idx_d->start_layout_index_shared = 0;
    file->idx_d->end_layout_index_shared = 0;

    //time->populate_idx_end_time = PIDX_get_time();
    //return PIDX_success;
  }

  ret = populate_idx_dataset_shared(file, start_var_index, end_var_index, hz_from, hz_to);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int hz_from_non_shared = total_partiton_level;//file_zero_level;//total_partiton_level;
  int hz_to_non_shared =  file->idx_d->maxh;//file_zero_level;// + 1;//total_partiton_level;//file->idx_d->maxh;
  if (hz_from_non_shared == hz_to_non_shared)
  {
    //time->populate_idx_end_time = PIDX_get_time();
    //return PIDX_success;
    file->idx_d->start_layout_index_non_shared = 0;
    file->idx_d->end_layout_index_non_shared = 0;
  }

  ret = populate_idx_dataset_non_shared(file, start_var_index, end_var_index, hz_from_non_shared, hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->populate_idx_end_time = PIDX_get_time();

  time->hz_s_time = PIDX_get_time();
  ret = create_hz_buffers(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->hz_e_time = PIDX_get_time();

  time->write_init_start[time->header_counter] = PIDX_get_time();
  if (file->idx_d->start_layout_index_shared != file->idx_d->end_layout_index_shared)
  {
    ret = write_headers_shared(file, start_var_index, end_var_index, 0);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  time->write_init_end[time->header_counter] = PIDX_get_time();
  time->header_counter++;

  time->write_init_start[time->header_counter] = PIDX_get_time();
  if (file->idx_d->start_layout_index_non_shared != file->idx_d->end_layout_index_non_shared)
  {
    ret = write_headers_non_shared(file, start_var_index, end_var_index, 0);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }
  time->write_init_end[time->header_counter] = PIDX_get_time();
  time->header_counter++;

  int start_index = 0;
  int i = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (/*file->idx_d->var_pipe_length + */1))
  {
    int agg_io_level_non_shared = 0, no_of_aggregators = 0;
    if (file->idx_d->agg_type == 1)
    {
      for (i = file->idx_d->start_layout_index_non_shared; i < file->idx_d->end_layout_index_non_shared ; i++)
      {
        no_of_aggregators = file->idx->variable[start_var_index]->block_layout_by_level_files[1][i - file->idx_d->start_layout_index_non_shared]->existing_file_count;
        if (no_of_aggregators <= nprocs)
          agg_io_level_non_shared = i + 1;
      }
      //agg_io_level_non_shared = agg_io_level_non_shared + 1;
    }

    if (file->idx->enable_agg == 0)
      agg_io_level_non_shared = file->idx_d->start_layout_index_non_shared;//0;

    //printf("[N] %d %d %d\n", file->idx_d->start_layout_index_non_shared, file->idx_d->end_layout_index_non_shared, agg_io_level_non_shared);
    ret = PIDX_global_aggregate(file, start_var_index, start_index, 1, file->idx_d->start_layout_index_non_shared, file->idx_d->end_layout_index_non_shared, file->idx_d->layout_count_non_shared, agg_io_level_non_shared);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    MPI_Status *status;
    file->idx_d->request = malloc(sizeof(*(file->idx_d->request)) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));
    memset(file->idx_d->request, 0, sizeof(*(file->idx_d->request)) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));

    file->idx_d->fp = malloc(sizeof(*(file->idx_d->fp)) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));
    memset(file->idx_d->fp, 0, sizeof(*(file->idx_d->fp)) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));

    status = malloc(sizeof(*status) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));
    memset(status, 0, sizeof(*status) * (agg_io_level_non_shared - file->idx_d->start_layout_index_non_shared));

    ret = PIDX_global_async_io(file, start_var_index, start_index, 1, file->idx_d->start_layout_index_non_shared, file->idx_d->end_layout_index_non_shared, file->idx_d->layout_count_non_shared, agg_io_level_non_shared);

    //ret = PIDX_global_io(file, start_var_index, start_index, 1, file->idx_d->start_layout_index_non_shared, file->idx_d->end_layout_index_non_shared, file->idx_d->layout_count_non_shared, agg_io_level_non_shared);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    int agg_io_level_shared = 0;
    if (file->idx_d->agg_type == 1)
    {
      for (i = file->idx_d->start_layout_index_shared; i < file->idx_d->end_layout_index_shared ; i++)
      {
        no_of_aggregators = file->idx->variable[start_var_index]->block_layout_by_level_files[0][i - file->idx_d->start_layout_index_shared]->existing_file_count;
        if (no_of_aggregators <= nprocs)
          agg_io_level_shared = i + 1;
      }
      //agg_io_level_shared = agg_io_level_shared + 1;
    }

    if (file->idx->enable_agg == 0)
      agg_io_level_shared = file->idx_d->start_layout_index_shared;

    //printf("[S] %d %d %d\n", file->idx_d->start_layout_index_shared, file->idx_d->end_layout_index_shared, agg_io_level_shared);
    ret = PIDX_global_aggregate(file, start_var_index, start_index, 0, file->idx_d->start_layout_index_shared, file->idx_d->end_layout_index_shared, file->idx_d->layout_count_shared, agg_io_level_shared);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret = PIDX_global_io(file, start_var_index, start_index, 0, file->idx_d->start_layout_index_shared, file->idx_d->end_layout_index_shared, file->idx_d->layout_count_shared, agg_io_level_shared);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    for (i = file->idx_d->start_layout_index_non_shared; i < (agg_io_level_non_shared); i++)
    {
      if (file->idx_d->request[i - file->idx_d->start_layout_index_non_shared] != 0)
      {
        //printf("[A] %d %p\n", &(file->idx_d->request[i]), &(file->idx_d->request[i]));
        ret = MPI_Wait(&(file->idx_d->request[i - file->idx_d->start_layout_index_non_shared]), &(status[i - file->idx_d->start_layout_index_non_shared]));
        if (ret != MPI_SUCCESS)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
        }

        MPI_File_close(&(file->idx_d->fp[i - file->idx_d->start_layout_index_non_shared]));
      }

      ret = PIDX_agg_buf_destroy(file->tagg_id[start_index][i], file->idx_d->agg_buffer[start_index][i]);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_agg;
      }

      free(file->idx_d->agg_buffer[start_index][i]);
      PIDX_agg_finalize(file->tagg_id[start_index][i]);
    }
    for(i = file->idx_d->start_layout_index_non_shared ; i < file->idx_d->end_layout_index_non_shared; i++)
      PIDX_file_io_finalize(file->tio_id[start_index][i]);


    free(status);
    free(file->idx_d->request);
    free(file->idx_d->fp);
  }

  delete_idx_dataset_shared(file, start_var_index, end_var_index, hz_from, hz_to);
  delete_idx_dataset_non_shared(file, start_var_index, end_var_index, hz_from_non_shared, hz_to_non_shared);
#endif

  free(file->idx_d->rank_buffer);
  ret = destroy_hz_buffers(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = finalize_agg_io(file, start_var_index);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = partition_destroy(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

#endif
  return PIDX_success;
}


PIDX_return_code PIDX_global_idx_io_finalize(PIDX_global_idx_io file)
{
#if PIDX_HAVE_MPI
  //if (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2] != 1)
  //  MPI_Comm_free(&(file->comm));
#endif

  free(file);
  file = 0;

  return PIDX_success;
}
