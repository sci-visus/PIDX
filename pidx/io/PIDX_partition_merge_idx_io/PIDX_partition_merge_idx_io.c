#include "../PIDX_io.h"

#if PIDX_HAVE_MPI

static int regular_bounds[PIDX_MAX_DIMENSIONS] = {512, 512, 512, 1, 1};
static PIDX_return_code populate_idx_layout(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level);
static PIDX_return_code delete_idx_dataset(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int hz_level_from, int hz_level_to);
static PIDX_return_code populate_idx_dataset(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int hz_level_from, int hz_level_to);
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_partition_merge_idx_io file, char* file_name, int current_time_step);
//static PIDX_return_code PIDX_parameter_validate(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index);
static PIDX_return_code partition_setup(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int pipe_len);
static PIDX_return_code write_headers(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int layout_type);
static PIDX_return_code initialize_once_per_idx(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index);
static PIDX_return_code partition(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index);
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);

struct PIDX_partition_merge_idx_io_descriptor
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


/*
static PIDX_return_code PIDX_parameter_validate(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index)
{
#if PIDX_HAVE_MPI
  int nprocs = 1;
  int local_patch_count = 0, total_patch_count = 0;
  int ret;
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm, &nprocs);
    PIDX_variable var0 = file->idx->variable[start_var_index];
    if (file->idx->enable_rst == 1)
    {
      /// Calculate the total number of patches across all processes.
      local_patch_count = var0->sim_patch_count;
      ret = MPI_Allreduce(&local_patch_count, &total_patch_count, 1, MPI_INT, MPI_SUM, file->comm);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI error\n", __FILE__, __LINE__);
        return PIDX_err_mpi;
      }

      /// If total patch count != total number of processes
      /// Then NO restructuring
      if (total_patch_count != nprocs)
        file->idx->enable_rst = 0;
    }
  }
  else
  {
    file->idx->enable_rst = 0;
    file->idx->enable_agg = 0;
  }
#else
  file->idx->enable_rst = 0;
  file->idx->enable_agg = 0;
#endif

  if (file->idx->compression_type == PIDX_CHUNKING_ONLY || file->idx->compression_type == PIDX_CHUNKING_ZFP)
  {
    // No Chunking and compression without restructuring
    if (file->idx->enable_rst != 1)
    {
      file->idx->compression_type = PIDX_NO_COMPRESSION;
      int d = 0;
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        file->idx->chunk_size[d] = 1;
      file->idx->compression_bit_rate = 64;
    }
  }

  return PIDX_success;
}
*/


static PIDX_return_code populate_idx_file_structure(PIDX_partition_merge_idx_io file)
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

  int64_t* cb = file->idx->chunked_bounds;
  bounds_point.x = (int) cb[0];
  bounds_point.y = (int) cb[1];
  bounds_point.z = (int) cb[2];
  bounds_point.u = (int) cb[3];
  bounds_point.v = (int) cb[4];
  GuessBitmaskPattern(file->idx->bitSequence, bounds_point);
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


static PIDX_return_code populate_idx_layout(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level)
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



static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}



static PIDX_return_code delete_idx_dataset(PIDX_partition_merge_idx_io file, int start_index, int end_index, int start_hz_level, int end_hz_level)
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


/// TODO: get rid of this function
PIDX_return_code PIDX_file_initialize_time_step(PIDX_partition_merge_idx_io file, char* filename, int current_time_step)
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



static PIDX_return_code populate_idx_dataset(PIDX_partition_merge_idx_io file, int start_index, int end_index, int hz_level_from, int hz_level_to)
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


    //if (rank == 0)
    //  printf("LI %d %d\n", file->idx_d->start_layout_index, file->idx_d->end_layout_index);
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


  if (rank == 0 && nprocs == 2)
  {
    printf("[B] Final Block Bitmap\n");
    PIDX_blocks_print_layout(block_layout);
  }


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


PIDX_partition_merge_idx_io PIDX_partition_merge_idx_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_partition_merge_idx_io partition_merge_idx_io;
  partition_merge_idx_io = malloc(sizeof (*partition_merge_idx_io));
  memset(partition_merge_idx_io, 0, sizeof (*partition_merge_idx_io));

  partition_merge_idx_io->idx = idx_meta_data;
  partition_merge_idx_io->idx_d = idx_derived_ptr;
  partition_merge_idx_io->idx_dbg = idx_dbg;

  return (partition_merge_idx_io);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_partition_merge_idx_io_set_communicator(PIDX_partition_merge_idx_io id, MPI_Comm comm)
{
  if (id == NULL)
    return PIDX_err_id;

  id->global_comm = comm;

  return PIDX_success;
}
#endif



static PIDX_return_code partition(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index)
{
  int *colors;
  int index_i = 0, index_j = 0, index_k = 0;
  int i = 0, j = 0, k = 0, d = 0;
  //PIDX_return_code ret;

  int rank = 0;
  MPI_Comm_rank(file->global_comm, &rank);

  //printf("Count %d %d %d\n", file->idx_d->idx_count[0], file->idx_d->idx_count[1], file->idx_d->idx_count[2]);

  colors = malloc(sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  memset(colors, 0, sizeof(*colors) * file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  file->idx_d->color = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]) + 1;//MPI_UNDEFINED;

  for (k = 0; k < file->idx_d->idx_count[2]; k++)
    for (j = 0; j < file->idx_d->idx_count[1]; j++)
      for (i = 0; i < file->idx_d->idx_count[0]; i++)
      {
        colors[(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i] = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i;
        //if (rank == 0)
        //  printf("COLORS [%d %d %d] [%d] --> %d\n", k, j, i, (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i, colors[(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * k) + (file->idx_d->idx_count[0] * j) + i]);
      }

  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));

  PIDX_variable var0 = file->idx->variable[start_var_index];
  if (var0->patch_group_count == 1)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = file->idx->variable[start_var_index]->rst_patch_group[0]->reg_patch->offset[d];
      local_proc_patch->size[d] = file->idx->variable[start_var_index]->rst_patch_group[0]->reg_patch->size[d];
    }
    //printf("[R %d] : %d %d %d - %d %d %d\n", rank, local_proc_patch->offset[0], local_proc_patch->offset[1], local_proc_patch->offset[2], local_proc_patch->size[0], local_proc_patch->size[1], local_proc_patch->size[2]);

    Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));
    int distance_x = 0, distance_y = 0, distance_z = 0;

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

    if ((regular_bounds[0]) > file->idx->bounds[0])
      regular_bounds[0] = file->idx->bounds[0];
    if ((regular_bounds[1]) > file->idx->bounds[1])
      regular_bounds[1] = file->idx->bounds[1];
    if ((regular_bounds[2]) > file->idx->bounds[2])
      regular_bounds[2] = file->idx->bounds[2];

    for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + regular_bounds[0], index_i++)
    {
      for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + regular_bounds[1], index_j++)
      {
        for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + regular_bounds[2], index_k++)
        {
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = regular_bounds[0];
          reg_patch->size[1] = regular_bounds[1];
          reg_patch->size[2] = regular_bounds[2];

          //Edge regular patches
          /*
          if ((i + regular_bounds[0]) > file->idx->bounds[0])
            reg_patch->size[0] = file->idx->bounds[0] - i;
          if ((j + regular_bounds[1]) > file->idx->bounds[1])
            reg_patch->size[1] = file->idx->bounds[1] - j;
          if ((k + regular_bounds[2]) > file->idx->bounds[2])
            reg_patch->size[2] = file->idx->bounds[2] - k;
            */

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

            //int clr = (file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * index_k) + (file->idx_d->idx_count[0] * index_j) + (index_i);
            //printf("[R %d] %d ----> %d Color %d\n", rank, clr, z_order, colors[z_order]);

            //file->idx_d->color = colors[clr];
            file->idx_d->color = colors[z_order];

            distance_x = index_i * regular_bounds[0];
            distance_y = index_j * regular_bounds[1];
            distance_z = index_k * regular_bounds[2];

            int start_index = 0;
            //printf("start_var_index = %d end_var_index = %d pipe %d\n", start_var_index, end_var_index, file->idx_d->var_pipe_length);
            for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + 1)
            {
              PIDX_variable var = file->idx->variable[start_index];

              var->rst_patch_group[0]->reg_patch->offset[0] = var->rst_patch_group[0]->reg_patch->offset[0] - distance_x;
              var->rst_patch_group[0]->reg_patch->offset[1] = var->rst_patch_group[0]->reg_patch->offset[1] - distance_y;
              var->rst_patch_group[0]->reg_patch->offset[2] = var->rst_patch_group[0]->reg_patch->offset[2] - distance_z;


              Ndim_patch_group out_patch = var->chunk_patch_group[0];
              out_patch->count = 1;

              out_patch->reg_patch->offset[0] = var->rst_patch_group[0]->reg_patch->offset[0];
              out_patch->reg_patch->offset[1] = var->rst_patch_group[0]->reg_patch->offset[1];
              out_patch->reg_patch->offset[2] = var->rst_patch_group[0]->reg_patch->offset[2];
              out_patch->reg_patch->offset[3] = 0;
              out_patch->reg_patch->offset[4] = 0;

              out_patch->reg_patch->size[0] = var->rst_patch_group[0]->reg_patch->size[0];
              out_patch->reg_patch->size[1] = var->rst_patch_group[0]->reg_patch->size[1];
              out_patch->reg_patch->size[2] = var->rst_patch_group[0]->reg_patch->size[2];
              out_patch->reg_patch->size[3] = 1;
              out_patch->reg_patch->size[4] = 1;

              memcpy(out_patch->patch[0]->size, out_patch->reg_patch->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
              memcpy(out_patch->patch[0]->offset, out_patch->reg_patch->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

              //var->chunk_patch_group[p]->reg_patch->size

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
            }

            file->idx->bounds[0] = reg_patch->size[0];
            file->idx->bounds[1] = reg_patch->size[1];
            file->idx->bounds[2] = reg_patch->size[2];

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

    int start_index = 0;
    for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
    {
      PIDX_variable var = file->idx->variable[start_index];
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
  }

  MPI_Comm_split(file->global_comm, file->idx_d->color, rank, &(file->comm));
  free(colors);


  int nprocs;
  MPI_Comm_rank(file->comm, &rank);
  MPI_Comm_size(file->comm, &nprocs);
  //printf("NP [%d] ------ R %d\n", nprocs, rank);

  /*
  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';

  if (file->idx_d->idx_count[0] != 1 || file->idx_d->idx_count[1] != 1 || file->idx_d->idx_count[2] != 1)
    sprintf(file->idx->filename, "%s_%d.idx", file_name_skeleton, file->idx_d->color);
  */

  //printf("[%d %d] File name = %s %d\n", rank, nprocs, file->idx->filename, file->idx_d->color);
  free(local_proc_patch);
  local_proc_patch = 0;

  return PIDX_success;
}


static PIDX_return_code partition_setup(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int pipe_len)
{

  int d = 0, ret = 0, nprocs = 1;
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

  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (pipe_len + 1))
  {
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the chunking ID */
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the compression ID */
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

#if PIDX_HAVE_MPI
    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    /* Attaching the communicator to the chunking phase */
    //
    ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
    //
#endif

    int reg_box_size = 32;
    int64_t reg_patch_size[PIDX_MAX_DIMENSIONS] = {1,1,1,1,1};

    restructure_loop:

    for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
      reg_patch_size[d] = reg_box_size;

    memcpy(file->idx->reg_patch_size, reg_patch_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));

    ret = PIDX_rst_meta_data_create(file->rst_id);
    if (ret != PIDX_success)
      return PIDX_err_rst;

    int mpgc;
    int pgc = file->idx->variable[start_index]->patch_group_count;
    MPI_Allreduce(&pgc, &mpgc, 1, MPI_INT, MPI_MAX, file->comm);
    if (mpgc > 1)
    {
      ret = PIDX_rst_meta_data_destroy(file->rst_id);
      if (ret != PIDX_success)
        return PIDX_err_rst;

      reg_box_size = reg_box_size * 2;
      goto restructure_loop;
    }


    //printf("[%d] Reg Box Count %d\n", reg_box_size, file->idx->reg_patch_size[0]);

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

    //
    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }
    //
#if 1
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
      /* Attaching the communicator to the compression phase */
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }

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
#endif
  }
  return PIDX_success;
}

static PIDX_return_code initialize_once_per_idx(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index)
{
  //PIDX_return_code ret;
  int total_header_size;

  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, file->idx_d->perm_layout_count);

    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;

    //ret = PIDX_parameter_validate(file, start_var_index, end_var_index);
    //if (ret != PIDX_success)
    //  return PIDX_err_file;

    file->one_time_initializations = 1;
  }
  return PIDX_success;
}

static PIDX_return_code write_headers(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int layout_type)
{

  /*  STEP 1: create files and folders based on the extents of the variable group
   *  STEP 2: if flush used
   *            when variable count is met, then write header information
   *          else
   *            pass the header buffer to the agg phase when no var pipelining is done (else if pipe, then go to if)
   *  STEP 3: at the end of all IO, write the .idx file
   */

  PIDX_return_code ret;
#if !SIMULATE_IO
  if (file->idx_dbg->debug_do_io == 1 && file->idx_d->maxh != 0)
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
      /*
      // Write the header
      if (file->flush_used == 1)
      {
        ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout, 1);
        if (ret != PIDX_success)
          return PIDX_err_header;
      }
      else if ( caching_state == 0)
      {
        ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout, 1);
        if (ret != PIDX_success)
          return PIDX_err_header;
      }
      */
      ret = PIDX_header_io_file_write(file->header_io_id, file->idx->variable[start_var_index]->global_block_layout, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    if (layout_type == 0)
    {
      ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename, file->idx->current_time_step);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }

#endif

  return PIDX_success;
}


static PIDX_return_code PIDX_partition_merge_write_io(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index, int start_res_level, int end_res_level, int partition_index)
{
  PIDX_return_code ret;
  int start_index = 0, end_index = 0;
  int i = 0, j = 0;
  int nprocs, rank;
  int gnprocs, grank;

  //PIDX_time time = file->idx_d->time;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);

    MPI_Comm_size(file->global_comm,  &gnprocs);
    MPI_Comm_rank(file->global_comm,  &grank);
  }
#endif

  PIDX_time time = file->idx_d->time;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    time->startup_start[time->variable_counter] = PIDX_get_time();
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    int agg_io_level = 0, no_of_aggregators = 0;
    if (file->idx_d->agg_type != 0)
    {
      //for (i = 0; i < file->idx_d->layout_count ; i++)
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
      //printf("[%d] file->idx->variable[start_var_index]->global_block_layout->existing_file_count = %d %d\n", rank, file->idx->variable[start_var_index]->global_block_layout->existing_file_count, nprocs);
      no_of_aggregators = file->idx->variable[start_var_index]->global_block_layout->existing_file_count;
      if (no_of_aggregators <= nprocs)
        agg_io_level = file->idx_d->end_layout_index;//file->idx_d->layout_count;
      else
        agg_io_level = file->idx_d->start_layout_index;//0;
    }

    if (file->idx->enable_agg == 0)
      agg_io_level = file->idx_d->start_layout_index;//0;


    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/

    /* Create the HZ encoding ID */
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the chunking ID */
    //file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);


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
      {
        //printf("LC %d j %d index %d\n", file->idx_d->layout_count, j, (j - file->idx_d->start_layout_index));
        file->tio_id[i][j - file->idx_d->start_layout_index] = PIDX_file_io_init(file->idx, file->idx_d, start_var_index, i, agg_end_index);
      }
    }

    file->idx_d->agg_buffer = malloc(sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);
    memset(file->idx_d->agg_buffer, 0, sizeof(*(file->idx_d->agg_buffer)) * file->idx->variable_count);

    //printf("(agg_io_level - file->idx_d->start_layout_index) = [%d - %d] %d\n", agg_io_level, file->idx_d->start_layout_index, (agg_io_level - file->idx_d->start_layout_index));
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      file->idx_d->agg_buffer[i] = malloc(sizeof(*(file->idx_d->agg_buffer[i])) * (agg_io_level - file->idx_d->start_layout_index));
      memset(file->idx_d->agg_buffer[i], 0, sizeof(*(file->idx_d->agg_buffer[i])) * (agg_io_level - file->idx_d->start_layout_index));

      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      //for(j = 0 ; j < agg_io_level; j++)
      {
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index] = malloc(sizeof(*(file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index])) );
        memset(file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index], 0, sizeof(*(file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index])) );

        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->file_number = -1;
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->var_number = -1;
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->sample_number = -1;

        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->no_of_aggregators = 0;
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->aggregator_interval = 0;
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]->aggregation_factor = 1;
      }

      /*
      for (j = 1 ; j < agg_io_level; j++)
        file->idx_d->agg_buffer[i][j]->aggregation_factor = 1;//(int)pow(2, (agg_io_level - j));
      */
    }
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/




    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {

      /* Attaching the communicator to the HZ encodig phase phase */
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      ret = PIDX_hz_encode_set_global_communicator(file->hz_id, file->global_comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }

      /* Attaching the communicator to the aggregation phase */
      /* Attaching the communicator to the I/O phase */
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        //for(j = 0 ; j < file->idx_d->layout_count; j++)
        {
          ret = PIDX_agg_set_communicator(file->tagg_id[i][j - file->idx_d->start_layout_index], file->comm);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_agg;
          }
        }

        for(j = file->idx_d->start_layout_index ; j < file->idx_d->end_layout_index; j++)
        {
          ret = PIDX_file_io_set_communicator(file->tio_id[i][j - file->idx_d->start_layout_index], file->comm);
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

    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      {
        ret = PIDX_agg_meta_data_create(file->tagg_id[i][j - file->idx_d->start_layout_index], file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->global_block_layout, file->idx->variable[start_var_index]->block_layout_by_level[j - file->idx_d->start_layout_index]);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }
    }
    time->startup_end[time->variable_counter] = PIDX_get_time();


    /*---------------------------------------------HZ [start]------------------------------------------------*/
    time->hz_start[time->variable_counter] = PIDX_get_time();

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
      //printf("RN %d %d ----> %d %d\n", grank, gnprocs, start_res_level, end_res_level);
      ret = PIDX_hz_encode_write(file->hz_id);
      //ret = PIDX_hz_encode_write_inverse(file->hz_id, 0, file->idx_d->maxh);
      //ret = PIDX_hz_encode_write_inverse(file->hz_id, start_res_level, end_res_level);
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

    time->hz_end[time->variable_counter] = PIDX_get_time();
    /*----------------------------------------------HZ [end]-------------------------------------------------*/



    /*----------------------------------------------Agg [staart]-----------------------------------------------*/


    /* Creating the buffers required for Aggregation */
    //static_var_counter = 0;
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      //for (j = 0 ; j < agg_io_level; j++)
      {
        time->agg_buf_start[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();

        /* Creating the buffers required for Aggregation */
        ret = PIDX_agg_buf_create(file->tagg_id[i][j - file->idx_d->start_layout_index], file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->block_layout_by_level[j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->global_block_layout, i, j);
        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_buf_end[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();
      }
      //static_var_counter++;
    }

#if PIDX_DEBUG_OUTPUT
    l_agg_buf = 1;
    MPI_Allreduce(&l_agg_buf, &g_agg_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg_buf == nprocs)
      printf("[A] Aggregation Buffer Created\n");
#endif


    /* Perform Aggregation */
    if (file->idx_dbg->debug_do_agg == 1)
    {
      //static_var_counter = 0;
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        //for(j = 0 ; j < agg_io_level; j++)
        for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        {
           time->agg_start[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();

           ret = PIDX_agg(file->tagg_id[i][j - file->idx_d->start_layout_index], file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index], j - file->idx_d->start_layout_index, file->idx->variable[start_var_index]->block_layout_by_level[j - file->idx_d->start_layout_index], PIDX_WRITE);


           if (ret != PIDX_success)
           {
             fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
             return PIDX_err_rst;
           }

           time->agg_end[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();
        }
      }
    }

#if 1
    if (file->idx_dbg->debug_do_io == 1)
    {
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        for(j = agg_io_level ; j < file->idx_d->end_layout_index; j++)
        //for(j = agg_io_level; j < file->idx_d->layout_count; j++)
        {
          time->io_per_process_start[i + partition_index * file->idx->variable_count][j - agg_io_level] = PIDX_get_time();

          ret = PIDX_file_io_per_process(file->tio_id[i][j - agg_io_level], file->idx->variable[start_var_index]->block_layout_by_level[j - agg_io_level], PIDX_WRITE);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
          time->io_per_process_end[i + partition_index * file->idx->variable_count][j - agg_io_level] = PIDX_get_time();
        }
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_agg = 1;
    MPI_Allreduce(&l_agg, &g_agg, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_agg == nprocs)
      printf("[A] Aggregation Completed\n");
#endif

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    /*--------------------------------------------Agg [end]--------------------------------------------------*/




    /*--------------------------------------------IO [start]--------------------------------------------------*/
    if (file->idx_dbg->debug_do_io == 1)
    {
      for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
      {
        //for(j = agg_io_level - 1 ; j < agg_io_level; j++)
        for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        //for(j = 0 ; j < agg_io_level; j++)
        {
          time->io_start[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();

          ret = PIDX_aggregated_io(file->tio_id[i][j - file->idx_d->start_layout_index], file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->block_layout_by_level[j - file->idx_d->start_layout_index], PIDX_WRITE);
          if (ret != PIDX_success)
          {
            fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }

          time->io_end[i + partition_index * file->idx->variable_count][j - file->idx_d->start_layout_index] = PIDX_get_time();
        }
      }
    }

#if PIDX_DEBUG_OUTPUT
    l_io = 1;
    MPI_Allreduce(&l_io, &g_io, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
    if (rank == 0 && g_io == nprocs)
      printf("[I] I/O completed\n");
#endif

    time->finalize_start[time->variable_counter] = PIDX_get_time();
    /* Destroy buffers allocated during aggregation phase */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = 0 ; j < agg_io_level; j++)
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      {
        PIDX_agg_buf_destroy(file->tagg_id[i][j - file->idx_d->start_layout_index], file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_agg;
        }
      }
    }


    /*----------------------------------------------IO [end]--------------------------------------------------*/


    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = 0 ; j < agg_io_level; j++)
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      {
        free(file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index]);
        file->idx_d->agg_buffer[i][j - file->idx_d->start_layout_index] = 0;
      }
      free(file->idx_d->agg_buffer[i]);
      file->idx_d->agg_buffer[i] = 0;
    }
    free(file->idx_d->agg_buffer);\
    file->idx_d->agg_buffer = 0;


    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      //for (j = 0 ; j < agg_io_level; j++)
      {
        ret = PIDX_agg_meta_data_destroy(file->tagg_id[i][j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->block_layout_by_level[j - file->idx_d->start_layout_index], file->idx->variable[start_var_index]->global_block_layout);

        if (ret != PIDX_success)
        {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
      }
    }
 #if 1
    ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }


    /*-------------------------------------------finalize [start]---------------------------------------------*/

    //if (file->small_agg_comm == 1)
    //  PIDX_destroy_local_aggregation_comm(file->agg_id);

    /* Deleting the I/O ID */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = 0 ; j < file->idx_d->layout_count; j++)
      //for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
      for(j = file->idx_d->start_layout_index ; j < file->idx_d->end_layout_index; j++)
        PIDX_file_io_finalize(file->tio_id[i][j - file->idx_d->start_layout_index]);
    }

    /* Deleting the aggregation ID */
    for(i = start_index ; i < (end_index + 1) ; i = i + (agg_var_pipe + 1))
    {
      //for(j = 0 ; j < file->idx_d->layout_count; j++)
      for(j = file->idx_d->start_layout_index ; j < agg_io_level; j++)
        PIDX_agg_finalize(file->tagg_id[i][j - file->idx_d->start_layout_index]);
    }

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

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);

#endif

    time->finalize_end[time->variable_counter] = PIDX_get_time();
    time->variable_counter++;
    /*-----------------------------------------finalize [end]--------------------------------------------------*/
#endif

    //if (rank == 0)
    //  printf("Finished Writing %d variables\n", end_index - start_index + 1);
    //start_index++;
  }

  return PIDX_success;
}


PIDX_return_code PIDX_partition_merge_idx_write(PIDX_partition_merge_idx_io file, int start_var_index, int end_var_index)
{
  int d = 0;
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
  //int nprocs = 1;


  PIDX_time time = file->idx_d->time;

  ret = populate_idx_file_structure(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

#if 1
  file->idx_d->perm_layout_count = (file->idx_d->maxh - (file->idx->bits_per_block + log2(file->idx->blocks_per_file)));
  if (file->idx_d->perm_layout_count <= 0)
    file->idx_d->perm_layout_count = 1;
#else
  file->idx_d->perm_layout_count = 1;
#endif

  ret = initialize_once_per_idx(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
    return PIDX_err_file;

  time->populate_idx_start_time = PIDX_get_time();

  ret = partition_setup(file, start_var_index, end_var_index, file->idx_d->var_pipe_length);
  if (ret != PIDX_success)
    return PIDX_err_file;


  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    file->idx_d->idx_count[d] = file->idx->bounds[d] / regular_bounds[d];
    if (file->idx->bounds[d] % regular_bounds[d] != 0)
      file->idx_d->idx_count[d]++;

    //printf("DD %d ----> %d (%d / %d)\n", d, file->idx_d->idx_count[d], file->idx->bounds[d], regular_bounds[d]);
    file->idx_d->idx_count[d] = pow(2, (int)ceil(log2(file->idx_d->idx_count[d])));
  }

  int partion_level = (int) log2(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2]);
  //partion_level = (int) pow(2, ((int)log2(file->idx_d->idx_count[0] * file->idx_d->idx_count[1] * file->idx_d->idx_count[2])));
  //printf("PL [%d %d %d]: %d\n", file->idx_d->idx_count[0], file->idx_d->idx_count[1], file->idx_d->idx_count[2], partion_level);
  //printf("file->idx->blocks_per_file = %d %d\n", file->idx->blocks_per_file, log2(file->idx->blocks_per_file));


  int total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (total_partiton_level > file->idx_d->maxh)
    total_partiton_level = file->idx_d->maxh;

  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    file->idx_d->idx_count[d] = 1;
  file->idx_d->color = 0;

  ret = populate_idx_dataset(file, start_var_index, end_var_index, 0, total_partiton_level /*file->idx_d->maxh*/);
  if (ret != PIDX_success)
    return PIDX_err_file;

  int remainder_level = file->idx_d->maxh - total_partiton_level;

  //printf("To level [%d + %d (%d) + %d] : %d RL %d \n", file->idx->bits_per_block,
  //                                          (int)log2(file->idx->blocks_per_file) + 1,
  //                                          file->idx->blocks_per_file,
  //                                          partion_level,
  //                                          file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level,
  //                                          remainder_level);

  //printf("[1] Maxh = %d %d [%d %d]\n", file->idx_d->maxh, file->idx_d->layout_count, file->idx_d->start_layout_index, file->idx_d->end_layout_index);


  ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
  if (ret != PIDX_success)
    return PIDX_err_file;

  time->populate_idx_end_time = PIDX_get_time();

  time->write_init_start[time->header_counter] = PIDX_get_time();
  ret = write_headers(file, start_var_index, end_var_index, 0);
  if (ret != PIDX_success)
    return PIDX_err_file;
  time->write_init_end[time->header_counter] = PIDX_get_time();
  time->header_counter++;

  ret = PIDX_partition_merge_write_io(file, start_var_index, end_var_index, 0, total_partiton_level, 0);
  if (ret != PIDX_success)
    return PIDX_err_file;


  delete_idx_dataset(file, start_var_index, end_var_index, 0, total_partiton_level);
#if 1
  time->partition_start_time = PIDX_get_time();
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    file->idx_d->idx_count[d] = file->idx->bounds[d] / regular_bounds[d];
    if (file->idx->bounds[d] % regular_bounds[d] != 0)
      file->idx_d->idx_count[d]++;

    file->idx_d->idx_count[d] = pow(2, (int)ceil(log2(file->idx_d->idx_count[d])));
  }

  ret = partition(file, start_var_index, end_var_index);
  if (ret != PIDX_success)
    return PIDX_err_file;


  ret = populate_idx_file_structure(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in populate_idx_file_structure\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  //printf("[%d %d %d %d %d] max H = %d\n", file->idx->chunked_bounds[0], file->idx->chunked_bounds[1], file->idx->chunked_bounds[2], file->idx->chunked_bounds[3], file->idx->chunked_bounds[4], file->idx_d->maxh);

  ret = populate_idx_dataset(file, start_var_index, end_var_index, file->idx_d->maxh - remainder_level, file->idx_d->maxh);
  if (ret != PIDX_success)
    return PIDX_err_file;
  time->partition_end_time = PIDX_get_time();

  if (file->idx_d->maxh == 0 || remainder_level == 0)
    goto cleanup;

  //printf("[B] Maxh = [%d %d] %d [%d %d]\n", file->idx_d->maxh, file->idx_d->maxh - remainder_level, file->idx_d->layout_count, file->idx_d->start_layout_index, file->idx_d->end_layout_index);


  //ret = PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->current_time_step);
  //if (ret != PIDX_success)
  //  return PIDX_err_file;

  time->write_init_start[time->header_counter] = PIDX_get_time();
  ret = write_headers(file, start_var_index, end_var_index, 1);
  if (ret != PIDX_success)
    return PIDX_err_file;
  time->write_init_end[time->header_counter] = PIDX_get_time();
  time->header_counter++;


  ret = PIDX_partition_merge_write_io(file, start_var_index, end_var_index, file->idx_d->maxh - remainder_level, file->idx_d->maxh, 1);
  if (ret != PIDX_success)
    return PIDX_err_file;


  /* Destroy buffers allocated during restructuring phase */
  ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  delete_idx_dataset(file, start_var_index, end_var_index, 0, total_partiton_level);

  cleanup:
  ret = PIDX_chunk_buf_destroy(file->chunk_id);
  if (ret != PIDX_success)
    return PIDX_err_chunk;

  ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
  if (ret != PIDX_success)
    return PIDX_err_rst;

  ret = PIDX_rst_meta_data_destroy(file->rst_id);
  if (ret != PIDX_success)
    return PIDX_err_rst;

  /* Deleting the compression ID */
  PIDX_compression_finalize(file->comp_id);

  /* Deleting the chunking ID */
  PIDX_chunk_finalize(file->chunk_id);

  /* Deleting the restructuring ID */
  PIDX_rst_finalize(file->rst_id);


#endif


#if PIDX_DEBUG_OUTPUT
  l_pidx = 1;
  MPI_Allreduce(&l_pidx, &g_pidx, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, file->comm);
  if (rank == 0 && g_pidx == nprocs)
    printf("PIDX closing file\n");

#endif

  free(file->idx_d->rank_r_offset);
  file->idx_d->rank_r_offset = 0;

  free(file->idx_d->rank_r_count);
  file->idx_d->rank_r_count = 0;



  return PIDX_success;
}



PIDX_return_code PIDX_partition_merge_idx_io_finalize(PIDX_partition_merge_idx_io file)
{

  free(file);
  file = 0;

  return PIDX_success;
}
#endif
