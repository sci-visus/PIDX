#include "../PIDX_io.h"

static PIDX_return_code write_headers_shared(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, int layout_type);

static PIDX_return_code write_headers_non_shared(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, int layout_type);

static PIDX_return_code write_headers_file_zero(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, int layout_type);

static PIDX_return_code write_headers_layout(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl);

static PIDX_return_code PIDX_file_initialize_time_step(PIDX_hybrid_idx_io file, char* filename, char* filename_template, int current_time_step);

PIDX_return_code write_headers(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, int mode)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  //printf("F: %d %d\n", var_grp->f0_start_layout_index, var_grp->f0_end_layout_index);
  //printf("S: %d %d\n", var_grp->shared_start_layout_index, var_grp->shared_end_layout_index);
  //printf("N: %d %d\n", var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index);

  if (var_grp->shared_start_layout_index != var_grp->shared_end_layout_index)
  {
    ret = write_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_partition, file->idx->filename_template_partition, var_grp->shared_block_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  if (var_grp->nshared_start_layout_index != var_grp->nshared_end_layout_index)
  {
    if (mode == PIDX_IDX_IO)
    {
      ret = write_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename, file->idx->filename_template, var_grp->nshared_block_layout);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
    else
    {
      ret = write_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_global, file->idx->filename_template_global, var_grp->nshared_block_layout);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
  }

  if (var_grp->f0_start_layout_index != var_grp->f0_end_layout_index)
  {
    ret = write_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_file_zero, file->idx->filename_template_file_zero, var_grp->f0_block_layout);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  return PIDX_success;
}




static PIDX_return_code write_headers_layout(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

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
    ret = PIDX_header_io_filename_create(file->header_io_id, bl, filename_template);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_filename_write(file->header_io_id, bl, filename, filename_template,  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_filename_write(file->header_io_id, bl, filename, filename_template, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
  return PIDX_success;
}



static PIDX_return_code write_headers_shared(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, int layout_type)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

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
    ret = PIDX_header_io_filename_create(file->header_io_id, var_grp->shared_block_layout, file->idx->filename_template_partition);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_filename_write(file->header_io_id, var_grp->shared_block_layout, filename, filename_template,  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_filename_write(file->header_io_id, var_grp->shared_block_layout, filename, filename_template, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename_partition, /*file->idx->filename_template_partition,*/ file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }
  return PIDX_success;
}


static PIDX_return_code write_headers_file_zero(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template,  int layout_type)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

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
    ret = PIDX_header_io_filename_create(file->header_io_id, var_grp->f0_block_layout, file->idx->filename_template_file_zero);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      //ret = PIDX_header_io_file_write(file->header_io_id, var_grp->f0_block_layout,  0);
      ret = PIDX_header_io_filename_write(file->header_io_id, var_grp->f0_block_layout, filename, filename_template, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_filename_write(file->header_io_id, var_grp->f0_block_layout, filename, filename_template, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename_file_zero, /*file->idx->filename_template_global,*/ file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }

  return PIDX_success;
}


static PIDX_return_code write_headers_non_shared(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template,  int layout_type)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

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
    ret = PIDX_header_io_filename_create(file->header_io_id, var_grp->nshared_block_layout, file->idx->filename_template_global);
    if (ret != PIDX_success)
      return PIDX_err_header;

    /* STEP 2 */
    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_file_write(file->header_io_id, var_grp->nshared_block_layout,  0);
      if (ret != PIDX_success)
        return PIDX_err_header;
      //file->flush_used = 1;
    }

    if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_filename_write(file->header_io_id, var_grp->nshared_block_layout, filename, filename_template, 1);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    /* STEP 3 */
    ret = PIDX_header_io_write_idx (file->header_io_id, file->idx->filename_global, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }

  return PIDX_success;
}


PIDX_return_code one_time_initialize(PIDX_hybrid_idx_io file, int mode, int io_type)
{
  int total_header_size;
  /// Initialization ONLY ONCE per IDX file
  if (file->one_time_initializations == 0)
  {
    PIDX_init_timming_buffers1(file->idx_d->time, file->idx->variable_count);
    PIDX_init_timming_buffers2(file->idx_d->time, file->idx->variable_count, file->idx_d->perm_layout_count);

    PIDX_file_initialize_time_step(file, file->idx->filename, file->idx->filename_template, file->idx->current_time_step);

    if (mode == PIDX_WRITE && io_type != PIDX_IDX_IO)
    {
      PIDX_file_initialize_time_step(file, file->idx->filename_global, file->idx->filename_template_global, file->idx->current_time_step);

      PIDX_file_initialize_time_step(file, file->idx->filename_partition, file->idx->filename_template_partition, file->idx->current_time_step);

      PIDX_file_initialize_time_step(file, file->idx->filename_file_zero, file->idx->filename_template_file_zero, file->idx->current_time_step);

    }

    total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
    file->idx_d->start_fs_block = total_header_size / file->idx_d->fs_block_size;
    if (total_header_size % file->idx_d->fs_block_size)
      file->idx_d->start_fs_block++;

    file->one_time_initializations = 1;
  }

  return PIDX_success;
}


/// TODO: get rid of this function
static PIDX_return_code PIDX_file_initialize_time_step(PIDX_hybrid_idx_io file, char* filename, char* filename_template, int current_time_step)
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
