#include "../../PIDX_inc.h"
#include "timming.h"


static PIDX_return_code write_idx_headers_layout(PIDX_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl);

PIDX_return_code write_global_idx(PIDX_io file, int start_var_index, int end_var_index, int mode)
{
  int ret = 0;

  generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename, file->idx->current_time_step, file->idx->filename_template);

  if (mode == PIDX_READ)
    return PIDX_success;

  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, file->idx_c, start_var_index, end_var_index);

    ret = PIDX_header_io_global_idx_write (file->header_io_id, file->idx->filename);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    PIDX_header_io_finalize(file->header_io_id);
  }

  return PIDX_success;
}


PIDX_return_code write_headers(PIDX_io file, int group_index, int start_var_index, int end_var_index, int mode)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  // When using non-partitioned idx IO there is no need to create different IDX file per partitions
  // so we use one single .idx file.
  // Note: This strcpy is a little hacky, but works.
  // The .idx file will be written exactly as a single partition file
  if(file->idx->io_type == PIDX_IDX_IO)
    strcpy(file->idx->filename_partition,file->idx->filename);
  
  generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename_partition, file->idx->current_time_step, file->idx->filename_template_partition);

  //fprintf(stderr, "FN %s FNT %s\n", file->idx->filename_partition, file->idx->filename_template_partition);

  if (mode == PIDX_READ)
    return PIDX_success;

  ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_partition, file->idx->filename_template_partition, var_grp->block_layout);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



static PIDX_return_code write_idx_headers_layout(PIDX_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  if (file->idx_dbg->debug_do_io == 1)
  {
    /* STEP 1 */
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, file->idx_c, start_var_index, end_var_index);

    ret = PIDX_header_io_idx_file_create(file->header_io_id, bl, filename_template);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    /* STEP 2 */
    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_idx_file_write(file->header_io_id, bl, filename_template,  0);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_header;
      }
    }

    if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_idx_file_write(file->header_io_id, bl, filename_template, 1);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_header;
      }
    }

    /* STEP 3 */
    ret = PIDX_header_io_partition_idx_write(file->header_io_id, filename);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }
  }

  return PIDX_success;
}



PIDX_return_code init_raw_headers_layout(PIDX_io file, int group_index, int start_var_index, int end_var_index, char* filename)
{
  int ret = 0;

  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, file->idx_c, start_var_index, end_var_index);

    ret = PIDX_header_io_raw_dir_create(file->header_io_id, filename);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    ret = PIDX_header_io_raw_idx_write (file->header_io_id, filename);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_header;
    }
  }

  return PIDX_success;
}
