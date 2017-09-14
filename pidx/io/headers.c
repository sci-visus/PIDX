#include "../PIDX_inc.h"
#include "timming.h"


static PIDX_return_code write_idx_headers_layout(PIDX_io file, int group_index, int start_var_index, int end_var_index, char* filename, char* filename_template, PIDX_block_layout bl);


PIDX_return_code write_headers(PIDX_io file, int group_index, int start_var_index, int end_var_index, int mode)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename, file->idx->current_time_step, file->idx->filename_template);
  generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename_global, file->idx->current_time_step, file->idx->filename_template_global);
  generate_file_name_template(file->idx_d->maxh, file->idx->bits_per_block, file->idx->filename_partition, file->idx->current_time_step, file->idx->filename_template_partition);

  //fprintf(stderr, "F0: %d %d\n", var_grp->f0_start_layout_index, var_grp->f0_end_layout_index);
  //fprintf(stderr, "SL: %d %d\n", var_grp->shared_start_layout_index, var_grp->shared_end_layout_index);
  //fprintf(stderr, "NS: %d %d\n", var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index);

  ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename, file->idx->filename_template, var_grp->block_layout);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

#if 0
  if (var_grp->shared_start_layout_index != var_grp->shared_end_layout_index)
  {
    ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_partition, file->idx->filename_template_partition, var_grp->shared_block_layout);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  if (var_grp->nshared_start_layout_index != var_grp->nshared_end_layout_index)
  {
    if (mode == PIDX_IDX_IO)
    {
      ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename, file->idx->filename_template, var_grp->nshared_block_layout);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
    else
    {
      if (mode == PIDX_LOCAL_PARTITION_IDX_IO)
      {
        ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_partition, file->idx->filename_template_partition, var_grp->nshared_block_layout);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
        }
      }
      else
      {
        //fprintf(stderr, "Filename %s\n", file->idx->filename_global);
        ret = write_idx_headers_layout(file, group_index, start_var_index, end_var_index, file->idx->filename_global, file->idx->filename_template_global, var_grp->nshared_block_layout);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
        }
      }
    }
  }
#endif

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



PIDX_return_code init_raw_headers_layout(PIDX_io file, int group_index, int start_var_index, int end_var_index, char* filename)
{
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[group_index];

  if (file->idx_dbg->debug_do_io == 1)
  {
    file->header_io_id = PIDX_header_io_init(file->idx, file->idx_d, file->idx_c, start_var_index, end_var_index);

    ret = PIDX_header_io_enable_raw_dump (file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;

    if (var_grp->variable_index_tracker < file->idx->variable_count )
    {
      // Create the header
      ret = PIDX_header_io_raw_file_write(file->header_io_id, filename);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }
    else if (var_grp->variable_index_tracker == file->idx->variable_count)
    {
      ret = PIDX_header_io_raw_file_write(file->header_io_id, filename);
      if (ret != PIDX_success)
        return PIDX_err_header;
    }

    ret = PIDX_header_io_write_idx (file->header_io_id, filename, file->idx->current_time_step);
    if (ret != PIDX_success)
      return PIDX_err_header;

    ret = PIDX_header_io_finalize(file->header_io_id);
    if (ret != PIDX_success)
      return PIDX_err_header;
  }

  return PIDX_success;
}
