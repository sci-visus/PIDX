#include "../PIDX_inc.h"
static int lgi = 0;

PIDX_return_code data_io(PIDX_io file, int gi, int lvi, int svi, int end_index, int mode)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  int j;
  int ret = 0;
  svi = svi - lvi;

  assert(var_grp->shared_start_layout_index == 0);
  for(j = var_grp->shared_start_layout_index; j < var_grp->agg_level; j++)
  {
    Agg_buffer temp_agg = idx->agg_buffer[svi][j];
    PIDX_block_layout temp_layout = var_grp->block_layout_by_level[j];

    file->io_id[svi][j] = PIDX_file_io_init(file->idx, file->idx_d, file->idx_c, svi + lvi, svi + lvi);

    if (file->idx_dbg->debug_do_io == 1)
    {
      ret = PIDX_async_aggregated_io(file->io_id[svi][j], temp_agg, temp_layout, &(idx->request1[j]), &(idx->fp1[j]), file->idx->filename_template_partition, mode);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
  }

  return PIDX_success;
}



PIDX_return_code data_aggregate(PIDX_io file, int gi, int lvi, int svi, int evi, int agg_mode, int mode )
{
  lgi = gi;
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  int j;
  PIDX_time time = file->idx_d->time;
  svi = svi - lvi;
  evi = evi - lvi;

  assert(var_grp->shared_start_layout_index == 0);
  for (j = var_grp->shared_start_layout_index; j < var_grp->agg_level; j++)
  {
    if (agg_mode == AGG_SETUP_AND_PERFORM || agg_mode == AGG_SETUP)
    {
      time->agg_init_start[lgi][svi][j] = PIDX_get_time();

      file->agg_id[svi][j] = PIDX_agg_init(file->idx, file->idx_d, file->idx_c, svi + lvi, evi + lvi, lvi);
      idx->agg_buffer[svi][j] = malloc(sizeof(*(idx->agg_buffer[svi][j])));
      memset(idx->agg_buffer[svi][j], 0, sizeof(*(idx->agg_buffer[svi][j])));

      idx->agg_buffer[svi][j]->file_number = -1;
      idx->agg_buffer[svi][j]->var_number = -1;
      idx->agg_buffer[svi][j]->sample_number = -1;

      idx->agg_buffer[svi][j]->no_of_aggregators = 0;
      idx->agg_buffer[svi][j]->aggregator_interval = 0;
      idx->agg_buffer[svi][j]->agg_f = 1;
      time->agg_init_end[lgi][svi][j] = PIDX_get_time();

      time->agg_meta_start[lgi][svi][j] = PIDX_get_time();
      ret = PIDX_agg_meta_data_create(file->agg_id[svi][j], idx->agg_buffer[svi][j], var_grp->block_layout_by_level[j]);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_meta_end[lgi][svi][j] = PIDX_get_time();

      time->agg_buf_start[lgi][svi][j] = PIDX_get_time();

      //ret = PIDX_agg_random_buf_create_multiple_level(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j], j, svi, file_status);

      ret = PIDX_agg_localized_aggregation(file->agg_id[svi][j], idx->agg_buffer[svi][j], var_grp->block_layout_by_level[j], j, svi);

      //ret = PIDX_agg_buf_create_local_uniform_dist(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j], j, svi, file_status);

      // working
      //ret = PIDX_agg_buf_create_global_uniform_dist(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j], j, svi, file_status);
      //ret = PIDX_agg_buf_create_multiple_level(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j], j, svi, file_status);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_buf_end[lgi][svi][j] = PIDX_get_time();
    }

    if (agg_mode == AGG_SETUP_AND_PERFORM || agg_mode == AGG_PERFORM)
    {
      if (file->idx_dbg->debug_do_agg == 1)
      {
        time->agg_start[lgi][svi][j] = PIDX_get_time();
        ret = PIDX_agg_global_and_local(file->agg_id[svi][j], idx->agg_buffer[svi][j], j, var_grp->block_layout_by_level[j], mode);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_end[lgi][svi][j] = PIDX_get_time();

        time->agg_compress_start[lgi][svi][j] = PIDX_get_time();
        ret = PIDX_agg_buffer_compress(file->agg_id[svi][j], idx->agg_buffer[svi][j], j, var_grp->block_layout_by_level[j], mode);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_compress_end[lgi][svi][j] = PIDX_get_time();

      }


      time->agg_meta_cleanup_start[lgi][svi][j] = PIDX_get_time();
      ret = PIDX_agg_meta_data_destroy(file->agg_id[svi][j], var_grp->block_layout_by_level[j]);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_meta_cleanup_end[lgi][svi][j] = PIDX_get_time();
    }
  }


  return PIDX_success;
}
