#include "../../PIDX_inc.h"


PIDX_return_code create_agg_io_buffer(PIDX_io file, int gi, int svi, int evi)
{
  int lc = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  assert (var_grp->shared_start_layout_index == 0);

  int vc =  file->idx->variable_count;// (evi - svi);
  if (vc <= 0)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  idx_dataset_derived_metadata idx = file->idx_d;

  file->agg_id = malloc(sizeof(*(file->agg_id)) * vc);
  file->io_id = malloc(sizeof(*(file->io_id)) * vc);
  memset(file->agg_id, 0, sizeof(*(file->agg_id)) * vc);
  memset(file->io_id, 0, sizeof(*(file->io_id)) * vc);

  idx->agg_buffer = malloc(sizeof(*(idx->agg_buffer)) * vc);
  memset(idx->agg_buffer, 0, sizeof(*(idx->agg_buffer)) * vc);

  int v = 0;
  for (v = 0; v < vc; v++)
  {
    lc = (var_grp->agg_level - var_grp->shared_start_layout_index);
    file->agg_id[v] = malloc(sizeof(*(file->agg_id[v])) * lc);
    file->io_id[v] = malloc(sizeof(*(file->io_id[v])) * lc);
    memset(file->agg_id[v], 0, sizeof(*(file->agg_id[v])) * lc);
    memset(file->io_id[v], 0, sizeof(*(file->io_id[v])) * lc);

    idx->agg_buffer[v] = malloc(sizeof(*(idx->agg_buffer[v])) * lc);
    memset(idx->agg_buffer[v], 0, sizeof(*(idx->agg_buffer[v])) * lc);
  }

  return PIDX_success;
}


PIDX_return_code destroy_agg_io_buffer(PIDX_io file, int svi, int evi)
{
  //int vc = evi - svi;
  int vc =  file->idx->variable_count;// (evi - svi);
  idx_dataset_derived_metadata idx = file->idx_d;

  int v = 0;
  for (v = 0; v < vc; v++)
  {
    free(file->agg_id[v]);
    free(file->io_id[v]);
    free(idx->agg_buffer[v]);
  }

  free(file->agg_id);
  free(file->io_id);
  free(idx->agg_buffer);

  return PIDX_success;
}



PIDX_return_code finalize_aggregation(PIDX_io file, int gi, int start_index)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  int ret;
  int i = 0;
  int i_1 = 0;
  //start_index = start_index - local_var_index;

  int sli = var_grp->shared_start_layout_index;
  int agg_i = var_grp->agg_level;

  //fprintf(stderr, "[%d] sli and agg_i %d %d si %d\n", file->idx_c->grank, sli, agg_i, start_index);
  for (i = sli; i < agg_i; i++)
  {
    i_1 = i - sli;
    ret = PIDX_agg_buf_destroy(file->idx_d->agg_buffer[start_index][i_1]);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx_d->agg_buffer[start_index][i_1]);
    PIDX_agg_finalize(file->agg_id[start_index][i_1]);
    PIDX_file_io_finalize(file->io_id[start_index][i_1]);
  }

  return PIDX_success;
}
