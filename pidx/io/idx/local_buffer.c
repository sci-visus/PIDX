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


PIDX_return_code create_async_buffers(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  file->idx_d->status1 = malloc(sizeof(*(file->idx_d->status1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));
  memset(file->idx_d->status1, 0, sizeof(*(file->idx_d->status1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));

  file->idx_d->request1 = malloc(sizeof(*(file->idx_d->request1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));
  memset(file->idx_d->request1, 0, sizeof(*(file->idx_d->request1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));

  file->idx_d->fp1 = malloc(sizeof(*(file->idx_d->fp1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));
  memset(file->idx_d->fp1, 0, sizeof(*(file->idx_d->fp1)) * (var_grp->agg_level - var_grp->shared_start_layout_index));

  return PIDX_success;
}


PIDX_return_code wait_and_destroy_async_buffers(PIDX_io file, int gi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  int sli = var_grp->shared_start_layout_index;
  int agg_i = var_grp->agg_level;

  assert (sli == 0);
  int i = 0;
  int ret;
  for (i = sli; i < (agg_i); i++)
  {
    if (file->idx_d->request1[i - sli] != 0)
    {
      ret = MPI_Wait(&(file->idx_d->request1[i - sli]), &(file->idx_d->status1[i - sli]));
      if (ret != MPI_SUCCESS)
      {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
      }

      MPI_File_close(&(file->idx_d->fp1[i - sli]));
    }
  }

  free(file->idx_d->status1);
  free(file->idx_d->request1);
  free(file->idx_d->fp1);

  return PIDX_success;
}


PIDX_return_code finalize_aggregation(PIDX_io file, int gi, int local_var_index, int start_index)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  int ret;
  int i = 0;
  int i_1 = 0;
  //start_index = start_index - local_var_index;

  int sli = var_grp->shared_start_layout_index;
  int agg_i = var_grp->agg_level;

  //printf("sli and agg_i %d %d\n", sli, agg_i);
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
