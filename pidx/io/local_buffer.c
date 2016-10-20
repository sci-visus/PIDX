#include "../PIDX_inc.h"

static PIDX_return_code create_non_shared_async_buffers(PIDX_io file, int start_layout_index_non_shared, int agg_io_level_non_shared);
static PIDX_return_code create_shared_async_buffers(PIDX_io file, int start_layout_index_shared, int agg_io_level_shared);
static PIDX_return_code create_file_zero_async_buffers(PIDX_io file, int start_layout_index_file_zero, int agg_io_level_file_zero);

static PIDX_return_code wait_and_destroy_non_shared_async_buffers(PIDX_io file, int start_layout_index_non_shared, int agg_io_level_non_shared);
static PIDX_return_code wait_and_destroy_shared_async_buffers(PIDX_io file, int start_layout_index_shared, int agg_io_level_shared);
static PIDX_return_code wait_and_destroy_file_zero_async_buffers(PIDX_io file, int start_layout_index_file_zero, int agg_io_level_file_zero);

static PIDX_return_code destroy_non_shared_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_non_shared, int end_layout_index_non_shared, int agg_io_level_non_shared);
static PIDX_return_code destroy_shared_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_shared, int end_layout_index_shared, int agg_io_level_shared);
static PIDX_return_code destroy_file_zero_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_file_zero, int end_layout_index_file_zero, int agg_io_level_file_zero);

static PIDX_return_code create_non_shared_async_buffers(PIDX_io file, int start_layout_index_non_shared, int agg_io_level_non_shared)
{
  file->idx_d->status_non_shared = malloc(sizeof(*(file->idx_d->status_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));
  memset(file->idx_d->status_non_shared, 0, sizeof(*(file->idx_d->status_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));

  file->idx_d->request_non_shared = malloc(sizeof(*(file->idx_d->request_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));
  memset(file->idx_d->request_non_shared, 0, sizeof(*(file->idx_d->request_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));

  file->idx_d->fp_non_shared = malloc(sizeof(*(file->idx_d->fp_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));
  memset(file->idx_d->fp_non_shared, 0, sizeof(*(file->idx_d->fp_non_shared)) * (agg_io_level_non_shared - start_layout_index_non_shared));

  return PIDX_success;
}


static PIDX_return_code create_shared_async_buffers(PIDX_io file, int start_layout_index_shared, int agg_io_level_shared)
{
  file->idx_d->status_shared = malloc(sizeof(*(file->idx_d->status_shared)) * (agg_io_level_shared - start_layout_index_shared));
  memset(file->idx_d->status_shared, 0, sizeof(*(file->idx_d->status_shared)) * (agg_io_level_shared - start_layout_index_shared));

  file->idx_d->request_shared = malloc(sizeof(*(file->idx_d->request_shared)) * (agg_io_level_shared - start_layout_index_shared));
  memset(file->idx_d->request_shared, 0, sizeof(*(file->idx_d->request_shared)) * (agg_io_level_shared - start_layout_index_shared));

  file->idx_d->fp_shared = malloc(sizeof(*(file->idx_d->fp_shared)) * (agg_io_level_shared - start_layout_index_shared));
  memset(file->idx_d->fp_shared, 0, sizeof(*(file->idx_d->fp_shared)) * (agg_io_level_shared - start_layout_index_shared));

  return PIDX_success;
}


static PIDX_return_code create_file_zero_async_buffers(PIDX_io file, int start_layout_index_file_zero, int agg_io_level_file_zero)
{
  file->idx_d->status_file_zero = malloc(sizeof(*(file->idx_d->status_shared)) * (agg_io_level_file_zero - start_layout_index_file_zero));
  memset(file->idx_d->status_file_zero, 0, sizeof(*(file->idx_d->status_shared)) * (agg_io_level_file_zero - start_layout_index_file_zero));

  file->idx_d->request_file_zero = malloc(sizeof(*(file->idx_d->request_file_zero)) * (agg_io_level_file_zero - start_layout_index_file_zero));
  memset(file->idx_d->request_file_zero, 0, sizeof(*(file->idx_d->request_file_zero)) * (agg_io_level_file_zero - start_layout_index_file_zero));

  file->idx_d->fp_file_zero = malloc(sizeof(*(file->idx_d->fp_file_zero)) * (agg_io_level_file_zero - start_layout_index_file_zero));
  memset(file->idx_d->fp_file_zero, 0, sizeof(*(file->idx_d->fp_file_zero)) * (agg_io_level_file_zero - start_layout_index_file_zero));

  return PIDX_success;
}


static PIDX_return_code wait_and_destroy_non_shared_async_buffers(PIDX_io file, int start_layout_index_non_shared, int agg_io_level_non_shared)
{
  int i = 0;
  int ret;
  for (i = start_layout_index_non_shared; i < (agg_io_level_non_shared); i++)
  {
    if (file->idx_d->request_non_shared[i - start_layout_index_non_shared] != 0)
    {
      ret = MPI_Wait(&(file->idx_d->request_non_shared[i - start_layout_index_non_shared]), &(file->idx_d->status_non_shared[i - start_layout_index_non_shared]));
      if (ret != MPI_SUCCESS)
      {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
      }

      MPI_File_close(&(file->idx_d->fp_non_shared[i - start_layout_index_non_shared]));
    }
  }

  free(file->idx_d->status_non_shared);
  free(file->idx_d->request_non_shared);
  free(file->idx_d->fp_non_shared);

  return PIDX_success;
}


static PIDX_return_code wait_and_destroy_shared_async_buffers(PIDX_io file, int start_layout_index_shared, int agg_io_level_shared)
{
  int i = 0;
  int ret;
  for (i = start_layout_index_shared; i < (agg_io_level_shared); i++)
  {
    if (file->idx_d->request_shared[i - start_layout_index_shared] != 0)
    {
      ret = MPI_Wait(&(file->idx_d->request_shared[i - start_layout_index_shared]), &(file->idx_d->status_shared[i - start_layout_index_shared]));
      if (ret != MPI_SUCCESS)
      {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
      }
      MPI_File_close(&(file->idx_d->fp_shared[i - start_layout_index_shared]));
    }
  }

  free(file->idx_d->status_shared);
  free(file->idx_d->request_shared);
  free(file->idx_d->fp_shared);

  return PIDX_success;
}


static PIDX_return_code wait_and_destroy_file_zero_async_buffers(PIDX_io file, int start_layout_index_file_zero, int agg_io_level_file_zero)
{
  int i = 0;
  int ret;
  for (i = start_layout_index_file_zero; i < (agg_io_level_file_zero); i++)
  {
    if (file->idx_d->request_file_zero[i - start_layout_index_file_zero] != 0)
    {
      ret = MPI_Wait(&(file->idx_d->request_file_zero[i - start_layout_index_file_zero]), &(file->idx_d->status_file_zero[i - start_layout_index_file_zero]));
      if (ret != MPI_SUCCESS)
      {
          fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_file;
      }
      MPI_File_close(&(file->idx_d->fp_file_zero[i - start_layout_index_file_zero]));
    }
  }

  free(file->idx_d->status_file_zero);
  free(file->idx_d->request_file_zero);
  free(file->idx_d->fp_file_zero);

  return PIDX_success;
}

static PIDX_return_code destroy_non_shared_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_non_shared, int end_layout_index_non_shared, int agg_io_level_non_shared)
{
  int ret;
  int i = 0;
  for (i = start_layout_index_non_shared; i < (agg_io_level_non_shared); i++)
  {
    ret = PIDX_agg_buf_destroy(file->idx_d->nshared_agg_buffer[start_index][i]);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx_d->nshared_agg_buffer[start_index][i]);
    PIDX_agg_finalize(file->nshared_agg_id[start_index][i]);
  }

  for(i = start_layout_index_non_shared ; i < end_layout_index_non_shared; i++)
    PIDX_file_io_finalize(file->nshared_io_id[start_index][i]);

  return PIDX_success;
}


static PIDX_return_code destroy_shared_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_shared, int end_layout_index_shared, int agg_io_level_shared)
{
  int i = 0;
  int ret;
  for (i = start_layout_index_shared; i < (agg_io_level_shared); i++)
  {
    ret = PIDX_agg_buf_destroy(file->idx_d->shared_agg_buffer[start_index][i]);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx_d->shared_agg_buffer[start_index][i]);
    PIDX_agg_finalize(file->shared_agg_id[start_index][i]);
  }

  for(i = start_layout_index_shared ; i < end_layout_index_shared; i++)
    PIDX_file_io_finalize(file->shared_io_id[start_index][i]);

  return PIDX_success;
}


static PIDX_return_code destroy_file_zero_ids_and_buffers(PIDX_io file, int start_index, int start_layout_index_file_zero, int end_layout_index_file_zero, int agg_io_level_file_zero)
{
  int i = 0;
  int ret;
  for (i = start_layout_index_file_zero; i < (agg_io_level_file_zero); i++)
  {
    ret = PIDX_agg_buf_destroy(file->idx_d->f0_agg_buffer[start_index][i]);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }

    free(file->idx_d->f0_agg_buffer[start_index][i]);
    PIDX_agg_finalize(file->f0_agg_id[start_index][i]);
  }

  for(i = start_layout_index_file_zero ; i < end_layout_index_file_zero; i++)
    PIDX_file_io_finalize(file->f0_io_id[start_index][i]);

  return PIDX_success;
}


PIDX_return_code create_agg_io_buffer(PIDX_io file, int group_index)
{
  int vc = file->idx->variable_count;
  if (vc <= 0)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int lc = file->idx_d->perm_layout_count;
  if (lc <= 0)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  idx_dataset_derived_metadata idx = file->idx_d;

  file->f0_agg_id = malloc(sizeof(*(file->f0_agg_id)) * vc);
  file->f0_io_id = malloc(sizeof(*(file->f0_io_id)) * vc);
  memset(file->f0_agg_id, 0, sizeof(*(file->f0_agg_id)) * vc);
  memset(file->f0_io_id, 0, sizeof(*(file->f0_io_id)) * vc);

  file->shared_agg_id = malloc(sizeof(*(file->shared_agg_id)) * vc);
  file->shared_io_id = malloc(sizeof(*(file->shared_io_id)) * vc);
  memset(file->shared_agg_id, 0, sizeof(*(file->shared_agg_id)) * vc);
  memset(file->shared_io_id, 0, sizeof(*(file->shared_io_id)) * vc);

  file->nshared_agg_id = malloc(sizeof(*(file->nshared_agg_id)) * vc);
  file->nshared_io_id = malloc(sizeof(*(file->nshared_io_id)) * vc);
  memset(file->nshared_agg_id, 0, sizeof(*(file->nshared_agg_id)) * vc);
  memset(file->nshared_io_id, 0, sizeof(*(file->nshared_io_id)) * vc);

  idx->f0_agg_buffer = malloc(sizeof(*(idx->f0_agg_buffer)) * vc);
  memset(idx->f0_agg_buffer, 0, sizeof(*(idx->f0_agg_buffer)) * vc);

  idx->shared_agg_buffer = malloc(sizeof(*(idx->shared_agg_buffer)) * vc);
  memset(idx->shared_agg_buffer, 0, sizeof(*(idx->shared_agg_buffer)) * vc);

  idx->nshared_agg_buffer = malloc(sizeof(*(idx->nshared_agg_buffer)) * vc);
  memset(idx->nshared_agg_buffer, 0, sizeof(*(idx->nshared_agg_buffer)) * vc);

  int v = 0;
  for (v = 0; v < vc; v++)
  {
    file->f0_agg_id[v] = malloc(sizeof(*(file->f0_agg_id[v])) * lc);
    file->f0_io_id[v] = malloc(sizeof(*(file->f0_io_id[v])) * lc);
    memset(file->f0_agg_id[v], 0, sizeof(*(file->f0_agg_id[v])) * lc);
    memset(file->f0_io_id[v], 0, sizeof(*(file->f0_io_id[v])) * lc);

    file->shared_agg_id[v] = malloc(sizeof(*(file->shared_agg_id[v])) * lc);
    file->shared_io_id[v] = malloc(sizeof(*(file->shared_io_id[v])) * lc);
    memset(file->shared_agg_id[v], 0, sizeof(*(file->shared_agg_id[v])) * lc);
    memset(file->shared_io_id[v], 0, sizeof(*(file->shared_io_id[v])) * lc);

    file->nshared_agg_id[v] = malloc(sizeof(*(file->nshared_agg_id[v])) * lc);
    file->nshared_io_id[v] = malloc(sizeof(*(file->nshared_io_id[v])) * lc);
    memset(file->nshared_agg_id[v], 0, sizeof(*(file->nshared_agg_id[v])) * lc);
    memset(file->nshared_io_id[v], 0, sizeof(*(file->nshared_io_id[v])) * lc);

    idx->f0_agg_buffer[v] = malloc(sizeof(*(idx->f0_agg_buffer[v])) * lc);
    memset(idx->f0_agg_buffer[v], 0, sizeof(*(idx->f0_agg_buffer[v])) * lc);

    idx->shared_agg_buffer[v] = malloc(sizeof(*(idx->shared_agg_buffer[v])) * lc);
    memset(idx->shared_agg_buffer[v], 0, sizeof(*(idx->shared_agg_buffer[v])) * lc);

    idx->nshared_agg_buffer[v] = malloc(sizeof(*(idx->nshared_agg_buffer[v])) * lc);
    memset(idx->nshared_agg_buffer[v], 0, sizeof(*(idx->nshared_agg_buffer[v])) * lc);
  }

  return PIDX_success;
}


PIDX_return_code destroy_agg_io_buffer(PIDX_io file)
{
  int vc = file->idx->variable_count;
  idx_dataset_derived_metadata idx = file->idx_d;

  int v = 0;
  for (v = 0; v < vc; v++)
  {
    free(file->f0_agg_id[v]);
    free(file->f0_io_id[v]);
    free(idx->f0_agg_buffer[v]);

    free(file->shared_agg_id[v]);
    free(file->shared_io_id[v]);
    free(idx->shared_agg_buffer[v]);

    free(file->nshared_agg_id[v]);
    free(file->nshared_io_id[v]);
    free(idx->nshared_agg_buffer[v]);
  }

  free(file->f0_agg_id);
  free(file->f0_io_id);

  free(file->shared_agg_id);
  free(file->shared_io_id);

  free(file->nshared_agg_id);
  free(file->nshared_io_id);

  free(idx->f0_agg_buffer);
  free(idx->shared_agg_buffer);
  free(idx->nshared_agg_buffer);

  return PIDX_success;
}


PIDX_return_code create_async_buffers(PIDX_io file, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  create_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);
  create_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);
  create_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);

  return PIDX_success;
}




PIDX_return_code wait_and_destroy_async_buffers(PIDX_io file, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  wait_and_destroy_file_zero_async_buffers(file, var_grp->f0_start_layout_index, agg_io_level_file_zero);
  wait_and_destroy_non_shared_async_buffers(file, var_grp->nshared_start_layout_index, agg_io_level_non_shared);
  wait_and_destroy_shared_async_buffers(file, var_grp->shared_start_layout_index, agg_io_level_shared);

  return PIDX_success;
}


PIDX_return_code finalize_aggregation(PIDX_io file, int start_index, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  destroy_file_zero_ids_and_buffers(file, start_index, var_grp->f0_start_layout_index, var_grp->f0_end_layout_index, agg_io_level_file_zero);

  destroy_non_shared_ids_and_buffers(file, start_index, var_grp->nshared_start_layout_index, var_grp->nshared_end_layout_index, agg_io_level_non_shared);

  destroy_shared_ids_and_buffers(file, start_index, var_grp->shared_start_layout_index, var_grp->shared_end_layout_index, agg_io_level_shared);

  return PIDX_success;
}
