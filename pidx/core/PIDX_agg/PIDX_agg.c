#include "PIDX_agg.h"


struct PIDX_agg_struct
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
  MPI_Comm global_comm;
  MPI_Comm local_comm;
  MPI_Win win;
#endif

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  PIDX_local_agg_id local_id;
  PIDX_global_agg_id global_id;

  int init_index;
  int first_index;
  int last_index;
};



PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  agg_id->idx = idx_meta_data;
  agg_id->idx_d = idx_d;

  agg_id->local_id = PIDX_local_agg_init(idx_meta_data, idx_d, init_index, first_index, last_index);
  agg_id->global_id = PIDX_global_agg_init(idx_meta_data, idx_d, init_index, first_index, last_index);

  agg_id->init_index = init_index;
  agg_id->first_index = first_index;
  agg_id->last_index = last_index;

  return agg_id;
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_agg_set_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->comm = comm;

  PIDX_local_agg_set_communicator(agg_id->local_id, comm);
  PIDX_global_agg_set_communicator(agg_id->global_id, comm);

  return PIDX_success;
}


PIDX_return_code PIDX_create_local_aggregation_comm(PIDX_agg_id agg_id)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  int color;
  int rank;
  int nprocs;

  MPI_Comm temp_comm = agg_id->comm;
  MPI_Comm_rank(temp_comm, &rank);
  MPI_Comm_size(agg_id->comm, &nprocs);
  //printf("[Before] Agg size = %d\n", nprocs);

  PIDX_variable var0 = agg_id->idx->variable[agg_id->first_index];
  if (var0->patch_group_count == 0)
    color = 0;
  else
    color = 1;

  MPI_Comm_split(temp_comm, color, rank, &(agg_id->comm));
  MPI_Comm_size(agg_id->comm, &nprocs);

  return PIDX_success;
}

PIDX_return_code PIDX_destroy_local_aggregation_comm(PIDX_agg_id agg_id)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  MPI_Comm_free(&(agg_id->comm));

  return PIDX_success;
}


PIDX_return_code PIDX_agg_set_global_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  if (agg_id == NULL)
    return PIDX_err_id;

  agg_id->global_comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout global_block_layout, PIDX_block_layout local_block_layout)
{
  PIDX_return_code ret = PIDX_success;
  if (agg_id->idx_d->agg_type == 0)
    ret = PIDX_global_agg_meta_data_create(agg_id->global_id, agg_buffer, global_block_layout);

  else if (agg_id->idx_d->agg_type == 1)
    ret = PIDX_local_agg_meta_data_create(agg_id->local_id, agg_buffer, local_block_layout);

  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;

}


PIDX_return_code PIDX_agg_buf_create(PIDX_agg_id agg_id, Agg_buffer agg_buffer, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout, int i1, int j1)
{
  PIDX_return_code ret = PIDX_success;
  int nprocs = 1;
#if PIDX_HAVE_MPI
  if (agg_id->idx_d->parallel_mode == 1)
    MPI_Comm_size(agg_id->comm,  &nprocs);
#endif

  if (agg_id->idx_d->agg_type == 0)
  {
    if (global_block_layout->existing_file_count * agg_id->idx->variable_count <= nprocs)
      ret = PIDX_global_agg_buf_create(agg_id->global_id, agg_buffer, local_block_layout, global_block_layout, i1, j1 - agg_id->idx_d->start_layout_index);
    else
      ret = PIDX_global_agg_buf_create(agg_id->global_id, agg_buffer, local_block_layout, global_block_layout, 0, j1 - agg_id->idx_d->start_layout_index);
  }

  else if (agg_id->idx_d->agg_type == 1)
  {
    if (local_block_layout->existing_file_count * agg_id->idx->variable_count <= nprocs)
      ret = PIDX_local_agg_buf_create(agg_id->local_id, agg_buffer, local_block_layout, j1 - agg_id->idx_d->start_layout_index);
    else
      ret = PIDX_local_agg_buf_create(agg_id->local_id, agg_buffer, local_block_layout, 0);
  }

  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int layout_id, PIDX_block_layout local_block_layout,  int PIDX_MODE, int vi, int bi)
{
  PIDX_return_code ret = PIDX_success;

  if (agg_id->idx_d->agg_type == 0)
    ret = PIDX_global_agg(agg_id->global_id, agg_buffer, layout_id, local_block_layout, PIDX_MODE, vi, bi);

  else if (agg_id->idx_d->agg_type == 1)
    ret = PIDX_local_agg(agg_id->local_id, agg_buffer, layout_id, local_block_layout, PIDX_MODE);

  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_destroy(PIDX_agg_id agg_id, Agg_buffer agg_buffer)
{
  PIDX_return_code ret = PIDX_success;

  if (agg_id->idx_d->agg_type == 0)
    ret = PIDX_global_agg_buf_destroy(agg_id->global_id, agg_buffer);

  if (agg_id->idx_d->agg_type == 1)
    ret = PIDX_local_agg_buf_destroy(agg_id->local_id, agg_buffer);

  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;

}


PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id agg_id, PIDX_block_layout local_block_layout, PIDX_block_layout global_block_layout)
{
  PIDX_return_code ret = PIDX_success;

  if (agg_id->idx_d->agg_type == 0)
    ret = PIDX_global_agg_meta_data_destroy(agg_id->global_id, global_block_layout);

  else if (agg_id->idx_d->agg_type == 1)
    ret = PIDX_local_agg_meta_data_destroy(agg_id->local_id, local_block_layout);

  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_rst;
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id)
{

  PIDX_global_agg_finalize(agg_id->global_id);
  PIDX_local_agg_finalize(agg_id->local_id);

  free(agg_id);
  agg_id = 0;

  return PIDX_success;
}
