#include "../../PIDX_inc.h"


PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, int init_index, int first_index, int last_index)
{
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  agg_id->idx = idx_meta_data;
  agg_id->idx_d = idx_d;

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

PIDX_return_code PIDX_agg_finalize(PIDX_agg_id agg_id)
{

  free(agg_id);
  agg_id = 0;

  return PIDX_success;
}
