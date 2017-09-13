#include "../PIDX_inc.h"


PIDX_return_code idx_init(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  var_grp->rank_buffer = malloc(file->idx_c->gnprocs * sizeof(*var_grp->rank_buffer));
  memset(var_grp->rank_buffer, 0, file->idx_c->gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&(file->idx_c->grank), 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->idx_c->global_comm);

#if 1
  //PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  file->idx->all_offset = malloc(sizeof(*(file->idx->all_offset)) * PIDX_MAX_DIMENSIONS * file->idx_c->gnprocs);
  file->idx->all_size = malloc(sizeof(*(file->idx->all_size)) * PIDX_MAX_DIMENSIONS * file->idx_c->gnprocs);

  MPI_Allgather(var0->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx->all_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->global_comm);
  MPI_Allgather(var0->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx->all_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->global_comm);
#endif
  return PIDX_success;
}


PIDX_return_code idx_finalize(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  free(var_grp->rank_buffer);

  free(file->idx->all_offset);
  free(file->idx->all_size);

  return PIDX_success;
}
