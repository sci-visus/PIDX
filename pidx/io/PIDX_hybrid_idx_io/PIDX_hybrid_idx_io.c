/*
 * 1. One global dataset.
 * 2. Partitions in local index space.
 * 3. Partitions in global index space.
 * 4. Partition divided into non-shared, shared, file 0 (global index)
 * 5. Partition divided into non-shared, shared, file 0 (local index)
 */


#include "../PIDX_io.h"
#include "idx_io.h"
#include "local_partition_idx_io.h"
#include "global_partition_idx_io.h"
#include "raw_io.h"


PIDX_hybrid_idx_io PIDX_hybrid_idx_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_hybrid_idx_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_d = idx_derived_ptr;
  idx_io_id->idx_dbg = idx_dbg;

  return (idx_io_id);
}


PIDX_return_code PIDX_hybrid_idx_io_set_communicator(PIDX_hybrid_idx_io id, MPI_Comm comm)
{
  if (id == NULL)
    return PIDX_err_id;

  id->global_comm = comm;
  id->comm = comm;

  return PIDX_success;
}



PIDX_return_code PIDX_hybrid_idx_write(PIDX_hybrid_idx_io file, int gi, int svi, int evi, int MODE)
{
  PIDX_time time = file->idx_d->time;
  time->SX = PIDX_get_time();

  PIDX_return_code ret;

  if (MODE == PIDX_IDX_IO)
    ret = PIDX_idx_write(file, gi, svi, evi);

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
    ret = PIDX_local_partition_idx_write(file, gi, svi, evi);

  else if (MODE == PIDX_GLOBAL_PARTITION_IDX_IO)
    ret = PIDX_global_partition_idx_write(file, gi, svi, evi);

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_write(file, gi, svi, evi);

  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  time->EX = PIDX_get_time();
  return PIDX_success;
}



PIDX_return_code PIDX_hybrid_idx_read(PIDX_hybrid_idx_io file, int gi, int svi, int evi, int MODE)
{

  PIDX_time time = file->idx_d->time;
  time->SX = PIDX_get_time();

  PIDX_return_code ret;

  if (MODE == PIDX_IDX_IO)
    ret = PIDX_idx_read(file, gi, svi, evi);

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
    ret = PIDX_local_partition_idx_read(file, gi, svi, evi);

  else if (MODE == PIDX_GLOBAL_PARTITION_IDX_IO)
    ret = PIDX_global_partition_idx_read(file, gi, svi, evi);

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_read(file, gi, svi, evi);

  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}


  time->EX = PIDX_get_time();

  return PIDX_success;
}




PIDX_return_code PIDX_hybrid_idx_io_finalize(PIDX_hybrid_idx_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}

