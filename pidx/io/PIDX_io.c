/*
 * 1. One global dataset.
 * 2. Partitions in local index space.
 * 3. Partitions in global index space.
 * 4. Partition divided into non-shared, shared, file 0 (global index)
 * 5. Partition divided into non-shared, shared, file 0 (local index)
 */


#include "../PIDX_inc.h"



PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_io idx_io_id;
  idx_io_id = malloc(sizeof (*idx_io_id));
  memset(idx_io_id, 0, sizeof (*idx_io_id));

  idx_io_id->idx = idx_meta_data;
  idx_io_id->idx_d = idx_derived_ptr;
  idx_io_id->idx_dbg = idx_dbg;
  idx_io_id->idx_c = idx_c;

  return (idx_io_id);
}



PIDX_return_code PIDX_write(PIDX_io file, int gi, int svi, int evi, int MODE)
{
  PIDX_return_code ret = 0;

  file->idx_d->time->SX = PIDX_get_time();
  if (MODE == PIDX_IDX_IO)
  {
    file->idx_d->io_mode = 0;
    ret = PIDX_idx_write(file, gi, svi, evi);
  }

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
  {
    file->idx_d->io_mode = 1;
    ret = PIDX_local_partition_idx_write(file, gi, svi, evi);
  }

  else if (MODE == PIDX_GLOBAL_PARTITION_IDX_IO)
  {
    file->idx_d->io_mode = 0;
    ret = PIDX_global_partition_idx_write(file, gi, svi, evi);
  }

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_write(file, gi, svi, evi);

  else if (MODE == PIDX_MERGE_TREE_ANALYSIS)
    ret = PIDX_idx_insitu(file, gi, svi, evi);

  else if (MODE == PIDX_WAVELET_IO)
  {
    file->idx_d->io_mode = 0;
    ret = PIDX_wavelet_write(file, gi, svi, evi);
  }

  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  file->idx_d->time->EX = PIDX_get_time();

  return PIDX_success;
}



PIDX_return_code PIDX_read(PIDX_io file, int gi, int svi, int evi, int MODE)
{
  PIDX_return_code ret = 0;

  file->idx_d->time->SX = PIDX_get_time();
  if (MODE == PIDX_IDX_IO)
  {
    file->idx_d->io_mode = 0;
    ret = PIDX_idx_read(file, gi, svi, evi);
  }

  else if (MODE == PIDX_LOCAL_PARTITION_IDX_IO)
  {
    file->idx_d->io_mode = 1;
    ret = PIDX_local_partition_idx_read(file, gi, svi, evi);
  }

  else if (MODE == PIDX_GLOBAL_PARTITION_IDX_IO)
  {
    file->idx_d->io_mode = 0;
    ret = PIDX_global_partition_idx_read(file, gi, svi, evi);
  }

  else if (MODE == PIDX_RAW_IO)
    ret = PIDX_raw_read(file, gi, svi, evi);

  if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_file;}

  file->idx_d->time->EX = PIDX_get_time();

  return PIDX_success;
}




PIDX_return_code PIDX_io_finalize(PIDX_io file)
{
  free(file);
  file = 0;

  return PIDX_success;
}

