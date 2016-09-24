#ifndef __AGG_IO_H
#define __AGG_IO_H

static PIDX_return_code PIDX_global_aggregate(PIDX_hybrid_idx_io file, PIDX_agg_id** agg_id, Agg_buffer** agg_buffer, PIDX_block_layout* block_layout_by_level, PIDX_block_layout global_block_layout_files, MPI_Comm comm, int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level, int agg_factor, int file_status);

static PIDX_return_code PIDX_global_async_io(PIDX_hybrid_idx_io file, PIDX_file_io_id **io_id, Agg_buffer **agg_buffer, PIDX_block_layout* block_layout_by_level,  MPI_File *fp, MPI_Request *request,  int init_index, int var_index, int index, int layout_start, int layout_end, int layout_count, int agg_io_level, int file_zero, int async_status);


#endif
