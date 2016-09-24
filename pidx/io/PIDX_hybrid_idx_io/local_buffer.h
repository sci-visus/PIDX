#ifndef __PIDX_LOCAL_BUFFER_H
#define __PIDX_LOCAL_BUFFER_H

PIDX_return_code create_non_shared_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_non_shared, int agg_io_level_non_shared);

PIDX_return_code create_shared_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_shared, int agg_io_level_shared);

PIDX_return_code create_file_zero_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_file_zero, int agg_io_level_file_zero);

PIDX_return_code wait_and_destroy_non_shared_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_non_shared, int agg_io_level_non_shared);

PIDX_return_code wait_and_destroy_shared_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_shared, int agg_io_level_shared);

PIDX_return_code wait_and_destroy_file_zero_async_buffers(PIDX_hybrid_idx_io file, int start_layout_index_file_zero, int agg_io_level_file_zero);

PIDX_return_code destroy_non_shared_ids_and_buffers(PIDX_hybrid_idx_io file, int start_index, int start_layout_index_non_shared, int end_layout_index_non_shared, int agg_io_level_non_shared);

PIDX_return_code destroy_shared_ids_and_buffers(PIDX_hybrid_idx_io file, int start_index, int start_layout_index_shared, int end_layout_index_shared, int agg_io_level_shared);

PIDX_return_code destroy_file_zero_ids_and_buffers(PIDX_hybrid_idx_io file, int start_index, int start_layout_index_file_zero, int end_layout_index_file_zero, int agg_io_level_file_zero);

PIDX_return_code init_agg_io_buffer(PIDX_hybrid_idx_io file, int group_index);

PIDX_return_code destroy_agg_io_buffer(PIDX_hybrid_idx_io file, int group_index);

#endif
