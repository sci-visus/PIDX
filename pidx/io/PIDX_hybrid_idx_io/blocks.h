#ifndef __BLOCKS_H
#define __BLOCKS_H

PIDX_return_code create_file_zero_block_layout(PIDX_hybrid_idx_io file, int gi, int hz_level_from, int hz_level_to);

PIDX_return_code create_shared_block_layout(PIDX_hybrid_idx_io file, int gi, int hz_level_from, int hz_level_to);

PIDX_return_code create_non_shared_block_layout(PIDX_hybrid_idx_io file, int gi, int hz_level_from, int hz_level_to);

PIDX_return_code populate_bit_string(PIDX_hybrid_idx_io file);

PIDX_return_code populate_idx_block_layout(PIDX_hybrid_idx_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int group_index, int start_index, int end_index, int hz_level_from, int hz_level_to);

PIDX_return_code destroy_file_zero_block_layout(PIDX_hybrid_idx_io file, int gi);
PIDX_return_code destroy_shared_block_layout(PIDX_hybrid_idx_io file, int gi);
PIDX_return_code destroy_non_shared_block_layout(PIDX_hybrid_idx_io file, int gi);

#endif
