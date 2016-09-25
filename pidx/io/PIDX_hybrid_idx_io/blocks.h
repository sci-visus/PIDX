#ifndef __BLOCKS_H
#define __BLOCKS_H

PIDX_return_code populate_bit_string(PIDX_hybrid_idx_io file);

PIDX_return_code populate_idx_dataset(PIDX_hybrid_idx_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int* start_layout_index, int* end_layout_index, int* layout_count, int group_index, int start_index, int end_index, int hz_level_from, int hz_level_to);

#endif
