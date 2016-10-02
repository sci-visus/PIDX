#ifndef __BLOCKS_H
#define __BLOCKS_H

PIDX_return_code populate_bit_string(PIDX_hybrid_idx_io file, int mode);

PIDX_return_code delete_block_layout(PIDX_hybrid_idx_io file, int gi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

PIDX_return_code populate_block_layouts(PIDX_hybrid_idx_io file, int gi, int svi, int evi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

#endif
