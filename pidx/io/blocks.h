#ifndef __BLOCKS_H
#define __BLOCKS_H

PIDX_return_code populate_global_bit_string(PIDX_io file, int mode);

PIDX_return_code populate_local_bit_string(PIDX_io file, int mode);

PIDX_return_code delete_block_layout(PIDX_io file, int gi, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

PIDX_return_code populate_block_layouts(PIDX_io file, int gi, int svi, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

#endif
