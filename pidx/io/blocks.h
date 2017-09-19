#ifndef __BLOCKS_H
#define __BLOCKS_H


PIDX_return_code delete_block_layout(PIDX_io file, int gi);

PIDX_return_code populate_block_layouts(PIDX_io file, int gi, int svi, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

#endif
