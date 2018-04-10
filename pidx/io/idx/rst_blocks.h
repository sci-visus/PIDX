#ifndef __RST_BLOCKS_H
#define __RST_BLOCKS_H


PIDX_return_code delete_rst_block_layout(PIDX_io file, int gi);

PIDX_return_code populate_rst_block_layouts(PIDX_io file, int gi, int svi, int hz_from_shared, int hz_to_non_shared);

#endif
