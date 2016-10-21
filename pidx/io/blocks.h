#ifndef __BLOCKS_H
#define __BLOCKS_H

PIDX_return_code populate_bit_string(PIDX_io file, int mode, int io_type);

PIDX_return_code delete_block_layout(PIDX_io file, int gi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared);

PIDX_return_code populate_block_layouts(PIDX_io file, int gi, int svi, int evi, int hz_from_file_zero, int hz_to_file_zero, int hz_from_shared, int hz_to_shared, int hz_from_non_shared, int hz_to_non_shared, int io_type);

#endif
