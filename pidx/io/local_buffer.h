#ifndef __PIDX_LOCAL_BUFFER_H
#define __PIDX_LOCAL_BUFFER_H

PIDX_return_code create_async_buffers(PIDX_io file, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared);

PIDX_return_code wait_and_destroy_async_buffers(PIDX_io file, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared);

PIDX_return_code finalize_aggregation(PIDX_io file, int start_index, int gi, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared);

PIDX_return_code create_agg_io_buffer(PIDX_io file);

PIDX_return_code destroy_agg_io_buffer(PIDX_io file);

#endif
