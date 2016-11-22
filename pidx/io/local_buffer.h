#ifndef __PIDX_LOCAL_BUFFER_H
#define __PIDX_LOCAL_BUFFER_H

PIDX_return_code create_async_buffers(PIDX_io file, int gi);

PIDX_return_code wait_and_destroy_async_buffers(PIDX_io file, int gi);

PIDX_return_code finalize_aggregation(PIDX_io file, int start_index, int gi);

PIDX_return_code create_agg_io_buffer(PIDX_io file, int gi, int svi, int evi);

PIDX_return_code destroy_agg_io_buffer(PIDX_io file, int svi, int evi);

#endif
