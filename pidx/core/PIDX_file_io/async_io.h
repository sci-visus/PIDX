#ifndef __PIDX_ASYNC_IO_H
#define __PIDX_ASYNC_IO_H


PIDX_return_code PIDX_async_aggregated_write(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template);


PIDX_return_code PIDX_async_aggregated_read(PIDX_file_io_id io_id, Agg_buffer agg_buf, PIDX_block_layout block_layout, MPI_Request* request, MPI_File* fh, char* filename_template);

#endif
