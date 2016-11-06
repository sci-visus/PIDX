#ifndef __AGG_IO_H
#define __AGG_IO_H

PIDX_return_code data_io(PIDX_io file, int gi, int start_index, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared, int mode );

PIDX_return_code data_aggregate(PIDX_io file, int gi, int start_index, int agg_io_level_file_zero, int agg_io_level_shared, int agg_io_level_non_shared, int agg_mode, int mode );

#endif
