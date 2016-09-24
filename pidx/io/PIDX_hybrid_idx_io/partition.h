#ifndef __PIDX_PARTITION_IO_H
#define __PIDX_PARTITION_IO_H

PIDX_return_code partition_domain(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index);

PIDX_return_code partition_communicator(PIDX_hybrid_idx_io file);
PIDX_return_code partition_communicator_destroy(PIDX_hybrid_idx_io file);
PIDX_return_code partition_destroy(PIDX_hybrid_idx_io file, int group_index);

#endif
