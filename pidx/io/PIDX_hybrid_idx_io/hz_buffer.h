#ifndef __HZ_BUFFERS_H
#define __HZ_BUFFERS_H

PIDX_return_code setup_hz_buffers(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index);

PIDX_return_code populate_hz_buffers(PIDX_hybrid_idx_io file, int svi, int evi);

PIDX_return_code create_hz_buffers(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index);

PIDX_return_code destroy_hz_buffers(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index);

#endif
