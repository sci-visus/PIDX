#ifndef __HEADERS_H
#define __HEADERS_H


PIDX_return_code write_headers(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index, int layout_type);

PIDX_return_code one_time_initialize(PIDX_hybrid_idx_io file, int mode);

#endif
