#ifndef __HZ_BUFFERS_H
#define __HZ_BUFFERS_H


PIDX_return_code hz_encode_setup(PIDX_io file, int group_index, int start_var_index, int end_var_index);

PIDX_return_code hz_encode(PIDX_io file, int mode);

PIDX_return_code hz_io(PIDX_io file, int gi, int mode);

PIDX_return_code hz_encode_cleanup(PIDX_io file);

#endif
