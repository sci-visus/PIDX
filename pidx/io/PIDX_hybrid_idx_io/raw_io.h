#ifndef __RAW_IO_H
#define __RAW_IO_H

PIDX_return_code PIDX_raw_write(PIDX_hybrid_idx_io file, int gi, int svi, int evi);

PIDX_return_code PIDX_raw_read(PIDX_hybrid_idx_io file, int gi, int svi, int evi);

#endif
