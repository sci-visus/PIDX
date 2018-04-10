#ifndef __SERIAL_IDX_IO_H
#define __SERIAL_IDX_IO_H

///
/// \brief PIDX_idx_write
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code PIDX_serial_idx_write(PIDX_io file, int gi, int svi, int evi);



///
/// \brief PIDX_parallel_local_partition_idx_read
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code PIDX_parallel_local_partition_idx_read(PIDX_io file, int gi, int svi, int evi);

#endif
