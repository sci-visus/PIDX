#ifndef __PIDX_PARTITION_IO_H
#define __PIDX_PARTITION_IO_H


///
/// \brief partition
/// \param file
/// \param gi
/// \param svi
/// \param mode
/// \return
///
//PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode);


///
/// \brief partition_setup
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code partition_setup(PIDX_io file, int gi, int svi);



///
/// \brief create_local_comm
/// \param file
/// \return
///
PIDX_return_code create_local_comm(PIDX_io file);



///
/// \brief destroy_local_comm
/// \param file
/// \return
///
PIDX_return_code destroy_local_comm(PIDX_io file);



///
/// \brief find_partition_count
/// \param file
/// \return
///
PIDX_return_code find_partition_count(PIDX_io file);

#endif
