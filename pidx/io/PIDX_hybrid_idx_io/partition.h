#ifndef __PIDX_PARTITION_IO_H
#define __PIDX_PARTITION_IO_H


///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code restructure(PIDX_hybrid_idx_io file, int gi, int svi, int evi);



///
/// \brief partition_setup
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code partition_setup(PIDX_hybrid_idx_io file, int gi, int svi);



///
/// \brief create_local_comm
/// \param file
/// \return
///
PIDX_return_code create_local_comm(PIDX_hybrid_idx_io file, int gi);


///
/// \brief destroy_local_comm
/// \param file
/// \return
///
PIDX_return_code destroy_local_comm(PIDX_hybrid_idx_io file);

///
/// \brief restructure_cleanup
/// \param file
/// \param group_index
/// \return
///
PIDX_return_code restructure_cleanup(PIDX_hybrid_idx_io file, int group_index);

#endif
