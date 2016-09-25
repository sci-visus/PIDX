#ifndef __PIDX_PARTITION_IO_H
#define __PIDX_PARTITION_IO_H


///
/// \brief partition
/// \param file
/// \param group_index
/// \param start_var_index
/// \param end_var_index
/// \return
///
PIDX_return_code partition(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index);


///
/// \brief create_local_comm
/// \param file
/// \return
///
PIDX_return_code create_local_comm(PIDX_hybrid_idx_io file);


///
/// \brief destroy_local_comm
/// \param file
/// \return
///
PIDX_return_code destroy_local_comm(PIDX_hybrid_idx_io file);

///
/// \brief partition_cleanup
/// \param file
/// \param group_index
/// \return
///
PIDX_return_code partition_cleanup(PIDX_hybrid_idx_io file, int group_index);

#endif
