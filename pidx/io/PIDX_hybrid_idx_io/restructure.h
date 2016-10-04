#ifndef __PIDX_RESTRUCTURE_H
#define __PIDX_RESTRUCTURE_H


///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code restructure_init(PIDX_hybrid_idx_io file, int gi, int svi, int evi);



///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code restructure(PIDX_hybrid_idx_io file, int mode);


///
/// \brief restructure_cleanup
/// \param file
/// \param group_index
/// \return
///
PIDX_return_code restructure_cleanup(PIDX_hybrid_idx_io file, int group_index);

#endif
