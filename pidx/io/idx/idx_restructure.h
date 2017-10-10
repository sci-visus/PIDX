#ifndef __PIDX_IDX_RESTRUCTURE_H
#define __PIDX_IDX_RESTRUCTURE_H


///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code idx_restructure_setup(PIDX_io file, int gi, int svi, int evi, int mode);



///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code idx_restructure(PIDX_io file, int mode);



///
/// \brief restructure_io
/// \param file
/// \param mode
/// \return
///
PIDX_return_code idx_restructure_io(PIDX_io file, int mode);



///
/// \brief restructure_cleanup
/// \param file
/// \return
///
PIDX_return_code idx_restructure_cleanup(PIDX_io file);



///
/// \brief restructure_forced_read
/// \param file
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code idx_restructure_forced_read(PIDX_io file, int svi, int evi);



PIDX_return_code idx_restructure_comm_create(PIDX_io file, int gi, int svi);



PIDX_return_code free_restructured_communicators(PIDX_io file, int gi);


#endif
