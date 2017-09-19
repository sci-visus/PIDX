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
PIDX_return_code restructure_setup(PIDX_io file, int gi, int svi, int evi, int mode);



///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code restructure(PIDX_io file, int mode);



///
/// \brief restructure_io
/// \param file
/// \param mode
/// \return
///
PIDX_return_code restructure_io(PIDX_io file, int mode);



///
/// \brief restructure_cleanup
/// \param file
/// \return
///
PIDX_return_code restructure_cleanup(PIDX_io file);



///
/// \brief restructure_forced_read
/// \param file
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code restructure_forced_read(PIDX_io file, int svi, int evi);



///
/// \brief create_restructured_communictors
/// \param file
/// \param gi
/// \param svi
/// \return
///
PIDX_return_code create_restructured_communictors(PIDX_io file, int gi, int svi);


PIDX_return_code create_restructured_communicators(PIDX_io file, int gi, int svi);


PIDX_return_code free_rst_box(PIDX_io file);


PIDX_return_code free_restructured_communicators(PIDX_io file, int gi);


#endif
