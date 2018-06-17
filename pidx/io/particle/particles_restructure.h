#ifndef __PIDX_PARTICLES_RESTRUCTURE_H
#define __PIDX_PARTICLES_RESTRUCTURE_H


///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code particles_restructure_setup(PIDX_io file, int svi, int evi);



///
/// \brief restructure
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code particles_restructure(PIDX_io file);



///
/// \brief restructure_io
/// \param file
/// \param mode
/// \return
///
PIDX_return_code particles_restructure_io(PIDX_io file);



///
/// \brief restructure_cleanup
/// \param file
/// \return
///
PIDX_return_code particles_restructure_cleanup(PIDX_io file);


PIDX_return_code particles_set_rst_box_size_for_raw_write(PIDX_io file, int svi);

#endif
