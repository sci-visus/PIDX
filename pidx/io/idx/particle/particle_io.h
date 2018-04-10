#ifndef __PARTICLE_IO_H
#define __PARTICLE_IO_H

///
/// \brief PIDX_particle_write
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code PIDX_particle_write(PIDX_io file, int gi, int svi, int evi);


///
/// \brief PIDX_particle_read
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code PIDX_particle_read(PIDX_io file, int gi, int svi, int evi);


#endif
