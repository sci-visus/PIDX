#ifndef __TIMMING_H
#define __TIMMING_H

///
/// \brief PIDX_init_timming_buffers1
/// \param time
/// \param group_count
/// \param variable_count
/// \param layout_count
///
void PIDX_init_timming_buffers1(PIDX_time time, int group_count, int variable_count, int layout_count);


///
/// \brief PIDX_delete_timming_buffers1
/// \param time
/// \param group_count
/// \param variable_count
///
void PIDX_delete_timming_buffers1(PIDX_time time, int group_count, int variable_count);

#endif
