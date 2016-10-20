#ifndef __TIMMING_H
#define __TIMMING_H


void PIDX_init_timming_buffers1(PIDX_time time, int variable_count);


void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count);


void PIDX_delete_timming_buffers1(PIDX_time time);


void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count);


#endif
