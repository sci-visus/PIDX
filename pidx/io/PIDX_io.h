#include "../PIDX_inc.h"
#include "./PIDX_idx_io/PIDX_idx_io.h"
#include "./PIDX_raw_io/PIDX_raw_io.h"
#include "./PIDX_partitioned_idx_io/PIDX_partitioned_idx_io.h"
#include "./PIDX_partition_merge_idx_io/PIDX_partition_merge_idx_io.h"

#ifndef __PIDX_IO_H
#define __PIDX_IO_H

#define PIDX_RAW_IO                              1
#define PIDX_IDX_IO                              2
#define PIDX_PARTITIONED_IDX_IO                  3
#define PIDX_PARTITION_MERGE_IDX_IO              4

struct PIDX_io_descriptor;
typedef struct PIDX_io_descriptor* PIDX_io;

///
/// \brief PIDX_get_time
/// \return
///
double PIDX_get_time();


///
/// \brief PIDX_init_timming_buffers1
/// \param time
/// \param variable_count
///
void PIDX_init_timming_buffers1(PIDX_time time, int variable_count);


///
/// \brief PIDX_init_timming_buffers2
/// \param time
/// \param variable_count
/// \param layout_count
///
void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count);



///
/// \brief PIDX_delete_timming_buffers1
/// \param time
///
void PIDX_delete_timming_buffers1(PIDX_time time);


///
/// \brief PIDX_delete_timming_buffers2
/// \param time
/// \param variable_count
///
void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count);



///
/// \brief PIDX_io_init
/// \param idx_meta_data
/// \param idx_derived_ptr
/// \param idx_dbg
/// \return
///
PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg);


///
/// \brief PIDX_io
/// \param mode
/// \return
///
PIDX_return_code PIDX_io_io(PIDX_io file, int mode, int io_type, int start_var_index, int end_var_index);



#if PIDX_HAVE_MPI
PIDX_return_code PIDX_io_set_communicator(PIDX_io file, MPI_Comm comm);
#endif



///
/// \brief PIDX_io_finalize
/// \param file
/// \return
///
PIDX_return_code PIDX_io_finalize(PIDX_io file);



///
//PIDX_return_code PIDX_parameter_validate(idx_dataset file, int start_var_index, int end_var_index);

#endif
