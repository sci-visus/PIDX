#ifndef __PIDX_MULTIPATCH_RAW_IO_H
#define __PIDX_MULTIPATCH_RAW_IO_H

struct PIDX_multi_patch_raw_io_descriptor;
typedef struct PIDX_multi_patch_raw_io_descriptor* PIDX_multi_patch_raw_io;


///
PIDX_multi_patch_raw_io PIDX_multi_patch_raw_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg);


#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id restructuring id
/// \param comm the communicator
/// \return error code
PIDX_return_code PIDX_multi_patch_raw_io_set_communicator(PIDX_multi_patch_raw_io id, MPI_Comm comm);
#endif

PIDX_return_code PIDX_multi_patch_rst_staged_write(PIDX_multi_patch_rst_id rst_id);

///
PIDX_return_code PIDX_multi_patch_raw_write(PIDX_multi_patch_raw_io file, int start_var_index, int end_var_index);


///
PIDX_return_code PIDX_multi_patch_raw_read(PIDX_multi_patch_raw_io file, int start_var_index, int end_var_index);



///
PIDX_return_code PIDX_multi_patch_forced_raw_read(PIDX_multi_patch_raw_io file, int start_var_index, int end_var_index);


///
PIDX_return_code PIDX_multi_patch_raw_io_finalize(PIDX_multi_patch_raw_io file);

#endif
