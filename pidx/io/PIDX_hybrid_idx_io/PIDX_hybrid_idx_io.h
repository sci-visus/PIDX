#ifndef __PIDX_HYBRID_IDX_IO_H
#define __PIDX_HYBRID_IDX_IO_H


struct PIDX_hybrid_idx_io_descriptor
{

#if PIDX_HAVE_MPI
  MPI_Comm global_comm;                               ///< MPI sub-communicator (including all processes per IDX file)
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
#endif

  PIDX_header_io_id header_io_id;              ///< IDX metadata id
  PIDX_rst_id rst_id;                          ///< Restructuring phase id
  PIDX_chunk_id chunk_id;              ///< Block restructuring id (prepration for compression)
  PIDX_comp_id comp_id;          ///< Compression (lossy and lossless) id
  PIDX_hz_encode_id hz_id;                     ///< HZ encoding phase id

  PIDX_agg_id** f0_agg_id;                          ///< Aggregation phase id
  PIDX_agg_id** shared_agg_id;                          ///< Aggregation phase id
  PIDX_agg_id** nshared_agg_id;                          ///< Aggregation phase id

  PIDX_file_io_id** f0_io_id;                          ///< Aggregation phase id
  PIDX_file_io_id** shared_io_id;                          ///< Aggregation phase id
  PIDX_file_io_id** nshared_io_id;                          ///< Aggregation phase id

  int one_time_initializations;                ///<

  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated
  idx_debug idx_dbg;
};
typedef struct PIDX_hybrid_idx_io_descriptor* PIDX_hybrid_idx_io;


///
PIDX_hybrid_idx_io PIDX_hybrid_idx_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg);


#if PIDX_HAVE_MPI
/// Attach the communicator wit the ID.
/// \param id restructuring id
/// \param comm the communicator
/// \return error code
PIDX_return_code PIDX_hybrid_idx_io_set_communicator(PIDX_hybrid_idx_io id, MPI_Comm comm);
#endif


///
PIDX_return_code PIDX_hybrid_idx_write(PIDX_hybrid_idx_io file, int group_index, int start_var_index, int end_var_index);


///
//PIDX_return_code PIDX_idx_read(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index);


///
PIDX_return_code PIDX_hybrid_idx_io_finalize(PIDX_hybrid_idx_io file);

#endif
