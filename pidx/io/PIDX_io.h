#ifndef __PIDX_IO_H
#define __PIDX_IO_H

///
/// \brief The PIDX_io_descriptor struct
///
struct PIDX_io_descriptor
{
  PIDX_header_io_id header_io_id;                   ///< IDX metadata id

  PIDX_raw_rst_id raw_rst_id;                       ///< Multi-patch restructuring phase id for raw writes

  PIDX_idx_rst_id idx_rst_id;                       ///< Multi-patch restructuring phase id for idx writes

  PIDX_chunk_id chunk_id;                           ///< Block restructuring id (prepration for compression)
  PIDX_comp_id comp_id;                             ///< Compression (lossy and lossless) id
  PIDX_hz_encode_id hz_id;                          ///< HZ encoding phase id

  PIDX_wavelet_id wavelet_id;                       ///< HZ encoding phase id

  PIDX_agg_id** agg_id;                             ///< Aggregation phase id for all shared file
  PIDX_file_io_id** io_id;                          ///< Aggregation phase id for all nonshared file


  idx_comm idx_c;


  idx_dataset idx;                                  ///< Contains all relevant IDX file info
                                                    ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;               ///< Contains all derieved IDX file info
                                                    ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;                                ///<

  idx_metadata_cache idx_cache;

  int hz_from_non_shared;
  int hz_from_shared;
  int hz_to_non_shared;
  int hz_to_shared;
};
typedef struct PIDX_io_descriptor* PIDX_io;


///
/// \brief PIDX_io_init
/// \param idx_meta_data
/// \param idx_derived_ptr
/// \param idx_c
/// \param idx_dbg
/// \return
///
PIDX_io PIDX_io_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, idx_metadata_cache cache);



///
/// \brief PIDX_write
/// \param file
/// \param group_index
/// \param start_var_index
/// \param end_var_index
/// \return
///
PIDX_return_code PIDX_write(PIDX_io file, int group_index, int start_var_index, int end_var_index, int MODE);


///
/// \brief PIDX_read
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code PIDX_read(PIDX_io file, int gi, int svi, int evi, int MODE);



///
/// \brief PIDX_io_finalize
/// \param file
/// \return
///
PIDX_return_code PIDX_io_finalize(PIDX_io file);

#endif
