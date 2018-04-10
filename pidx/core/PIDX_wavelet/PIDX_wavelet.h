#ifndef __PIDX_WAVELET_IDX_STENCIL_H
#define __PIDX_WAVELET_IDX_STENCIL_H


//Struct for restructuring ID
struct PIDX_wavelet_struct
{
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;

  idx_debug idx_dbg;


  idx_comm idx_c;

  int first_index;
  int last_index;
  int group_index;
};
typedef struct PIDX_wavelet_struct* PIDX_wavelet_id;

PIDX_wavelet_id PIDX_wavelet_init( idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_comm idx_c, idx_debug idx_dbg, int var_start_index, int var_end_index);


///
/// \brief wavelet
/// \param file
/// \param gi
/// \param svi
/// \param evi
/// \return
///
PIDX_return_code idx_stencil_wavelet(PIDX_wavelet_id file, int gi, int svi, int evi, int mode);



PIDX_return_code PIDX_wavelet_finalize(PIDX_wavelet_id id);

#endif

