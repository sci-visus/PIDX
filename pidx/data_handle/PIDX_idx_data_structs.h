/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

/**
 * \file PIDX_idx_data_structs.h
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 *
 * IDX related data structs:
 * -- PIDX_variable
 * -- idx_dataset
 * -- idx_dataset_derived_metadata
 */

#ifndef __PIDX_IDX_DATA_STRUCTS_H
#define __PIDX_IDX_DATA_STRUCTS_H

#include "PIDX_memory_layout_data_structs.h"


struct PIDX_timming_struct
{
  //double a1, a2, a3, a4, a5;
  double SX, EX;
  double sim_start, sim_end;

  double init_start, init_end;
  double set_reg_box_start, set_reg_box_end;
  double bit_string_start, bit_string_end;
  double layout_start, layout_end;
  double header_io_start, header_io_end;
  double group_cleanup_start, group_cleanup_end;
  double partition_start, partition_end;
  double partition_cleanup_start, partition_cleanup_end;

  double ***w_stencil_comm_x_odd_start, ***w_stencil_comm_x_odd_end;
  double ***w_stencil_comm_y_odd_start, ***w_stencil_comm_y_odd_end;
  double ***w_stencil_comm_z_odd_start, ***w_stencil_comm_z_odd_end;

  double ***w_stencil_comm_x_even_start, ***w_stencil_comm_x_even_end;
  double ***w_stencil_comm_y_even_start, ***w_stencil_comm_y_even_end;
  double ***w_stencil_comm_z_even_start, ***w_stencil_comm_z_even_end;

  double ***w_stencil_comp_x_odd_start, ***w_stencil_comp_x_odd_end;
  double ***w_stencil_comp_y_odd_start, ***w_stencil_comp_y_odd_end;
  double ***w_stencil_comp_z_odd_start, ***w_stencil_comp_z_odd_end;

  double ***w_stencil_comp_x_even_start, ***w_stencil_comp_x_even_end;
  double ***w_stencil_comp_y_even_start, ***w_stencil_comp_y_even_end;
  double ***w_stencil_comp_z_even_start, ***w_stencil_comp_z_even_end;

  double ***w_rst_comp_x_start, ***w_rst_comp_x_end;
  double ***w_rst_comp_y_start, ***w_rst_comp_y_end;
  double ***w_rst_comp_z_start, ***w_rst_comp_z_end;

  double **rst_init_start, **rst_init_end;
  double **rst_meta_data_create_start, **rst_meta_data_create_end;
  double **rst_meta_data_io_start, **rst_meta_data_io_end;
  double **rst_buffer_start, **rst_buffer_end;
  double **rst_write_read_start, **rst_write_read_end;
  double **rst_buff_agg_start, **rst_buff_agg_end;
  double **rst_buff_agg_free_start, **rst_buff_agg_free_end;
  double **rst_buff_agg_io_start, **rst_buff_agg_io_end;
  double **rst_cleanup_start, **rst_cleanup_end;

  double **hz_init_start, **hz_init_end;
  double **hz_meta_start, **hz_meta_end;
  double **hz_buffer_start, **hz_buffer_end;
  double **hz_start, **hz_end;
  double **hz_compress_start, **hz_compress_end;
  double **hz_buffer_free_start, **hz_buffer_free_end;
  double **hz_cleanup_start, **hz_cleanup_end;
  double ***hz_io_start, ***hz_io_end;

  double **chunk_init_start, **chunk_init_end;
  double **chunk_meta_start, **chunk_meta_end;
  double **chunk_buffer_start, **chunk_buffer_end;
  double **chunk_start, **chunk_end;
  double **chunk_buffer_free_start, **chunk_buffer_free_end;
  double **chunk_cleanup_start, **chunk_cleanup_end;

  double **compression_init_start, **compression_init_end;
  double **compression_start, **compression_end;

  double ***agg_init_start, ***agg_init_end;
  double ***agg_meta_start, ***agg_meta_end;
  double ***agg_buf_start, ***agg_buf_end;
  double ***agg_start, ***agg_end;
  double ***agg_compress_start, ***agg_compress_end;
  double ***agg_meta_cleanup_start, ***agg_meta_cleanup_end;

  double **io_start, **io_end;
};
typedef struct PIDX_timming_struct* PIDX_time;


struct PIDX_variable_struct
{
  // General Info
  char var_name[1024];                                       ///< Variable name
  int vps;                                                   ///< values per sample, Vector(3), scalar(1), or n
  int bpv;                                                   ///< Number of bits each need
  PIDX_data_type type_name;                                  ///< Name of the type uint8, bob
  PIDX_data_layout data_layout;                              ///< Row major or column major

  // buffer (before, after HZ encoding phase)
  int sim_patch_count;                                       ///< The actual number of patches (application layout), most probably more than 1 in uintah
  PIDX_patch sim_patch[1024];                                ///< Pointer to the patches
  HZ_buffer hz_buffer;                                      ///< HZ encoded buffer of the patches

  // buffer before aggregation
  int patch_group_count;                                     ///< Number of groups of patches to be passed to aggregation phase
  PIDX_super_patch rst_patch_group;                         ///< Pointer to the patch groups
  PIDX_super_patch chunk_patch_group;                       ///< Pointer to the patch group after block restructuring

  //PIDX_patch rst_wavelet_patch;
};
typedef struct PIDX_variable_struct* PIDX_variable;



struct PIDX_variable_group_struct
{
  int variable_index_tracker;
  int variable_count;
  int variable_tracker[256];
  PIDX_variable variable[256];

  int local_variable_index;
  int local_variable_count;

  //int agg_l_shared;
  int shared_start_layout_index;
  int shared_end_layout_index;
  int shared_layout_count;

  //int agg_l_nshared;
  int nshared_start_layout_index;
  int nshared_end_layout_index;
  int nshared_layout_count;

  int *rank_buffer;

  int agg_level;
  PIDX_block_layout block_layout;
  PIDX_block_layout* block_layout_by_level;

  // Block level layout
  //PIDX_block_layout shared_block_layout;
  //PIDX_block_layout* shared_block_layout_by_level;

  //PIDX_block_layout nshared_block_layout;
  //PIDX_block_layout* nshared_block_layout_by_level;
};
typedef struct PIDX_variable_group_struct* PIDX_variable_group;


/// Communicator related struct
struct idx_comm_struct
{
  int lrank;
  int lnprocs;

  int grank;
  int grank_x;
  int grank_y;
  int grank_z;

  int gnprocs;
  int gnproc_x;
  int gnproc_y;
  int gnproc_z;

  /// Names
  MPI_Comm global_comm;
  MPI_Comm local_comm;
  MPI_Comm rst_comm;
};
typedef struct idx_comm_struct* idx_comm;


/// idx_file
struct idx_file_struct
{
  int io_type;
  int current_time_step;

  int variable_pipe_length;
  int variable_count;
  int variable_group_count;
  int group_index_tracker;
  PIDX_variable_group variable_grp[16];
  
  char agg_list_filename[1024];

  char filename[1024];
  char filename_global[1024];
  char filename_partition[1024];
  char filename_template_global[1024];
  char filename_template_partition[1024];
  char filename_template[1024];


  int first_tstep;
  int last_tstep;

  int bits_per_block;
  int blocks_per_file;
  unsigned long long bounds[PIDX_MAX_DIMENSIONS];
  unsigned long long box_bounds[PIDX_MAX_DIMENSIONS];
  double transform[16];
  char bitSequence[512];
  char bitPattern[512];

  double zfp_precisison;
  int compression_type;
  int enable_rst;

  /// 0 No aggregation
  /// 1 Only aggregation
  int enable_agg;

  int compression_start_level;
  int compression_factor;
  float compression_bit_rate;
  unsigned long long chunk_size[PIDX_MAX_DIMENSIONS];

  int file_zero_merge;

  /// 1 for little endian
  /// 0 for big endian
  int endian;

  /// 1 for flipping endian
  /// 0 for big endian
  int flip_endian;

  int agg_counter;

  int cached_ts;

  //unsigned long long* all_offset;
  //unsigned long long* all_size;

  //int random_agg_counter;
  //int *random_agg_list;
};
typedef struct idx_file_struct* idx_dataset;


/// idx_dataset_derived_metadata
struct idx_dataset_derived_metadata_struct
{
  int pidx_version;
  int io_mode;

  int w_nx;
  int w_px;
  int w_ny;
  int w_py;
  int w_nz;
  int w_pz;
  int wavelet_levels;
  int wavelet_imeplementation_type;

  PIDX_restructured_grid restructured_grid;

  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
  
  int fs_block_size;
  off_t start_fs_block;

  Agg_buffer **agg_buffer;

  int color;
  int partition_count[PIDX_MAX_DIMENSIONS];
  int partition_size[PIDX_MAX_DIMENSIONS];
  int partition_offset[PIDX_MAX_DIMENSIONS];

  int var_pipe_length;
  
  int start_layout_index;
  int end_layout_index;

  MPI_Status *status1;
  MPI_File *fp1;
  MPI_Request *request1;

  int layout_count;
  int reduced_res_from;
  int reduced_res_to;

  PIDX_time time;

  int raw_io_pipe_length;

  int aggregator_multiplier;
  int data_core_count;

  //int shared_block_level;
  int total_partiton_level;

  int **block_bitmap;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;



struct idx_debug_info_struct
{
  int simulate_rst_io;
  int simulate_rst;

  //FILE *rst_dump_fp;
  //int dump_rst_info;
  //char rst_dump_dir_name[512];

  //FILE *agg_dump_fp;
  //int dump_agg_info;
  //char agg_dump_dir_name[512];

  //FILE *io_dump_fp;
  //int dump_io_info;
  //char io_dump_dir_name[512];

  //FILE *process_size_and_offset_dump_fp;
  //int dump_process_state;
  //char process_state_dump_dir_name[512];

  FILE *local_dump_fp;
  FILE *mpi_dump_fp;
  int state_dump;

  int debug_rst;                               ///< Debug restructuring phase, works only on the test application
  int debug_hz;                                ///< Debug HZ encoding phase, works only on the test application

  /// Flags set by user
  int debug_do_rst;                            ///< User controlled flag to activate/deactivate restructuring phase
  int debug_do_chunk;                          ///< User controlled flag to activate/deactivate chunking phase
  int debug_do_compress;                       ///< User controlled flag to activate/deactivate compression
  int debug_do_hz;                             ///< User controlled flag to activate/deactivate hz encoding phase
  int debug_do_agg;                            ///< User controlled flag to activate/deactivate aggregation phase
  int debug_do_io;                             ///< User controlled flag to activate/deactivate I/O phase
};
typedef struct idx_debug_info_struct* idx_debug;

#endif
