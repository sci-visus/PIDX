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
   double SX, EX;
   double sim_start, sim_end;

   double file_create_time;

   double init_start, init_end;
   double set_reg_box_start, set_reg_box_end;
   double bit_string_start, bit_string_end;
   double layout_start, layout_end;
   double header_io_start, header_io_end;
   double group_cleanup_start, group_cleanup_end;
   double partition_start, partition_end;
   double partition_cleanup_start, partition_cleanup_end;

   double *rst_init_start, *rst_init_end;
   double *rst_meta_data_create_start, *rst_meta_data_create_end;
   double *rst_meta_data_io_start, *rst_meta_data_io_end;
   double *rst_buffer_start, *rst_buffer_end;
   double *rst_write_read_start, *rst_write_read_end;
   double *rst_buff_agg_start, *rst_buff_agg_end;
   double *rst_buff_agg_free_start, *rst_buff_agg_free_end;
   double *rst_buff_agg_io_start, *rst_buff_agg_io_end;
   double *rst_cleanup_start, *rst_cleanup_end;

   double *hz_init_start, *hz_init_end;
   double *hz_meta_start, *hz_meta_end;
   double *hz_buffer_start, *hz_buffer_end;
   double *hz_start, *hz_end;
   double *hz_buffer_free_start, *hz_buffer_free_end;
   double *hz_cleanup_start, *hz_cleanup_end;

   double *chunk_init_start, *chunk_init_end;
   double *chunk_meta_start, *chunk_meta_end;
   double *chunk_buffer_start, *chunk_buffer_end;
   double *chunk_start, *chunk_end;
   double *chunk_buffer_free_start, *chunk_buffer_free_end;
   double *chunk_cleanup_start, *chunk_cleanup_end;

   double *compression_init_start, *compression_init_end;
   double *compression_start, *compression_end;

   double **agg_init_start, **agg_init_end;
   double **agg_meta_start, **agg_meta_end;
   double **agg_buf_start, **agg_buf_end;
   double **agg_start, **agg_end;
   double **agg_meta_cleanup_start, **agg_meta_cleanup_end;

   double *io_start, *io_end;
};
typedef struct PIDX_timming_struct* PIDX_time;


///
struct PIDX_variable_struct
{
  // General Info
  char var_name[1024];                                                  ///< Variable name
  int vps;                                                              ///< values per sample, Vector(3), scalar(1), or n
  int bpv;                                                   ///< Number of bits each need
  PIDX_data_type type_name;                                                 ///< Name of the type uint8, bob
  PIDX_data_layout data_layout;                                         ///< Row major or column major

  // buffer (before, after HZ encoding phase)
  int sim_patch_count;                                                  ///< The actual number of patches (application layout), most probably more than 1 in uintah
  Ndim_patch sim_patch[1024];                                           ///< Pointer to the patches
  HZ_buffer* hz_buffer;                                                 ///< HZ encoded buffer of the patches

  // buffer before aggregation
  int patch_group_count;                                                ///< Number of groups of patches to be passed to aggregation phase
  Ndim_patch_group* rst_patch_group;                                    ///< Pointer to the patch groups
  Ndim_patch_group* chunk_patch_group;                                  ///< Pointer to the patch group after block restructuring
};
typedef struct PIDX_variable_struct* PIDX_variable;




struct PIDX_variable_group_struct
{
  int variable_index_tracker;
  int variable_count;
  int variable_tracker[128];
  PIDX_variable variable[128];

  int local_variable_index;
  int local_variable_count;

  int f0_start_layout_index;
  int f0_end_layout_index;
  int f0_layout_count;

  int shared_start_layout_index;
  int shared_end_layout_index;
  int shared_layout_count;

  int nshared_start_layout_index;
  int nshared_end_layout_index;
  int nshared_layout_count;

  //extents of meta-data
  int *rank_buffer;

  // Block level layout
  PIDX_block_layout f0_block_layout;
  PIDX_block_layout* f0_block_layout_by_level;

  PIDX_block_layout shared_block_layout;
  PIDX_block_layout* shared_block_layout_by_level;

  PIDX_block_layout nshared_block_layout;
  PIDX_block_layout* nshared_block_layout_by_level;
};
typedef struct PIDX_variable_group_struct* PIDX_variable_group;


/// idx_file
struct idx_file_struct
{
  int io_type;
  int current_time_step;                                                ///< Time step tracker

  int variable_pipe_length;
  int variable_count;
  int variable_group_count;
  int group_index_tracker;
  PIDX_variable_group variable_grp[16];
  
  char filename[1024];
  char filename_global[1024];
  char filename_partition[1024];
  char filename_file_zero[1024];
  char filename_template_global[1024];
  char filename_template_partition[1024];
  char filename_template[1024];
  char filename_template_file_zero[1024];

  int first_tstep;                                                     ///< The first timestep index in the IDX file
  int last_tstep;                                                      ///< The last timestep index in the IDX file

  int bits_per_block;
  int blocks_per_file;
  unsigned long long bounds[PIDX_MAX_DIMENSIONS];
  double transform[16];
  char bitSequence[512];
  char bitPattern[512];
  
  int reg_box_set;
  unsigned long long reg_patch_size[PIDX_MAX_DIMENSIONS];
  
  int compression_type;                                               ///< counter to enable/disable (1/0) compression
  int enable_rst;                                               ///< counter to enable/disable (1/0) compression

  /// 0 No aggregation
  /// 1 Hybrid aggregation
  /// 2 Only aggregation
  int enable_agg;                                               ///< counter to enable/disable (1/0) compression

  int compression_factor;
  int compression_bit_rate;
  unsigned long long chunk_size[PIDX_MAX_DIMENSIONS];                              ///< size of the block at which compression is applied eg. (4x4x4)
                                                                        ///< the current compression schemes only work in three dimensions
  unsigned long long chunked_bounds[PIDX_MAX_DIMENSIONS];                ///< Compressed global extents

  /// 1 for little endian
  /// 0 for big endian
  int endian;

  /// 1 for flipping endian
  /// 0 for big endian
  int flip_endian;

};
typedef struct idx_file_struct* idx_dataset;


/// idx_dataset_derived_metadata
struct idx_dataset_derived_metadata_struct
{
  //int io_mode;

  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
  
  int fs_block_size;
  off_t start_fs_block;

  Agg_buffer **f0_agg_buffer;
  Agg_buffer **shared_agg_buffer;
  Agg_buffer **nshared_agg_buffer;

  FILE *rst_dump_fp;
  int dump_rst_info;
  char rst_dump_dir_name[512];
  int dump_agg_info;
  char agg_dump_dir_name[512];
  int dump_io_info;
  char io_dump_dir_name[512];

  int color;
  int partition_count[PIDX_MAX_DIMENSIONS];          ///< Number of idx files in each dimensions
  int partition_size[PIDX_MAX_DIMENSIONS];          ///< Number of idx files in each dimensions

  int var_pipe_length;
  int parallel_mode;  
  
  //int staged_aggregation;
  int agg_type;

  int start_layout_index;
  int end_layout_index;


  MPI_Status *status_non_shared;
  MPI_File *fp_non_shared;
  MPI_Request *request_non_shared;

  MPI_Status *status_shared;
  MPI_File *fp_shared;
  MPI_Request *request_shared;


  MPI_Status *status_file_zero;
  MPI_File *fp_file_zero;
  MPI_Request *request_file_zero;

  MPI_File *fp;
  MPI_Request *request;

  int perm_layout_count;
  int layout_count;
  int reduced_res_from;
  int reduced_res_to;

  PIDX_time time;

  int raw_io_pipe_length;

  int aggregator_multiplier;
  int data_core_count;

  int shared_block_level;
  int total_partiton_level;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;



struct idx_debug_info_struct
{
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
