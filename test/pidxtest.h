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

#ifndef PIDXTEST_H
#define PIDXTEST_H

#include <PIDX.h>

/// Kind of test to run
enum Kind { DEFAULT = 0, SERIAL_READER, PARALLEL_READER, SERIAL_WRITER, PARALLEL_WRITER, PARALLEL_CONVERTER, PARALLEL_MULTI_PATCH_WRITER, HDF5_WRITER, HDF5_READER};

/// kindToStr
char* kindToStr(enum Kind k);

/// strToKind
enum Kind strToKind(const char *str);

/// Args
struct Args
{
  /// type of test to run
  enum Kind kind;

  /// global dimensions of 3D volume
  int64_t extents[5];

  /// per-process dimensions of each sub-block
  int count_local[5];

  /// per-process dimensions of each sub-block
  int restructured_box_size[5];

  /// Number of time-steps
  int time_step;

  /// Number of Variables
  int variable_count;

  /// Number of IDX files partitions in each direction (x, y, z)
  int idx_count[3];

  /// output IDX file Name Template
  char output_file_template[512];

  /// output IDX filename
  char *output_file_name;

  /// One of the parameters that controls the aggegator count
  int aggregation_factor;

  /// Number of IDX blocks for every IDX file
  int blocks_per_file;

  /// Controls the size of every IDX block
  int bits_per_block;

  /// Works only for this test program verifies the functioning of the HZ encoding phase
  /// Default is 0
  int debug_hz;

  /// Works only for this test program verifies the functioning of the restructing phase
  int debug_rst;

  /// Dumps all the reguired offsets and counts for the aggregation phase
  /// Default is 0
  int dump_agg;

  /// All the agg dumps will go here (file equal to number of processes will be ceated)
  char agg_dump_dir_name[512];

  ///
  ///
  int perform_brst;

  ///
  ///
  int perform_compression;

  /// 0 for no aggregation (test IO only) and 1 for aggregation
  /// Default is 1
  int perform_agg;

  /// 0 for no HZ encoding 1 for HZ encoding
  /// Default is 0
  int perform_hz;

  /// 0 for no IO (tests only aggregation phase) 1 for IO
  /// Default is 0
  int perform_io;

  ///
  int64_t compression_block_size[5];

  /// 1 for lossy 0 for lossless
  int compression_type;


  ///
  int compression_bit_rate;


  ///
  int is_rank_z_ordering;

  ///
  int is_global_indexing;

  ///
  int hz_from;

  ///
  int hz_to;
};

/// main
int main(int argc, char **argv);

/// parse_args
int parse_args(struct Args *args, int argc, char **argv);

/// print_error
int print_error(char *error_message, char* file, int line);

#endif
