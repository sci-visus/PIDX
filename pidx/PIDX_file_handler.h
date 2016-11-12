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


#include "PIDX.h"

///
/// \brief The PIDX_file_descriptor struct is the PIDX File descriptor
/// (equivalent to the descriptor returned by) POSIX or any other IO framework
///
struct PIDX_file_descriptor
{
  int flags;

  PIDX_access access;                          ///< serial or parallel access
  PIDX_io io;

  int local_group_index;                    ///<
  int local_group_count;                    ///<

  int flush_used;
  int write_on_close;                          ///< HPC Writes

  int ROI_writes;

  idx_comm idx_c;


  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;
};
