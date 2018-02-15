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
 * \file PIDX_meta_data_cache.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 *
 */
 
#ifndef __PIDX_META_DATA_CACHE_H
#define __PIDX_META_DATA_CACHE_H


struct PIDX_meta_data_cache_struct
{
  int is_set;               /// flag to specify if cache buffer is populated
  int element_count;        /// Number of elements (tuples) in the cache
  int *xyz_mapped_index;    /// The xyz index (application row-order index)
  int *hz_level;            /// Corresponding HZ index to the xyz index
  int *index_level;         /// The hz index level
};
typedef struct PIDX_meta_data_cache_struct* PIDX_meta_data_cache;


PIDX_return_code PIDX_create_meta_data_cache(PIDX_meta_data_cache* cache);


PIDX_return_code PIDX_free_meta_data_cache(PIDX_meta_data_cache cache);


#endif
