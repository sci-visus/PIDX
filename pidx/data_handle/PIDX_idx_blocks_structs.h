/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2010-2018 ViSUS L.L.C.,
 * Scientific Computing and Imaging Institute of the University of Utah
 *
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 *
 */

#ifndef __PIDX_IDX_BLOCKS_STRUCTS_H
#define __PIDX_IDX_BLOCKS_STRUCTS_H

/*

  IDX File with resolution 8 x 8
  Bits per block 3
  Blocks per file 2
  Two processes writing out this dataset

  Level     #elements   #file_id    #agg_group_id
  0         1           0           0
  1         1           0           0
  2         2           0           0
  3         4           0           0
  4         8           0           0
  5         16          1           1
  6         32          2,3         2

  Number of files (8 x 8) / ((2 ^ 3) x 2) = 4


  Case A:
  No partitioning (one global IDX file)

  // partition#     hz_file0_from      hz_file0_to       hz_n_file0_from       hz_n_file0_to
  //          1                 0                5                     5                   7


  0   16   4   18  |  1    24   6    26         0    5    3    5  |  1    5    3    5       0    1    0    1  |  0    1    0    1
                   |                                              |                                           |
  32  33   36  37  |  48   49   52   53         6    6    6    6  |  6    6    6    6       2    2    2    2  |  3    3    3    3
                   |                                              |                                           |
  8   17   9   19  |  12   25   13   27         4    5    4    5  |  4    5    4    5       0    1    0    1  |  0    1    0    1
                   |                                              |                                           |
  34  35   38  39  |  50   51   54   55         6    6    6    6  |  6    6    6    6       2    2    2    2  |  3    3    3    3
                   |                                              |                                           |
  -------------------------------------        --------------------------------------       --------------------------------------
                   |                                              |                                           |
  2   20   5   22  |  3    28   7    30         2    5    3    5  |  2    5    3    5       0    1    0    1  |  0    1    0    1
                   |                                              |                                           |
  40  41   44  45  |  56   57   60   61         6    6    6    6  |  6    6    6    6       2    2    2    2  |  3    3    3    3
                   |                                              |                                           |
  10  21   11  23  |  14   29   15   31         4    5    4    5  |  4    5    4    5       0    1    0    1  |  0    1    0    1
                   |                                              |                                           |
  42  43   46  47  |  58   59   62   63         6    6    6    6  |  6    6    6    6       2    2    2    2  |  3    3    3    3



  Case B:
  Four partitions

  // partition#     hz_file0_from      hz_file0_to       hz_n_file0_from       hz_n_file0_to
  //          4                 0                5                     5                   5

  0    4    1    6  |  0    4    1    6         0    3    1    3  |  0    3    1    3       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
  8    9   12   13  |  8    9   12   13         4    4    4    4  |  4    4    4    4       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
  2    5    3    7  |  2    5    3    7         2    3    2    3  |  2    3    2    3       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
 10   11   14   15  | 10   11   14   15         4    4    4    4  |  4    4    4    4       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
  -------------------------------------        --------------------------------------       --------------------------------------
                    |                                             |                                           |
  0    4    1    6  |  0    4    1    6         0    3    1    3  |  0    3    1    3       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
  8    9   12   13  |  8    9   12   13         4    4    4    4  |  4    4    4    4       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
  2    5    3    7  |  2    5    3    7         2    3    2    3  |  2    3    2    3       0    0    0    0  |  0    0    0    0
                    |                                             |                                           |
 10   11   14   15  | 10   11   14   15         4    4    4    4  |  4    4    4    4       0    0    0    0  |  0    0    0    0



 */

/// IDX Blocks related struct
struct idx_blocks_struct
{

  // The reason we seperate file 0 from other files is because. File 0 and File 1 (from the non-file-zero group)
  // share the same extents w.r.t. the number of processes that contributes to those files.
  // By treating file 0 seperately, we avoid clash between file 0 and file 1 while assigning aggregators.

  // partition#      hz_file0_from       hz_file0_to          hz_n_file0_from          hz_n_file0_to
  //          1                  0                 5                        5                      7
  //          2                  0                 5                        5                      6
  //          4                  0                 5                        5                      5

  // Note that the total number of levels decreases with more partitioning, so hz_n_file0_to range also adjusts accordingly.

  // hz range for file 0 [hz_file0_from, hz_file0_to)
  int hz_file0_from;                            /// this is set to 0
  int hz_file0_to;                              /// this is bits_per_block + log2(blocks_per_file) + 1

  // hz range for all other file [hz_n_file0_from, hz_n_file0_to)
  int hz_n_file0_from;                          /// this is bits_per_block + log2(blocks_per_file) + 1
  int hz_n_file0_to;                            /// this is maxh

  // this is the aggregation group for file 0
  int file0_agg_group_from_index;               /// defaults to 0
  int file0_agg_group_to_index;                 /// defaults to 1
  int file0_agg_group_count;                    /// defaults to 1 - 0 = 1

  // all other aggregation groups
  int nfile0_agg_group_from_index;              /// defaults to 1
  int nfile0_agg_group_to_index;                /// is computed using maxh
  int nfile0_agg_group_count;

  int agg_level;
  PIDX_block_layout block_layout;               /// block layout for the entire idx file
  PIDX_block_layout* block_layout_by_agg_group; /// block layout for every agg group (file0_agg_group_count + nfile0_agg_group_count)

  // for reduced resolution io
  int reduced_res_from;                         /// skip hz levels from 0 to rediced_res_from
  int reduced_res_to;                           /// write only upto reduced_res_to

  // the bitmap stores if a block b exists within a file n; if (block_bitmap[n][b] == 1)
  int **block_bitmap;

  // you can ignore this, it is only used in one case (serial run), it stores offset of every block, for every variable in every idx file
  int ***block_offset_bitmap;
};
typedef struct idx_blocks_struct* idx_blocks;

#endif
