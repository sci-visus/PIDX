/*
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
 

/*
Block size = 1024
Bits per block = 10
Blocks per file = 4
Dataset size = 32 x 32 x 32 (32768)
maxh = 16
#of agg levels = maxh - (bits_per_block + log2(blocks per file))
Level #Elements [#Block_id] (HZ_start HZ_end) (HZ_start_form HZ_end_form)
0 1 [0] (0,0) [0]
1 1 [0] (1,1) [0] (pow(2, level-1), pow(2, level) - 1)
2 2 [0] (2,3) [0]
3 4 [0] (4,7) [0]
4 8 [0] (8,15) [0]
5 16 [0] (16, 31) [0]
6 32 [0] (32, 63) [0]
7 64 [0] (64, 127) [0]
8 128 [0] (128, 255) [0]
9 256 [0] (256, 511) [0]
10 512 [0] (512, 1023) [0]
11 1024 [1] (1024, 2047) [0]
12 2048 [2, 3] (2048, 4095) [0] ----> 1

13 4096 [4, 5, 6, 7] (4096, 8191) [1]  ----> 2
14 8192 [8, 9, 10, 11, 12, 13, 14, 15] (8192, 16383) [2, 3]  ----> 3, 4
15 16384 [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31] (16384, 32767) [4, 5, 6, 7]  ----> 5, 6, 7, 8
*/

#include "../PIDX_inc.h"

//const int PIDX_default_bits_per_block       = 15;
//const int PIDX_default_blocks_per_file      = 256;

int PIDX_blocks_initialize_layout (PIDX_block_layout layout, int resolution_from, int resolution_to, int maxh, int bits_per_block)
{
  int ctr = 1, j;

  layout->resolution_from = resolution_from;
  layout->resolution_to = resolution_to;

  if (maxh == 0)
    return PIDX_success;

  layout->hz_block_number_array = malloc(sizeof(int*) * maxh);
  memset(layout->hz_block_number_array, 0, sizeof(int*) * maxh);

  for (j = 0 ; j < PIDX_MIN(maxh, bits_per_block + 1) ; j++)
  {
    layout->hz_block_number_array[j] = malloc(sizeof(int));
    memset(layout->hz_block_number_array[j], 0, sizeof(int));
  }

  ctr = 1;
  for (j = bits_per_block + 1 ; j < maxh ; j++)
  {
    layout->hz_block_number_array[j] = malloc(sizeof(int) * ctr);
    memset(layout->hz_block_number_array[j], 0, sizeof(int) * ctr);
    ctr = ctr * 2;
  }

  return PIDX_success;
}


///
int PIDX_blocks_create_layout (int bounding_box[2][5], int maxH, int bits_per_block, const char* bitPattern, PIDX_block_layout layout, int res_from, int res_to)
{
  int m = 0, n_blocks = 1, t = 0, block_number = 1;
  unsigned long long hz_from = 0, hz_to = 0;
  unsigned long long *ZYX_from, *ZYX_to;

  int adjusted_res_to = layout->resolution_to;
  if (layout->resolution_to > maxH - res_to)
    adjusted_res_to = maxH - res_to;//layout->resolution_to - res_to;

  ZYX_to = (unsigned long long*) malloc(sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  ZYX_from = (unsigned long long*) malloc(sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  //fprintf(stderr, "res from to to : %d %d bpb %d\n", layout->resolution_from, layout->resolution_to, bits_per_block);
  if (layout->resolution_from <= bits_per_block)
  {
    for (m = layout->resolution_from ; m <= bits_per_block; m++)
    {
      if (m == 0)
      {
        hz_from = 0;
        hz_to = 0;
      }
      else
      {
        hz_from = (unsigned long long)pow(2, m-1);
        hz_to = (unsigned long long)pow(2, m) - 1;
      }

      memset(ZYX_to, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
      memset(ZYX_from, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

      Hz_to_xyz(bitPattern, maxH - 1, hz_from, ZYX_from);
      Hz_to_xyz(bitPattern, maxH - 1, hz_to, ZYX_to);

      if (ZYX_to[0] >= bounding_box[0][0] && ZYX_from[0] < bounding_box[1][0] && ZYX_to[1] >= bounding_box[0][1] && ZYX_from[1] < bounding_box[1][1] && ZYX_to[2] >= bounding_box[0][2] && ZYX_from[2] < bounding_box[1][2])
        layout->hz_block_number_array[m][0] = 0;
    }

    for (m = bits_per_block + 1 ; m < adjusted_res_to; m++)
    {
      n_blocks = pow(2, (m - (bits_per_block + 1)));
      for (t = 0 ; t < n_blocks ; t++)
      {
        hz_from = (unsigned long long)(block_number) * pow(2, bits_per_block);
        hz_to = (unsigned long long)((block_number + 1) * pow(2, bits_per_block)) - 1;

        memset(ZYX_to, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
        memset(ZYX_from, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

        Hz_to_xyz(bitPattern, maxH - 1, hz_from, ZYX_from);
        Hz_to_xyz(bitPattern, maxH - 1, hz_to, ZYX_to);
      
        if (ZYX_to[0] >= bounding_box[0][0] && ZYX_from[0] < bounding_box[1][0] && ZYX_to[1] >= bounding_box[0][1] && ZYX_from[1] < bounding_box[1][1] && ZYX_to[2] >= bounding_box[0][2] && ZYX_from[2] < bounding_box[1][2])
          layout->hz_block_number_array[m][t] = block_number;
      
        block_number++;
      }
    }
  }
  else
  {
    for (m = bits_per_block + 1 ; m < adjusted_res_to; m++)
    {
      n_blocks = pow(2, (m - (bits_per_block + 1)));
      if (m >= layout->resolution_from)
      {
        for (t = 0 ; t < n_blocks ; t++)
        {
          hz_from = (unsigned long long)(block_number) * pow(2, bits_per_block);
          hz_to = (unsigned long long)((block_number + 1) * pow(2, bits_per_block)) - 1;

          memset(ZYX_to, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
          memset(ZYX_from, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

          Hz_to_xyz(bitPattern, maxH - 1, hz_from, ZYX_from);
          Hz_to_xyz(bitPattern, maxH - 1, hz_to, ZYX_to);

          if (ZYX_to[0] >= bounding_box[0][0] && ZYX_from[0] < bounding_box[1][0] && ZYX_to[1] >= bounding_box[0][1] && ZYX_from[1] < bounding_box[1][1] && ZYX_to[2] >= bounding_box[0][2] && ZYX_from[2] < bounding_box[1][2])
            layout->hz_block_number_array[m][t] = block_number;

          block_number++;
        }
      }
      else
        block_number = block_number + n_blocks;
    }
  }
  free(ZYX_from);
  ZYX_from = 0;
  free(ZYX_to);
  ZYX_to = 0;

  return 0;
}


void PIDX_blocks_print_layout(PIDX_block_layout layout, int bits_per_block)
{
  int ctr = 1, i, j;
  int res_level = 0;
  int res_index = 0;

  if (layout->resolution_from <= bits_per_block)
  {
    for (i = layout->resolution_from; i <= bits_per_block; i++)
    {
      fprintf(stderr, "A Number of blocks at level %d = %d :: ", i, ctr);
      for (j = 0 ; j <  ctr; j++)
      {
        if (layout->hz_block_number_array[i][j] == 0)
          res_level = 0;
        else
        {
          res_level = log2 (layout->hz_block_number_array[i][j]) + 1 + bits_per_block;
          res_index = layout->hz_block_number_array[i][j] % ((int) pow(2, (res_level - 1 - bits_per_block)));
        }
        fprintf(stderr, "[%d (%d %d)] ", layout->hz_block_number_array[i][j], res_level, res_index);
      }
      fprintf(stderr, "\n");
    }

    for (i = bits_per_block + 1; i < layout->resolution_to; i++)
    {
      fprintf(stderr, "B Number of blocks at level %d = %d :: ", i, ctr);
      for (j = 0 ; j <  ctr; j++)
      {
        if (layout->hz_block_number_array[i][j] == 0)
          res_level = 0;
        else
        {
          res_level = log2 (layout->hz_block_number_array[i][j]) + 1 + bits_per_block;
          res_index = layout->hz_block_number_array[i][j] % ((int) pow(2, (res_level - 1 - bits_per_block)));
        }
        fprintf(stderr, "[%d (%d %d)] ", layout->hz_block_number_array[i][j], res_level, res_index);
      }

      ctr = ctr * 2;
      fprintf(stderr, "\n");
    }
  }
  else
  {
    for (i = bits_per_block + 1; i < layout->resolution_to; i++)
    {
      if (i >= layout->resolution_from)
      {
        fprintf(stderr, "C Number of blocks at level %d = %d :: ", i, ctr);
        for (j = 0 ; j <  ctr; j++)
        {
          if (layout->hz_block_number_array[i][j] == 0)
            res_level = 0;
          else
          {
            res_level = log2 (layout->hz_block_number_array[i][j]) + 1 + bits_per_block;
            res_index = layout->hz_block_number_array[i][j] % ((int) pow(2, (res_level - 1 - bits_per_block)));
          }
          fprintf(stderr, "[%d (%d %d)] ", layout->hz_block_number_array[i][j], res_level, res_index);
        }
      }
      ctr = ctr * 2;
      if (i >= layout->resolution_from)
        fprintf(stderr, "\n");
    }
  }
}


int PIDX_blocks_is_block_present(int block_number, int bits_per_block, PIDX_block_layout layout)
{
  int res_level = 0, res_index = 0;

  if (block_number == 0)
  {
    if (layout->resolution_from <= bits_per_block)
      return 1;

    return 0;
    //res_level = bits_per_block;
  }
  else
  {
    res_level = log2 (block_number) + 1 + bits_per_block;
    res_index = block_number % ((int) pow(2, (res_level - 1 - bits_per_block)));
  }

  if (res_level < layout->resolution_from || res_level >= layout->resolution_to)
    return 0;

  if (layout->hz_block_number_array[res_level][res_index] == block_number)
    return 1;

  return 0;
}


int PIDX_blocks_find_negative_offset(int blocks_per_file, int bits_per_block, int block_number, PIDX_block_layout layout)
{
  //PIDX_blocks_print_layout(layout);
  int b;
  int file_no = block_number/blocks_per_file;
  int block_offset = 0;
  for (b = file_no * blocks_per_file ; b < block_number; b++)
    if (!PIDX_blocks_is_block_present(b, bits_per_block, layout))
      block_offset++;
  return block_offset;
}


void PIDX_blocks_free_layout(int bits_per_block, int maxh, PIDX_block_layout layout)
{
  int j = 0;
  for (j = 0 ; j < bits_per_block + 1 ; j++)
  {
    free(layout->hz_block_number_array[j]);
    layout->hz_block_number_array[j] = 0;
  }

  for (j = bits_per_block + 1 ; j < maxh ; j++)
  {
    free(layout->hz_block_number_array[j]);
    layout->hz_block_number_array[j] = 0;
  }

  free(layout->hz_block_number_array);
  layout->hz_block_number_array = 0;

  return;
}

void PIDX_free_layout(PIDX_block_layout layout)
{
  free(layout->file_bitmap);
  layout->file_bitmap = 0;

  free(layout->bcpf);
  layout->bcpf = 0;

  free(layout->lbi);
  layout->lbi = 0;

  free(layout->existing_file_index);
  layout->existing_file_index = 0;

  free(layout->inverse_existing_file_index);
  layout->inverse_existing_file_index = 0;

  free(layout->file_index);
  layout->file_index = 0;
}
