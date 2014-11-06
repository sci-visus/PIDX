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
 
#include "PIDX_inc.h"

const int PIDX_default_bits_per_block       = 15;
const int PIDX_default_blocks_per_file      = 256;

int is_block_present(int block_number, block_layout* layout)
{
  long i, j;
  
  if( block_number == 0)
    return 1;
  
  for(i = 1 ; i < layout->levels ; i++)
    for(j = 0 ; j < layout->hz_block_count_array[i] ; j++)
      if ( layout->hz_block_number_array[i][j] == block_number)
	return 1;
    
  return(0);
}

int find_block_negative_offset(int blocks_per_file, int block_number, block_layout* layout)
{
  int b;
  int file_no = block_number/blocks_per_file;
  int block_offset = 0;
  for(b = file_no * blocks_per_file ; b < block_number; b++)
    if (!is_block_present(b, layout))
      block_offset++;
  return block_offset;
}

void destroyBlockBitmap(block_layout* layout)
{
  int j = 0;
  free(layout->hz_block_count_array);
  for(j = 0 ; j < (layout->levels) ; j++)
  {     
    free(layout->hz_block_number_array[j]);
    layout->hz_block_number_array[j] = 0;
  }
  free(layout->hz_block_number_array);
  layout->hz_block_number_array = 0;
  layout->levels = 0;
}

void print_block_layout(block_layout* layout)
{
  int ctr = 1, i, j;
  printf("Number of levels: %d\n", layout->levels);
  for (i = 1; i < (layout->levels); i++) 
  {
    printf("Number of blocks at level %d = %d::", i, layout->hz_block_count_array[i]);
    for(j = 0 ; j <  ctr; j++)
      if(layout->hz_block_number_array[i][j] != 0)
	printf("%d ", layout->hz_block_number_array[i][j]);
    
    ctr = ctr * 2;
    printf("\n");
  }
}

int initialize_block_layout(block_layout* layout, int maxh, int bits_per_block)
{
  int ctr = 1, i, j;
  if(maxh < bits_per_block)
    layout->levels = 1;
  else
    layout->levels = maxh - bits_per_block;

  layout->hz_block_count_array = (int*)malloc(sizeof(int) * (layout->levels));
  layout->hz_block_count_array[0] = 1;
  layout->hz_block_number_array = (int**)malloc(sizeof(int*) * (layout->levels));
  layout->hz_block_number_array[0] = (int*)malloc(sizeof(int) * ctr);
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    layout->hz_block_count_array[j] = 0;
    layout->hz_block_number_array[j] = (int*)malloc(sizeof(int) * ctr);
    ctr = ctr * 2;
  }

  ctr = 1;
  layout->hz_block_number_array[0][0] = 0;
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    for(i = 0 ; i < ctr ; i++)
    {
      layout->hz_block_number_array[j][i] = 0;
    }
    ctr = ctr * 2;
  }
  layout->hz_block_count_array[0] = 1;			//This block contains data upto level "bits_per_block"
  return 1;
}

int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, const char* bitPattern, block_layout* layout)
{
  long long hz_from = 0, hz_to = 0, block_number = 1;
  int i, j, m, n_blocks = 1, ctr = 1;
  long long *ZYX_from, *ZYX_to;
  
  if(maxH < bits_per_block)
    layout->levels = 1;
  else
    layout->levels = maxH - bits_per_block;

  layout->hz_block_count_array = (int*)malloc(sizeof(int) * (layout->levels));
  layout->hz_block_count_array[0] = 1;
  layout->hz_block_number_array = (int**)malloc(sizeof(int*) * (layout->levels));
  layout->hz_block_number_array[0] = (int*)malloc(sizeof(int) * ctr);
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    layout->hz_block_count_array[j] = 0;
    layout->hz_block_number_array[j] = (int*)malloc(sizeof(int) * ctr);
    ctr = ctr * 2;
  }
  ctr = 1;
  layout->hz_block_number_array[0][0] = 0;
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    for(i = 0 ; i < ctr ; i++)
    {
      layout->hz_block_number_array[j][i] = 0;
    }
    ctr = ctr * 2;
  }
  layout->hz_block_count_array[0] = 1;			//This block contains data upto level "bits_per_block"
  
  hz_from = (long long)(block_number - 1) * pow(2, bits_per_block);
  hz_to = (long long)(block_number * pow(2, bits_per_block)) - 1;
  
  for(m = 1 ; m < (maxH - bits_per_block); m++)
  {
    n_blocks = pow(2, (m - 1));
    int t = 0;
    for(t = 0 ; t < n_blocks ; t++)
    {
      block_number = block_number + 1;
      
      hz_from = (long long)(block_number - 1) * pow(2, bits_per_block);
      hz_to = (long long)(block_number * pow(2, bits_per_block)) - 1;
      
      ZYX_to = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
      ZYX_from = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
      
      Hz_to_xyz(bitPattern, maxH - 1, hz_from, ZYX_from);
      Hz_to_xyz(bitPattern, maxH - 1, hz_to, ZYX_to);
      
      if(ZYX_to[0] >= bounding_box[0][0] && ZYX_from[0] < bounding_box[1][0] && 
	 ZYX_to[1] >= bounding_box[0][1] && ZYX_from[1] < bounding_box[1][1] && 
	 ZYX_to[2] >= bounding_box[0][2] && ZYX_from[2] < bounding_box[1][2] &&
	 ZYX_to[3] >= bounding_box[0][3] && ZYX_from[3] < bounding_box[1][3] &&
	 ZYX_to[4] >= bounding_box[0][4] && ZYX_from[4] < bounding_box[1][4])
      {
	layout->hz_block_count_array[m] = layout->hz_block_count_array[m] + 1;
	layout->hz_block_number_array[m][t] = block_number - 1;
      }
      free(ZYX_from);
      ZYX_from = 0;
      free(ZYX_to);
      ZYX_to = 0;
    }
  }
  return 0;
}