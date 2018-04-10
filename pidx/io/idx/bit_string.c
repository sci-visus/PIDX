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


#include "../../PIDX_inc.h"


PIDX_return_code populate_bit_string(PIDX_io file, int mode)
{
  size_t* cs = file->idx->chunk_size;

  if (mode == PIDX_WRITE)
  {
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];
    char partition_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) ceil((float)file->idx_d->restructured_grid->patch_size[0] / cs[0]);
    rpp.y = (int) ceil((float)file->idx_d->restructured_grid->patch_size[1] / cs[1]);
    rpp.z = (int) ceil((float)file->idx_d->restructured_grid->patch_size[2] / cs[2]);
    if (rpp.x == 0)  rpp.x = 1;
    if (rpp.y == 0)  rpp.y = 1;
    if (rpp.z == 0)  rpp.z = 1;
    guess_bit_string_ZYX(reg_patch_bs, rpp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[1] %s : %d %d %d\n", reg_patch_bs, rpp.x, rpp.y, rpp.z);
#endif

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) file->idx_d->partition_size[0] / file->idx_d->restructured_grid->patch_size[0];
    prcp.y = (int) file->idx_d->partition_size[1] / file->idx_d->restructured_grid->patch_size[1];
    prcp.z = (int) file->idx_d->partition_size[2] / file->idx_d->restructured_grid->patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[2] %s : %d %d %d\n", process_bs, prcp.x, prcp.y, prcp.z);
#endif

    // Last part of the bitstring
    Point3D pcp;
    pcp.x = (int) file->idx_d->partition_count[0];
    pcp.y = (int) file->idx_d->partition_count[1];
    pcp.z = (int) file->idx_d->partition_count[2];
    guess_bit_string(partition_bs, pcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[3] %s : %d %d %d\n", partition_bs, pcp.x, pcp.y, pcp.z);
#endif

    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    strcpy(file->idx->bitSequence, partition_bs);
    strcat(file->idx->bitSequence, temp_bs + 1);
  }

  return PIDX_success;
}


PIDX_return_code populate_global_bit_string(PIDX_io file, int mode)
{
  int i = 0;
  unsigned long long cb[PIDX_MAX_DIMENSIONS];
  size_t* cs = file->idx->chunk_size;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx->bounds[i] % file->idx->chunk_size[i] == 0)
      cb[i] = (int) file->idx->bounds[i] / file->idx->chunk_size[i];
    else
      cb[i] = (int) (file->idx->bounds[i] / file->idx->chunk_size[i]) + 1;
  }

  if (mode == PIDX_WRITE)
  {
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];
    char partition_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx_d->restructured_grid->patch_size[0] / cs[0];
    rpp.y = (int) file->idx_d->restructured_grid->patch_size[1] / cs[1];
    rpp.z = (int) file->idx_d->restructured_grid->patch_size[2] / cs[2];
    if (rpp.x == 0)  rpp.x = 1;
    if (rpp.y == 0)  rpp.y = 1;
    if (rpp.z == 0)  rpp.z = 1;
    guess_bit_string_ZYX(reg_patch_bs, rpp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[1] %s : %d %d %d\n", reg_patch_bs, rpp.x, rpp.y, rpp.z);
#endif

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) file->idx_d->partition_size[0] / file->idx_d->restructured_grid->patch_size[0];
    prcp.y = (int) file->idx_d->partition_size[1] / file->idx_d->restructured_grid->patch_size[1];
    prcp.z = (int) file->idx_d->partition_size[2] / file->idx_d->restructured_grid->patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[2] %s : %d %d %d\n", process_bs, prcp.x, prcp.y, prcp.z);
#endif

    // Last part of the bitstring
    Point3D pcp;
    pcp.x = (int) file->idx_d->partition_count[0];
    pcp.y = (int) file->idx_d->partition_count[1];
    pcp.z = (int) file->idx_d->partition_count[2];
    guess_bit_string(partition_bs, pcp);

#if DETAIL_OUTPUT
    if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
      fprintf(stderr, "[3] %s : %d %d %d\n", partition_bs, pcp.x, pcp.y, pcp.z);
#endif

    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    strcpy(file->idx->bitSequence, partition_bs);
    strcat(file->idx->bitSequence, temp_bs + 1);
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

#if DETAIL_OUTPUT
  if (file->idx_c->grank == 0 && file->idx->cached_ts == file->idx->current_time_step)
    fprintf(stderr, "Bitstring %s maxh %d\n", file->idx->bitSequence, file->idx_d->maxh);
#endif

  unsigned long long total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  unsigned long long max_sample_per_file = (unsigned long long) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

  int partion_level = (int) log2(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  return PIDX_success;
}



PIDX_return_code populate_local_bit_string(PIDX_io file, int mode)
{
  int i = 0;
  unsigned long long cb[PIDX_MAX_DIMENSIONS];

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx->box_bounds[i] % file->idx->chunk_size[i] == 0)
      cb[i] = (int) file->idx->box_bounds[i] / file->idx->chunk_size[i];
    else
      cb[i] = (int) (file->idx->box_bounds[i] / file->idx->chunk_size[i]) + 1;
  }

  //if (mode == PIDX_WRITE)
  //{
    char temp_bs[512];
    char reg_patch_bs[512];
    char process_bs[512];

    // First part of the bitstring
    Point3D rpp;
    rpp.x = (int) file->idx_d->restructured_grid->patch_size[0];
    rpp.y = (int) file->idx_d->restructured_grid->patch_size[1];
    rpp.z = (int) file->idx_d->restructured_grid->patch_size[2];
    guess_bit_string_ZYX(reg_patch_bs, rpp);
    //if (file->idx_c->lrank == 0)
    //  fprintf(stderr, "[1X %d] %s : %d %d %d\n", file->idx_d->color, reg_patch_bs, rpp.x, rpp.y, rpp.z);

    // Middle part of the bitstring
    Point3D prcp;
    prcp.x = (int) getPowerOf2(file->idx->box_bounds[0]) / file->idx_d->restructured_grid->patch_size[0];
    prcp.y = (int) getPowerOf2(file->idx->box_bounds[1]) / file->idx_d->restructured_grid->patch_size[1];
    prcp.z = (int) getPowerOf2(file->idx->box_bounds[2]) / file->idx_d->restructured_grid->patch_size[2];
    //prcp.x = (int) file->idx_d->partition_size[0] / file->idx->reg_patch_size[0];
    //prcp.y = (int) file->idx_d->partition_size[1] / file->idx->reg_patch_size[1];
    //prcp.z = (int) file->idx_d->partition_size[2] / file->idx->reg_patch_size[2];
    if (prcp.x == 0)  prcp.x = 1;
    if (prcp.y == 0)  prcp.y = 1;
    if (prcp.z == 0)  prcp.z = 1;
    guess_bit_string_Z(process_bs, prcp);
    //if (file->idx_d->color == 1)
    //  fprintf(stderr, "[2Y %d] %s : %d %d %d - %d %d %d\n", file->idx_d->color, process_bs, prcp.x, prcp.y, prcp.z, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2]);


    // Concatenating the three components to get the final bit string
    strcpy(temp_bs, process_bs);
    strcat(temp_bs, reg_patch_bs + 1);
    //strcpy(file->idx->bitSequence, partition_bs);
    //strcat(file->idx->bitSequence, temp_bs + 1);
    strcpy(file->idx->bitSequence, temp_bs);
    //fprintf(stderr, "[%d] - bit sequence %s\n", file->idx_d->color, file->idx->bitSequence);
  //}

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

  if (file->idx->bits_per_block >= file->idx_d->maxh)
  {
    file->idx->bits_per_block = file->idx_d->maxh - 1;
    file->idx_d->samples_per_block = pow(2, file->idx->bits_per_block);
  }

  //if (file->idx_c->lrank == 0)
  //  fprintf(stderr, "%d Bitstring %s maxh %d\n", file->idx_d->color, file->idx->bitSequence, file->idx_d->maxh);

  unsigned long long total_reg_sample_count = (getPowerOf2(cb[0]) * getPowerOf2(cb[1]) * getPowerOf2(cb[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  unsigned long long max_sample_per_file = (unsigned long long) file->idx_d->samples_per_block * file->idx->blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, file->idx_d->samples_per_block, file->idx->blocks_per_file);
    return PIDX_err_file;
  }
  //fprintf(stderr, "total reg sample count %d %d\n", total_reg_sample_count, max_sample_per_file);

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;

#if 1
  int partion_level = (int) log2(/*file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]*/1);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  if (cb[0] == 0 && cb[1] == 0 && cb[2] == 0)
  {
    file->idx_d->maxh = 0;
    file->idx_d->max_file_count = 0;
  }
#endif
  return PIDX_success;
}
