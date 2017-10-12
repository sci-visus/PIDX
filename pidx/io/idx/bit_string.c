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


#include "../../PIDX_inc.h"


PIDX_return_code populate_bit_string(PIDX_io file, int mode)
{
  int i = 0;
  unsigned long long cb[PIDX_MAX_DIMENSIONS];
  unsigned long long* cs = file->idx->chunk_size;

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
  unsigned long long* cs = file->idx->chunk_size;

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

  if (mode == PIDX_WRITE)
  {
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
  }

  // maxh calculation
  file->idx_d->maxh = strlen(file->idx->bitSequence);
  for (i = 0; i <= file->idx_d->maxh; i++)
    file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);

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

  file->idx_d->max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    file->idx_d->max_file_count++;


  int partion_level = (int) log2(/*file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]*/1);
  file->idx_d->total_partiton_level = file->idx->bits_per_block + (int)log2(file->idx->blocks_per_file) + 1 + partion_level;
  if (file->idx_d->total_partiton_level >= file->idx_d->maxh)
    file->idx_d->total_partiton_level = file->idx_d->maxh;

  if (cb[0] == 0 && cb[1] == 0 && cb[2] == 0)
  {
    file->idx_d->maxh = 0;
    file->idx_d->max_file_count = 0;
  }

  return PIDX_success;
}
