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
#include "../../../PIDX_inc.h"


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int svi, int evi);
static PIDX_return_code create_block_layout_and_buffers(PIDX_io file, int svi, int evi, int mode);
static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout global_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int start_index, int hz_level_from, int hz_level_to);
static PIDX_return_code destroy_block_layout(PIDX_io file);
static PIDX_return_code delete_sim_block_layout(PIDX_io file);
static PIDX_return_code populate_sim_block_layouts(PIDX_io file, int svi, int hz_file0_from, int hz_n_file0_to);


PIDX_return_code PIDX_serial_idx_write(PIDX_io file, int svi, int evi)
{
  int bytes_for_datatype;
  int i = 0, j = 0, k = 0;
  int si = 0;
  PIDX_return_code ret;
  char file_name[PATH_MAX];
  MPI_File fp = 0;
  MPI_Status status;
  set_rst_box_size_for_write(file, svi);

  if (create_block_layout_and_buffers(file, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (si = svi; si < evi; si++)
  {
    bytes_for_datatype = ((file->idx->variable[si]->bpv / 8) * file->idx->variable[si]->vps);
    file->idx->variable_tracker[si] = 1;

    uint64_t index = 0;
    uint64_t hz;
    uint64_t xyz[PIDX_MAX_DIMENSIONS];
    unsigned char* block_buffer = malloc(file->idx->samples_per_block * bytes_for_datatype);

    for (i = 0; i < file->idx->max_file_count; i++)
    {
      if (generate_file_name(file->idx->blocks_per_file, file->idx->filename_template_partition, i, file_name, PATH_MAX) == 1)
      {
        fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      if (MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() filename %s failed.\n", __FILE__, __LINE__, file_name);
        return PIDX_err_io;
      }

      for (j = 0; j < file->idx->blocks_per_file; j++)
      {
        if (file->idx_b->block_bitmap[i][j] != 0)
        {
          //fprintf(stderr, "File number %d Block number %d\n", i, j);
          for (k = 0; k < file->idx->samples_per_block; k++)
          {
            hz = (i * file->idx->blocks_per_file * file->idx->samples_per_block) + (j * file->idx->samples_per_block) + k;
            Hz_to_xyz(file->idx->bitPattern, file->idx->maxh, hz, xyz);

            index = (file->idx->variable[si]->sim_patch[0]->size[0] * file->idx->variable[si]->sim_patch[0]->size[1] * xyz[2])
                + (file->idx->variable[si]->sim_patch[0]->size[0] * xyz[1])
                + xyz[0];

            if (xyz[0] >= file->idx->bounds[0] || xyz[1] >= file->idx->bounds[1] || xyz[2] >= file->idx->bounds[2])
              continue;

            memcpy(block_buffer + (k * bytes_for_datatype),
                 file->idx->variable[si]->sim_patch[0]->buffer + (index * bytes_for_datatype),
                bytes_for_datatype);
          }

          if (MPI_File_write_at(fp, file->idx_b->block_offset_bitmap[si][i][j], block_buffer, file->idx->samples_per_block * bytes_for_datatype, MPI_BYTE, &status) != MPI_SUCCESS)
          {
            fprintf(stderr, "[%s] [%d] MPI_File_open() failed.\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }
        }
      }
      MPI_File_close(&fp);
    }
    free(block_buffer);
  }

  // Step 9
  ret = group_meta_data_finalize(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



static PIDX_return_code create_block_layout_and_buffers(PIDX_io file, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->time;

  time->bit_string_start = PIDX_get_time();

  // calculates maxh and bitstring
  ret = populate_global_bit_string(file, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // selects layout levels based on maxh
  select_io_mode(file);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layoutven this is pure IDX only non-share block layout is populated
  ret = populate_sim_block_layouts(file, svi, file->idx_b->hz_file0_from, file->idx_b->hz_n_file0_to);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->layout_end = PIDX_get_time();



  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  ret = write_headers(file, svi, evi, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();


  return PIDX_success;
}


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int svi, int evi)
{
  int ret;
  PIDX_time time = file->time;

  time->group_cleanup_start = PIDX_get_time();
  ret = delete_sim_block_layout(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}



static PIDX_return_code populate_sim_block_layouts(PIDX_io file, int svi, int hz_file0_from, int hz_n_file0_to)
{
  int i = 0, v = 0, ret;

  file->idx_b->block_bitmap = malloc(file->idx->max_file_count * sizeof (*file->idx_b->block_bitmap));
  memset(file->idx_b->block_bitmap, 0, file->idx->max_file_count * sizeof (*file->idx_b->block_bitmap));
  for (i = 0; i < file->idx->max_file_count; i++)
  {
    file->idx_b->block_bitmap[i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_b->block_bitmap[i]));
    memset(file->idx_b->block_bitmap[i], 0, file->idx->blocks_per_file * sizeof (*file->idx_b->block_bitmap[i]));
  }

  file->idx_b->block_offset_bitmap = malloc(file->idx->variable_count * sizeof (*file->idx_b->block_offset_bitmap));
  memset(file->idx_b->block_offset_bitmap, 0, file->idx->variable_count * sizeof (*file->idx_b->block_offset_bitmap));
  for (v = 0; v < file->idx->variable_count; v++)
  {
    file->idx_b->block_offset_bitmap[v] = malloc(file->idx->max_file_count * sizeof (*(file->idx_b->block_offset_bitmap[v])));
    memset(file->idx_b->block_offset_bitmap[v], 0, file->idx->max_file_count * sizeof (*(file->idx_b->block_offset_bitmap[v])));

    for (i = 0; i < file->idx->max_file_count; i++)
    {
      file->idx_b->block_offset_bitmap[v][i] = malloc(file->idx->blocks_per_file * sizeof (*file->idx_b->block_offset_bitmap[v][i]));
      memset(file->idx_b->block_offset_bitmap[v][i], 0, file->idx->blocks_per_file * sizeof (*file->idx_b->block_offset_bitmap[v][i]));
    }
  }

  int total_layout_count = file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count;
  file->idx_b->block_layout = malloc(sizeof (*file->idx_b->block_layout));
  memset(file->idx_b->block_layout, 0, sizeof (*file->idx_b->block_layout));

  file->idx_b->block_layout_by_agg_group = malloc(sizeof (*file->idx_b->block_layout_by_agg_group) * total_layout_count);
  memset(file->idx_b->block_layout_by_agg_group, 0, sizeof (*file->idx_b->block_layout_by_agg_group) * total_layout_count);

  for (i = 0; i < total_layout_count ; i++)
  {
    file->idx_b->block_layout_by_agg_group[i] = malloc(sizeof(*(file->idx_b->block_layout_by_agg_group[i])));
    memset(file->idx_b->block_layout_by_agg_group[i], 0, sizeof(*(file->idx_b->block_layout_by_agg_group[i])));
  }

  ret = populate_idx_block_layout(file,
                  file->idx_b->block_layout, file->idx_b->block_layout_by_agg_group,
                  file->idx_b->file0_agg_group_from_index, file->idx_b->nfile0_agg_group_to_index,
                  file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count,
                  svi,
                  hz_file0_from, hz_n_file0_to);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}

static PIDX_return_code populate_idx_layout_serial(PIDX_io file, int start_var_index, PIDX_block_layout block_layout, int lower_hz_level, int higher_hz_level)
{
  int i, j;
  int p = 0, ctr = 1;
  PIDX_return_code ret_code;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  int lvi = start_var_index;


  PIDX_block_layout all_patch_local_block_layout = malloc(sizeof (*all_patch_local_block_layout));
  memset(all_patch_local_block_layout, 0, sizeof (*all_patch_local_block_layout));
  ret_code = PIDX_blocks_initialize_layout(all_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable var = file->idx->variable[lvi];


    for (p = 0 ; p < var->sim_patch_count ; p++)
    {
      for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      {
        bounding_box[0][i] = var->sim_patch[p]->offset[i];
        bounding_box[1][i] = var->sim_patch[p]->size[i] + var->sim_patch[p]->offset[i];

        bounding_box[0][i] = (bounding_box[0][i] / file->idx->chunk_size[i]);

        if (bounding_box[1][i] % file->idx->chunk_size[i] == 0)
          bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]);
        else
          bounding_box[1][i] = (bounding_box[1][i] / file->idx->chunk_size[i]) + 1;
      }

      PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
      memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
      ret_code = PIDX_blocks_initialize_layout(per_patch_local_block_layout, lower_hz_level, higher_hz_level, file->idx->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = PIDX_blocks_create_layout (bounding_box, file->idx->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_b->reduced_resolution_factor);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      if (all_patch_local_block_layout->resolution_from <= file->idx->bits_per_block)
      {
        for (i = all_patch_local_block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
        {
          if (per_patch_local_block_layout->hz_block_number_array[i][0] == 0)
          {
            all_patch_local_block_layout->hz_block_number_array[i][0] = per_patch_local_block_layout->hz_block_number_array[i][0];
            break;
          }
        }

        ctr = 1;
        for (i = file->idx->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
        {
          for (j = 0 ; j < ctr ; j++)
          {
            if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
              all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
          }
          ctr = ctr * 2;
        }
      }
      else
      {
        ctr = 1;
        for (i = file->idx->bits_per_block + 1 ; i < all_patch_local_block_layout->resolution_to ; i++)
        {
          if (i >= all_patch_local_block_layout->resolution_from)
          {
            for (j = 0 ; j < ctr ; j++)
            {
              if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
                all_patch_local_block_layout->hz_block_number_array[i][j] = per_patch_local_block_layout->hz_block_number_array[i][j];
            }
          }
          ctr = ctr * 2;
        }
      }

      PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, per_patch_local_block_layout);
      free(per_patch_local_block_layout);
      per_patch_local_block_layout = 0;
    }


  if (block_layout->resolution_from <=   file->idx->bits_per_block)
  {
    int level_count = 1;
    for (i = block_layout->resolution_from; i <=   file->idx->bits_per_block; i++)
    {
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->partition_comm);
    }

    for (i =   file->idx->bits_per_block + 1; i < (block_layout->resolution_to); i++)
    {
      MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->partition_comm);
      level_count = level_count * 2;
    }
  }
  else
  {
    int level_count = 1;
    for (i =   file->idx->bits_per_block + 1; i < (block_layout->resolution_to); i++)
    {
      if (i >= block_layout->resolution_from)
      {
        MPI_Allreduce(all_patch_local_block_layout->hz_block_number_array[i], block_layout->hz_block_number_array[i], level_count, MPI_INT, MPI_BOR, file->idx_c->partition_comm);

      }
      level_count = level_count * 2;
    }
  }


  PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, all_patch_local_block_layout);
  free(all_patch_local_block_layout);
  all_patch_local_block_layout = 0;

  block_layout->file_bitmap = malloc(file->idx->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx->max_file_count));

  block_layout->bcpf = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->bcpf, 0, sizeof(int) * (file->idx->max_file_count));

  block_layout->lbi = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->lbi, 0, sizeof(int) * (file->idx->max_file_count));

  int file_number = 0;
  if (block_layout->resolution_from <=   file->idx->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->bcpf[file_number]++;
        block_layout->lbi[file_number] = 0;
        break;
      }
    }

    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->bcpf[file_number]++;
          block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
        }
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->bcpf[file_number]++;
            block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
          }
        }
      }
      ctr = ctr * 2;
    }
  }
  block_layout->efc = 0;
  for (i = 0; i < file->idx->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->efc++;

  block_layout->existing_file_index = (int*) malloc(block_layout->efc * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->efc * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;

      count++;
    }
  }

  return PIDX_success;
}


static PIDX_return_code populate_idx_block_layout(PIDX_io file, PIDX_block_layout block_layout, PIDX_block_layout* layout_by_level, int start_layout_index, int end_layout_index, int layout_count, int si, int hz_level_from, int hz_level_to)
{
  if (hz_level_from == 0 && hz_level_to == 0)
    return PIDX_success;

  PIDX_return_code ret_code;

  int i = 0, j = 0, ctr, file_number = 0;

  int lower_hz_level = 0, higher_hz_level = 0;
  int lower_level_low_layout = 0, higher_level_low_layout = 0;
  int lower_level_higher_layout = 0, higher_level_higher_layout = 0;

  lower_hz_level = hz_level_from;
  higher_hz_level = hz_level_to;
  ret_code = PIDX_blocks_initialize_layout(block_layout, lower_hz_level, higher_hz_level, file->idx->maxh, file->idx->bits_per_block);
  if (ret_code != PIDX_success)
  {
    fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (start_layout_index == 0)
  {
    lower_level_low_layout = 0;
    higher_level_low_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;

    if (higher_level_low_layout >= higher_hz_level)
      higher_level_low_layout = higher_hz_level;

    ret_code = PIDX_blocks_initialize_layout(layout_by_level[0], lower_level_low_layout, higher_level_low_layout, file->idx->maxh, file->idx->bits_per_block);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    ret_code = populate_idx_layout_serial(file, si, layout_by_level[0], lower_level_low_layout, higher_level_low_layout);
    if (ret_code != PIDX_success)
    {
      fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (j = lower_hz_level ; j < file->idx->bits_per_block + 1 ; j++)
      memcpy(block_layout->hz_block_number_array[j], layout_by_level[0]->hz_block_number_array[j], sizeof(int));

    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
    {
      memcpy(block_layout->hz_block_number_array[j], layout_by_level[0]->hz_block_number_array[j], sizeof(int) * ctr);
      ctr = ctr * 2;
    }

    for (i = 1; i < layout_count; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout, file->idx->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = populate_idx_layout_serial(file, si, layout_by_level[i], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], layout_by_level[i]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    int temp_level = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1;
    if (temp_level >= higher_hz_level)
      temp_level = higher_hz_level;
    for (j = file->idx->bits_per_block + 1 ; j < temp_level ; j++)
      ctr = ctr * 2;

    ctr = (int)pow(2, start_layout_index - 1) * file->idx->blocks_per_file;
    for (i = start_layout_index; i < end_layout_index; i++)
    {
      lower_level_higher_layout = file->idx->bits_per_block + log2(file->idx->blocks_per_file) + 1 + (i - 1);
      higher_level_higher_layout = lower_level_higher_layout + 1;

      ret_code = PIDX_blocks_initialize_layout(layout_by_level[i - start_layout_index], lower_level_higher_layout, higher_level_higher_layout, file->idx->maxh, file->idx->bits_per_block);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      ret_code = populate_idx_layout_serial(file, si, layout_by_level[i - start_layout_index], lower_level_higher_layout, higher_level_higher_layout);
      if (ret_code != PIDX_success)
      {
        fprintf(stderr, "[%s] [%d ]Error in populate_idx_layout\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      memcpy(block_layout->hz_block_number_array[lower_level_higher_layout], layout_by_level[i - start_layout_index]->hz_block_number_array[lower_level_higher_layout], sizeof(int) * ctr);
      ctr = ctr * 2;
    }
  }

  block_layout->file_bitmap = malloc(file->idx->max_file_count * sizeof (int));
  memset(block_layout->file_bitmap, 0, file->idx->max_file_count * sizeof (int));

  block_layout->file_index = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->file_index, 0, sizeof(int) * (file->idx->max_file_count));

  block_layout->bcpf = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->bcpf, 0, sizeof(int) * (file->idx->max_file_count));

  block_layout->lbi = malloc(sizeof(int) * (file->idx->max_file_count));
  memset(block_layout->lbi, 0, sizeof(int) * (file->idx->max_file_count));

  if (block_layout->resolution_from <=   file->idx->bits_per_block)
  {
    for (i = block_layout->resolution_from ; i <= file->idx->bits_per_block ; i++)
    {
      if (block_layout->hz_block_number_array[i][0] == 0)
      {
        file_number = block_layout->hz_block_number_array[i][0] / file->idx->blocks_per_file;
        block_layout->file_bitmap[file_number] = 1;
        file->idx_b->block_bitmap[file_number][block_layout->hz_block_number_array[i][0] % file->idx->blocks_per_file] = 1;
        block_layout->file_index[file_number] = 1;
        block_layout->bcpf[file_number]++;
        break;
      }
    }

    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      for (j = 0; j < ctr; j++)
      {
        if (block_layout->hz_block_number_array[i][j] != 0)
        {
          file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
          block_layout->file_bitmap[file_number] = 1;
          file->idx_b->block_bitmap[file_number][block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file] = 1;
          block_layout->file_index[file_number] = 1;
          block_layout->bcpf[file_number]++;
          block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
        }
      }
      ctr = ctr * 2;
    }
  }
  else
  {
    ctr = 1;
    for (i =   file->idx->bits_per_block + 1 ; i < block_layout->resolution_to ; i++)
    {
      if (i >= block_layout->resolution_from)
      {
        for (j = 0; j < ctr; j++)
        {
          if (block_layout->hz_block_number_array[i][j] != 0)
          {
            file_number = block_layout->hz_block_number_array[i][j] / file->idx->blocks_per_file;
            block_layout->file_bitmap[file_number] = 1;
            file->idx_b->block_bitmap[file_number][block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file] = 1;
            block_layout->file_index[file_number] = 1;
            block_layout->bcpf[file_number]++;
            block_layout->lbi[file_number] = block_layout->hz_block_number_array[i][j] % file->idx->blocks_per_file;
          }
        }
      }
      ctr = ctr * 2;
    }
  }

  block_layout->efc = 0;
  for (i = 0; i < file->idx->max_file_count; i++)
    if (block_layout->file_index[i] == 1)
      block_layout->efc++;

  block_layout->existing_file_index = (int*) malloc(block_layout->efc * sizeof (int));
  memset(block_layout->existing_file_index, 0, block_layout->efc * sizeof (int));

  block_layout->inverse_existing_file_index = (int*) malloc(file->idx->max_file_count * sizeof (int));
  memset(block_layout->inverse_existing_file_index, 0, file->idx->max_file_count * sizeof (int));

  int count = 0;
  for (i = 0; i < file->idx->max_file_count; i++)
  {
    if (block_layout->file_index[i] == 1)
    {
      block_layout->existing_file_index[count] = i;
      block_layout->inverse_existing_file_index[i] = count;
      count++;
    }
  }

  //if (file->idx_c->simulation_rank == 0)
  //  PIDX_blocks_print_layout(block_layout, file->idx->bits_per_block);

  return PIDX_success;
}


static PIDX_return_code delete_sim_block_layout(PIDX_io file)
{
  int i;

  PIDX_free_layout(file->idx_b->block_layout);
  PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, file->idx_b->block_layout);

  for (i = 0; i < file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count ; i++)
  {
    PIDX_free_layout(file->idx_b->block_layout_by_agg_group[i]);
    PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, file->idx_b->block_layout_by_agg_group[i]);
  }
  destroy_block_layout(file);

  int v = 0;
  for (v = 0; v < file->idx->variable_count; v++)
  {
    for (i = 0; i < file->idx->max_file_count; i++)
      free(file->idx_b->block_offset_bitmap[v][i]);
    free(file->idx_b->block_offset_bitmap[v]);
  }
  free(file->idx_b->block_offset_bitmap);

  for (i = 0; i < file->idx->max_file_count; i++)
    free(file->idx_b->block_bitmap[i]);
  free(file->idx_b->block_bitmap);

  return PIDX_success;
}

static PIDX_return_code destroy_block_layout(PIDX_io file)
{
  int i = 0;

  for (i = 0; i < file->idx_b->file0_agg_group_count + file->idx_b->nfile0_agg_group_count; i++)
    free(file->idx_b->block_layout_by_agg_group[i]);

  free(file->idx_b->block_layout);
  free(file->idx_b->block_layout_by_agg_group);

  return PIDX_success;
}
