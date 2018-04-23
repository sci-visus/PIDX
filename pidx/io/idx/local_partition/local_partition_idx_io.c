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
static PIDX_return_code partition(PIDX_io file, int svi, int mode);
static PIDX_return_code adjust_offsets(PIDX_io file, int svi, int evi);
static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int svi, int evi, int mode);
static PIDX_return_code create_midx(PIDX_io file, int svi);
static PIDX_return_code parse_local_partition_idx_file(PIDX_io file, int partition_index);
static PIDX_return_code read_block(PIDX_io file, int vi, int p, int block_number, uint64_t *patch_offset, uint64_t *patch_size, unsigned char* patch_buffer);

/// local Partitioned IDX Write Steps
/*********************************************************
*  Step 0:  group and IDX related meta data        *
*                            *
*  Step 1:  Restrucure setup               *
*  Step 2:  Restrucure                   *
*  Step 3:  Partition                  *
*                            *
*  Step 4:  Adjust offset                *
*                            *
*  Step 5:  Post partition group meta data         *
*                            *
*  Step 6:  Setup HZ encoding Phase            *
*  Step 7:  Perform HZ encoding              *
*  Step 8:  Setup aggregation buffers          *
*  Step 9:  Perform data aggregation           *
*  Step 10: Perform actual file IO             *
*  Step 11: cleanup for Steps 6              *
*                            *
*  Step 12: Cleanup the group and IDX related meta-data  *
*                            *
*  Step 13: Partition cleanup              *
*  Step 14: Restructuring cleanup            *
**********************************************************/

PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->time;

  // Step 1:  Restrucure setup
  set_rst_box_size_for_write(file, svi);

  if (idx_restructure_setup(file, svi, evi - 1, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2:  Restrucure
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable var0 = file->idx->variable[svi];

  // TODO make a utility print error function for all the error prints
  
  if (var0->restructured_super_patch_count == 1)
  {
    
    // Step 3:  Partition
    if (partition(file, svi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"Error [File %s Line %d]: partition\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    
    ret = populate_bit_string(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 5:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, svi, evi, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if 1
    for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 6:  Setup HZ encoding Phase
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 7: Perform HZ encoding
      if (hz_encode(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      if (hz_io(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 8: Setup aggregation buffers
      ret = data_aggregate(file, si, ei, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 9: Performs data aggregation
      ret = data_aggregate(file, si, ei, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 10: Performs actual file io
      time->io_start[li] = PIDX_get_time();

      ret = data_io(file, si, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      finalize_aggregation(file, si);
      time->io_end[li] = PIDX_get_time();


      // Step 11: Cleanup for step 6
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 12: Cleanup the group and IDX related meta-data
    ret = group_meta_data_finalize(file, svi, evi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 13: Partition cleanup
    time->partition_cleanup_start = MPI_Wtime();
    if (destroy_partition_comm(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->partition_cleanup_end = MPI_Wtime();
#endif
  }

  //free_restructured_communicators(file);

  // Step 14: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}


/// local Partitioned IDX Read Steps
/*********************************************************
*  Step 0:  group and IDX related meta data              *
*                                                        *
*  Step 1:  Restrucure setup                             *
*  Step 2:  Partition                                    *
*                                                        *
*  Step 3:  Adjust offsets                               *
*                                                        *
*  Step 4:  Post partition group meta data               *
*                                                        *
*  Step 5:  Setup HZ encoding Phase                      *
*  Step 6:  Setup aggregation buffers                    *
*  Step 7:  Perform actual file IO                       *
*  Step 8:  Perform data aggregation                     *
*  Step 9:  Perform HZ encoding                          *
*  Step 10:  cleanup for Steps 6                         *
*                                                        *
*  Step 11: Cleanup the group and IDX related meta-data  *
*                                                        *
*  Step 12: Partition cleanup                            *
*  Step 13: Restrucure                                   *
*  Step 14: Restructuring cleanup                        *
**********************************************************/

PIDX_return_code PIDX_local_partition_idx_read(PIDX_io file, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->time;

  // Step 1:  Restrucure setup

  set_rst_box_size_for_write(file, svi);

  if (idx_restructure_setup(file, svi, evi - 1, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable var0 = file->idx->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    // Step 2:  Partition
    ret = partition(file, svi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, svi, evi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 5:  Setup HZ encoding Phase
      ret = hz_encode_setup(file, si, ei);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      if (hz_io(file, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 6: Setup aggregation buffers
      ret = data_aggregate(file, si, ei, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }


      // Setup 7: Performs actual file io
      time->io_start[li] = PIDX_get_time();

      ret = data_io(file, si, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[li] = PIDX_get_time();

#if 1
      //
      // Setup 8: Performs data aggregation
      ret = data_aggregate(file, si, ei, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      finalize_aggregation(file, si);

      //

      // Step 9: Perform HZ encoding
      ret = hz_encode(file, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 10: Cleanup for step 6
      ret = hz_encode_cleanup(file);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      #endif
    }

    // Step 11: Cleanup the group and IDX related meta-data
    ret = group_meta_data_finalize(file, svi, evi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 12: Partition cleanup
    time->partition_cleanup_start = MPI_Wtime();
    if (destroy_partition_comm(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->partition_cleanup_end = MPI_Wtime();

    int i = 0;
    PIDX_variable var = file->idx->variable[svi];
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] + file->idx->partition_offset[i];

  }

#if 1
  // Step 13:  Restrucure
  if (idx_restructure(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 14: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_parallel_local_partition_idx_read(PIDX_io file, int svi, int evi)
{
  int i = 0, j = 0, p = 0;
  int si = 0;
  int ret = 0;

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  ret = populate_global_bit_string(file, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (si = svi; si < evi; si++)
  {
    PIDX_variable var = file->idx->variable[si];

    for (p = 0; p < file->idx->variable[si]->sim_patch_count; p++)
    {
      int par = 0;
      for (par = 0; par < file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]; par++)
      {
        if (parse_local_partition_idx_file(file, par) == -1) continue;
        char dirname[1024], basename[1024];
        VisusSplitFilename(file->idx->filename_template_partition, dirname, basename);
        sprintf(file->idx->filename_template_partition, "%s/time%09d/%s", dirname, file->idx->current_time_step, basename );

        //fprintf(stderr, "Partition %d ----> offset %d %d %d size %d %d %d -> %s\n", par, file->idx->partition_offset[0], file->idx->partition_offset[1], file->idx->partition_offset[2], file->idx->partition_size[0], file->idx->partition_size[1], file->idx->partition_size[2], file->idx->filename_template_partition);

        int d = 0, check_bit = 0;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          check_bit = check_bit || (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) < file->idx->partition_offset[d] || (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1) < var->sim_patch[p]->offset[d];

        if (!check_bit)
        {
          uint64_t intersected_box_offset[PIDX_MAX_DIMENSIONS];
          uint64_t intersected_box_size[PIDX_MAX_DIMENSIONS];

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            if (var->sim_patch[p]->offset[d] <= file->idx->partition_offset[d] && (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) <= (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1))
            {
              intersected_box_offset[d] = file->idx->partition_offset[d];
              intersected_box_size[d] = (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) - file->idx->partition_offset[d] + 1;
            }
            else if (file->idx->partition_offset[d] <= var->sim_patch[p]->offset[d] && (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) >= (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1))
            {
              intersected_box_offset[d] = var->sim_patch[p]->offset[d];
              intersected_box_size[d] = (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1) - var->sim_patch[p]->offset[d] + 1;
            }
            else if (( file->idx->partition_offset[d] + file->idx->partition_size[d] - 1) <= (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) && file->idx->partition_offset[d] >= var->sim_patch[p]->offset[d])
            {
              intersected_box_offset[d] = file->idx->partition_offset[d];
              intersected_box_size[d] = file->idx->partition_size[d];
            }
            else if (( var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) <= (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1) && var->sim_patch[p]->offset[d] >= file->idx->partition_offset[d])
            {
              intersected_box_offset[d] = var->sim_patch[p]->offset[d];
              intersected_box_size[d] = var->sim_patch[p]->size[d];
            }

            intersected_box_offset[d] = intersected_box_offset[d] - file->idx->partition_offset[d];
          }

          for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
          {
            bounding_box[0][i] = intersected_box_offset[i];
            bounding_box[1][i] = intersected_box_offset[i] + intersected_box_size[i];
          }
          int bytes_for_datatype = ((file->idx->variable[si]->bpv / 8) * file->idx->variable[si]->vps);
          unsigned char* intersected_box_buffer = malloc(intersected_box_size[0] * intersected_box_size[1] * intersected_box_size[2] * bytes_for_datatype);

          PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
          memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
          ret = PIDX_blocks_initialize_layout(per_patch_local_block_layout, 0, file->idx->maxh, file->idx->maxh, file->idx->bits_per_block);
          if (ret != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          ret = PIDX_blocks_create_layout (bounding_box, file->idx->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_b->reduced_res_from, file->idx_b->reduced_res_to);
          if (ret != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          //PIDX_blocks_print_layout(per_patch_local_block_layout, file->idx->bits_per_block);
          int ctr = 1;
          int block_number = 0;
          read_block(file, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer);
          for (i = file->idx->bits_per_block + 1 ; i < per_patch_local_block_layout->resolution_to ; i++)
          {
            for (j = 0 ; j < ctr ; j++)
            {
              if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
              {
                block_number = per_patch_local_block_layout->hz_block_number_array[i][j];
                read_block(file, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer);
              }
            }
            ctr = ctr * 2;
          }

          PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, per_patch_local_block_layout);
          free(per_patch_local_block_layout);

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            intersected_box_offset[d] = intersected_box_offset[d] + file->idx->partition_offset[d];

          int k1, j1, i1, index = 0, recv_o = 0, send_o = 0, send_c = 0;
          for (k1 = intersected_box_offset[2]; k1 < intersected_box_offset[2] + intersected_box_size[2]; k1++)
          {
            for (j1 = intersected_box_offset[1]; j1 < intersected_box_offset[1] + intersected_box_size[1]; j1++)
            {
              for (i1 = intersected_box_offset[0]; i1 < intersected_box_offset[0] + intersected_box_size[0]; i1 = i1 + intersected_box_size[0])
              {
                index = ((intersected_box_size[0])* (intersected_box_size[1]) * (k1 - intersected_box_offset[2])) + ((intersected_box_size[0]) * (j1 - intersected_box_offset[1])) + (i1 - intersected_box_offset[0]);
                send_o = index * var->vps * (var->bpv/8);
                send_c = (intersected_box_size[0]);
                recv_o = (file->idx->variable[si]->sim_patch[p]->size[0] * file->idx->variable[si]->sim_patch[p]->size[1] * (k1 - file->idx->variable[si]->sim_patch[p]->offset[2])) + (file->idx->variable[si]->sim_patch[p]->size[0] * (j1 - file->idx->variable[si]->sim_patch[p]->offset[1])) + (i1 - file->idx->variable[si]->sim_patch[p]->offset[0]);

                memcpy(file->idx->variable[si]->sim_patch[p]->buffer + (recv_o * var->vps * (var->bpv/8)), intersected_box_buffer + send_o, send_c * var->vps * (var->bpv/8));
              }
            }
          }
          free(intersected_box_buffer);
        }
      }
    }
  }

  return PIDX_success;
}


static PIDX_return_code parse_local_partition_idx_file(PIDX_io file, int partition_index)
{

    char file_name_skeleton[1024];
    strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
    file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';
    sprintf(file->idx->filename_partition, "%s_%d.idx", file_name_skeleton, partition_index);

    char *pch;
    int count = 0;
    char line [ 512 ];

    FILE *fp = fopen(file->idx->filename_partition, "r");
    if (fp == NULL)  return -1;


    while (fgets(line, sizeof (line), fp) != NULL)
    {
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "(partition index)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        file->idx_c->color = atoi(line);
      }

      if (strcmp(line, "(partition size)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          file->idx->partition_size[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
      }

      if (strcmp(line, "(bitsperblock)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        file->idx->bits_per_block = atoi(line);
        file->idx->samples_per_block = (int)pow(2, file->idx->bits_per_block);
      }

      if (strcmp(line, "(bits)") == 0)
      {
        int i = 0;
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        strcpy(file->idx->bitSequence, line);
        file->idx->maxh = strlen(file->idx->bitSequence);
        for (i = 0; i <= file->idx->maxh; i++)
          file->idx->bitPattern[i] = RegExBitmaskBit(file->idx->bitSequence, i);
        //fprintf(stderr, "BS %s MH %d\n", file->idx->bitSequence, file->idx->maxh);
      }

      if (strcmp(line, "(filename_template)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        memset(file->idx->filename_template_partition, 0, 1024 * sizeof(char));
        strcpy(file->idx->filename_template_partition, line);
      }

      if (strcmp(line, "(partition offset)") == 0)
      {
        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          file->idx->partition_offset[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
      }
    }
    fclose(fp);

    return PIDX_success;
}


static PIDX_return_code read_block(PIDX_io file, int vi, int p, int block_number, uint64_t* patch_offset, uint64_t* patch_size, unsigned char* patch_buffer)
{
  int k = 0, ret = 0;
  MPI_File fp = 0;
  MPI_Status status;
  uint64_t index = 0;
  uint64_t hz;
  uint64_t xyz[PIDX_MAX_DIMENSIONS];
  char file_name[PATH_MAX];

  int file_number = block_number / file->idx->blocks_per_file;

  int bytes_for_datatype = ((file->idx->variable[vi]->bpv / 8) * file->idx->variable[vi]->vps);

  unsigned char* block_buffer = malloc(file->idx->samples_per_block * bytes_for_datatype);

  //generate_file_name_template(file->idx->maxh, file->idx->bits_per_block, file->idx->filename_partition, file->idx->current_time_step, file->idx->filename_template_partition);

  if (generate_file_name(file->idx->blocks_per_file, file->idx->filename_template_partition, file_number, file_name, PATH_MAX) == 1)
  {
    fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  if (MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_open() block number %d file number %d filename %s failed.\n", __FILE__, __LINE__, block_number, file_number, file_name);
    return PIDX_err_io;
  }

  uint32_t *headers;
  int total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  headers = malloc(total_header_size);
  memset(headers, 0, total_header_size);

  ret = MPI_File_read_at(fp, 0, headers, total_header_size , MPI_BYTE, &status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Data offset = [%s] [%d] MPI_File_write_at() failed for filename %s.\n", __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }
  int read_count = 0;
  MPI_Get_count(&status, MPI_BYTE, &read_count);
  if (read_count != total_header_size)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed. %d != %dd\n", __FILE__, __LINE__, read_count, total_header_size);
    return PIDX_err_io;
  }

  uint64_t data_offset = htonl(headers[12 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);
  uint64_t data_size = htonl(headers[14 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);

  //fprintf(stderr, "File number %d Block number %d\n", file_number, block_number);
  assert (data_size != 0);

  ret = MPI_File_read_at(fp, data_offset, block_buffer, data_size, MPI_BYTE, &status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }

  for (k = 0; k < file->idx->samples_per_block; k++)
  {
    hz = (block_number * file->idx->samples_per_block) + k;
    Hz_to_xyz(file->idx->bitPattern, file->idx->maxh, hz, xyz);

    if ( ((xyz[0] < patch_offset[0] || xyz[0] >= patch_offset[0] + patch_size[0]) ||
        (xyz[1] < patch_offset[1] || xyz[1] >= patch_offset[1] + patch_size[1]) ||
        (xyz[2] < patch_offset[2] || xyz[2] >= patch_offset[2] + patch_size[2]) ) )
      continue;

    index = (patch_size[0] * patch_size[1] * (xyz[2] - patch_offset[2]))
        + (patch_size[0] * (xyz[1] - patch_offset[1]))
        + (xyz[0] - patch_offset[0]);

    memcpy(patch_buffer + (index * bytes_for_datatype),
         block_buffer + (k * bytes_for_datatype),
         bytes_for_datatype);

  }

  MPI_File_close(&fp);

  free(headers);
  free(block_buffer);

  return PIDX_success;
}


// TODO move this function into some utils file, it is used in too many files already
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];
  
  return !(check_bit);
}

static PIDX_return_code create_midx(PIDX_io file, int svi)
{
  int *colors;
  int i = 0, j = 0, k = 0, d = 0;
  
  int num_parts = file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2];
  
  colors = malloc(sizeof(*colors) * num_parts);
  memset(colors, 0, sizeof(*colors) * num_parts);

    for (k = 0; k < file->idx->partition_count[2]; k++)
      for (j = 0; j < file->idx->partition_count[1]; j++)
        for (i = 0; i < file->idx->partition_count[0]; i++)
        {
          colors[(file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i] = (file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i;
          
        }
  
    PIDX_variable var = file->idx->variable[svi];

    PIDX_patch local_p = (PIDX_patch)malloc(sizeof (*local_p));
    memset(local_p, 0, sizeof (*local_p));
    
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_p->offset[d] = var->restructured_super_patch->restructured_patch->offset[d];
      local_p->size[d] = var->restructured_super_patch->restructured_patch->size[d];
    }
    
    int wc = file->idx_c->partition_nprocs * (PIDX_MAX_DIMENSIONS);
  
    uint64_t* global_patch_offset = malloc(wc * sizeof(*global_patch_offset));
    memset(global_patch_offset, 0, wc * sizeof(*global_patch_offset));
    
    uint64_t* global_patch_size = malloc(wc * sizeof(*global_patch_size));
    memset(global_patch_size, 0, wc * sizeof(*global_patch_size));
    
    MPI_Allgather(&local_p->offset[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->partition_comm);
    
    MPI_Allgather(&local_p->size[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->partition_comm);
  
    //free(local_p);
  
  if (file->idx_c->simulation_rank == 0)
  {
    // TODO make utility functions to retrieve filenames
    // and use those in all the code for all filename template related things
    char file_name_skeleton[PIDX_STRING_SIZE];
    strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
    file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';
    
    char* partition_filenames = malloc(PIDX_STRING_SIZE * num_parts);
    memset(partition_filenames, 0, PIDX_STRING_SIZE * num_parts);
    int* offsets = malloc(sizeof(*offsets) * PIDX_MAX_DIMENSIONS * num_parts);
    memset(offsets, 0, sizeof(*offsets) * PIDX_MAX_DIMENSIONS * num_parts);
    
    int proc=0;
    for (proc = 0; proc < file->idx_c->partition_nprocs; proc++)
    {
      PIDX_patch curr_local_p = (PIDX_patch)malloc(sizeof (*curr_local_p));
      
      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        curr_local_p->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * proc + d];
        curr_local_p->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * proc + d];
      }
      
      PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
      memset(reg_patch, 0, sizeof (*reg_patch));
      
      Point3D bounds_point;
      int maxH = 0;
      bounds_point.x = (int) file->idx->partition_count[0];
      bounds_point.y = (int) file->idx->partition_count[1];
      bounds_point.z = (int) file->idx->partition_count[2];
      char bitSequence[512];
      char bitPattern[512];
      GuessBitmaskPattern(bitSequence, bounds_point);
      maxH = strlen(bitSequence);
      
      for (i = 0; i <= maxH; i++)
        bitPattern[i] = RegExBitmaskBit(bitSequence, i);
      
      int z_order = 0;
      int number_levels = maxH - 1;
      int index_i = 0, index_j = 0, index_k = 0;
      for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx->partition_size[0], index_i++)
      {
        for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx->partition_size[1], index_j++)
        {
          for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx->partition_size[2], index_k++)
          {
            reg_patch->offset[0] = i;
            reg_patch->offset[1] = j;
            reg_patch->offset[2] = k;
            reg_patch->size[0] = file->idx->partition_size[0];
            reg_patch->size[1] = file->idx->partition_size[1];
            reg_patch->size[2] = file->idx->partition_size[2];
            
            if (intersectNDChunk(reg_patch, curr_local_p))
            {
              Point3D xyzuv_Index;
              xyzuv_Index.x = index_i;
              xyzuv_Index.y = index_j;
              xyzuv_Index.z = index_k;
              
              z_order = 0;
              Point3D zero;
              zero.x = 0;
              zero.y = 0;
              zero.z = 0;
              memset(&zero, 0, sizeof (Point3D));
              
              int cnt = 0;
              for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
              {
                int bit = bitPattern[number_levels];
                z_order |= ((uint64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                PGET(xyzuv_Index, bit) >>= 1;
              }
     
              //printf("found color %d offset %d %d %d\n", colors[z_order], i, j, k);
            
              sprintf(&partition_filenames[colors[z_order]*PIDX_STRING_SIZE], "%s_%d.idx", file_name_skeleton, colors[z_order]);
              
              uint64_t curr_off = colors[z_order]*PIDX_MAX_DIMENSIONS;
              offsets[curr_off+0] = i;
              offsets[curr_off+1] = j;
              offsets[curr_off+2] = k;
              
              break;
            }
          }
        }
      }
      
      free(curr_local_p);
      free(reg_patch);
      
    }
    
    char midx_filename[PIDX_STRING_SIZE];
    sprintf(midx_filename,"%s.midx", file_name_skeleton);
    
    FILE* midx_file = fopen(midx_filename, "w");
    if (!midx_file)
    {
      fprintf(stderr, "Error [%s] [%d]: Cannot create filename %s.\n", __FILE__, __LINE__, midx_filename);
      return -1;
    }
    
    fprintf(midx_file, "<dataset typename='IdxMultipleDataset'>\n");
    
    for (i = 0; i < num_parts; i++)
    {
      uint64_t curr_off = i*PIDX_MAX_DIMENSIONS;
      
      //printf("filename %s offset  %d %d %d\n", partition_filenames+(PIDX_STRING_SIZE*i), offsets[curr_off+0],offsets[curr_off+1], offsets[curr_off+2]);
      
      fprintf(midx_file, "\t<dataset url=\"file://$(CurrentFileDirectory)/%s\" name=\"%s_%d\" offset=\"%d %d %d\"/>\n",
              partition_filenames+(PIDX_STRING_SIZE*i), file_name_skeleton, i, offsets[curr_off+0],offsets[curr_off+1], offsets[curr_off+2]);
    }
    
    fprintf(midx_file,"</dataset>\n");
    
    fclose(midx_file);
    
    free(colors);
    free(global_patch_size);
    free(global_patch_offset);
  
  }
  
  return PIDX_success;
}


static PIDX_return_code partition(PIDX_io file, int svi, int mode)
{
  int ret;
  PIDX_time time = file->time;

  time->partition_start = MPI_Wtime();
  // Calculates the number of partititons
  ret = find_partition_count(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  
  // assign same color to processes within the same partition
  ret = partition_setup(file, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Create midx file that will point to different .idx partition files
  if (create_midx(file, svi) != PIDX_success)
  {
    fprintf(stderr,"Error [File %s Line %d]: create_midx\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Splits the local communicator into local communicators
  ret = create_partition_comm(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_end = MPI_Wtime();

  return PIDX_success;
}



static PIDX_return_code adjust_offsets(PIDX_io file, int svi, int evi)
{
  int i = 0, v = 0;

  //fprintf(stderr, "SE %d %d\n", svi, evi);
  for (v = svi; v < evi; v++)
  {
  PIDX_variable var = file->idx->variable[v];

  if (var->restructured_super_patch_count != 1)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] - file->idx->partition_offset[i];
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx->partition_offset[i] + file->idx->partition_size[i] <= file->idx->box_bounds[i])
      file->idx->box_bounds[i] = file->idx->partition_size[i];
    else
      file->idx->box_bounds[i] = file->idx->box_bounds[i] - file->idx->partition_offset[i];

    if (file->restructured_grid->patch_size[i] > file->idx->box_bounds[i])
      file->restructured_grid->patch_size[i] = getPowerOf2(file->idx->box_bounds[i]);
  }

  memcpy(file->idx->bounds, file->idx->box_bounds, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  //fprintf(stderr, "%d - %d %d %d -- PO %d %d %d\n", file->idx_c->simulation_rank, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], file->idx->partition_offset[0], file->idx->partition_offset[1], file->idx->partition_offset[2]);

  return PIDX_success;
}



static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->time;

  //fprintf(stderr, "%d %d %d\n", file->restructured_grid->patch_size[0], file->restructured_grid->patch_size[1], file->restructured_grid->patch_size[2]);

  time->bit_string_start = PIDX_get_time();
  // calculates maxh and bitstring
  ret = populate_local_bit_string(file, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // selects levels based on mode and maxh
  select_io_mode(file);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layoutven this is pure IDX only non-share block layout is populated
  //if (file->idx_d->color == 2)
  //{
  ret = populate_rst_block_layouts(file, svi, file->idx_b->hz_file0_from, file->idx_b->hz_n_file0_to);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  //}
#if 1
  // Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Creates the agg and io ids
  ret = create_agg_io_buffer(file);
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
#endif
  return PIDX_success;
}


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int svi, int evi)
{
  int ret;
  PIDX_time time = file->time;

  time->group_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = delete_rst_block_layout(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
