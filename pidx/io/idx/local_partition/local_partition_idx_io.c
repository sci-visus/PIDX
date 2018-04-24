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

static PIDX_return_code parse_local_partition_idx_file(PIDX_io file, int partition_index);
static PIDX_return_code read_block(PIDX_io file, int vi, int p, int block_number, uint64_t *patch_offset, uint64_t *patch_size, unsigned char* patch_buffer);

/// local Partitioned IDX Write Steps
PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int svi, int evi)
{
  PIDX_return_code ret;
  PIDX_time time = file->time;

  // Step 1: Compute the restructuring box (grid) size
  // We need to identify the restructuring box size (super patch first because we directly use that to populate the
  // bitstring, which is first written out to the .idx file and used in almost every phase of IO
  if (set_rst_box_size_for_write(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 2: Setting the stage for restructuring (meta data)
  if (idx_restructure_setup(file, svi, evi - 1) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 3: Perform data restructuring
  // restructuring is done keeping partitioning in mind
  // restructuring phase ensures that a process is not holding more that one super patch
  if (idx_restructure(file, PIDX_WRITE) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 4: Group the processes holding the super patch into new communicator rst_comm
  // This step is intended to reduce MPI overhead by having a comm with fewer processes
  if (idx_restructure_rst_comm_create(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  // proceed only if a process holds a superpatch, others just wait at completion of io
  // essential for partitioning
  PIDX_variable var0 = file->idx->variable[svi];
  if (var0->restructured_super_patch_count == 1)
  {
    // Step 5:  Perform partitioning
    // Create new communicator
    if (partition(file, svi) != PIDX_success)
    {
      fprintf(stderr,"Error [File %s Line %d]: partition\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    
    // Step 6: Populate the bit string of the global dataset
    ret = populate_bit_string(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 7: write the global .idx file
    if (write_global_idx(file, svi, evi, PIDX_WRITE) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 8: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 9:  Popukate block layout for the partitions (this is a local step, from now onwards the processes work within its partition only)
    if (populate_block_layout_and_buffers(file, svi, evi, PIDX_WRITE, PIDX_LOCAL_PARTITION_IDX_IO) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Iterate through the variables, variable_pipe_length at a point
    // variable_pipe_length is computed based on the configuration of the run. If there are enough number of processes
    // then all the variables are aggregated at once and then variable_pipe_length is (variable_count - 1), otherwise
    // the variables are worked in smaller packs.
    for (uint32_t si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      uint32_t ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_tracker[si] = 1;

      // Step 10:  Setup HZ encoding Phase
      if (hz_encode_setup(file, si, ei) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 11: Perform HZ encoding
      if (hz_encode(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 12: This is not performed by default, it only happens when aggregation is
      // turned off or when there are a limited number of aggregators
      if (hz_io(file, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 13: Setup aggregation buffers
      ret = aggregation_setup(file, si, ei);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 14: Performs data aggregation
      ret = aggregation(file, si, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 15: Performs actual file io
      ret = file_io(file, si, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 16: free aggregation buffers
      if (aggregation_cleanup(file, si) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 17: Cleanup HZ buffers
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 18: free block layout and agg related buffer
    ret = destroy_block_layout_and_buffers(file, svi, evi);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 19: free partitioning buffer
    time->partition_cleanup_start = MPI_Wtime();
    if (partition_cleanup(file) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->partition_cleanup_end = MPI_Wtime();
  }

  // Step 20: free the restructured communicator
  if (free_restructured_communicators(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Step 21: Restructuring cleanup
  if (idx_restructure_cleanup(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_parallel_local_partition_idx_read(PIDX_io file, int svi, int evi)
{
  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  if (populate_global_bit_string(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (uint32_t si = svi; si < evi; si++)
  {
    PIDX_variable var = file->idx->variable[si];
    for (uint32_t p = 0; p < file->idx->variable[si]->sim_patch_count; p++)
    {
      for (uint32_t par = 0; par < file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]; par++)
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

          for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
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

          for (uint32_t i = 0; i < PIDX_MAX_DIMENSIONS; i++)
          {
            bounding_box[0][i] = intersected_box_offset[i];
            bounding_box[1][i] = intersected_box_offset[i] + intersected_box_size[i];
          }
          int bytes_for_datatype = ((file->idx->variable[si]->bpv / 8) * file->idx->variable[si]->vps);
          unsigned char* intersected_box_buffer = malloc(intersected_box_size[0] * intersected_box_size[1] * intersected_box_size[2] * bytes_for_datatype);

          PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
          memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
          if (PIDX_blocks_initialize_layout(per_patch_local_block_layout, 0, file->idx->maxh, file->idx->maxh, file->idx->bits_per_block) != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          if (PIDX_blocks_create_layout (bounding_box, file->idx->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_b->reduced_res_from, file->idx_b->reduced_res_to) != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          //PIDX_blocks_print_layout(per_patch_local_block_layout, file->idx->bits_per_block);
          int ctr = 1;
          int block_number = 0;
          read_block(file, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer);
          for (uint32_t i = file->idx->bits_per_block + 1 ; i < per_patch_local_block_layout->resolution_to ; i++)
          {
            for (uint32_t j = 0 ; j < ctr ; j++)
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

          for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            intersected_box_offset[d] = intersected_box_offset[d] + file->idx->partition_offset[d];

          uint32_t index = 0, recv_o = 0, send_o = 0, send_c = 0;
          for (uint32_t k1 = intersected_box_offset[2]; k1 < intersected_box_offset[2] + intersected_box_size[2]; k1++)
          {
            for (uint32_t j1 = intersected_box_offset[1]; j1 < intersected_box_offset[1] + intersected_box_size[1]; j1++)
            {
              for (uint32_t i1 = intersected_box_offset[0]; i1 < intersected_box_offset[0] + intersected_box_size[0]; i1 = i1 + intersected_box_size[0])
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





