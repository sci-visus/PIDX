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


PIDX_return_code PIDX_local_partition_idx_generic_read(PIDX_io file, int svi, int evi)
{
  // We iterate through all the patches, for each of the patches
  // iterate through all the partitions and see if there is an intersection
  // If an intersection is found, find the intersection bounding box
  // and make box query in that particular partition.
  // The box query is performed by finding out which idx block are present
  // corresponding to the intersection box, make those reads, and
  // iterate through the elements of the block and if the sample is in the
  // box, copy it to the patch buffer.


  // Use this function to compute maxh and also the maximum number of files
  if (populate_global_bit_string(file, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (uint32_t si = svi; si < evi; si++)
  {
    PIDX_variable var = file->idx->variable[si];

    // Iterate through all the patches
    for (uint32_t p = 0; p < file->idx->variable[si]->sim_patch_count; p++)
    {
      // for every patch iterate through all the partitions
      for (uint32_t par = 0; par < file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]; par++)
      {
        // If the particular partition does not exist then move to the next partition
        if (parse_local_partition_idx_file(file, par) != PIDX_success)
            continue;

        // populate the local partition template using the current time step index
        char dirname[1024], basename[1024];
        VisusSplitFilename(file->idx->filename_template_partition, dirname, basename);
        char time_template[512];
        sprintf(time_template, "%%s/%s/%%s", file->idx->filename_time_template);
        sprintf(file->idx->filename_template_partition, time_template, dirname, file->idx->current_time_step, basename );

        // Checking if the simulation patch intersects with the partition
        int d = 0, check_bit = 0;
        for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          check_bit = check_bit || (var->sim_patch[p]->offset[d] + var->sim_patch[p]->size[d] - 1) < file->idx->partition_offset[d] || (file->idx->partition_offset[d] + file->idx->partition_size[d] - 1) < var->sim_patch[p]->offset[d];

        // if the patch intersects with the partition
        if (!check_bit)
        {
          // find intersection bounding box
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

          // intersection bounding box
          int bounding_box[2][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
          for (uint32_t i = 0; i < PIDX_MAX_DIMENSIONS; i++)
          {
            bounding_box[0][i] = intersected_box_offset[i];
            bounding_box[1][i] = intersected_box_offset[i] + intersected_box_size[i];
          }

          // allocate buffer to hold data corresponding to the intersection bounding box
          int bytes_for_datatype = ((file->idx->variable[si]->bpv / 8) * file->idx->variable[si]->vps);
          unsigned char* intersected_box_buffer = malloc(intersected_box_size[0] * intersected_box_size[1] * intersected_box_size[2] * bytes_for_datatype);

          // For the intersection bounding box, find out what all idx box to query
          PIDX_block_layout per_patch_local_block_layout = malloc(sizeof (*per_patch_local_block_layout));
          memset(per_patch_local_block_layout, 0, sizeof (*per_patch_local_block_layout));
          if (PIDX_blocks_initialize_layout(per_patch_local_block_layout, 0, file->idx->maxh, file->idx->maxh, file->idx->bits_per_block) != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_initialize_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }

          if (PIDX_blocks_create_layout (bounding_box, file->idx->maxh, file->idx->bits_per_block,  file->idx->bitPattern, per_patch_local_block_layout, file->idx_b->reduced_resolution_factor) != PIDX_success)
          {
            fprintf(stderr, "[%s] [%d ]Error in PIDX_blocks_create_layout", __FILE__, __LINE__);
            return PIDX_err_file;
          }


          uint32_t ctr = 1;
          uint32_t block_number = 0;

          // read the first block, the first block contains data from a lot of hz levels, so just make one accumulated read
          if (read_block(file, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer) != PIDX_success)
          {
            fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
            return PIDX_err_io;
          }

          // Iterate through the blocks from HZ level file->idx->bits_per_block + 1 to per_patch_local_block_layout->resolution_to
          for (uint32_t i = file->idx->bits_per_block + 1 ; i < per_patch_local_block_layout->resolution_to ; i++)
          {
            for (uint32_t j = 0 ; j < ctr ; j++)
            {
              if (per_patch_local_block_layout->hz_block_number_array[i][j] != 0)
              {
                block_number = per_patch_local_block_layout->hz_block_number_array[i][j];

                // Iterating through all other blocks
                if (read_block(file, si, p, block_number, intersected_box_offset, intersected_box_size, intersected_box_buffer) != PIDX_success)
                {
                  fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
                  return PIDX_err_io;
                }
              }
            }
            ctr = ctr * 2;
          }

          // With the reads completed, free the block bitmap
          PIDX_blocks_free_layout(file->idx->bits_per_block, file->idx->maxh, per_patch_local_block_layout);
          free(per_patch_local_block_layout);

          // Adjust the intersection box to the global index space
          for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            intersected_box_offset[d] = intersected_box_offset[d] + file->idx->partition_offset[d];

          // copy the data from the intersected bounding box the the application patch
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
  // Parse the partition idx file to get partition specific meta data

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
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      strcpy(file->idx->bitSequence, line);
      file->idx->maxh = strlen(file->idx->bitSequence);
      for (uint32_t i = 0; i <= file->idx->maxh; i++)
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
  MPI_File fp = 0;
  MPI_Status status;

  // populate the name of the binary file to read
  char file_name[PATH_MAX];
  int file_number = block_number / file->idx->blocks_per_file;
  if (generate_file_name(file->idx->blocks_per_file, file->idx->filename_template_partition, file_number, file_name, PATH_MAX) == 1)
  {
    fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }
  
  char *directory_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  char full_path_file_name[PATH_MAX];
  char *lastdir = strrchr(file->idx->filename, '/');
  if (lastdir != NULL && file_name[0] == '.') { // if using relative paths use absolute path
    strncpy(directory_path, file->idx->filename, lastdir - file->idx->filename + 1);
    sprintf(full_path_file_name, "%s/%s", directory_path,file_name);
  }
  else{ 
    sprintf(full_path_file_name, "%s", file_name);
  }

  // open the binary file
  if (MPI_File_open(MPI_COMM_SELF, full_path_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
  {
    fprintf(stderr, "[%s] [%d] MPI_File_open() block number %d file number %d filename %s failed.\n", __FILE__, __LINE__, block_number, file_number, full_path_file_name);
    return PIDX_err_io;
  }

  // read the header
  uint32_t *headers;
  int total_header_size = (10 + (10 * file->idx->blocks_per_file)) * sizeof (uint32_t) * file->idx->variable_count;
  headers = malloc(total_header_size);
  memset(headers, 0, total_header_size);

  if (MPI_File_read_at(fp, 0, headers, total_header_size , MPI_BYTE, &status) != MPI_SUCCESS)
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

  // use the header to find the offset in file that contains data corresponding to block_number
  uint64_t data_offset = htonl(headers[12 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);
  uint64_t data_size = htonl(headers[14 + (( (block_number % file->idx->blocks_per_file) + (file->idx->blocks_per_file * vi))*10 )]);
  assert (data_size != 0);

  // read the data
  int bytes_for_datatype = ((file->idx->variable[vi]->bpv / 8) * file->idx->variable[vi]->vps);
  unsigned char* block_buffer = malloc(file->idx->samples_per_block * bytes_for_datatype);
  if (MPI_File_read_at(fp, data_offset, block_buffer, data_size, MPI_BYTE, &status) != MPI_SUCCESS)
  {
    fprintf(stderr, "Data offset = %lld [%s] [%d] MPI_File_write_at() failed for filename %s.\n", (long long)  data_offset, __FILE__, __LINE__, file_name);
    return PIDX_err_io;
  }

  // copy the data from the block space to the box space
  uint64_t xyz[PIDX_MAX_DIMENSIONS];
  for (uint64_t k = 0; k < file->idx->samples_per_block; k++)
  {
    uint64_t hz = (block_number * file->idx->samples_per_block) + k;
    Hz_to_xyz(file->idx->bitPattern, file->idx->maxh, hz, xyz);

    // check if the sample in the block is within the box query
    if ( ((xyz[0] < patch_offset[0] || xyz[0] >= patch_offset[0] + patch_size[0]) || (xyz[1] < patch_offset[1] || xyz[1] >= patch_offset[1] + patch_size[1]) || (xyz[2] < patch_offset[2] || xyz[2] >= patch_offset[2] + patch_size[2]) ) )
      continue;

    uint64_t index = (patch_size[0] * patch_size[1] * (xyz[2] - patch_offset[2])) + (patch_size[0] * (xyz[1] - patch_offset[1])) + (xyz[0] - patch_offset[0]);
    memcpy(patch_buffer + (index * bytes_for_datatype), block_buffer + (k * bytes_for_datatype), bytes_for_datatype);

  }

  MPI_File_close(&fp);

  // free buffers
  free(headers);
  free(block_buffer);
  free(directory_path);

  return PIDX_success;
}
