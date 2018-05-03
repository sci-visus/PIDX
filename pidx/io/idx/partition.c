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

static int intersectNDChunk(PIDX_patch A, PIDX_patch B);
static PIDX_return_code find_partition_size(PIDX_io file);
static PIDX_return_code partition_setup(PIDX_io file, int svi);
static PIDX_return_code create_partition_comm(PIDX_io file);
static PIDX_return_code create_midx(PIDX_io file, int svi);

PIDX_return_code partition(PIDX_io file, int svi)
{
  // 1. First find the size of the partitions using global bounds and partition count
  // 2. Associate color (unique identifier) to processes within same partition
  // 3. Create midx file (make it compatible with visus)
  // 4. Create the actual partiitoned communicator

  PIDX_time time = file->time;

  time->partition_start = MPI_Wtime();
  // Calculates the size of the partititons
  if (find_partition_size(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // assign same color to processes within the same partition
  if (partition_setup(file, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Create midx file that will point to different .idx partition files
  if (file->idx->current_time_step == 0)
  {
    if (create_midx(file, svi) != PIDX_success)
    {
      fprintf(stderr,"Error [File %s Line %d]: create_midx\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  // Splits the local communicator into local communicators
  if (create_partition_comm(file) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_end = MPI_Wtime();

  return PIDX_success;
}



PIDX_return_code adjust_offsets(PIDX_io file, int svi, int evi)
{
  // Shifting the coordinates of super patches from the global index space to local index space
  for (uint32_t v = svi; v < evi; v++)
  {
    PIDX_variable var = file->idx->variable[v];
    if (var->restructured_super_patch_count != 1)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    for (uint32_t i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] - file->idx->partition_offset[i];
  }

  // Adjusting the size of the bounding box (query), making it equal to the partition size
  for (uint32_t i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    //if (file->idx->partition_offset[i] + file->idx->partition_size[i] <= file->idx->box_bounds[i])
      file->idx->box_bounds[i] = file->idx->partition_size[i];
    //else
    //  file->idx->box_bounds[i] = file->idx->box_bounds[i] - file->idx->partition_offset[i];

    // Edge case when after partitioning the edge partitions have super patch size greater that the entire partition
    //if (file->restructured_grid->patch_size[i] > file->idx->box_bounds[i])
    //  file->restructured_grid->patch_size[i] = getPowerOf2(file->idx->box_bounds[i]);
  }

  // Also adjusting the size of entire bounding box
  memcpy(file->idx->bounds, file->idx->box_bounds, PIDX_MAX_DIMENSIONS * sizeof(uint64_t));

  return PIDX_success;
}



PIDX_return_code re_adjust_offsets(PIDX_io file, int svi, int evi)
{
  for (uint32_t v = svi; v < evi; v++)
  {
    PIDX_variable var = file->idx->variable[v];
    for (uint32_t i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] + file->idx->partition_offset[i];
  }

  return PIDX_success;
}

static PIDX_return_code partition_setup(PIDX_io file, int svi)
{
  // Function to associate color to every process
  // The color is aunique identifier donating which partition it belongs to.

  int *colors;
  colors = malloc(sizeof(*colors) * file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]);
  memset(colors, 0, sizeof(*colors) * file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]);
  file->idx_c->color = (file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2]) + 1;

  for (uint32_t k = 0; k < file->idx->partition_count[2]; k++)
    for (uint32_t j = 0; j < file->idx->partition_count[1]; j++)
      for (uint32_t i = 0; i < file->idx->partition_count[0]; i++)
      {
        colors[(file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i] = (file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i;
      }

  PIDX_patch local_p = (PIDX_patch)malloc(sizeof (*local_p));
  memset(local_p, 0, sizeof (*local_p));

  PIDX_variable var = file->idx->variable[svi];
  for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_p->offset[d] = var->restructured_super_patch->restructured_patch->offset[d];
    local_p->size[d] = var->restructured_super_patch->restructured_patch->size[d];
  }
  PIDX_patch partition_patch = (PIDX_patch)malloc(sizeof (*partition_patch));
  memset(partition_patch, 0, sizeof (*partition_patch));

  Point3D bounds_point;
  int maxH = 0;
  bounds_point.x = (int) file->idx->partition_count[0];
  bounds_point.y = (int) file->idx->partition_count[1];
  bounds_point.z = (int) file->idx->partition_count[2];
  char bitSequence[512];
  char bitPattern[512];
  GuessBitmaskPattern(bitSequence, bounds_point);
  maxH = strlen(bitSequence);

  for (uint32_t i = 0; i <= maxH; i++)
    bitPattern[i] = RegExBitmaskBit(bitSequence, i);

  int z_order = 0;
  int number_levels = maxH - 1;

  // iterate through all the partitions and finding out which partition a process belongs to
  for (uint32_t i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx->partition_size[0], index_i++)
  {
    for (uint32_t j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx->partition_size[1], index_j++)
    {
      for (uint32_t k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx->partition_size[2], index_k++)
      {
        partition_patch->offset[0] = i;
        partition_patch->offset[1] = j;
        partition_patch->offset[2] = k;
        partition_patch->size[0] = file->idx->partition_size[0];
        partition_patch->size[1] = file->idx->partition_size[1];
        partition_patch->size[2] = file->idx->partition_size[2];


        if (intersectNDChunk(partition_patch, local_p))
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

          // the partitions are arranged in z order
          file->idx_c->color = colors[z_order];

          file->idx->partition_offset[0] = i;
          file->idx->partition_offset[1] = j;
          file->idx->partition_offset[2] = k;

          break;
        }
      }
    }
  }
  free(partition_patch);

  free(colors);

  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx->filename_partition, "%s_%d.idx", file_name_skeleton, file->idx_c->color);

  free(local_p);

  return PIDX_success;
}



static PIDX_return_code create_partition_comm(PIDX_io file)
{
  // Grouping processes with same color and creating a partition for them
  if (MPI_Comm_split(file->idx_c->rst_comm, file->idx_c->color, file->idx_c->rrank, &(file->idx_c->partition_comm)) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  MPI_Comm_size(file->idx_c->partition_comm, &(file->idx_c->partition_nprocs));
  MPI_Comm_rank(file->idx_c->partition_comm, &(file->idx_c->partition_rank));

  return PIDX_success;
}



PIDX_return_code partition_cleanup(PIDX_io file)
{
  PIDX_time time = file->time;

  time->partition_cleanup_start = MPI_Wtime();
  // freeing the partitined comm
  if (MPI_Comm_free(&(file->idx_c->partition_comm)) != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->partition_cleanup_end = MPI_Wtime();

  return PIDX_success;
}



static PIDX_return_code find_partition_size(PIDX_io file)
{
  // calculate size of the partition partitions
  for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    file->idx->partition_size[d] = file->idx->bounds[d] / file->idx->partition_count[d];
    if (file->idx->bounds[d] % file->idx->partition_count[d] != 0)
      file->idx->partition_size[d]++;

    file->idx->partition_size[d] = pow(2, (int)ceil(log2(file->idx->partition_size[d])));

    if (file->idx->partition_size[d] < file->restructured_grid->patch_size[d])
    {
      fprintf(stderr,"Need to decrease the number of partition in %d dimension\n", (int)d);
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  uint64_t dims;
  uint64_t adjusted_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_bounds[0] = file->idx->partition_size[0] / file->idx->chunk_size[0];
  adjusted_bounds[1] = file->idx->partition_size[1] / file->idx->chunk_size[1];
  adjusted_bounds[2] = file->idx->partition_size[2] / file->idx->chunk_size[2];

  if (PIDX_inner_product(&dims, adjusted_bounds))
    return PIDX_err_size;

  if (dims < file->idx->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx->bits_per_block = getNumBits(file->idx->samples_per_block) - 1;
  }

  if (file->idx->bits_per_block == 0)
  {
    file->idx->bits_per_block = 0;
    file->idx->samples_per_block = 1;
  }

  return PIDX_success;
}


static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int check_bit = 0;
  for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}


static PIDX_return_code create_midx(PIDX_io file, int svi)
{

  int N;
  char dirname[1024], basename[1024];
  char file_temp[1024];

  int nbits_blocknumber = (file->idx->maxh - file->idx->bits_per_block - 1);
  VisusSplitFilename(file->idx->filename, dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--)
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

  //pidx does not do path remapping
  strcpy(file_temp, file->idx->filename);
  for (N = strlen(file_temp) - 1; N >= 0; N--)
  {
    int ch = file_temp[N];
    file_temp[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0)
    strcat(file_temp, "/%01x.bin");

  else
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4)
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8)
      strcat(file_temp, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12)
      strcat(file_temp, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16)
      strcat(file_temp, "/%04x.bin"); //no directories, 65536  files
    else
    {
      while (nbits_blocknumber > 16)
      {
        strcat(file_temp, "/%02x"); //256 subdirectories
        nbits_blocknumber -= 8;
      }
      strcat(file_temp, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }

  // making sure that the folder containing the .idx file is populated (the .idx file might be created other that ./)
  if (file->idx_c->rrank == 0)
  {
    int ret = 0;
    char last_path[PATH_MAX] = {0};
    char this_path[PATH_MAX] = {0};
    char tmp_path[PATH_MAX] = {0};
    char* pos;
    strcpy(this_path, file->idx->filename);
    if ((pos = strrchr(this_path, '/')))
    {
      pos[1] = '\0';
      if (!(strcmp(this_path, last_path) == 0))
      {
        //this file is in a previous directory than the last
        //one; we need to make sure that it exists and create
        //it if not.
        strcpy(last_path, this_path);
        memset(tmp_path, 0, PATH_MAX * sizeof (char));
        //walk up path and mkdir each segment
        for (int j = 0; j < (int)strlen(this_path); j++)
        {
          if (j > 0 && this_path[j] == '/')
          {
            //printf("mkdir %s\n", tmp_path);
            ret = mkdir(tmp_path, S_IRWXU | S_IRWXG | S_IRWXO);
            if (ret != 0 && errno != EEXIST)
            {
              //perror("mkdir");
              fprintf(stderr, "Error: failed to mkdir %s\n", tmp_path);
              return PIDX_err_file;
            }
          }
          tmp_path[j] = this_path[j];
        }
      }
    }
  }
  // Making sure that all processes wait for the folder containing the .idx file (which will also contain data) is populated
  MPI_Barrier(file->idx_c->rst_comm);

  int *colors;
  int num_parts = file->idx->partition_count[0] * file->idx->partition_count[1] * file->idx->partition_count[2];

  colors = malloc(sizeof(*colors) * num_parts);
  memset(colors, 0, sizeof(*colors) * num_parts);

  for (uint32_t k = 0; k < file->idx->partition_count[2]; k++)
    for (uint32_t j = 0; j < file->idx->partition_count[1]; j++)
      for (uint32_t i = 0; i < file->idx->partition_count[0]; i++)
        colors[(file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i] = (file->idx->partition_count[0] * file->idx->partition_count[1] * k) + (file->idx->partition_count[0] * j) + i;

  PIDX_variable var = file->idx->variable[svi];

  PIDX_patch local_p = (PIDX_patch)malloc(sizeof (*local_p));
  memset(local_p, 0, sizeof (*local_p));

  for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
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

    for (uint32_t proc = 0; proc < file->idx_c->partition_nprocs; proc++)
    {
      PIDX_patch curr_local_p = (PIDX_patch)malloc(sizeof (*curr_local_p));

      for (uint32_t d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      {
        curr_local_p->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * proc + d];
        curr_local_p->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * proc + d];
      }

      PIDX_patch partition_patch = (PIDX_patch)malloc(sizeof (*partition_patch));
      memset(partition_patch, 0, sizeof (*partition_patch));

      Point3D bounds_point;
      int maxH = 0;
      bounds_point.x = (int) file->idx->partition_count[0];
      bounds_point.y = (int) file->idx->partition_count[1];
      bounds_point.z = (int) file->idx->partition_count[2];
      char bitSequence[512];
      char bitPattern[512];
      GuessBitmaskPattern(bitSequence, bounds_point);
      maxH = strlen(bitSequence);

      for (uint32_t i = 0; i <= maxH; i++)
        bitPattern[i] = RegExBitmaskBit(bitSequence, i);

      int z_order = 0;
      int number_levels = maxH - 1;
      for (uint32_t i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx->partition_size[0], index_i++)
      {
        for (uint32_t j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx->partition_size[1], index_j++)
        {
          for (uint32_t k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx->partition_size[2], index_k++)
          {
            partition_patch->offset[0] = i;
            partition_patch->offset[1] = j;
            partition_patch->offset[2] = k;
            partition_patch->size[0] = file->idx->partition_size[0];
            partition_patch->size[1] = file->idx->partition_size[1];
            partition_patch->size[2] = file->idx->partition_size[2];

            if (intersectNDChunk(partition_patch, curr_local_p))
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

              for (uint32_t cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
              {
                int bit = bitPattern[number_levels];
                z_order |= ((uint64_t) PGET(xyzuv_Index, bit) & 1) << cnt;
                PGET(xyzuv_Index, bit) >>= 1;
              }
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
      free(partition_patch);

    }

    char midx_filename[PIDX_STRING_SIZE];
    sprintf(midx_filename,"%s.midx", file_name_skeleton);

    FILE* midx_file = fopen(midx_filename, "w");
    if (!midx_file)
    {
      fprintf(stderr, "Error [%s] [%d]: Cannot create filename %s.\n", __FILE__, __LINE__, midx_filename);
      return -1;
    }

    fprintf(midx_file, "<dataset typename='IdxMultipleDataset' mosaic='true'>\n");
    for (uint32_t i = 0; i < num_parts; i++)
    {
      if(strcmp(file_name_skeleton,"") == 0) // skip empty partitions
        continue;
      
      uint64_t curr_off = i*PIDX_MAX_DIMENSIONS;
      fprintf(midx_file, "\t<dataset url=\"file://$(CurrentFileDirectory)/%s\" name=\"%s_%d\" offset=\"%d %d %d\"/>\n",
              partition_filenames+(PIDX_STRING_SIZE*i), file_name_skeleton, i, offsets[curr_off+0],offsets[curr_off+1], offsets[curr_off+2]);
    }
    fprintf(midx_file,"</dataset>\n");

    fclose(midx_file);
    free(partition_filenames);
    free(offsets);
  }
  free(local_p);
  free(global_patch_size);
  free(global_patch_offset);
  free(colors);

  return PIDX_success;
}
