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

static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode);
static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);
static PIDX_return_code create_midx(PIDX_io file, int gi, int svi);


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

PIDX_return_code PIDX_local_partition_idx_write(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 1:  Restrucure setup
  set_rst_box_size_for_write(file, gi, svi);

  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE) != PIDX_success)
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

  if (idx_restructure_rst_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  // TODO make a utility print error function for all the error prints
  
  if (var0->restructured_super_patch_count == 1)
  {
    
    // Step 3:  Partition
    if (partition(file, gi, svi, PIDX_WRITE) != PIDX_success)
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
    if (adjust_offsets(file, gi, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 5:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
#if 1
    for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
    {
      ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      // Step 6:  Setup HZ encoding Phase
      if (hz_encode_setup(file, gi, si, ei) != PIDX_success)
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

      if (hz_io(file, gi, PIDX_WRITE) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 8: Setup aggregation buffers
      ret = data_aggregate(file, gi, si, ei, AGG_SETUP, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 9: Performs data aggregation
      ret = data_aggregate(file, gi, si, ei, AGG_PERFORM, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 10: Performs actual file io
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, si, ei, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      finalize_aggregation(file, gi, si);
      time->io_end[gi][li] = PIDX_get_time();


      // Step 11: Cleanup for step 6
      if (hz_encode_cleanup(file) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }

    // Step 12: Cleanup the group and IDX related meta-data
    ret = group_meta_data_finalize(file, gi, svi, evi);
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

  //free_restructured_communicators(file, gi);

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

PIDX_return_code PIDX_local_partition_idx_read(PIDX_io file, int gi, int svi, int evi)
{
  int li = 0;
  int si = 0, ei = 0;
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 1:  Restrucure setup

  set_rst_box_size_for_write(file, gi, svi);

  if (idx_restructure_setup(file, gi, svi, evi - 1, PIDX_READ) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  if (idx_restructure_rst_comm_create(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (var0->restructured_super_patch_count == 1)
  {
    // Step 2:  Partition
    ret = partition(file, gi, svi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 3: Adjust per process offsets and global bounds as per the partition
    if (adjust_offsets(file, gi, svi, evi) != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    // Step 4:  Post partition group meta data
    ret = post_partition_group_meta_data_init(file, gi, svi, evi, PIDX_READ);
    if (ret != PIDX_success)
    {
      fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }

    for (si = svi; si < evi; si = si + (file->idx_d->variable_pipe_length + 1))
    {
      ei = ((si + file->idx_d->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->variable_pipe_length);
      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      // Step 5:  Setup HZ encoding Phase
      ret = hz_encode_setup(file, gi, si, ei);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      if (hz_io(file, gi, PIDX_READ) != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Setup 6: Setup aggregation buffers
      ret = data_aggregate(file, gi, si, ei, AGG_SETUP, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }


      // Setup 7: Performs actual file io
      time->io_start[gi][li] = PIDX_get_time();

      ret = data_io(file, gi, si, ei, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      time->io_end[gi][li] = PIDX_get_time();

#if 1
      //
      // Setup 8: Performs data aggregation
      ret = data_aggregate(file, gi, si, ei, AGG_PERFORM, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
      finalize_aggregation(file, gi, si);

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
    ret = group_meta_data_finalize(file, gi, svi, evi);
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
    PIDX_variable var = var_grp->variable[svi];
    for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
      var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] + file->idx_d->partition_offset[i];

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


PIDX_return_code PIDX_local_partition_mapped_idx_read(PIDX_io file, int gi, int svi, int evi)
{

}

// TODO move this function into some utils file, it is used in too many files already
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];
  
  return !(check_bit);
}

static PIDX_return_code create_midx(PIDX_io file, int gi, int svi)
{
  int *colors;
  int i = 0, j = 0, k = 0, d = 0;
  
  int num_parts = file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2];
  
  colors = malloc(sizeof(*colors) * num_parts);
  memset(colors, 0, sizeof(*colors) * num_parts);

    for (k = 0; k < file->idx_d->partition_count[2]; k++)
      for (j = 0; j < file->idx_d->partition_count[1]; j++)
        for (i = 0; i < file->idx_d->partition_count[0]; i++)
        {
          colors[(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * k) + (file->idx_d->partition_count[0] * j) + i] = (file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * k) + (file->idx_d->partition_count[0] * j) + i;
          
        }
  
    PIDX_variable_group var_grp = file->idx->variable_grp[gi];
    PIDX_variable var = var_grp->variable[svi];

    PIDX_patch local_p = (PIDX_patch)malloc(sizeof (*local_p));
    memset(local_p, 0, sizeof (*local_p));
    
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_p->offset[d] = var->restructured_super_patch->restructured_patch->offset[d];
      local_p->size[d] = var->restructured_super_patch->restructured_patch->size[d];
    }
    
    int wc = file->idx_c->partition_nprocs * (PIDX_MAX_DIMENSIONS);
  
    off_t* global_patch_offset = malloc(wc * sizeof(*global_patch_offset));
    memset(global_patch_offset, 0, wc * sizeof(*global_patch_offset));
    
    size_t* global_patch_size = malloc(wc * sizeof(*global_patch_size));
    memset(global_patch_size, 0, wc * sizeof(*global_patch_size));
    
    MPI_Allgather(&local_p->offset[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->partition_comm);
    
    MPI_Allgather(&local_p->size[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->partition_comm);
  
    //free(local_p);
  
  if(file->idx_c->simulation_rank == 0)
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
      bounds_point.x = (int) file->idx_d->partition_count[0];
      bounds_point.y = (int) file->idx_d->partition_count[1];
      bounds_point.z = (int) file->idx_d->partition_count[2];
      char bitSequence[512];
      char bitPattern[512];
      GuessBitmaskPattern(bitSequence, bounds_point);
      maxH = strlen(bitSequence);
      
      for (i = 0; i <= maxH; i++)
        bitPattern[i] = RegExBitmaskBit(bitSequence, i);
      
      int z_order = 0;
      int number_levels = maxH - 1;
      int index_i = 0, index_j = 0, index_k = 0;
      for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx_d->partition_size[0], index_i++)
      {
        for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx_d->partition_size[1], index_j++)
        {
          for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx_d->partition_size[2], index_k++)
          {
            reg_patch->offset[0] = i;
            reg_patch->offset[1] = j;
            reg_patch->offset[2] = k;
            reg_patch->size[0] = file->idx_d->partition_size[0];
            reg_patch->size[1] = file->idx_d->partition_size[1];
            reg_patch->size[2] = file->idx_d->partition_size[2];
            
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
                z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
                PGET(xyzuv_Index, bit) >>= 1;
              }
     
              //printf("found color %d offset %d %d %d\n", colors[z_order], i, j, k);
            
              sprintf(&partition_filenames[colors[z_order]*PIDX_STRING_SIZE], "%s_%d.idx", file_name_skeleton, colors[z_order]);
              
              off_t curr_off = colors[z_order]*PIDX_MAX_DIMENSIONS;
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
    
    for(i = 0; i < num_parts; i++)
    {
      off_t curr_off = i*PIDX_MAX_DIMENSIONS;
      
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


static PIDX_return_code partition(PIDX_io file, int gi, int svi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->partition_start = MPI_Wtime();
  // Calculates the number of partititons
  ret = find_partition_count(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  
  // assign same color to processes within the same partition
  ret = partition_setup(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Create midx file that will point to different .idx partition files
  if (create_midx(file, gi, svi) != PIDX_success)
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



static PIDX_return_code adjust_offsets(PIDX_io file, int gi, int svi, int evi)
{
  int i = 0, v = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  //fprintf(stderr, "SE %d %d\n", svi, evi);
  for (v = svi; v < evi; v++)
  {
  PIDX_variable var = var_grp->variable[v];

  if (var->restructured_super_patch_count != 1)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    var->restructured_super_patch->restructured_patch->offset[i] = var->restructured_super_patch->restructured_patch->offset[i] - file->idx_d->partition_offset[i];
  }

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    if (file->idx_d->partition_offset[i] + file->idx_d->partition_size[i] <= file->idx->box_bounds[i])
      file->idx->box_bounds[i] = file->idx_d->partition_size[i];
    else
      file->idx->box_bounds[i] = file->idx->box_bounds[i] - file->idx_d->partition_offset[i];

    if (file->idx_d->restructured_grid->patch_size[i] > file->idx->box_bounds[i])
      file->idx_d->restructured_grid->patch_size[i] = getPowerOf2(file->idx->box_bounds[i]);
  }

  memcpy(file->idx->bounds, file->idx->box_bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  //fprintf(stderr, "%d - %d %d %d -- PO %d %d %d\n", file->idx_c->simulation_rank, file->idx->box_bounds[0], file->idx->box_bounds[1], file->idx->box_bounds[2], file->idx_d->partition_offset[0], file->idx_d->partition_offset[1], file->idx_d->partition_offset[2]);

  return PIDX_success;
}



static PIDX_return_code post_partition_group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  //fprintf(stderr, "%d %d %d\n", file->idx_d->restructured_grid->patch_size[0], file->idx_d->restructured_grid->patch_size[1], file->idx_d->restructured_grid->patch_size[2]);

  time->bit_string_start = PIDX_get_time();
  // calculates maxh and bitstring
  ret = populate_local_bit_string(file, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // selects levels based on mode and maxh
  select_io_mode(file, gi);
  time->bit_string_end = PIDX_get_time();

  time->layout_start = PIDX_get_time();
  // calculates the block layout, given this is pure IDX only non-share block layout is populated
  //if (file->idx_d->color == 2)
  //{
  ret = populate_rst_block_layouts(file, gi, svi, file->hz_from_shared, file->hz_to_non_shared);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  //}
#if 1
  // Calculate the hz level upto which aggregation is possible
  ret = find_agg_level(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  // Creates the agg and io ids
  ret = create_agg_io_buffer(file, gi, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->layout_end = PIDX_get_time();

  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  ret = write_headers(file, gi, svi, evi, mode);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();
#endif
  return PIDX_success;
}


static PIDX_return_code group_meta_data_finalize(PIDX_io file, int gi, int svi, int evi)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->group_cleanup_start = PIDX_get_time();
  ret = destroy_agg_io_buffer(file, svi, evi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  ret = delete_rst_block_layout(file, gi);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->group_cleanup_end = PIDX_get_time();

  return PIDX_success;
}
