/*****************************************************
 **  PIDX Parallel I/O Library            **
 **  Copyright (c) 2010-2014 University of Utah   **
 **  Scientific Computing and Imaging Institute   **
 **  72 S Central Campus Drive, Room 3750       **
 **  Salt Lake City, UT 84112             **
 **                         **
 **  PIDX is licensed under the Creative Commons  **
 **  Attribution-NonCommercial-NoDerivatives 4.0  **
 **  International License. See LICENSE.md.     **
 **                         **
 **  For information about this project see:    **
 **  http://www.cedmav.com/pidx           **
 **  or contact: pascucci@sci.utah.edu        **
 **  For support: PIDX-support@visus.net      **
 **                         **
 *****************************************************/

/**
 * \file PIDX_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_raw_rst.h
 *
 */

#include "../../PIDX_inc.h"


///
/// \brief intersectNDChunk
/// \param A
/// \param B
/// \return
///
static int intersectNDChunk(PIDX_patch A, PIDX_patch B);


///
/// \brief contains_patch
/// \param reg_patch
/// \param patches
/// \param count
/// \return
///
static int contains_patch(PIDX_patch reg_patch, PIDX_patch* patches, int count);


PIDX_return_code PIDX_particles_rst_meta_data_create(PIDX_particles_rst_id rst_id)
{

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  PIDX_variable var0 = var_grp->variable[rst_id->first_index];
  int p = 0, v = 0;

  int r, d, c;
  unsigned long long max_vol, patch_count, pc;
  size_t reg_patch_count = 0;

  int start_var_index = rst_id->first_index;

  MPI_Allreduce(&var_grp->variable[start_var_index]->sim_patch_count, &rst_id->sim_max_patch_group_count, 1, MPI_INT, MPI_MAX, rst_id->idx_c->global_comm);

  rst_id->sim_raw_r_count = malloc(sizeof (double) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count);
  memset(rst_id->sim_raw_r_count, 0, (sizeof (double) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count));
  rst_id->sim_raw_r_offset = malloc(sizeof (double) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count);
  memset(rst_id->sim_raw_r_offset, 0, (sizeof (double) * rst_id->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count));

  for(pc=0; pc < var_grp->variable[start_var_index]->sim_patch_count; pc++)
  {
    double* tempoff = var_grp->variable[start_var_index]->sim_patch[pc]->physical_offset;
    double* tempsize = var_grp->variable[start_var_index]->sim_patch[pc]->physical_size;

    unsigned long long index = rst_id->idx_c->grank * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
    double* curr_patch_offset = &rst_id->sim_raw_r_offset[index];
    double* curr_patch_size = &rst_id->sim_raw_r_count[index];

    memcpy(curr_patch_offset, tempoff,sizeof(double) * PIDX_MAX_DIMENSIONS);
    memcpy(curr_patch_size, tempsize,sizeof(double) * PIDX_MAX_DIMENSIONS);
  }

  size_t* count_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*count_buffer_copy));
  memset(count_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*count_buffer_copy));

  memcpy(count_buffer_copy, &rst_id->sim_raw_r_count[rst_id->idx_c->grank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*count_buffer_copy));

  MPI_Allgather(count_buffer_copy, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count*sizeof(double), MPI_BYTE, rst_id->sim_raw_r_count, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(double), MPI_BYTE, rst_id->idx_c->global_comm);
  free(count_buffer_copy);

  off_t* offset_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*offset_buffer_copy));
  memset(offset_buffer_copy, 0, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*offset_buffer_copy));

  memcpy(offset_buffer_copy, &rst_id->sim_raw_r_offset[rst_id->idx_c->grank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count * sizeof(*offset_buffer_copy));

  MPI_Allgather(offset_buffer_copy, PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count * sizeof(double), MPI_BYTE, rst_id->sim_raw_r_offset, PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count * sizeof(double), MPI_BYTE, rst_id->idx_c->global_comm);
  free(offset_buffer_copy);

  var0->patch_group_count = 0;

  /// extents for the local process(rst_id->idx_c->grank)
  double adjusted_bounds[PIDX_MAX_DIMENSIONS];
  memcpy(adjusted_bounds, rst_id->idx->physical_bounds, PIDX_MAX_DIMENSIONS * sizeof(size_t));

  size_t max_found_reg_patches = 1;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    max_found_reg_patches *= ceil((float)rst_id->idx->physical_box_bounds[d]/(float)rst_id->physical_reg_patch_size[d]);


  //fprintf(stderr, "max_found_reg_patches = %d [%d %d %d]\n", max_found_reg_patches, rst_id->physical_reg_patch_size[0], rst_id->physical_reg_patch_size[1], rst_id->physical_reg_patch_size[2]);
  rst_id->reg_raw_grp_count = 0;

  unsigned long long pc0 = 0, d0 = 0;
  PIDX_patch* found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  int found_reg_patches_count = 0;
  for(pc0 = 0; pc0 < var_grp->variable[start_var_index]->sim_patch_count; pc0++)
  {
    PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
    {
      local_proc_patch->physical_offset[d0] = var_grp->variable[start_var_index]->sim_patch[pc0]->physical_offset[d0];
      local_proc_patch->physical_size[d0] = var_grp->variable[start_var_index]->sim_patch[pc0]->physical_size[d0];
    }

    for (double i = 0; i < adjusted_bounds[0]; i = i + rst_id->physical_reg_patch_size[0])
      for (double j = 0; j < adjusted_bounds[1]; j = j + rst_id->physical_reg_patch_size[1])
        for (double k = 0; k < adjusted_bounds[2]; k = k + rst_id->physical_reg_patch_size[2])
        {
          PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
          memset(reg_patch, 0, sizeof (*reg_patch));

          //Interior regular patches
          reg_patch->physical_offset[0] = i;
          reg_patch->physical_offset[1] = j;
          reg_patch->physical_offset[2] = k;
          reg_patch->physical_size[0] = rst_id->physical_reg_patch_size[0];
          reg_patch->physical_size[1] = rst_id->physical_reg_patch_size[1];
          reg_patch->physical_size[2] = rst_id->physical_reg_patch_size[2];

          if (intersectNDChunk(reg_patch, local_proc_patch))
          {
            if(!contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
            {
              found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
              memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));

              found_reg_patches_count++;
              rst_id->reg_raw_grp_count++;
            }
          }
          free(reg_patch);
        }
    free(local_proc_patch);
  }

  for(int i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  found_reg_patches = malloc(sizeof(PIDX_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(PIDX_patch*)*max_found_reg_patches);

  found_reg_patches_count = 0;

  rst_id->reg_raw_grp = malloc(sizeof(*rst_id->reg_raw_grp) * rst_id->reg_raw_grp_count);
  memset(rst_id->reg_raw_grp, 0, sizeof(*rst_id->reg_raw_grp) * rst_id->reg_raw_grp_count);

  reg_patch_count = 0;

  /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with

  for (double i = 0; i < adjusted_bounds[0]; i = i + rst_id->physical_reg_patch_size[0])
    for (double j = 0; j < adjusted_bounds[1]; j = j + rst_id->physical_reg_patch_size[1])
      for (double k = 0; k < adjusted_bounds[2]; k = k + rst_id->physical_reg_patch_size[2])
      {
        PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));

        //Interior regular patches
        reg_patch->physical_offset[0] = i;
        reg_patch->physical_offset[1] = j;
        reg_patch->physical_offset[2] = k;
        reg_patch->physical_size[0] = rst_id->physical_reg_patch_size[0];
        reg_patch->physical_size[1] = rst_id->physical_reg_patch_size[1];
        reg_patch->physical_size[2] = rst_id->physical_reg_patch_size[2];

        for(pc0 = 0; pc0 < var_grp->variable[start_var_index]->sim_patch_count; pc0++)
        {
          PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
          memset(local_proc_patch, 0, sizeof (*local_proc_patch));
          for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
          {
            local_proc_patch->physical_offset[d0] = var_grp->variable[start_var_index]->sim_patch[pc0]->physical_offset[d0];
            local_proc_patch->physical_size[d0] = var_grp->variable[start_var_index]->sim_patch[pc0]->physical_size[d0];
          }

          /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
          if (intersectNDChunk(reg_patch, local_proc_patch) && !contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
          {
            found_reg_patches[found_reg_patches_count] = (PIDX_patch)malloc(sizeof (*reg_patch));
            memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));
            found_reg_patches_count++;

            rst_id->reg_raw_grp[reg_patch_count] = malloc(sizeof(*(rst_id->reg_raw_grp[reg_patch_count])));
            memset(rst_id->reg_raw_grp[reg_patch_count], 0, sizeof(*(rst_id->reg_raw_grp[reg_patch_count])));

            PIDX_super_patch patch_grp = rst_id->reg_raw_grp[reg_patch_count];

            patch_grp->source_patch = (PIDX_source_patch_index*)malloc(sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
            patch_grp->patch = malloc(sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
            patch_grp->restructured_patch = malloc(sizeof(*patch_grp->restructured_patch));
            memset(patch_grp->source_patch, 0, sizeof(PIDX_source_patch_index) * rst_id->maximum_neighbor_count);
            memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * rst_id->maximum_neighbor_count);
            memset(patch_grp->restructured_patch, 0, sizeof(*patch_grp->restructured_patch));

            patch_count = 0;
            patch_grp->patch_count = 0;

            //Iterate through all processes
            for (r = 0; r < rst_id->idx_c->gnprocs; r++)
            {
              for(pc = 0; pc < rst_id->sim_max_patch_group_count; pc++)
              {
                //Extent of process with rst_id->idx_c->grank r
                PIDX_patch curr_patch = malloc(sizeof (*curr_patch));
                memset(curr_patch, 0, sizeof (*curr_patch));

                for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                {
                  unsigned long long index = r * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
                  curr_patch->physical_offset[d] = rst_id->sim_raw_r_offset[index+d];
                  curr_patch->physical_size[d] = rst_id->sim_raw_r_count[index+d];
                }

                if(curr_patch->size[0] == 0)
                {
                  // not existing patch for current rank, skip
                  //fprintf(stderr, "%d: skipping not existing patch\n", rank);
                  free(curr_patch);
                  continue;
                }

                if (intersectNDChunk(reg_patch, curr_patch))
                {
                  patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
                  memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

                  for (int c = 0; c < curr_patch->particle_count; c++)
                  {
                      // populate patch_grp->patch[patch_count]->particle_count
                      //patch_grp->patch[patch_count]->buffer
                  }

                  patch_grp->source_patch[patch_count].rank = r;
                  patch_grp->source_patch[patch_count].index = pc;
                  patch_count++;

                  if (patch_count >= rst_id->maximum_neighbor_count)
                  {
                    rst_id->maximum_neighbor_count = rst_id->maximum_neighbor_count * 2;

                    PIDX_source_patch_index *temp_buffer2 = realloc(patch_grp->source_patch, rst_id->maximum_neighbor_count * sizeof(PIDX_source_patch_index));
                    if (temp_buffer2 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->source_patch = temp_buffer2;

                    PIDX_patch *temp_buffer3 = realloc(patch_grp->patch, rst_id->maximum_neighbor_count * sizeof(*patch_grp->patch));
                    if (temp_buffer3 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->patch = temp_buffer3;

                    if (rst_id->idx_c->grank == 0)
                      fprintf(stderr, "[ERROR] rst_id->maximum_neighbor_count needs to be increased\n");
                  }

                  patch_grp->patch_count = patch_count;
                }
                free(curr_patch);
              }
            }

            patch_grp->max_patch_rank = patch_grp->source_patch[0].rank;
            max_vol = 1;
            for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              max_vol = max_vol * patch_grp->patch[0]->size[d];

            unsigned long long c_vol = 1;
            for(c = 1; c < patch_grp->patch_count; c++)
            {
              c_vol = 1;
              for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                c_vol = c_vol * patch_grp->patch[c]->size[d];
              if(c_vol > max_vol)
              {
                max_vol = c_vol;
                patch_grp->max_patch_rank = patch_grp->source_patch[c].rank;
              }
            }

            if(rst_id->idx_c->grank == patch_grp->max_patch_rank)
              var0->patch_group_count = var0->patch_group_count + 1;

            reg_patch_count++;
          }

          free(local_proc_patch);
        }
        free(reg_patch);
      }

  for(int i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    var->patch_group_count = var0->patch_group_count;
    var->patch_group_count = var_grp->variable[rst_id->first_index]->patch_group_count;

    var->rst_patch_group = malloc(var->patch_group_count * sizeof(*(var->rst_patch_group)));
    memset(var->rst_patch_group, 0, var->patch_group_count * sizeof(*(var->rst_patch_group)));
    for (p = 0; p < var->patch_group_count; p++)
    {
      var->rst_patch_group[p] = malloc(sizeof(*(var->rst_patch_group[p])));
      memset(var->rst_patch_group[p], 0, sizeof(*(var->rst_patch_group[p])));
    }
  }

  v = 0;
  p = 0;

  int cnt = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    cnt = 0;

    for (int i = 0; i < rst_id->reg_raw_grp_count; i++)
    {
      if (rst_id->idx_c->grank == rst_id->reg_raw_grp[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->rst_patch_group[cnt]; // here use patch_group
        patch_group->patch_count = rst_id->reg_raw_grp[i]->patch_count;
        patch_group->is_boundary_patch = rst_id->reg_raw_grp[i]->is_boundary_patch;
        patch_group->patch = malloc(sizeof(*(patch_group->patch)) * rst_id->reg_raw_grp[i]->patch_count);
        memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * rst_id->reg_raw_grp[i]->patch_count);

        patch_group->restructured_patch = malloc(sizeof(*(patch_group->restructured_patch)));
        memset(patch_group->restructured_patch, 0, sizeof(*(patch_group->restructured_patch)));

        for(int j = 0; j < rst_id->reg_raw_grp[i]->patch_count; j++)
        {
          patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
          memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));

          memcpy(patch_group->patch[j]->physical_offset, rst_id->reg_raw_grp[i]->patch[j]->physical_offset, PIDX_MAX_DIMENSIONS * sizeof(double));
          memcpy(patch_group->patch[j]->size, rst_id->reg_raw_grp[i]->patch[j]->size, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
        }
        memcpy(patch_group->restructured_patch->physical_offset, rst_id->reg_raw_grp[i]->restructured_patch->physical_offset, sizeof(double) * PIDX_MAX_DIMENSIONS);
        memcpy(patch_group->restructured_patch->size, rst_id->reg_raw_grp[i]->restructured_patch->size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
        cnt++;
      }
    }

    if (cnt != var->patch_group_count) //TODO CHECK THIS
      return PIDX_err_rst;
  }

  free(rst_id->sim_raw_r_offset);
  free(rst_id->sim_raw_r_count);

  return PIDX_success;
}



PIDX_return_code PIDX_particles_rst_meta_data_write(PIDX_particles_rst_id file)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[file->group_index];

  double *global_patch;
  PIDX_variable var0 = var_grp->variable[file->first_index];
  int max_patch_count;
  int patch_count = var0->patch_group_count;
  // TODO WILL: Why the max patch count across all ranks? Ranks could have differing numbers of patches
  // right? So would we leave some unused space in the file here?
  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);

  double *local_patch = malloc(sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));
  memset(local_patch, 0, sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));

  int pcounter = 0;
  local_patch[0] = (double)patch_count;
  for (int i = 0; i < patch_count; i++)
  {
    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * PIDX_MAX_DIMENSIONS + d + 1] = var0->rst_patch_group[i]->restructured_patch->physical_offset[d];

    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * PIDX_MAX_DIMENSIONS + PIDX_MAX_DIMENSIONS + d + 1] = var0->rst_patch_group[i]->restructured_patch->physical_size[d];

    local_patch[i * PIDX_MAX_DIMENSIONS + 2*PIDX_MAX_DIMENSIONS + 1] = var0->rst_patch_group[i]->restructured_patch->particle_count;

    pcounter++;
  }

  global_patch = malloc((file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));
  memset(global_patch, 0,(file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));

  MPI_Allgather(local_patch, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, global_patch + 2, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, file->idx_c->global_comm);

  global_patch[0] = (double)file->idx_c->gnprocs;
  global_patch[1] = (double)max_patch_count;

  char file_path[PATH_MAX];

  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(file_path, "%s_OFFSET_SIZE", directory_path);
  free(directory_path);
  if (file->idx_c->grank == 1 || file->idx_c->gnprocs == 1)
  {
    int fp = open(file_path, O_CREAT | O_WRONLY, 0664);
    ssize_t write_count = pwrite(fp, global_patch, (file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double), 0);
    if (write_count != (file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);
  }

  free(local_patch);
  free(global_patch);

  return PIDX_success;
}



PIDX_return_code PIDX_particles_rst_meta_data_destroy(PIDX_particles_rst_id rst_id)
{
  int i, j, v;
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for(i = 0; i < var_grp->variable[v]->patch_group_count; i++)
    {
      for(j = 0; j < var_grp->variable[v]->rst_patch_group[i]->patch_count; j++)
        free(var->rst_patch_group[i]->patch[j]);

      free(var->rst_patch_group[i]->restructured_patch);
      free(var->rst_patch_group[i]->patch);
      free(var->rst_patch_group[i]);

    }
    free(var->rst_patch_group);
  }


  for (i = 0; i < rst_id->reg_raw_grp_count; i++)
  {
    for (j = 0; j < rst_id->reg_raw_grp[i]->patch_count ; j++ )
      free(rst_id->reg_raw_grp[i]->patch[j]);

    free(rst_id->reg_raw_grp[i]->source_patch);
    free(rst_id->reg_raw_grp[i]->patch);
    free(rst_id->reg_raw_grp[i]->restructured_patch);
    free(rst_id->reg_raw_grp[i]);
  }

  free(rst_id->reg_raw_grp);
  rst_id->reg_raw_grp = 0;

  return PIDX_success;
}

/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}


static int contains_patch(PIDX_patch reg_patch, PIDX_patch* patches, int count)
{
  int i=0;

  for(i=0; i<count; i++)
  {
    int d=0;
    int matches = 0;
    for(d=0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      if(reg_patch->offset[d] == patches[i]->offset[d] && reg_patch->size[d] == patches[i]->size[d])
        matches++;
    }

    if(matches == PIDX_MAX_DIMENSIONS)
      return 1;
  }

  return 0;
}
