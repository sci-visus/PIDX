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
#define INVERT_ENDIANESS 0
#include "../../PIDX_inc.h"

#if INVERT_ENDIANESS
static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf);
static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf);
#endif

static int maximum_neighbor_count = 1024;
static int intersectNDChunk(PIDX_patch A, PIDX_patch B);

PIDX_return_code PIDX_raw_rst_forced_raw_read(PIDX_raw_rst_id rst_id)
{
  int temp_max_dim = 3;

  if (rst_id->idx->pidx_version == 0)
    temp_max_dim = 5;

  int svi = rst_id->first_index;
  int evi = rst_id->last_index;

  char *idx_directory_path;
  char offset_path[PATH_MAX];
  char size_path[PATH_MAX];

  idx_directory_path = malloc(sizeof(*idx_directory_path) * PATH_MAX);
  memset(idx_directory_path, 0, sizeof(*idx_directory_path) * PATH_MAX);
  strncpy(idx_directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  sprintf(offset_path, "%s_OFFSET", idx_directory_path);
  sprintf(size_path, "%s_SIZE", idx_directory_path);
  free(idx_directory_path);

  uint32_t number_cores = 0;
  int fp = open(size_path, O_RDONLY);
  if (fp < 0)
  {
    fprintf(stderr, "Error opening file %s Error code %d\n", size_path, errno);
    return PIDX_err_io;
  }

  uint64_t read_count = pread(fp, &number_cores, sizeof(uint32_t), 0);
  if (read_count != sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

#if INVERT_ENDIANESS
  uint32_t temp_number_cores = 0;
  if (rst_id->idx->flip_endian == 1)
  {
    bit32_reverse_endian((unsigned char*)&number_cores, (unsigned char*)&temp_number_cores);
    number_cores = temp_number_cores;
  }
#endif

  uint32_t max_patch_count = 0;
  read_count = pread(fp, &max_patch_count, sizeof(uint32_t), sizeof(uint32_t));
  if (read_count != sizeof(uint32_t))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }


#if INVERT_ENDIANESS
  uint32_t temp_max_patch_count = 0;
  if (rst_id->idx->flip_endian == 1)
  {
    bit32_reverse_endian((unsigned char*)&max_patch_count, (unsigned char*)&temp_max_patch_count);
    max_patch_count = temp_max_patch_count;
  }
#endif

  int buffer_read_size = (number_cores * (max_patch_count * temp_max_dim + 1)) * sizeof(uint32_t);

  uint32_t *size_buffer = malloc(buffer_read_size);
  memset(size_buffer, 0, buffer_read_size);

  read_count = pread(fp, size_buffer, buffer_read_size, 2 * sizeof(uint32_t));
  if (read_count != buffer_read_size)
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }


  close(fp);

#if INVERT_ENDIANESS
  int buff_i;
  for (buff_i = 0; buff_i < (number_cores * (max_patch_count * temp_max_dim + 1)); buff_i++){
    uint32_t temp_value = 0;
    if (rst_id->idx->flip_endian == 1)
    {
      bit32_reverse_endian((unsigned char*)&size_buffer[buff_i], (unsigned char*)&temp_value);
      size_buffer[buff_i] = temp_value;
    }
  }
#endif

  uint32_t *offset_buffer = malloc(buffer_read_size);
  memset(offset_buffer, 0, buffer_read_size);

  int fp1 = open(offset_path, O_RDONLY);
  read_count = pread(fp1, offset_buffer, buffer_read_size, 2 * sizeof(uint32_t));
  if (read_count != buffer_read_size)
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }
  close(fp1);

#if INVERT_ENDIANESS
  for (buff_i = 0; buff_i<  (number_cores * (max_patch_count * temp_max_dim + 1)); buff_i++){
    uint32_t temp_value = 0;
    if (rst_id->idx->flip_endian == 1)
    {
      bit32_reverse_endian((unsigned char*)&offset_buffer[buff_i], (unsigned char*)&temp_value);
      offset_buffer[buff_i] = temp_value;
    }
  }
#endif

  char *file_name;
  file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, rst_id->idx->current_time_step);

  for (int pc1 = 0; pc1 < rst_id->idx->variable[svi]->sim_patch_count; pc1++)
  {
    PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = rst_id->idx->variable[svi]->sim_patch[pc1]->offset[d];
      local_proc_patch->size[d] = rst_id->idx->variable[svi]->sim_patch[pc1]->size[d];
    }

    // PC - - - - - - - - - -   PC - - - - - - - - - -  PC - - - - -
    // 0  1 2 3 4 5 6 7 8 9 10   11 12 13 14 15 16 17 18 19 20 21  22
    PIDX_patch n_proc_patch = (PIDX_patch)malloc(sizeof (*n_proc_patch));
    memset(n_proc_patch, 0, sizeof (*n_proc_patch));
    int p_counter = 1;
    int pc = 0, pc_index = 0;

    PIDX_super_patch patch_grp;
    patch_grp = malloc(sizeof(*(patch_grp)));
    memset(patch_grp, 0, sizeof(*(patch_grp)));

    int *source_patch_id = malloc(sizeof(int) * maximum_neighbor_count);
    memset(source_patch_id, 0, sizeof(int) * maximum_neighbor_count);

    patch_grp->source_patch = malloc(sizeof(*patch_grp->source_patch) * maximum_neighbor_count);
    patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
    patch_grp->restructured_patch = malloc(sizeof(*patch_grp->restructured_patch));
    memset(patch_grp->source_patch, 0, sizeof(*patch_grp->source_patch) * maximum_neighbor_count);
    memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
    memset(patch_grp->restructured_patch, 0, sizeof(*patch_grp->restructured_patch));

    int patch_count = 0;
    for (int n = 0; n < number_cores; n++)
    {
      pc = (int)offset_buffer[n * (max_patch_count * temp_max_dim + 1)];
      pc_index = n * (max_patch_count * temp_max_dim + 1);
      for (int m = 0; m < pc; m++)
      {
        for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          n_proc_patch->offset[d] = (uint64_t)offset_buffer[pc_index + m * temp_max_dim + d + 1];
          n_proc_patch->size[d] = (uint64_t)size_buffer[pc_index + m * temp_max_dim + d + 1];
        }

        if (intersectNDChunk(local_proc_patch, n_proc_patch))
        {
          //sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, n, m);
          patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
          memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));
          for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
            if (n_proc_patch->offset[d] <= local_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) - local_proc_patch->offset[d] + 1;
            }
            else if (local_proc_patch->offset[d] <= n_proc_patch->offset[d] && (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) >= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1))
            {
              patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) - n_proc_patch->offset[d] + 1;
            }
            else if (( local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) <= (n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) && local_proc_patch->offset[d] >= n_proc_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = local_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = local_proc_patch->size[d];
            }
            else if (( n_proc_patch->offset[d] + n_proc_patch->size[d] - 1) <= (local_proc_patch->offset[d] + local_proc_patch->size[d] - 1) && n_proc_patch->offset[d] >= local_proc_patch->offset[d])
            {
              patch_grp->patch[patch_count]->offset[d] = n_proc_patch->offset[d];
              patch_grp->patch[patch_count]->size[d] = n_proc_patch->size[d];
            }
          }
          patch_grp->source_patch[patch_count].rank = n;
          source_patch_id[patch_count] = m;
          //fprintf(stderr, "patch_count %d -> %d %d %d\n", patch_count, patch_grp->patch[patch_count]->size[0], patch_grp->patch[patch_count]->size[1], patch_grp->patch[patch_count]->size[2]);

          patch_count++;
          if (patch_count >= maximum_neighbor_count)
          {
            maximum_neighbor_count = maximum_neighbor_count * 2;
            //fprintf(stderr, "patch_count = %d maximum_neighbor_count = %d\n", patch_count, maximum_neighbor_count);

            int *temp_buffer1 = realloc(source_patch_id, sizeof(int) * maximum_neighbor_count);
            if (temp_buffer1 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              source_patch_id = temp_buffer1;

            PIDX_source_patch_index *temp_buffer2 = realloc(patch_grp->source_patch, maximum_neighbor_count * sizeof(*(patch_grp->source_patch)));
            if (temp_buffer2 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->source_patch = temp_buffer2;

            PIDX_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
            if (temp_buffer3 == NULL)
            {
              fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
              return PIDX_err_rst;
            }
            else
              patch_grp->patch = temp_buffer3;
          }
        }
        p_counter++;
      }
    }

    free(local_proc_patch);
    local_proc_patch = 0;
    free(n_proc_patch);
    n_proc_patch = 0;

    uint64_t send_c = 0;
    uint64_t sim_patch_offsetx[PIDX_MAX_DIMENSIONS];
    uint64_t sim_patch_countx[PIDX_MAX_DIMENSIONS];

    unsigned char **temp_patch_buffer2;
    temp_patch_buffer2 = malloc(sizeof(*temp_patch_buffer2) * (evi - svi + 1));
    memset(temp_patch_buffer2, 0, sizeof(*temp_patch_buffer2) * (evi - svi + 1));

    int other_offset = 0;
    for (int start_index = svi; start_index <= evi; start_index = start_index + 1)
    {
      int index = start_index - svi;
      PIDX_variable var = rst_id->idx->variable[start_index];

      for (uint64_t i = 0; i < patch_count; i++)
      {
        pc_index = patch_grp->source_patch[i].rank * (max_patch_count * temp_max_dim + 1);
        uint64_t total_sample_count = 1;
        for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          sim_patch_offsetx[d] = (uint64_t)offset_buffer[pc_index + source_patch_id[i] * temp_max_dim + d + 1];
          sim_patch_countx[d] = (uint64_t)size_buffer[pc_index + source_patch_id[i] * temp_max_dim + d + 1];
          total_sample_count = total_sample_count * (uint64_t)size_buffer[pc_index + source_patch_id[i] * temp_max_dim + d + 1];
        }

        temp_patch_buffer2[index] = malloc((uint64_t) sizeof(*(temp_patch_buffer2[index])) * total_sample_count * var->bpv/8 * var->vps);
        memset(temp_patch_buffer2[index], 0, (uint64_t) sizeof(*(temp_patch_buffer2[index])) * total_sample_count * var->bpv/8 * var->vps);

        sprintf(file_name, "%s/time%09d/%d_%d", directory_path, rst_id->idx->current_time_step, patch_grp->source_patch[i].rank, source_patch_id[i]);

        int fpx = open(file_name, O_RDONLY);
        other_offset = 0;
        for (int v1 = 0; v1 < start_index; v1++)
        {
          PIDX_variable var1 = rst_id->idx->variable[v1];
          other_offset = other_offset + ((var1->bpv/8) * var1->vps * sim_patch_countx[0] * sim_patch_countx[1] * sim_patch_countx[2]);
        }

        PIDX_variable var = rst_id->idx->variable[start_index];
        uint64_t preadc = pread(fpx, temp_patch_buffer2[index], total_sample_count * var->vps * var->bpv/8, other_offset);
        if (preadc != total_sample_count * var->vps * var->bpv/8)
        {
          fprintf(stderr, "[%s] [%d] Error in pread [%d %d]\n", __FILE__, __LINE__, (int)preadc, (int)send_c * var->bpv/8);
          return PIDX_err_rst;
        }
        close(fpx);

        uint64_t recv_o;
        for (uint64_t k1 = patch_grp->patch[i]->offset[2]; k1 < patch_grp->patch[i]->offset[2] + patch_grp->patch[i]->size[2]; k1++)
        {
          for (uint64_t j1 = patch_grp->patch[i]->offset[1]; j1 < patch_grp->patch[i]->offset[1] + patch_grp->patch[i]->size[1]; j1++)
          {
            for (uint64_t i1 = patch_grp->patch[i]->offset[0]; i1 < patch_grp->patch[i]->offset[0] + patch_grp->patch[i]->size[0]; i1 = i1 + patch_grp->patch[i]->size[0])
            {
              uint64_t send_index = (sim_patch_countx[0] * sim_patch_countx[1] * (k1 - sim_patch_offsetx[2])) +
                  (sim_patch_countx[0] * (j1 - sim_patch_offsetx[1])) +
                  (i1 - sim_patch_offsetx[0]);


              PIDX_variable var = rst_id->idx->variable[start_index];

              send_c = (patch_grp->patch[i]->size[0]);
              recv_o = (rst_id->idx->variable[start_index]->sim_patch[pc1]->size[0] * rst_id->idx->variable[start_index]->sim_patch[pc1]->size[1] * (k1 - rst_id->idx->variable[start_index]->sim_patch[pc1]->offset[2])) + (rst_id->idx->variable[start_index]->sim_patch[pc1]->size[0] * (j1 - rst_id->idx->variable[start_index]->sim_patch[pc1]->offset[1])) + (i1 - rst_id->idx->variable[start_index]->sim_patch[pc1]->offset[0]);

              memcpy(rst_id->idx->variable[start_index]->sim_patch[pc1]->buffer + ((uint64_t) recv_o * var->vps * (var->bpv/8)), temp_patch_buffer2[index] + send_index * var->vps * (var->bpv/8), send_c * var->vps * (var->bpv/8));

#if INVERT_ENDIANESS
              if (rst_id->idx->flip_endian == 1)
              {
                if (var->bpv/8 == 4 || var->bpv/8 == 12)
                {
                  float temp;
                  float temp2;

                  for (int y = 0; y < send_c * var->vps * (var->bpv/8) / sizeof(float); y++)
                  {
                    memcpy(&temp, rst_id->idx->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(float)), sizeof(float));
                    bit32_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
                    memcpy(rst_id->idx->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(float)), &temp2, sizeof(float));
                  }
                }
                else if (var->bpv/8 == 8 || var->bpv/8 == 24)
                {
                  double temp;
                  double temp2;

                  for (int y = 0; y < send_c * var->vps * (var->bpv/8) / sizeof(double); y++)
                  {
                    memcpy(&temp, rst_id->idx->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(double)), sizeof(double));
                    bit64_reverse_endian((unsigned char*)&temp, (unsigned char*)&temp2);
                    memcpy(rst_id->idx->variable[start_index]->sim_patch[pc1]->buffer + (recv_o * var->vps * (var->bpv/8)) + (y * sizeof(double)), &temp2, sizeof(double));
                  }
                }
              }
#endif

            }
          }
        }
        free(patch_grp->patch[i]);
        free(temp_patch_buffer2[index]);
      }
    }

    free(temp_patch_buffer2);

    free(patch_grp->source_patch);
    free(patch_grp->patch);
    free(patch_grp->restructured_patch);
    free(patch_grp);
    free(source_patch_id);
  }

  free(file_name);
  free(data_set_path);
  free(directory_path);

  free(offset_buffer);
  free(size_buffer);

  return PIDX_success;
}




static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int check_bit = 0;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}



#if INVERT_ENDIANESS

static void bit32_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
  unsigned char *data = ((unsigned char *)val) + 3;
  unsigned char *out = (unsigned char *)outbuf;

  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out = *data;

  return;
}

static void bit64_reverse_endian(unsigned char* val, unsigned char *outbuf)
{
  unsigned char *data = ((unsigned char *)val) + 7;
  unsigned char *out = (unsigned char *)outbuf;

  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out++ = *data--;
  *out = *data;

  return;
}

#endif
