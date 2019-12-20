/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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
#include "../PIDX_file_handler.h"

/// Function to populate file descriptor when opening an existing IDX file version 6
PIDX_return_code PIDX_metadata_parse_v6_0(FILE *fp, PIDX_file* file)
{
  int variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];

  while (fgets(line, sizeof (line), fp) != NULL)
  {
    line[strcspn(line, "\r\n")] = 0;

    // Note: assuming version is the first info in the .idx
    if (strcmp(line, "(metadata version)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      
      strncpy((*file)->idx->metadata_version, line, 8);
    }
    
    if (strcmp(line, "(box)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        if (count % 2 == 1 && count / 2 < PIDX_MAX_DIMENSIONS)
        {
          (*file)->idx->bounds[count / 2] = atoi(pch) + 1;
          (*file)->idx->box_bounds[count / 2] = atoi(pch) + 1;
        }
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(physical box)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        if (count % 2 == 1 && count / 2 < PIDX_MAX_DIMENSIONS)
        {
          (*file)->idx->physical_bounds[count / 2] = atof(pch);
          (*file)->idx->physical_box_bounds[count / 2] = atof(pch);
        }
        count++;
        pch = strtok(NULL, " ");
      }
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
        (*file)->idx->partition_size[count] = atoi(pch);
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(partition count)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        (*file)->idx->partition_count[count] = atoi(pch);
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(restructure box size)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        (*file)->restructured_grid->patch_size[count] = atoi(pch);
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(raw_dump)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        (*file)->restructured_grid->patch_size[count] = atoi(pch);
        count++;
        pch = strtok(NULL, " ");
      }

      (*file)->idx->io_type = PIDX_RAW_IO;
      (*file)->idx->pidx_version = 0;
    }

    if (strcmp(line, "(io mode)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      int mode = atoi(line);
        
      if (mode == 0)
        (*file)->idx->io_type = PIDX_IDX_IO;
      else if (mode == 1)
        (*file)->idx->io_type = PIDX_LOCAL_PARTITION_IDX_IO;
      else if (mode == 2)
        (*file)->idx->io_type = PIDX_RAW_IO;
      else if (mode == 3)
        (*file)->idx->io_type = PIDX_PARTICLE_IO;

    }

    if (strcmp(line, "(endian)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      
      (*file)->idx->endian = atoi(line);
    }

    if (strcmp(line, "(fields)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      count = 0;
      variable_counter = 0;

      while (line[0] != '(')
      {
        (*file)->idx->variable[variable_counter] = malloc(sizeof (*((*file)->idx->variable[variable_counter])));
        if ((*file)->idx->variable[variable_counter] == NULL)
          return PIDX_err_file;

        memset((*file)->idx->variable[variable_counter], 0, sizeof (*((*file)->idx->variable[variable_counter])));

        pch1 = strtok(line, " +");
        while (pch1 != NULL)
        {
          if (count == 0)
          {
            char* temp_name = strdup(pch1);
            strcpy((*file)->idx->variable[variable_counter]->var_name, /*strdup(pch1)*/temp_name);
            free(temp_name);
          }

          if (count == 1)
          {
            len = strlen(pch1) - 1;
            if (pch1[len] == '\n')
              pch1[len] = 0;

            strcpy((*file)->idx->variable[variable_counter]->type_name, pch1);
            int ret;
            int bits_per_sample = 0;
            ret = PIDX_default_bits_per_datatype((*file)->idx->variable[variable_counter]->type_name, &bits_per_sample);
            if (ret != PIDX_success)  return PIDX_err_file;

            (*file)->idx->variable[variable_counter]->bpv = bits_per_sample;
            (*file)->idx->variable[variable_counter]->vps = 1;
          }
          count++;
          pch1 = strtok(NULL, " +");
        }
        count = 0;

        if ( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        variable_counter++;
      }
      (*file)->idx->variable_count = variable_counter;
    }

    if (strcmp(line, "(bits)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      strcat((*file)->idx->bitSequence, line);
    }

    if (strcmp(line, "(bitsperblock)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      (*file)->idx->bits_per_block = atoi(line);
      (*file)->idx->samples_per_block = (int)pow(2, (*file)->idx->bits_per_block);
    }

    if (strcmp(line, "(compression type)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      (*file)->idx->compression_type = atoi(line);
      if ((*file)->idx->compression_type != PIDX_NO_COMPRESSION)
      {
        int i1 = 0;
        for (i1 = 0; i1 < PIDX_MAX_DIMENSIONS; i1++)
          (*file)->idx->chunk_size[i1] = 4;
      }
    }

    if (strcmp(line, "(compressed box)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        (*file)->idx->chunk_size[count] = atoi(pch);
        count++;
        pch = strtok(NULL, " ");
      }

      //if ((*file)->idx->chunk_size[0] < 0 || (*file)->idx->chunk_size[1] < 0 || (*file)->idx->chunk_size[2] < 0)
      //  return PIDX_err_box;
    }

    if (strcmp(line, "(compression bit rate)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      (*file)->idx->compression_bit_rate = atof(line);
    }

    if (strcmp(line, "(blocksperfile)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      (*file)->idx->blocks_per_file= atoi(line);
    }

    if (strcmp(line, "(file system block size)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
      (*file)->fs_block_size = atoi(line);
    }

    if (strcmp(line, "(filename_template)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;
      line[strcspn(line, "\r\n")] = 0;
    }

    if (strcmp(line, "(time)") == 0)
    {
      if ( fgets(line, sizeof line, fp) == NULL)
        return PIDX_err_file;

      line[strcspn(line, "\r\n")] = 0;

      pch = strtok(line, " ");

      if (pch != NULL)
        (*file)->idx->first_tstep = atoi(pch);
      else
        return PIDX_err_file;

      pch = strtok(NULL, " ");

      if (pch != NULL)
        (*file)->idx->last_tstep = atoi(pch);
      else
        return PIDX_err_file;

      pch = strtok(NULL, " ");

      if (pch != NULL) {
        strcpy((*file)->idx->filename_time_template, pch);
        replace_str((*file)->idx->filename_time_template, "%","%%");
      }
      else
        return PIDX_err_file;
    }
  }

  return PIDX_success;
}
