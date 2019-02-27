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
 
#include "../PIDX_inc.h"

void mira_create_folder_name(char* bin_file, char* folder_name)
{
  char *ch_dot;
  ch_dot=strrchr(bin_file,'.');
  //fprintf(stderr, "Last occurence of '.' found at %d \n", (ch_dot - bin_file + 1) );
  
  char *ch_slash;
  ch_slash=strrchr(bin_file,'/');
  //fprintf(stderr, "Last occurence of '/' found at %d \n", (ch_slash - bin_file + 1) );
  
  //fprintf(stderr, "Name length = %d\n", (ch_dot - ch_slash));
  
  char file_name[PIDX_FILE_PATH_LENGTH];
  memset(file_name, 0, (ch_dot - ch_slash));
  strncpy(file_name, bin_file + (ch_slash - bin_file + 1), (ch_dot - ch_slash - 1));
  //fprintf(stderr, "File name = %s\n", file_name);
  
  char bin_file_first_part[PIDX_FILE_PATH_LENGTH];
  strncpy(bin_file_first_part, bin_file, (ch_slash - bin_file));
  //fprintf(stderr, "File name copy %s\n", bin_file_first_part);
  
  sprintf(folder_name, "%s/%s", bin_file_first_part, file_name);
  //fprintf(stderr, "Final File name = %s\n", folder_name);
    
  return;
}

void adjust_file_name(char* bin_file, char* adjusted_name)
{
  char *ch_dot;
  ch_dot=strrchr(bin_file,'.');
  //fprintf(stderr, "Last occurence of '.' found at %d \n", (ch_dot - bin_file + 1) );
  
  char *ch_slash;
  ch_slash=strrchr(bin_file,'/');
  //fprintf(stderr, "Last occurence of '/' found at %d \n", (ch_slash - bin_file + 1) );
  
  //fprintf(stderr, "Name length = %d\n", (ch_dot - ch_slash));
  
  char file_name[PIDX_FILE_PATH_LENGTH];
  memset(file_name, 0, (ch_dot - ch_slash));
  strncpy(file_name, bin_file + (ch_slash - bin_file + 1), (ch_dot - ch_slash - 1));
  //fprintf(stderr, "File name = %s\n", file_name);
  
  char bin_file_first_part[PIDX_FILE_PATH_LENGTH];
  strncpy(bin_file_first_part, bin_file, (ch_slash - bin_file));
  //fprintf(stderr, "File name copy %s\n", bin_file_first_part);
  
  sprintf(adjusted_name, "%s/%s/%s.bin", bin_file_first_part, file_name, file_name);
  //fprintf(stderr, "Final File name = %s\n", adjusted_name);
    
  return;
}

int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen)
{
  uint64_t address = 0;
  unsigned int segs[PIDX_MAX_TEMPLATE_DEPTH] = {0};
  int seg_count = 0;
  char* pos;
  int ret;

  //fprintf(stderr, "[generate_file_name]: %d %s %d :: %s\n", file_number, filename, maxlen, filename_template);
  // determine the first HZ address for the file in question 
  address = file_number * blocks_per_file;

  // walk backwards through the file name template to find places where we need to substitute strings
  for (pos = &filename_template[strlen(filename_template) - 1];
          pos != filename_template;
          pos--) 
  {
    // be careful not to lo0 past the end of the array 
    if (pos - filename_template > (strlen(filename_template) - 3))
      continue;

    if (pos[0] == '%' && pos[1] == '0' && pos[3] == 'x') 
    {
      // TODO: for now we have a hard coded max depth 
      if (seg_count >= PIDX_MAX_TEMPLATE_DEPTH)
      {
        fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
        return 1;
      }

      // found an occurance of %0 in the template; check the next character to see how many bits to use here 

      switch (pos[2]) 
      {
        case '1':
            segs[seg_count] += address & 0xf;
            address = address >> 4;
            break;
        case '2':
            segs[seg_count] += address & 0xff;
            address = address >> 8;
            break;
        case '3':
            segs[seg_count] += address & 0xfff;
            address = address >> 12;
            break;
        case '4':
            segs[seg_count] += address & 0xffff;
            address = address >> 16;
            break;
        case '5':
            segs[seg_count] += address & 0xfffff;
            address = address >> 20;
            break;
        default:
            // TODO: generalize this to work for any value 
            fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
            return 1;
      }
      seg_count++;
    }
  }
  switch (seg_count) 
  {
    case 0:
        ret = strlen(filename_template);
        if (ret < maxlen) {
            strcpy(filename, filename_template);
        }
        break;
    case 1:
        ret = snprintf(filename, maxlen, filename_template, segs[0]);
        break;
    case 2:
        ret = snprintf(filename, maxlen, filename_template,
                segs[1], segs[0]);
        break;
    case 3:
        ret = snprintf(filename, maxlen, filename_template,
                segs[2], segs[1], segs[0]);
        break;
    case 4:
        ret = snprintf(filename, maxlen, filename_template,
                segs[3], segs[2], segs[1], segs[0]);
        break;
    case 5:
        ret = snprintf(filename, maxlen, filename_template,
                segs[4], segs[3], segs[2], segs[1], segs[0]);
        break;
    case 6:
        ret = snprintf(filename, maxlen, filename_template,
                segs[5], segs[4], segs[3], segs[2], segs[1], segs[0]);
        break;
    default:
        // TODO: generalize this 
        fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
        return 1;
        break;
  }
  // make sure that the resulting string fit into the buffer ok 
  if (ret >= maxlen - 1) 
  {
    fprintf(stderr, "Error: filename too short in generate_filename()\n");
    return 1;
  }
  return 0;
}


/////////////////////////////////////////////////
int generate_file_name_template(int maxh, int bits_per_block, char* filename, char* time_template, int current_time_step, char* filename_template)
{
  int N;
  char dirname[PIDX_FILE_PATH_LENGTH], basename[PIDX_FILE_PATH_LENGTH];
  int nbits_blocknumber;
  char* directory_path;
  char* data_set_path;
  
  data_set_path = malloc(sizeof(*data_set_path) * 1024);
  memset(data_set_path, 0, sizeof(*data_set_path) * 1024);

  directory_path = malloc(sizeof(*directory_path) * 1024);
  memset(directory_path, 0, sizeof(*directory_path) * 1024);

  strncpy(directory_path, filename, strlen(filename) - 4);
  char this_time_template[512];
  sprintf(this_time_template, "%%s/%s.idx", time_template);

  sprintf(data_set_path, this_time_template, directory_path, current_time_step);
  free(directory_path);

  nbits_blocknumber = (maxh - bits_per_block - 1);
  VisusSplitFilename(data_set_path, dirname, basename);
  //fprintf(stderr, "dirname %s basename %s\n", dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--) 
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

#if 0
  //if i put . as the first character, if I move files VisusOpen can do path remapping
  sprintf(pidx->filename_template, "./%s", basename);
#endif
  //pidx does not do path remapping 
  strcpy(filename_template, data_set_path);
  for (N = strlen(filename_template) - 1; N >= 0; N--) 
  {
    int ch = filename_template[N];
    filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
    strcat(filename_template, "/%01x.bin");
   
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
      strcat(filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12) 
      strcat(filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16) 
      strcat(filename_template, "/%04x.bin"); //no directories, 65536  files
    else 
    {
      while (nbits_blocknumber > 16) 
      {
        strcat(filename_template, "/%02x"); //256 subdirectories
        nbits_blocknumber -= 8;
      }
      strcat(filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }
  
  free(data_set_path);
  return 0;
}
