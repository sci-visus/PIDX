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

/*
             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        IDX Format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*
*/

#include <stdarg.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>

#define MAX_DIMENSIONS 5
#define MAX_TEMPLATE_DEPTH 6

//#include <PIDX.h>

enum { X, Y, Z, NUM_DIMS };
//static int process_count = 1, rank = 0;
//static int cores = 1;
static int rank = 0;
static char output_file_template[512] = "test";
static char output_file_name[512] = "test.idx";
static char partition_file_name[512];
static char partition_file_template[512] = "test";
static int idx_count[MAX_DIMENSIONS];

static int current_time_step = 0;
static int compression_type;
static int bits_per_block = 0;
static int samples_per_block;
static int blocks_per_file = 0;
static int bounds[MAX_DIMENSIONS];
static int chunked_bounds[MAX_DIMENSIONS];
static int idx_size[MAX_DIMENSIONS];
static int idx_count[MAX_DIMENSIONS];
static char bitPattern[512];
static char bitSequence[512];
static int compression_bit_rate = 64;
static int compression_factor = 1;
static int chunk_size[MAX_DIMENSIONS];
static int data_core_count;
static char var_name[128][512];
static char type_name[128][512];
static int bpv[128];
static int vps[128];
static int variable_count = 0;
static int fs_block_size = 0;
static int maxh = 0;
static int max_file_count = 0;
static int start_time_step = 0;
static int end_time_step = 0;

static int get_default_bits_per_datatype(char* type, int* bits);
static unsigned long long getPowerOf2(int x);
static int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen) ;
static double get_time();

static char *usage = "Serial Usage: ./checkpoint -g 32x32x32 -l 32x32x32 -v 3 -t 16 -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./checkpoint -g 32x32x32 -l 16x16x16 -f output_idx_file_name -v 3 -t 16\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -t: number of timesteps\n"
                     "  -v: number of variables\n";


//----------------------------------------------------------------
static void terminate()
{
#if 0
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
}

//----------------------------------------------------------------
static void terminate_with_error_msg(const char *format, ...)
{
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
  terminate();
}

//----------------------------------------------------------------
static void rank_0_print(const char *format, ...)
{
  if (rank != 0) return;
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
}

//----------------------------------------------------------------
static void init_mpi(int argc, char **argv)
{
#if 0
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
#endif
}

//----------------------------------------------------------------
static void shutdown_mpi()
{
#if 0
  MPI_Finalize();
#endif
}





///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  char flags[] = "f:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", usage);
      sprintf(output_file_name, "%s", output_file_template);
      break;

    //case('c'): // number of timesteps
    //  if (sscanf(optarg, "%d", &cores) < 0)
    //    terminate_with_error_msg("Invalid variable file\n%s", usage);
    //  break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

static int IDX_file_open(const char* filename)
{
  int i;
  //int ret;
  int rank = 0;

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return (-1);

  for (i = 0; i < MAX_DIMENSIONS; i++)
    idx_count[i] = 1;

  current_time_step = 0;
  compression_type = 0;

  bits_per_block = 0;
  blocks_per_file = 0;

  for (i=0;i<MAX_DIMENSIONS;i++)
    bounds[i]=65535;

  memset(bitPattern, 0, 512);
  memset(bitSequence, 0, 512);

  compression_bit_rate = 64;
  compression_factor = 1;
  for (i=0;i<MAX_DIMENSIONS;i++)
    chunk_size[i] = 1;

  samples_per_block = (int)pow(2, 0);

  int var = 0, variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];

  if (rank == 0)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
      fprintf(stderr, "Error Opening %s\n", filename);
      return (-1);
    }

    while (fgets(line, sizeof (line), fp) != NULL)
    {
      //fprintf(stderr, "%s", line);
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "(box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          if (count % 2 == 1)
            bounds[count / 2] = atoi(pch) + 1;
          count++;
          pch = strtok(NULL, " ");
        }
      }

      if (strcmp(line, "(partition size)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          if (count % 2 == 1)
            idx_size[count / 2] = atoi(pch) ;
          count++;
          pch = strtok(NULL, " ");
        }
      }

      if (strcmp(line, "(partition count)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          if (count % 2 == 1)
            idx_count[count / 2] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }
      }


      if (strcmp(line, "(cores)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        data_core_count = atoi(line);
      }

      if (strcmp(line, "(fields)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        count = 0;
        variable_counter = 0;

        while (line[0] != '(')
        {
          pch1 = strtok(line, " +");
          while (pch1 != NULL)
          {
            if (count == 0)
            {
              char* temp_name = strdup(pch1);
              strcpy(var_name[variable_counter], temp_name);
              free(temp_name);
            }

            if (count == 1)
            {
              len = strlen(pch1) - 1;
              if (pch1[len] == '\n')
                pch1[len] = 0;

              strcpy(type_name[variable_counter], pch1);
              int ret;
              int bits_per_sample = 0;
              ret = get_default_bits_per_datatype(type_name[variable_counter], &bits_per_sample);
              if (ret != 0)  return (-1);

              bpv[variable_counter] = bits_per_sample;
              vps[variable_counter] = 1;
            }
            count++;
            pch1 = strtok(NULL, " +");
          }
          count = 0;

          if( fgets(line, sizeof line, fp) == NULL)
            return (-1);
          line[strcspn(line, "\r\n")] = 0;
          variable_counter++;
        }
        variable_count = variable_counter;
      }

      if (strcmp(line, "(bits)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        strcpy(bitSequence, line);
      }

      if (strcmp(line, "(bitsperblock)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        bits_per_block = atoi(line);
        samples_per_block = (int)pow(2, bits_per_block);
      }

      if (strcmp(line, "(compression type)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        compression_type = atoi(line);
        if (compression_type != 0)
        {
          int i1 = 0;
          for (i1 = 0; i1 < MAX_DIMENSIONS; i1++)
            chunk_size[i1] = 4;
        }
      }

      if (strcmp(line, "(compressed box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          chunk_size[count] = atoi(pch);
          count++;
          pch = strtok(NULL, " ");
        }

        if(chunk_size[0] < 0 || chunk_size[1] < 0 || chunk_size[2] < 0)
          return (-1);
      }

      if (strcmp(line, "(compression bit rate)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        compression_bit_rate = atoi(line);
      }

      if (strcmp(line, "(blocksperfile)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
        blocks_per_file= atoi(line);
      }

      if (strcmp(line, "(filename_template)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;
      }

      if (strcmp(line, "(time)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return (-1);
        line[strcspn(line, "\r\n")] = 0;

        pch = strtok(line, " ");
        count = 0;
        while (pch != NULL)
        {
          if (count == 0)
            start_time_step = atoi(pch);

          if (count == 1)
            end_time_step = atoi(pch);

          count++;
          pch = strtok(NULL, " ");
        }
      }
    }
    fclose(fp);
  }

  fprintf(stderr, "Start time step %d End time step %d\n", start_time_step, end_time_step);

#if 0

  MPI_Bcast(bounds, 5, MPI_UNSIGNED_LONG_LONG, 0, comm);
  MPI_Bcast(idx_count, 5, MPI_UNSIGNED_LONG_LONG, 0, comm);
  MPI_Bcast(idx_size, 5, MPI_UNSIGNED_LONG_LONG, 0, comm);
  MPI_Bcast(chunk_size, 5, MPI_UNSIGNED_LONG_LONG, 0, comm);
  MPI_Bcast(&(blocks_per_file), 1, MPI_INT, 0, comm);
  MPI_Bcast(&(bits_per_block), 1, MPI_INT, 0, comm);
  MPI_Bcast(&(variable_count), 1, MPI_INT, 0, comm);
  MPI_Bcast(&(data_core_count), 1, MPI_INT, 0, comm);
  MPI_Bcast(&(compression_bit_rate), 1, MPI_INT, 0, comm);
  MPI_Bcast(&(compression_type), 1, MPI_INT, 0, comm);


  if (compression_type == PIDX_CHUNKING_ZFP)
  {
    if (compression_bit_rate == 32)
      compression_factor = 2;
    if (compression_bit_rate == 16)
      compression_factor = 4;
    if (compression_bit_rate == 8)
      compression_factor = 8;
    if (compression_bit_rate == 4)
      compression_factor = 16;
    if (compression_bit_rate == 2)
      compression_factor = 32;
    if (compression_bit_rate == 1)
      compression_factor = 64;
  }

  int d = 0;
  for (d = 0; d < MAX_DIMENSIONS; d++)
  {
    if (bounds[d] % chunk_size[d] == 0)
      chunked_bounds[d] = (int) bounds[d] / chunk_size[d];
    else
      chunked_bounds[d] = (int) (bounds[d] / chunk_size[d]) + 1;
  }

  samples_per_block = (int)pow(2, bits_per_block);

#endif

  for (var = 0; var < variable_count; var++)
  {
#if 0
    MPI_Bcast(&(bpv[var]), 1, MPI_INT, 0, comm);
    MPI_Bcast(&(vps[var]), 1, MPI_INT, 0, comm);
    MPI_Bcast(var_name[var], 512, MPI_CHAR, 0, comm);
    MPI_Bcast(type_name[var], 512, MPI_CHAR, 0, comm);
#endif
  }


  unsigned long long total_reg_sample_count = (getPowerOf2(chunked_bounds[0]) * getPowerOf2(chunked_bounds[1]) * getPowerOf2(chunked_bounds[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return (-1);
  }

  unsigned long long max_sample_per_file = (unsigned long long) samples_per_block * blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, samples_per_block, blocks_per_file);
    return (-1);
  }

  max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    max_file_count++;

  //fprintf(stderr, "%d %d %d\n", chunk_size[0], chunk_size[1], chunk_size[2]);

  if (rank == 0)
  {
    int ret;
    struct stat stat_buf;
    ret = stat(filename, &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return (-1);
    }
    fs_block_size = stat_buf.st_blksize;
  }

#if 0
  MPI_Bcast(&(fs_block_size), 1, MPI_INT, 0, comm);
#endif

  return 0;
}


int SplitFilename(const char* filename,char* dirname,char* basename)
{
  int i;
  int N=strlen(filename);

  if (!N)
    return 0;

  //find the last separator
  for (i=N-1;i>=0;i--)
  {
    if (filename[i]=='/' || filename[i]=='\\')
    {
      strncpy(dirname,filename,i);
      dirname[i]=0;
      strcpy(basename,filename+i+1);
      return 1;
    }
  }
  //assume is all filename (without directory name)
  dirname [0]=0;
  strcpy(basename,filename);
  return 1;
}


int file_initialize_time_step(int current_time_step, char* file_name, char* file_template)
{
  int N;
  char dirname[PIDX_FILE_PATH_LENGTH], basename[PIDX_FILE_PATH_LENGTH];
  int nbits_blocknumber;
  char *directory_path;
  char *data_set_path;

  data_set_path = malloc(sizeof(*data_set_path) * PIDX_FILE_PATH_LENGTH);
  memset(data_set_path, 0, sizeof(*data_set_path) * PIDX_FILE_PATH_LENGTH);

  directory_path = malloc(sizeof(*directory_path) * PIDX_FILE_PATH_LENGTH);
  memset(directory_path, 0, sizeof(*directory_path) * PIDX_FILE_PATH_LENGTH);

  strncpy(directory_path, file_name, strlen(file_name) - 4);

  char time_template[512];
  sprintf(time_template, "%%s/%s.idx", file->idx->filename_time_template);
  sprintf(data_set_path, time_template, directory_path, current_time_step);
  free(directory_path);

  nbits_blocknumber = (maxh - bits_per_block - 1);
  SplitFilename(data_set_path, dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--)
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

#if 0
  //if i put . as the first character, if I move files VisusOpen can do path remapping
  sprintf(file_template, "./%s", basename);
#endif
  //pidx does not do path remapping
  strcpy(file_template, data_set_path);
  for (N = strlen(file_template) - 1; N >= 0; N--)
  {
    int ch = file_template[N];
    file_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0)
    strcat(file_template, "/%01x.bin");

  else
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4)
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8)
      strcat(file_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12)
      strcat(file_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16)
      strcat(file_template, "/%04x.bin"); //no directories, 65536  files
    else
    {
      while (nbits_blocknumber > 16)
      {
        strcat(file_template, "/%02x"); //256 subdirectories
        nbits_blocknumber -= 8;
      }
      strcat(file_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }

  free(data_set_path);
  return 0;
}


int main(int argc, char **argv)
{
  double start_time, end_time;
  start_time = get_time();
  int ret = 0;
  init_mpi(argc, argv);
  parse_args(argc, argv);

  rank_0_print("Merge Program\n");

#if 0
  comm = MPI_COMM_WORLD;
#endif

  ret = IDX_file_open(output_file_name);
  if (ret != 0)  terminate_with_error_msg("PIDX_file_create");

  maxh = strlen(bitSequence);

  fprintf(stderr, "Partition size :: and count %d %d %d :: %d %d %d\n", idx_count[0], idx_count[1], idx_count[2], idx_size[0], idx_size[1], idx_size[2]);
  fprintf(stderr, "bitstring %s maxh = %d\n", bitSequence, maxh);

  // shared_block_level is the level upto which the idx blocks are shared
  int shared_block_level = (int)log2(idx_count[0] * idx_count[1] * idx_count[2]) + bits_per_block + 1;
  if (shared_block_level >= maxh)
    shared_block_level = maxh;

  int shared_block_count = pow(2, shared_block_level - 1) / samples_per_block;
  fprintf(stderr, "Shared block level = %d Shared block count = %d\n", shared_block_level, shared_block_count);

  int level = 0;
  int ts = 0;

  // Iteration through all the timesteps
  for (ts = start_time_step; ts <= end_time_step; ts++)
  {
    // Iteration through all the shared blocks
    //for (level = 0; level < shared_block_level; level = level + 1)
    {
      int hz_index = (int)pow(2, level - 1);
      int file_no = hz_index / (blocks_per_file * samples_per_block);
      int file_count;
      char existing_file_name[PIDX_FILE_PATH_LENGTH];
      char new_file_name[PIDX_FILE_PATH_LENGTH];
      int ic = 0;
      if (level <= bits_per_block + log2(blocks_per_file) + 1)
        file_count = 1;
      else
        file_count = (int)pow(2, level - (bits_per_block + log2(blocks_per_file) + 1));

      // file_no is the index of the file that needs to be opened to read from all the partitions
      // they contain the shared blocks
      fprintf(stderr, "Opening file %d\n", file_no);

#if 1
      // iterate throuh all the files that contains the shared blocks
      // most likely this will be only the first file of all the partitions
      // so fc = 1
      int fc = 0;
      for (fc = file_no; fc < file_no + file_count; fc++)
      {
        // malloc for the header for the outpur blocks, i.e. the merged blocks
        uint32_t* write_binheader;
        int write_binheader_count;
        write_binheader_count = 10 + 10 * blocks_per_file * variable_count;
        int write_binheader_length = write_binheader_count * sizeof (*write_binheader);
        write_binheader = malloc(write_binheader_length);
        memset(write_binheader, 0, write_binheader_length);

        //iterate through all the variables/fields
        int var = 0;
        off_t var_offset = 0;
        for (var = 0; var < 1; var++)
        {
          unsigned char *write_data_buffer = malloc(samples_per_block * shared_block_count * bpv[var]/8);
          memset(write_data_buffer, 0, samples_per_block * shared_block_count * bpv[var]/8);
          //fprintf(stderr, "Write bufer size = %d [%d x %d x %d]\n", samples_per_block * shared_block_count * bpv[var]/8, (int)pow(2, bits_per_block), shared_block_count, bpv[var]/8);

          // shared block data
          // doube pointer (number o fpartitions x number of shared blocks)
          unsigned char **read_data_buffer = malloc(idx_count[0] * idx_count[1] * idx_count[2] * sizeof(*read_data_buffer));
          memset(read_data_buffer, 0, idx_count[0] * idx_count[1] * idx_count[2] * sizeof(*read_data_buffer));

          // shared block header info

          uint32_t** read_binheader = malloc(idx_count[0] * idx_count[1] * idx_count[2] * sizeof(*read_binheader));
          memset(read_binheader, 0, idx_count[0] * idx_count[1] * idx_count[2] * sizeof(*read_binheader));

          file_initialize_time_step(ts, output_file_name, output_file_template);
          generate_file_name(blocks_per_file, output_file_template, fc, new_file_name, PATH_MAX);
          //fprintf(stderr, "Merged blocks to be written in %s\n", new_file_name);

          // iterate through all the parttions
          for (ic = 0; ic < idx_count[0] * idx_count[1] * idx_count[2]; ic++)
          {
            char file_name_skeleton[PIDX_FILE_PATH_LENGTH];
            strncpy(file_name_skeleton, output_file_name, strlen(output_file_name) - 4);
            file_name_skeleton[strlen(output_file_name) - 4] = '\0';

            if (idx_count[0] != 1 || idx_count[1] != 1 || idx_count[2] != 1)
              sprintf(partition_file_name, "%s_%d.idx", file_name_skeleton, ic);
            else
              strcpy(partition_file_name, output_file_name);

            file_initialize_time_step(ts, partition_file_name, partition_file_template);
            generate_file_name(blocks_per_file, partition_file_template, fc, existing_file_name, PATH_MAX);

            int read_binheader_count;
            read_binheader_count = 10 + 10 * blocks_per_file * variable_count;
            read_binheader[ic] = (uint32_t*) malloc(sizeof (*read_binheader[ic])*(read_binheader_count));
            memset(read_binheader[ic], 0, sizeof (*(read_binheader[ic]))*(read_binheader_count));

            fprintf(stderr, "[%d] Partition File name %s\n", ic, existing_file_name);
            // file exists
            if ( access( partition_file_name, F_OK ) != -1 )
            {
              // contins data from the shared blocks
              read_data_buffer[ic] = malloc(samples_per_block * shared_block_count * bpv[var]/8);
              memset(read_data_buffer[ic], 0, samples_per_block * shared_block_count * bpv[var]/8);

              int fd;
              fd = open(existing_file_name, O_RDONLY | O_BINARY);
              if (fd < 0)
              {
                fprintf(stderr, "[File : %s] [Line : %d] open\n", __FILE__, __LINE__);
                continue;
                return 0;
              }

              // reads the header infor from binary file of the partitions
              ret = read(fd, read_binheader[ic], (sizeof (*(read_binheader[ic])) * read_binheader_count));
              if (ret < 0)
              {
                fprintf(stderr, "[File : %s] [Line : %d] read\n", __FILE__, __LINE__);
                return 0;
              }
              //assert(ret == (sizeof (*(read_binheader[ic])) * read_binheader_count));

              // copy the header from the partition file to the merged output file
              // do it only for first partition (this gets all info other than block offset nd count)
              if (ic == 0)
                memcpy(write_binheader, read_binheader[ic], 10 * sizeof (*write_binheader));

              int bpf = 0;
              size_t data_size = 0;
              off_t data_offset = 0;
              for (bpf = 0; bpf < shared_block_count; bpf++)
              {
                data_offset = ntohl(read_binheader[ic][(bpf + var * blocks_per_file)*10 + 12]);
                data_size = ntohl(read_binheader[ic][(bpf + var * blocks_per_file)*10 + 14]);
                fprintf(stderr, "[%s] [Partition %d Block %d Variable %d] --> Offset %d Count %d\n", partition_file_name, ic, bpf, var, (int)data_offset, (int)data_size);

                if (data_offset != 0 && data_size != 0)
                {
                  pread(fd, read_data_buffer[ic] + (bpf * samples_per_block * (bpv[var] / 8)), data_size, data_offset);

                  write_binheader[((bpf + var * blocks_per_file)*10 + 12)] = htonl(write_binheader_length + (bpf * data_size) + var * shared_block_count);
                  write_binheader[((bpf + var * blocks_per_file)*10 + 14)] = htonl(data_size);

                  // Merge happening while the shared block is being read
                  // Hardcoded stupid merge
                  // checks if value is not zero then copies to the write block
                  int m = 0;
                  for (m = 0; m < data_size / (bpv[var] / 8) ; m++)
                  {
                    double temp;
                    memcpy(&temp, read_data_buffer[ic] + (bpf * samples_per_block + m) * sizeof(double), sizeof(double));
                    if (temp != 0)
                      memcpy(write_data_buffer + ((bpf * samples_per_block) + m) * sizeof(double), &temp, sizeof(double));
                  }
                }
              }

              close(fd);
            }
            else
              continue;
          }

          //Merge after all the reads
          for (ic = 0; ic < idx_count[0] * idx_count[1] * idx_count[2]; ic++)
          {
            //input is read_data_buffer**
            //output is write_data_buffer*
          }

          if ( access( new_file_name, F_OK ) != -1 )
          {
            // file exists
            int fd;
            fd = open(new_file_name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
            {
            }
            close(fd);
          }
          else
          {
            // file doesn't exist
            /*
            int r;
            for (r = 0; r < (shared_block_count * samples_per_block * bpv[var]/8) / sizeof(double); r++)
            {
              double dval;
              memcpy(&dval, write_data_buffer + r * sizeof(double), sizeof(double));
              fprintf(stderr, "value at %d = %f\n", r, dval);
            }
            */

            int fd;
            fd = open(new_file_name, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
            pwrite(fd, write_binheader, sizeof (*write_binheader)*(write_binheader_count), 0);
            pwrite(fd, write_data_buffer, shared_block_count * samples_per_block * bpv[var]/8, sizeof (*write_binheader)*(write_binheader_count));
            close(fd);
          }
        }
      }
      #endif
    }
  }


  shutdown_mpi();

  end_time = get_time();
  fprintf(stderr, "Total time taken = %f %f\n", end_time, start_time);

  return 0;
}

static int get_default_bits_per_datatype(char* type, int* bits)
{
  if (strcmp(type, "1*float64") == 0)
    *bits = 64;
  else if (strcmp(type, "2*float64") == 0)
    *bits = 128;
  else if (strcmp(type, "3*float64") == 0)
    *bits = 192;
  else if (strcmp(type, "4*float64") == 0)
    *bits = 256;
  else
    *bits = 0;

  return 0;
}

static unsigned long long getPowerOf2(int x)
{
  /*  find the power of 2 of an integer value (example 5->8) */
  int n = 1;
  while (n < x) n <<= 1;
  return n;
}


static int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen)
{
  unsigned long long address = 0;
  unsigned int segs[MAX_TEMPLATE_DEPTH] = {0};
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
      if (seg_count >= MAX_TEMPLATE_DEPTH)
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

static double get_time()
{
#if PIDX_HAVE_MPI
  //return MPI_Wtime();
#else
  struct timeval temp;
  gettimeofday(&temp, NULL);
  return (double)(temp.tv_sec) + (double)(temp.tv_usec)/1000000.0;
#endif
}
