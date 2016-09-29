/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

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

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <fcntl.h>
#include <errno.h>
#include <PIDX.h>

enum { X, Y, Z, NUM_DIMS };
static int process_count = 1, rank = 0;
static int cores = 1;
static char output_file_template[512] = "test";
static char output_file_name[512] = "test.idx";
static char partition_file_name[512];
static char partition_file_template[512] = "test";
static int idx_count[PIDX_MAX_DIMENSIONS];

static int current_time_step = 0;
static int compression_type = PIDX_NO_COMPRESSION;
static int bits_per_block = PIDX_default_bits_per_block;
static int blocks_per_file = PIDX_default_blocks_per_file;
static int bounds[PIDX_MAX_DIMENSIONS];
static int chunked_bounds[PIDX_MAX_DIMENSIONS];
static int idx_size[PIDX_MAX_DIMENSIONS];
static int idx_count[PIDX_MAX_DIMENSIONS];
static char bitPattern[512];
static char bitSequence[512];
static int compression_bit_rate = 64;
static int compression_factor = 1;
static int chunk_size[PIDX_MAX_DIMENSIONS];
static int samples_per_block;
static int data_core_count;
static char var_name[128][512];
static char type_name[128][512];
static int bits_per_value[128];
static int values_per_sample[128];
static int variable_count = 0;
static int fs_block_size = 0;
static int maxh = 0;
static int max_file_count = 0;
static int start_layout_index_shared = 0;
static int end_layout_index_shared = 0;
static int layout_count_shared = 0;
static int start_layout_index_file_zero = 0;
static int end_layout_index_file_zero = 0;
static int layout_count_file_zero = 0;
static int start_time_step = 0;
static int end_time_step = 0;
static MPI_Comm comm;
static PIDX_block_layout global_block_layout_shared_files;
static PIDX_block_layout* local_block_layout_shared_files;
static PIDX_block_layout global_block_layout_file_zero;
static PIDX_block_layout* local_block_layout_file_zero;


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
#if PIDX_HAVE_MPI
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
#if PIDX_HAVE_MPI
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
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}





///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  char flags[] = "f:c:";
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

    case('c'): // number of timesteps
      if (sscanf(optarg, "%d", &cores) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

static PIDX_return_code IDX_file_open(const char* filename)
{
  int i;
  //int ret;
  int rank = 0;

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
    idx_count[i] = 1;

  current_time_step = 0;
  compression_type = PIDX_NO_COMPRESSION;

  bits_per_block = PIDX_default_bits_per_block;
  blocks_per_file = PIDX_default_blocks_per_file;

  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    bounds[i]=65535;

  memset(bitPattern, 0, 512);
  memset(bitSequence, 0, 512);

  compression_bit_rate = 64;
  compression_factor = 1;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    chunk_size[i] = 1;

  samples_per_block = (int)pow(2, PIDX_default_bits_per_block);

  int var = 0, variable_counter = 0, count = 0, len = 0;
  char *pch, *pch1;
  char line [ 512 ];

  if (rank == 0)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
      fprintf(stdout, "Error Opening %s\n", filename);
      return PIDX_err_file;
    }

    while (fgets(line, sizeof (line), fp) != NULL)
    {
      //printf("%s", line);
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "(box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
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
          return PIDX_err_file;
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
          return PIDX_err_file;
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
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        data_core_count = atoi(line);
      }

      if (strcmp(line, "(fields)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
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
              ret = PIDX_default_bits_per_datatype(type_name[variable_counter], &bits_per_sample);
              if (ret != PIDX_success)  return PIDX_err_file;

              bits_per_value[variable_counter] = bits_per_sample;
              values_per_sample[variable_counter] = 1;
            }
            count++;
            pch1 = strtok(NULL, " +");
          }
          count = 0;

          if( fgets(line, sizeof line, fp) == NULL)
            return PIDX_err_file;
          line[strcspn(line, "\r\n")] = 0;
          variable_counter++;
        }
        variable_count = variable_counter;
      }

      if (strcmp(line, "(bits)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        strcpy(bitSequence, line);
      }

      if (strcmp(line, "(bitsperblock)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        bits_per_block = atoi(line);
        samples_per_block = (int)pow(2, bits_per_block);
      }

      if (strcmp(line, "(compression type)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        compression_type = atoi(line);
        if (compression_type != PIDX_NO_COMPRESSION)
        {
          int i1 = 0;
          for (i1 = 0; i1 < PIDX_MAX_DIMENSIONS; i1++)
            chunk_size[i1] = 4;
        }
      }

      if (strcmp(line, "(compressed box)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
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
          return PIDX_err_box;
      }

      if (strcmp(line, "(compression bit rate)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        compression_bit_rate = atoi(line);
      }

      if (strcmp(line, "(blocksperfile)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
        blocks_per_file= atoi(line);
      }

      if (strcmp(line, "(filename_template)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
        line[strcspn(line, "\r\n")] = 0;
      }

      if (strcmp(line, "(time)") == 0)
      {
        if( fgets(line, sizeof line, fp) == NULL)
          return PIDX_err_file;
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

  printf("Start time step %d End time step %d\n", start_time_step, end_time_step);

#if PIDX_HAVE_MPI

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
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
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
#if PIDX_HAVE_MPI
    MPI_Bcast(&(bits_per_value[var]), 1, MPI_INT, 0, comm);
    MPI_Bcast(&(values_per_sample[var]), 1, MPI_INT, 0, comm);
    MPI_Bcast(var_name[var], 512, MPI_CHAR, 0, comm);
    MPI_Bcast(type_name[var], 512, MPI_CHAR, 0, comm);
#endif
  }


  unsigned long long total_reg_sample_count = (getPowerOf2(chunked_bounds[0]) * getPowerOf2(chunked_bounds[1]) * getPowerOf2(chunked_bounds[2]));
  if (total_reg_sample_count <= 0)
  {
    fprintf(stderr, "[%s] [%d ]File dimensions are wrong\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  unsigned long long max_sample_per_file = (unsigned long long) samples_per_block * blocks_per_file;
  if (max_sample_per_file <= 0)
  {
    fprintf(stderr, "[%s] [%d ]IDX dimensions are wrong %d %d\n", __FILE__, __LINE__, samples_per_block, blocks_per_file);
    return PIDX_err_file;
  }

  max_file_count = total_reg_sample_count / max_sample_per_file;
  if (total_reg_sample_count % max_sample_per_file)
    max_file_count++;

  //printf("%d %d %d %d %d\n", chunk_size[0], chunk_size[1], chunk_size[2]);

  if (rank == 0)
  {
    int ret;
    struct stat stat_buf;
    ret = stat(filename, &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    fs_block_size = stat_buf.st_blksize;
  }

#if PIDX_HAVE_MPI
  MPI_Bcast(&(fs_block_size), 1, MPI_INT, 0, comm);
#endif

  return PIDX_success;
}



PIDX_return_code file_initialize_time_step(int current_time_step, char* file_name, char* file_template)
{
  int N;
  char dirname[1024], basename[1024];
  int nbits_blocknumber;
  char *directory_path;
  char *data_set_path;

  data_set_path = malloc(sizeof(*data_set_path) * 1024);
  memset(data_set_path, 0, sizeof(*data_set_path) * 1024);

  directory_path = malloc(sizeof(*directory_path) * 1024);
  memset(directory_path, 0, sizeof(*directory_path) * 1024);

  strncpy(directory_path, file_name, strlen(file_name) - 4);
  sprintf(data_set_path, "%s/time%09d.idx", directory_path, current_time_step);
  free(directory_path);

  nbits_blocknumber = (maxh - bits_per_block - 1);
  VisusSplitFilename(data_set_path, dirname, basename);

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
  return PIDX_success;
}


int main(int argc, char **argv)
{
  int ret = 0;
  init_mpi(argc, argv);
  parse_args(argc, argv);

  rank_0_print("Merge Program\n");

  comm = MPI_COMM_WORLD;
  ret = IDX_file_open(output_file_name);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

  maxh = strlen(bitSequence);

  printf("Partition size and count: %d %d %d :: %d %d %d maxh = %d\n", idx_count[0], idx_count[1], idx_count[2], idx_size[0], idx_size[1], idx_size[2], maxh);

  int shared_block_level = (int)log2(idx_count[0] * idx_count[1] * idx_count[2]) + bits_per_block + 1;
  if (shared_block_level >= maxh)
    shared_block_level = maxh;

  int partion_level = (int) log2(idx_count[0] * idx_count[1] * idx_count[2]);
  int total_partiton_level = bits_per_block + (int)log2(blocks_per_file) + 1 + partion_level;
  if (total_partiton_level >= maxh)
    total_partiton_level = maxh;

  int level = 0;
  int ts = 0;
  for (ts = start_time_step; ts <= end_time_step; ts++)
  {
    for (level = shared_block_level; level < total_partiton_level; level = level + 1)
    {
      int hz_index = (int)pow(2, level - 1);
      int file_no = hz_index / (blocks_per_file * (int)pow(2, bits_per_block));
      int file_count;
      char existing_file_name[1024];
      char new_file_name[1024];
      int ic = 0;
      if (level <= bits_per_block + log2(blocks_per_file) + 1)
        file_count = 1;
      else
        file_count = (int)pow(2, level - (bits_per_block + log2(blocks_per_file) + 1));

      printf("File Index and Count at %d = %d %d\n", level, file_no, file_count);

      int fc = 0;
      for (fc = file_no; fc < file_no + file_count; fc++)
      {
        uint32_t* write_binheader;
        int write_binheader_count;
        write_binheader_count = 10 + 10 * blocks_per_file * variable_count;
        write_binheader = (uint32_t*) malloc(sizeof (*write_binheader)*(write_binheader_count));

        int var = 0;
        for (var = 0; var < variable_count; var++)
        {
          unsigned char *write_data_buffer = malloc((int)pow(2, bits_per_block) * blocks_per_file * bits_per_value[var]);
          memset(write_data_buffer, 0, (int)pow(2, bits_per_block) * blocks_per_file * bits_per_value[var]);

          int write_block_counter = 0;

          file_initialize_time_step(ts, output_file_name, output_file_template);
          generate_file_name(blocks_per_file, output_file_template, fc, new_file_name, PATH_MAX);

          int* block_header_offset = malloc(idx_count[0] * idx_count[1] * idx_count[2] * sizeof(int));
          memset(block_header_offset, 0, idx_count[0] * idx_count[1] * idx_count[2] * sizeof(int));


          int** block_header_block_offset = malloc(idx_count[0] * idx_count[1] * idx_count[2] * sizeof(int*));
          memset(block_header_block_offset, 0, idx_count[0] * idx_count[1] * idx_count[2] * sizeof(int*));

          size_t totl_data_size = 0;
          off_t totl_data_offset = 0;
          for (ic = 0; ic < idx_count[0] * idx_count[1] * idx_count[2]; ic++)
          {
            char file_name_skeleton[1024];
            strncpy(file_name_skeleton, output_file_name, strlen(output_file_name) - 4);
            file_name_skeleton[strlen(output_file_name) - 4] = '\0';

            if (idx_count[0] != 1 || idx_count[1] != 1 || idx_count[2] != 1)
              sprintf(partition_file_name, "%s_%d.idx", file_name_skeleton, ic);
            else
              strcpy(partition_file_name, output_file_name);

            file_initialize_time_step(ts, partition_file_name, partition_file_template);
            generate_file_name(blocks_per_file, partition_file_template, fc, existing_file_name, PATH_MAX);

            printf("[%d] OLD File name %s NEW File name %s\n", ic, existing_file_name, new_file_name);

            if ( access( partition_file_name, F_OK ) != -1 )
            {
              // file exists
              uint32_t* read_binheader;
              uint32_t adjusted_offset;
              int read_binheader_count;

              block_header_block_offset[ic] = malloc(blocks_per_file * sizeof(int));
              memset(block_header_block_offset[ic], 0, blocks_per_file * sizeof(int));

              read_binheader_count = 10 + 10 * blocks_per_file * variable_count;
              read_binheader = (uint32_t*) malloc(sizeof (*read_binheader)*(read_binheader_count));
              memset(read_binheader, 0, sizeof (*read_binheader)*(read_binheader_count));

              unsigned char *read_data_buffer = malloc((int)pow(2, bits_per_block) * blocks_per_file * bits_per_value[var]);
              memset(read_data_buffer, 0, (int)pow(2, bits_per_block) * blocks_per_file * bits_per_value[var]);

              int fd;
              fd = open(existing_file_name, O_RDONLY);
              if (fd < 0)
              {
                fprintf(stderr, "[File : %s] [Line : %d] open\n", __FILE__, __LINE__);
                continue;
                return 0;
              }

              ret = read(fd, read_binheader, (sizeof (*read_binheader) * read_binheader_count));
              if (ret < 0)
              {
                fprintf(stderr, "[File : %s] [Line : %d] read\n", __FILE__, __LINE__);
                return 0;
              }
              assert(ret == (sizeof (*read_binheader) * read_binheader_count));

              memcpy(write_binheader, read_binheader, 10 * sizeof (*write_binheader));

              int bpf = 0;
              size_t data_size = 0;
              off_t data_offset = 0;
              int read_block_counter = 0;

              int prev_ic = 0;
              int previous_block_offset = 0;
              for (prev_ic = 0; prev_ic < ic; prev_ic++)
                previous_block_offset = previous_block_offset + block_header_offset[prev_ic];

              for (bpf = 0; bpf < blocks_per_file; bpf++)
              {
                data_offset = ntohl(read_binheader[(bpf + var * blocks_per_file)*10 + 12]);
                data_size = ntohl(read_binheader[(bpf + var * blocks_per_file)*10 + 14]);
                printf("[%s] [%d %d %d] --> %d %d\n", partition_file_name, bpf, var, blocks_per_file, (int)data_offset, (int)data_size);

                if (data_offset != 0 && data_size != 0)
                {
                  //printf("[%d] --> %d %d\n", bpf, (int)data_offset, (int)data_size);
                  pread(fd, read_data_buffer + (read_block_counter * (int)pow(2, bits_per_block) * (bits_per_value[var] / 8)), data_size, data_offset);

                  adjusted_offset = data_offset + previous_block_offset;

                  write_binheader[((bpf + var * blocks_per_file)*10 + 12)] = htonl(adjusted_offset);
                  write_binheader[((bpf + var * blocks_per_file)*10 + 14)] = htonl(data_size);

                  printf("IC %d RBC %d WBC %d RO %d WO %d AO %d [%d + %d]\n", ic, read_block_counter, write_block_counter, (read_block_counter * (int)pow(2, bits_per_block) * (bits_per_value[var] / 8)), (write_block_counter * (int)pow(2, bits_per_block) * (bits_per_value[var] / 8)), adjusted_offset, (int)data_offset, previous_block_offset);
                  memcpy(write_data_buffer + (write_block_counter * (int)pow(2, bits_per_block) * (bits_per_value[var] / 8)), read_data_buffer + (read_block_counter * (int)pow(2, bits_per_block) * (bits_per_value[var] / 8)), (int)pow(2, bits_per_block) * (bits_per_value[var] / 8));

                  read_block_counter++;
                  write_block_counter++;

                  totl_data_size = totl_data_size + data_size;

                  block_header_offset[ic] = block_header_offset[ic] + data_size;
                  block_header_block_offset[ic][bpf] = adjusted_offset;
                }
              }

              close(fd);

              free(read_binheader);
            }
            else
              continue;
          }

          printf("Write Filename %s\n", new_file_name);
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
            off_t first_offset = 0;
            int break_counter = 0;
            for (ic = 0; ic < idx_count[0] * idx_count[1] * idx_count[2]; ic++)
            {
              int bpf = 0;
              for (bpf = 0; bpf < blocks_per_file; bpf++)
              {
                if (block_header_block_offset[ic][bpf] != 0)
                {
                  first_offset = block_header_block_offset[ic][bpf];
                  break_counter = 1;
                  break;
                }
              }
              if (break_counter == 1)
                break;
            }

            // file doesn't exist
            int fd;
            fd = open(new_file_name, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
            printf("Write File %s Offset %d Count %d\n", new_file_name, (int)first_offset, (int)totl_data_size);

            pwrite(fd, write_binheader, sizeof (*write_binheader)*(write_binheader_count), 0);
            pwrite(fd, write_data_buffer, totl_data_size, first_offset);
            close(fd);
          }
        }
      }
    }
  }

  shutdown_mpi();
  return 0;
}
