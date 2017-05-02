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

#define _XOPEN_SOURCE 600

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <math.h>

//#define PIDX_HAVE_ZFP 1
//#include "PIDX_config.h"

#if PIDX_HAVE_ZFP
  #include "zfp.h"
#endif

#define PIDX_MAX_DIMENSIONS 5
#define MAX_VARIABLE_COUNT 1024
#define MAX_TEMPLATE_DEPTH 6

#define Min2ab(a,b)      (((a)<=(b))?(a):(b))
#define Max2ab(a,b)      (((a)> (b))?(a):(b))
#define OffsetFor(_D_,_From_,_Off_) for ((_D_)=(_From_);(_D_)<(PIDX_MAX_DIMENSIONS+(_Off_));(_D_)++)
#define For(_D_) for ((_D_)=0;(_D_)<PIDX_MAX_DIMENSIONS;(_D_)++)
#define PGET(_Point_,_Coordinate_) ((&((_Point_).x))[(_Coordinate_)])

typedef struct {int x,y,z,u,v;} PointND;
struct block_layout_struct
{
  int levels;                          // Total number of Levels
  int *hz_block_count_array;           // Number of filled blocks per level
  int ** hz_block_number_array;        // Indices of filled blocks
};
typedef struct block_layout_struct* block_layout;

static int compression_block_size[PIDX_MAX_DIMENSIONS] = {1, 1, 1, 1, 1};
static int compression_bit_rate = 0;
static int compression_type = 0;

//static void revstr(char* str);
//static void GuessBitmaskPattern(char* _bits, PointND dims);
static int generate_file_name_template(int maxh, int bits_per_block, char* filename, int current_time_step, char* filename_template);
static int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen);
static int is_block_present(int block_number, block_layout layout);
static void destroyBlockBitmap(block_layout layout);
static int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, int res, const char* bitPattern, block_layout layout);
static int VisusSplitFilename(const char* filename,char* dirname,char* basename);
static void Hz_to_xyz(const char* bitmask,  int maxh, unsigned long long hzaddress, unsigned long long* xyz);
static int RegExBitmaskBit(const char* bitmask_pattern,int N);
static unsigned long long getPowerOf2(int x);
static int decompress(float* input_buffer, float* output_buffer, size_t buffer_size);

double swap(double d)
{
   double a;
   unsigned char *dst = (unsigned char *)&a;
   unsigned char *src = (unsigned char *)&d;

   dst[0] = src[7];
   dst[1] = src[6];
   dst[2] = src[5];
   dst[3] = src[4];
   dst[4] = src[3];
   dst[5] = src[2];
   dst[6] = src[1];
   dst[7] = src[0];

   return a;
}

float swap_float(float d)
{
   float a;
   unsigned char *dst = (unsigned char *)&a;
   unsigned char *src = (unsigned char *)&d;

   dst[0] = src[3];
   dst[1] = src[2];
   dst[2] = src[1];
   dst[3] = src[0];

   return a;
}

static int decompress(float* input_buffer, float* output_buffer, size_t buffer_size)
{
#if PIDX_HAVE_ZFP
  int i;
  zfp_params params;
  size_t typesize, outsize;

  typesize = sizeof(float);
  params.type = ZFP_TYPE_FLOAT;
  params.nx = compression_block_size[0];
  params.ny = compression_block_size[1];
  params.nz = compression_block_size[2];
  zfp_set_rate(&params, compression_bit_rate);
  
  outsize = compression_block_size[0] * compression_block_size[1] * compression_block_size[2] * typesize;
  //fprintf(stderr, "(Buffer size %ld) (outsize %ld) (compression_bit_rate %d)\n", buffer_size, outsize, compression_bit_rate);

  unsigned char* zip = (unsigned char*)malloc(outsize);
  for (i = 0; i < buffer_size; i = i + (compression_block_size[0] * compression_block_size[1] * compression_block_size[2] * typesize)/(32/compression_bit_rate))
  {
    //fprintf(stderr, "%d \n", (i*(32/compression_bit_rate)));
    zfp_decompress(&params, zip, (unsigned char*)input_buffer + (i), outsize/(32/compression_bit_rate));
    memcpy((unsigned char*)output_buffer + (i*(32/compression_bit_rate)), zip, outsize);
  }
  free(zip);
  
#endif
  return 0;
}

int main(int argc, char **argv)
{
  if ( argv[1] == NULL )
  {
    fprintf(stderr, "Missing arguments (provide a .idx file as command line input argument)\n");
    return 0;
  }

  int j, k, t, var, ret, counter = 0;
  char *pch, *pch1;
  char line [ 512 ];
  int count = 0, len = 0;
  int max_files;
  int maxh;
  char bitSequence[512];
  char bitPattern[512];
  int bits_per_block = 0;
  int samples_per_block = 0;
  int blocks_per_file = 0;
  int start_time_step = 0, end_time_step = 1;
  int global_bounds[PIDX_MAX_DIMENSIONS];
  int compressed_global_bounds[PIDX_MAX_DIMENSIONS];
  int vps[MAX_VARIABLE_COUNT];
  char variable_name[MAX_VARIABLE_COUNT][1024];
  char variable_type[MAX_VARIABLE_COUNT][1024];
  char filename_template[1024];
  int variable_count = 0;
  int resolution = 0;

  FILE *fp = fopen(argv[1], "r");
  while (fgets(line, sizeof (line), fp) != NULL)
  {
    //fprintf(stderr, "%s", line);
    len = strlen(line) - 1;
    if (line[len] == '\n')
      line[len] = 0;

    if (strcmp(line, "(version)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;
    }

    if (strcmp(line, "(box)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
	if (count % 2 == 1)
	  global_bounds[count / 2] = atoi(pch) + 1;
	count++;
	pch = strtok(NULL, " ");
      }
    }
    if (strcmp(line, "(fields)") == 0)
    {
      if( fgets(line, sizeof (line), fp) == NULL)
        exit(0);
      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      count = 0;
      variable_count = 0;

      while (strcmp(line, "(logic_to_physic)") != 0 && strcmp(line, "(version)") != 0 && strcmp(line, "(box)") != 0 && strcmp(line, "(bits)") && strcmp(line, "(bitsperblock)") != 0 && strcmp(line, "(blocksperfile)") != 0 && strcmp(line, "(filename_template)") != 0 && strcmp(line, "(time)") != 0 && strcmp(line, "(compression bit rate)") != 0)
      {
        pch1 = strtok(line, " *+");
        while (pch1 != NULL)
        {
          if (count == 0)
            strcpy(variable_name[variable_count], strdup(pch1));

          if (count == 1)
          {
            vps[variable_count] = atoi(pch1);
            fprintf(stderr, "vps[variable_count] = %d\n", vps[variable_count]);
          }

          if (count == 2)
          {
            len = strlen(pch1) - 1;
            if (pch1[len] == '\n')
              pch1[len] = 0;

            if (strcmp(pch1, "float64") == 0)
              strcpy (variable_type[variable_count], "float64");
            else if (strcmp(pch1, "uint64") == 0)
              strcpy (variable_type[variable_count], "uint64");
            else if (strcmp(pch1, "float32") == 0)
              strcpy (variable_type[variable_count], "float32");
            else if (strcmp(pch1, "int32") == 0)
              strcpy (variable_type[variable_count], "int32");
            else if (strcmp(pch1, "3*float64") == 0)
            {
              strcpy (variable_type[variable_count], "3*float64");
            }
            else
            {
              fprintf(stderr, "Currently supporting only float64 types\n");
              return 0;
            }
          }
          count++;
          pch1 = strtok(NULL, " *+");
        }
        count = 0;

        if( fgets(line, sizeof (line), fp) == NULL)
          exit(0);
        len = strlen(line) - 1;
        if (line[len] == '\n')
          line[len] = 0;
        variable_count++;
      }
    }

    if (strcmp(line, "(compressed box)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;

      pch = strtok(line, " ");
      count = 0;
      while (pch != NULL)
      {
        //if (count % 2 == 1)
        compression_block_size[count]/*global_bounds[count / 2]*/ = atoi(pch);
        //fprintf(stderr, "compressed box size %d %d\n", comp, count);
        count++;
        pch = strtok(NULL, " ");
      }
    }

    if (strcmp(line, "(bits)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;

      strcpy(bitSequence, line);

    }

    if (strcmp(line, "(bitsperblock)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      bits_per_block = atoi(line);
      samples_per_block = pow(2, bits_per_block);
    }

    if (strcmp(line, "(resolution)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      resolution = atoi(line);
    }
    
    if (strcmp(line, "(compression bit rate)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      compression_bit_rate = atoi(line);
    }

    if (strcmp(line, "(compression type)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      compression_type = atoi(line);
    }

    if (strcmp(line, "(blocksperfile)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      blocks_per_file= atoi(line);
    }

    if (strcmp(line, "(filename_template)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
      line[len] = 0;
    }

    if (strcmp(line, "(time)") == 0)
    {
      if( fgets(line, sizeof line, fp) == NULL)
        return 0;

      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;

      pch1 = strtok(line, " ");
      count = 0;
      while (pch1 != NULL)
      {
        //fprintf(stderr, "");
        if (count == 0)
          start_time_step = atoi(pch1);

        if (count == 1)
          end_time_step = atoi(pch1);

        count++;
        pch1 = strtok(NULL, " ");
      }
      count = 0;
    }
  }
  fclose(fp);

  fprintf(stderr, "Finished Parsing %s\n", argv[1]);
  fprintf(stderr, "Starting time step %d and Ending time step %d\n", start_time_step, end_time_step);
  fprintf(stderr, "(box)\n0 %d 0 %d 0 %d 0 %d 0 %d\n", global_bounds[0], global_bounds[1], global_bounds[2], global_bounds[3], global_bounds[4]);
  fprintf(stderr, "compressed box size %d %d %d %d %d\n", compression_block_size[0], compression_block_size[1], compression_block_size[2], compression_block_size[3], compression_block_size[4]);
  fprintf(stderr, "compressed bit rate %d\n", compression_bit_rate);
  fprintf(stderr, "(fields)\n");
  for(var = 0 ; var < variable_count ; var++)
  {
      fprintf(stderr, "%s %d*", variable_name[var], vps[var]);
      fprintf(stderr, "float64 ");
      if(var != variable_count - 1)
        fprintf(stderr, " + \n");
  }
  fprintf(stderr, "\n(bitsperblock)\n%d\n(blocksperfile)\n%d\n", bits_per_block, blocks_per_file);
  fprintf(stderr, "(filename)\n%s\n", argv[1]);
  
  if (global_bounds[0] % compression_block_size[0] == 0)
    compressed_global_bounds[0] = (int) global_bounds[0] / compression_block_size[0];
  else
    compressed_global_bounds[0] = (int) (global_bounds[0] / compression_block_size[0]) + 1;

  if (global_bounds[1] % compression_block_size[1] == 0)
    compressed_global_bounds[1] = (int) global_bounds[1] / compression_block_size[1];
  else
    compressed_global_bounds[1] = (int) (global_bounds[1] / compression_block_size[1]) + 1;

  if (global_bounds[2] % compression_block_size[2] == 0)
    compressed_global_bounds[2] = (int) global_bounds[2] / compression_block_size[2];
  else
    compressed_global_bounds[2] = (int) (global_bounds[2] / compression_block_size[2]) + 1;

  if (global_bounds[3] % compression_block_size[3] == 0)
    compressed_global_bounds[3] = (int) global_bounds[3] / compression_block_size[3];
  else
    compressed_global_bounds[3] = (int) (global_bounds[3] / compression_block_size[3]) + 1;

  if (global_bounds[4] % compression_block_size[4] == 0)
    compressed_global_bounds[4] = (int) global_bounds[4] / compression_block_size[4];
  else
    compressed_global_bounds[4] = (int) (global_bounds[4] / compression_block_size[4]) + 1;

  unsigned long long total_compression_block_size = compression_block_size[0] * compression_block_size[1] * compression_block_size[2] * compression_block_size[3] * compression_block_size[4];

  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  int i1 = 0;
  for (i1 = 0; i1 < PIDX_MAX_DIMENSIONS; i1++)
  {
    bounding_box[0][i1] = 0;
    bounding_box[1][i1] = compressed_global_bounds[i1];
  }


  /*
  PointND extents;
  extents.x = compressed_global_bounds[0];
  extents.y = compressed_global_bounds[1];
  extents.z = compressed_global_bounds[2];
  extents.u = compressed_global_bounds[3];
  extents.v = compressed_global_bounds[4];
  GuessBitmaskPattern(bitSequence, extents);
  */

  maxh = strlen(bitSequence);

  fprintf(stderr, "bitsequence = %s [%d]\n", bitSequence, maxh);
  for (i1 = 0; i1 <= maxh; i1++)
    bitPattern[i1] = RegExBitmaskBit(bitSequence, i1);

  block_layout global_block_layout =  (block_layout)malloc(sizeof (*global_block_layout));
  memset(global_block_layout, 0, sizeof (*global_block_layout));

  createBlockBitmap(bounding_box, blocks_per_file, bits_per_block, maxh, resolution, bitPattern, global_block_layout);

  k = 1;
  for (i1 = 1; i1 < (global_block_layout->levels); i1++)
  {
    counter = 0;
    for(j = 0 ; j < k ; j++)
    {
      if(global_block_layout->hz_block_number_array[i1][j] != 0)
      {
        global_block_layout->hz_block_number_array[i1][counter] = global_block_layout->hz_block_number_array[i1][j];
        counter++;
      }
    }
    k = k * 2;
  }

  /// maximum number of files possible
  max_files = (getPowerOf2(compressed_global_bounds[0]) * getPowerOf2(compressed_global_bounds[1]) * getPowerOf2(compressed_global_bounds[2]) * getPowerOf2(compressed_global_bounds[3]) * getPowerOf2(compressed_global_bounds[4])) / ((unsigned long long) pow(2, bits_per_block) * (unsigned long long) blocks_per_file);
  if ((getPowerOf2(compressed_global_bounds[0]) * getPowerOf2(compressed_global_bounds[1]) * getPowerOf2(compressed_global_bounds[2]) * getPowerOf2(compressed_global_bounds[3]) * getPowerOf2(compressed_global_bounds[4])) % ((unsigned long long) pow(2, bits_per_block) * (unsigned long long) blocks_per_file))
    max_files++;
  assert(max_files != 0);

  int file_no;
  int *existing_file_index;
  existing_file_index = (int*) malloc(max_files * sizeof (int));
  memset(existing_file_index, 0, max_files * sizeof (int));
  existing_file_index[0] = 1;

  for (i1 = 1; i1 < global_block_layout->levels; i1++)
  {
    for (j = 0; j < global_block_layout->hz_block_count_array[i1]; j++)
    {
      file_no = global_block_layout->hz_block_number_array[i1][j] / blocks_per_file;
      existing_file_index[file_no] = 1;
    }
  }

  uint32_t* binheader;
  int binheader_count;
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader = (uint32_t*) malloc(sizeof (*binheader)*(binheader_count));

  int fd;
  unsigned long long ZYX[PIDX_MAX_DIMENSIONS];

  for (t = start_time_step; t <= end_time_step; t++)
  {
    unsigned long long lost_element_count = 0, element_count = 0, element_count1 = 0, lost_element_count1= 0;
    int fi = 0;
    for (fi = 0; fi < max_files; fi++)
    //for (i = 3; i < 4; i++)
    {
      if (existing_file_index[fi] == 1)
      {
        generate_file_name_template(maxh, bits_per_block, argv[1], t, filename_template);

        char bin_file[PATH_MAX];
        ret = generate_file_name(blocks_per_file, filename_template, fi, bin_file, PATH_MAX);
        if (ret == -1)
        {
          fprintf(stderr, "[File : %s] [Line : %d] generate_file_name\n", __FILE__, __LINE__);
          return 0;
        }

        fprintf(stderr, "Parsing File : %s\n", bin_file);

        fd = open(bin_file, O_RDONLY);
        if (fd < 0)
        {
          fprintf(stderr, "[File : %s] [Line : %d] open %s\n", __FILE__, __LINE__, bin_file);
          continue;
          return 0;

        }

        ret = read(fd, binheader, (sizeof (*binheader) * binheader_count));
        if (ret < 0)
        {
          fprintf(stderr, "[File : %s] [Line : %d] read\n", __FILE__, __LINE__);
          return 0;
        }

        if (ret < (int) (sizeof (*binheader) * binheader_count))
        {
          fprintf(stderr, "[File : %s] [Line : %d] read : size\n", __FILE__, __LINE__);
          return 0;
        }


        int bpf = 0;
        //unsigned long long* long_long_buffer = NULL;
        double* double_buffer = NULL;
        float* float_buffer = NULL;
        unsigned long long* ulong_buffer = NULL;
        int32_t* iint_buffer = NULL;
        double* decompressed_double_buffer = NULL;
        float* decompressed_float_buffer = NULL;
        int check_bit = 1, s = 0;
        unsigned long long hz_index, hz_val;

        size_t data_size;
        off_t data_offset;
        for (var = 0; var < variable_count; var++)
        {
          for (bpf = 0; bpf < blocks_per_file; bpf++)
          {
            if (is_block_present((bpf + (blocks_per_file * fi)), global_block_layout))
            {
              data_offset = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 12]);
              data_size = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 14]);

              fprintf(stderr, "[F %d %s] [V %d] [B %d] Offset %lld Count %ld\n", fi, bin_file, var, bpf, (long long)data_offset, (long)data_size);

              int sample_size = 0;
              if (strcmp(variable_type[var], "float64") == 0)
              {
                double_buffer = (double*)malloc(data_size);
                memset(double_buffer, 0, data_size);

                ret = pread(fd, double_buffer, data_size, data_offset);
                //assert(ret == data_size);
              
                //if (compression_type == 2)
                //{
                //  decompressed_double_buffer = (double*) malloc(data_size * (64/compression_bit_rate));
                //  decompress(double_buffer, decompressed_double_buffer, data_size);
                //  free(double_buffer);
                //}
                sample_size = sizeof(double);
              }
              else if (strcmp(variable_type[var], "uint64") == 0)
              {
                ulong_buffer = (unsigned long long*)malloc(data_size);
                memset(ulong_buffer, 0, data_size);

                ret = pread(fd, ulong_buffer, data_size, data_offset);
                assert(ret == data_size);
                sample_size = sizeof(unsigned long long);
              }
              else if (strcmp(variable_type[var], "int32") == 0)
              {
                iint_buffer = (int32_t*)malloc(data_size);
                memset(iint_buffer, 0, data_size);

                ret = pread(fd, iint_buffer, data_size, data_offset);
                assert(ret == data_size);
                sample_size = sizeof(int32_t);
              }
              else if (strcmp(variable_type[var], "float32") == 0)
              {
                float_buffer = (float*)malloc(data_size);
                memset(float_buffer, 0, data_size);

                ret = pread(fd, float_buffer, data_size, data_offset);
                //assert(ret == data_size);
                sample_size = sizeof(float);

                if (compression_type == 2)
                {
                  //fprintf(stderr, "buffer size = %d x %d - %d\n", data_size, (32/compression_bit_rate), data_size * (32/compression_bit_rate));
                  decompressed_float_buffer = (float*) malloc(data_size * (32/compression_bit_rate));
                  if (decompressed_float_buffer == NULL)
                    fprintf(stderr, "Error!!!!\n");

                  decompress(float_buffer, decompressed_float_buffer, data_size);
                  free(float_buffer);
                  sample_size = compression_bit_rate/8;
                }

              }


              //fprintf(stderr, "Block %d - %d\n", bpf, data_size/(sample_size * vps[var] * total_compression_block_size));
              for (hz_val = 0; hz_val < data_size/(sample_size * vps[var] * total_compression_block_size); hz_val++)
              {
                hz_index = (blocks_per_file * fi * samples_per_block) + (bpf * samples_per_block) + hz_val;
                Hz_to_xyz(bitPattern, maxh - 1, hz_index, ZYX);

                if (ZYX[2] >= compressed_global_bounds[2] || ZYX[1] >= compressed_global_bounds[1] || ZYX[0] >= compressed_global_bounds[0])
                  continue;

                unsigned long long llhs, lrhs;
                double dlhs, drhs;
                float flhs, frhs;
                int iilhs, iirhs;
                check_bit = 1, s = 0;
                int i = 0, j = 0, k = 0, index = 0;
                for (k = 0; k < compression_block_size[2]; k++)
                {
                  for (j = 0; j < compression_block_size[1]; j++)
                  {
                    for (i = 0; i < compression_block_size[0]; i++)
                    {
                      index = (compression_block_size[0] * compression_block_size[1] * k) + (compression_block_size[0] * j) + i;
                      for (s = 0; s < vps[var]; s++)
                      {
                        int value_index = index + (ZYX[2] * compressed_global_bounds[0] * compressed_global_bounds[1] + ZYX[1] * compressed_global_bounds[0] + ZYX[0]) * total_compression_block_size;
                        int block_index = value_index / total_compression_block_size;
                        int block_index_x = block_index % compressed_global_bounds[0];
                        int block_index_y = ((block_index - block_index_x) % (compressed_global_bounds[0] * compressed_global_bounds[1])) / compressed_global_bounds[0];
                        int block_index_z = (block_index - block_index_x - block_index_y * compressed_global_bounds[0]) / (compressed_global_bounds[0] * compressed_global_bounds[1]);
                        int local_index = value_index % total_compression_block_size; // IMPORTANT: only works if the total size % compressed_block_size == 0
                        int local_index_x = local_index % compression_block_size[0];
                        int local_index_y = ((local_index - local_index_x) % (compression_block_size[0] * compression_block_size[1])) / compression_block_size[0];
                        int local_index_z = (local_index - local_index_x - local_index_y * compression_block_size[0]) / (compression_block_size[0] * compression_block_size[1]);
                        int index_x = block_index_x * compression_block_size[0] + local_index_x;
                        int index_y = block_index_y * compression_block_size[1] + local_index_y;
                        int index_z = block_index_z * compression_block_size[2] + local_index_z;

                        if (index_x >= global_bounds[0] || index_y >= global_bounds[1] || index_z >= global_bounds[2])
                          continue;

                        //fprintf(stderr, "compression_type = %d\n", compression_type);
                        if (strcmp(variable_type[var], "float64") == 0)
                        {

                          drhs =  100 + var + ((global_bounds[0] * global_bounds[1] * index_z)+(global_bounds[0]*index_y) + index_x);
                          if (compression_type == 0 || compression_type == 1)
                            dlhs = double_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s];
                          else
                            dlhs = decompressed_double_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s];

                          check_bit = check_bit && (dlhs == drhs);

                          //fprintf(stderr, "%d:  %f %f\n", var, dlhs, drhs);
                          //fprintf(stderr, "[value at %d %d %d] is %f\n", (int)ZYX[0], (int)ZYX[1], (int)ZYX[2], dlhs);
                          if (dlhs == drhs)
                          {
                            //fprintf(stderr, "[C] Expected %f Found %f\n", drhs, dlhs);
                            element_count1++;
                          }
                          //else
                            //fprintf(stderr, "[W] [F %d %s] [B %d] Expected %f Found %f\n", fi, bin_file, bpf, drhs, dlhs);
                        }
                        else if (strcmp(variable_type[var], "uint64") == 0)
                        {
                          lrhs = 100 + var + ((global_bounds[0] * global_bounds[1] * index_z)+(global_bounds[0]*index_y) + index_x);
                          llhs = ulong_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s];

                          check_bit = check_bit && (llhs == lrhs);

                          if (llhs == lrhs)
                            element_count1++;
                        }
                        else if (strcmp(variable_type[var], "int32") == 0)
                        {
                          iirhs = 100 + var + ((global_bounds[0] * global_bounds[1] * index_z)+(global_bounds[0]*index_y) + index_x);
                          iilhs = iint_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s];
                          check_bit = check_bit && (iilhs == iirhs);

                          //fprintf(stderr, "%d:  %d %d\n", var, iilhs, iirhs);
                          if (iilhs == iirhs)
                            element_count1++;
                        }
                        else if (strcmp(variable_type[var], "float32") == 0)
                        {
                          frhs = 100 + var + ((global_bounds[0] * global_bounds[1] * index_z)+(global_bounds[0]*index_y) + index_x);
                          if (compression_type == 2)
                            flhs = ((decompressed_float_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s]));
                          else
                            flhs = ((float_buffer[((hz_val * total_compression_block_size) + index) * vps[var] + s]));

                          check_bit = check_bit && (flhs == frhs);

                          if (flhs == frhs)
                          {
                            element_count1++;
                            //fprintf(stderr, "A [%d]: %f %f\n", bpf, flhs, frhs);
                          }
                          else
                          {
                            lost_element_count1++;
                            //fprintf(stderr, "B [%d]: %f %f\n", bpf, flhs, frhs);
                          }
                        }
                      }
                    }
                  }
                }

                if (check_bit == 0)
                {
                  //
                  //if (strcmp(variable_type[var], "float64") == 0)
                  //  fprintf(stderr, "%f %f\n", dlhs, drhs);
                  //else if (strcmp(variable_type[var], "uint64") == 0)
                  //    fprintf(stderr, "%lld %lld\n", (unsigned long long)llhs, (unsigned long long)lrhs);
                   //
                  //if (strcmp(variable_type[var], "float32") == 0)
                  //  fprintf(stderr, "%f %f\n", flhs, frhs);
                  lost_element_count++;
                  //break;
                }
                else
                {
                  element_count++;
                }

              }


              if (strcmp(variable_type[var], "float64") == 0)
              {
                if (compression_type == 1)
                {
                  free(decompressed_double_buffer);
                  decompressed_double_buffer = 0;
                }
                else
                {
                  free(double_buffer);
                  double_buffer = 0;
                }
              }
              else if (strcmp(variable_type[var], "float64") == 0)
              {
                free(ulong_buffer);
                ulong_buffer = 0;
              }
              else if (strcmp(variable_type[var], "float32") == 0)
              {
                if (compression_type == 2)
                {
                  free(decompressed_float_buffer);
                  decompressed_float_buffer = 0;
                }
                else
                {
                  free(float_buffer);
                  float_buffer = 0;
                }
              }
            }
          }
        }
      }
    }

    fprintf(stderr, "[=]%lld (%lld) + [!=]%lld (%lld) [%lld : %lld]\n", (long long) (element_count), (long long) (element_count1), (long long)lost_element_count, (long long)lost_element_count1, (long long) element_count + lost_element_count, (long long) compressed_global_bounds[0] * compressed_global_bounds[1] * compressed_global_bounds[2] * compressed_global_bounds[3] * compressed_global_bounds[4] * variable_count / (long long)pow(2, resolution));

    assert(element_count == (unsigned long long) compressed_global_bounds[0] * compressed_global_bounds[1] * compressed_global_bounds[2] * compressed_global_bounds[3] * compressed_global_bounds[4] * variable_count / pow(2, resolution));

  }

  destroyBlockBitmap(global_block_layout);
  free(global_block_layout);
  global_block_layout = 0;

  return 0;
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

/////////////////////////////////////////////////
static int generate_file_name_template(int maxh, int bits_per_block, char* filename, int current_time_step, char* filename_template)
{
  int N;
  char dirname[1024], basename[1024];
  int nbits_blocknumber;
  char* directory_path;
  char* data_set_path;

  directory_path = (char*) malloc(sizeof (char) * 1024);
  assert(directory_path);
  memset(directory_path, 0, sizeof (char) * 1024);

  data_set_path = (char*) malloc(sizeof (char) * 1024);
  assert(data_set_path);
  memset(data_set_path, 0, sizeof (char) * 1024);

  strncpy(directory_path, filename, strlen(filename) - 4);
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
      assert(!(nbits_blocknumber % 4));
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
      assert(nbits_blocknumber <= 0);
    }
  }

  free(data_set_path);
  return 0;
}

/*
static void revstr(char* str)
{
  unsigned long long i;
  char cpstr[strlen(str)+1];
  for(i=0; i < strlen(str); i++)
    cpstr[i] = str[strlen(str)-i-1];

  cpstr[i] = '\0';
  strcpy(str, cpstr);
}


static void GuessBitmaskPattern(char* _bits, PointND dims)
{
  int D,N,ordered;
  int dim = 1;
  char* p=_bits;

  PointND id,sorted_id;

  *p++='V';

  For(D)
  {
    PGET(dims,D)=( int)getPowerOf2(PGET(dims,D));
    PGET(id,D)=D;
  }

  //order is ASC order (from smaller dimension to bigger)
  for (ordered=0;!ordered;)
  {
    ordered=1;
    OffsetFor(D,0,-1)
    {
      int ref0=PGET(id,D  ),dim0=PGET(dims,ref0);
      int ref1=PGET(id,D+1),dim1=PGET(dims,ref1);
      if (!(dim0<dim1 || (dim0==dim1 && ref0<ref1)))
      {
        int _temp=PGET(id,D);
        PGET(id,D)=PGET(id,D+1);
        PGET(id,D+1)=_temp;
        ordered=0;
      }
    }
  }

  For(D)
  {
    //order in DESC order
    for (ordered=0,sorted_id=id;!ordered;)
    {
      ordered=1;
      OffsetFor(N,D,-1)
      {
        if (PGET(sorted_id,N+0)<PGET(sorted_id,N+1))
        {
          int _temp=PGET(sorted_id,N);
          PGET(sorted_id,N)=PGET(sorted_id,N+1);
          PGET(sorted_id,N+1)=_temp;
          ordered=0;
        }
      }
    }
    //while dim is not consumed
    for (;dim<PGET(dims,PGET(id,D));dim<<=1)
    {
      OffsetFor(N,D,0)
      {
        *p++='0'+PGET(sorted_id,N);
      }
    }
  }
  *p++=0;
  revstr(_bits+1);
  //strrev(_bits+1)
}
*/

static int is_block_present(int block_number, block_layout layout)
{
  long i, j;

  if( block_number == 0)
    return 1;

  for(i = 1 ; i < layout->levels ; i++)
    for(j = 0 ; j < layout->hz_block_count_array[i] ; j++)
      if ( layout->hz_block_number_array[i][j] == block_number)
        return 1;

  return(0);
}


static void destroyBlockBitmap(block_layout layout)
{
  int j = 0;
  free(layout->hz_block_count_array);
  for(j = 0 ; j < (layout->levels) ; j++)
  {
    free(layout->hz_block_number_array[j]);
    layout->hz_block_number_array[j] = 0;
  }
  free(layout->hz_block_number_array);
  layout->hz_block_number_array = 0;
  layout->levels = 0;
}

static int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, int res, const char* bitPattern, block_layout layout)
{
  unsigned long long hz_from = 0, hz_to = 0, block_number = 1;
  int i, j, m, n_blocks = 1, ctr = 1;
  unsigned long long *ZYX_from, *ZYX_to;

  if(maxH < bits_per_block)
    layout->levels = 1;
  else
    layout->levels = maxH - bits_per_block;

  layout->hz_block_count_array = (int*)malloc(sizeof(int) * (layout->levels));
  layout->hz_block_count_array[0] = 1;
  layout->hz_block_number_array = (int**)malloc(sizeof(int*) * (layout->levels));
  layout->hz_block_number_array[0] = (int*)malloc(sizeof(int) * ctr);
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    layout->hz_block_count_array[j] = 0;
    layout->hz_block_number_array[j] = (int*)malloc(sizeof(int) * ctr);
    ctr = ctr * 2;
  }
  ctr = 1;
  layout->hz_block_number_array[0][0] = 0;
  for(j = 1 ; j < (layout->levels) ; j++)
  {
    for(i = 0 ; i < ctr ; i++)
    {
      layout->hz_block_number_array[j][i] = 0;
    }
    ctr = ctr * 2;
  }
  layout->hz_block_count_array[0] = 1;                  //This block contains data upto level "bits_per_block"

  hz_from = (unsigned long long)(block_number - 1) * pow(2, bits_per_block);
  hz_to = (unsigned long long)(block_number * pow(2, bits_per_block)) - 1;

  for(m = 1 ; m < (maxH - bits_per_block - res); m++)
  {
    n_blocks = pow(2, (m - 1));
    int t = 0;
    for(t = 0 ; t < n_blocks ; t++)
    {
      block_number = block_number + 1;

      hz_from = (unsigned long long)(block_number - 1) * pow(2, bits_per_block);
      hz_to = (unsigned long long)(block_number * pow(2, bits_per_block)) - 1;

      ZYX_to = (unsigned long long*) malloc(sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
      ZYX_from = (unsigned long long*) malloc(sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

      Hz_to_xyz(bitPattern, maxH - 1, hz_from, ZYX_from);
      Hz_to_xyz(bitPattern, maxH - 1, hz_to, ZYX_to);

      if(ZYX_to[0] >= bounding_box[0][0] && ZYX_from[0] < bounding_box[1][0] &&
         ZYX_to[1] >= bounding_box[0][1] && ZYX_from[1] < bounding_box[1][1] &&
         ZYX_to[2] >= bounding_box[0][2] && ZYX_from[2] < bounding_box[1][2] &&
         ZYX_to[3] >= bounding_box[0][3] && ZYX_from[3] < bounding_box[1][3] &&
         ZYX_to[4] >= bounding_box[0][4] && ZYX_from[4] < bounding_box[1][4])
      {
        layout->hz_block_count_array[m] = layout->hz_block_count_array[m] + 1;
        layout->hz_block_number_array[m][t] = block_number - 1;
      }
      free(ZYX_from);
      ZYX_from = 0;
      free(ZYX_to);
      ZYX_to = 0;
    }
  }
  return 0;
}

static int RegExBitmaskBit(const char* bitmask_pattern,int N)
{
  const char *OpenRegEx;
  int S, L;
  assert(bitmask_pattern[0]=='V');

  if (!N)
    return bitmask_pattern[0];

  if ((OpenRegEx=strchr(bitmask_pattern,'{')))
  {
    S = 1+OpenRegEx-bitmask_pattern;
    L = strchr(bitmask_pattern,'}')-bitmask_pattern-S;

    if ((N+1)<S)
      return bitmask_pattern[N]-'0';
    else
      return bitmask_pattern[S+((N+1-S)%L)]-'0';
  }
  return bitmask_pattern[N]-'0';
}

static void Hz_to_xyz(const char* bitmask,  int maxh, unsigned long long hzaddress, unsigned long long* xyz)
{
  unsigned long long lastbitmask=((unsigned long long)1)<<maxh;

  hzaddress <<= 1;
  hzaddress  |= 1;
  while ((lastbitmask & hzaddress) == 0) hzaddress <<= 1;
    hzaddress &= lastbitmask - 1;

  PointND cnt;
  PointND p  ;
  int n = 0;

  memset(&cnt,0,sizeof(PointND));
  memset(&p  ,0,sizeof(PointND));

  for (;hzaddress; hzaddress >>= 1,++n, maxh--)
  {
    int bit= bitmask[maxh];
    PGET(p,bit) |= (hzaddress & 1) << PGET(cnt,bit);
    ++PGET(cnt,bit);
  }
  xyz[0] = p.x;
  xyz[1] = p.y;
  xyz[2] = p.z;
  xyz[3] = p.u;
  xyz[4] = p.v;
}

static int VisusSplitFilename(const char* filename,char* dirname,char* basename)
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

static unsigned long long getPowerOf2(int x)
{
  /*  find the power of 2 of an integer value (example 5->8) */
  int n = 1;
  while (n < x) n <<= 1;
  return n;
}
