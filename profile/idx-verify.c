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

#define PIDX_MAX_DIMENSIONS 5
#define MAX_VARIABLE_COUNT 1024
#define MAX_TEMPLATE_DEPTH 6

#define Min2ab(a,b)      (((a)<=(b))?(a):(b))
#define Max2ab(a,b)      (((a)> (b))?(a):(b))
#define OffsetFor(_D_,_From_,_Off_) for ((_D_)=(_From_);(_D_)<(PIDX_MAX_DIMENSIONS+(_Off_));(_D_)++)
#define For(_D_) for ((_D_)=0;(_D_)<PIDX_MAX_DIMENSIONS;(_D_)++)
#define PGET(_Point_,_Coordinate_) ((&((_Point_).x))[(_Coordinate_)])

typedef struct {int x,y,z,u,v;} PointND;
struct block_layout
{
  int    levels;                          // Total number of Levels
  int    *hz_block_count_array;           // Number of filled blocks per level
  int    ** hz_block_number_array;        // Indices of filled blocks
};
typedef struct block_layout block_layout;

static void revstr(char* str);
static void GuessBitmaskPattern(char* _bits, PointND dims);
static int generate_file_name_template(int maxh, int bits_per_block, char* filename, int current_time_step, char* filename_template);
static int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen);
static int is_block_present(int block_number, block_layout* layout);
static void destroyBlockBitmap(block_layout* layout);
static int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, const char* bitPattern, block_layout* layout);
static int VisusSplitFilename(const char* filename,char* dirname,char* basename);
static void Hz_to_xyz(const char* bitmask,  int maxh, long long hzaddress, long long* xyz);
static int RegExBitmaskBit(const char* bitmask_pattern,int N);
static unsigned long long getPowerOf2(int x);

int main(int argc, char **argv) 
{  
  if ( argv[1] == NULL )
  {
    fprintf(stderr, "Missing arguments (provide a .idx file as command line input argument)\n");
    return 0;
  }
  
  int i, j, k, t, var, ret, counter = 0;
  char *pch, *pch1;
  char line [ 512 ];
  int count = 0, len = 0;
  
  int max_files;
  int maxh;
  char bitSequence[512];
  char bitPattern[512];
  int bits_per_block;
  int samples_per_block;
  int blocks_per_file;
  int start_time_step, end_time_step;
  
  int global_bounds[PIDX_MAX_DIMENSIONS];
  int values_per_sample[MAX_VARIABLE_COUNT];
  char variable_name[MAX_VARIABLE_COUNT][1024];
  char filename_template[1024];
  int variable_count = 0;
  
  FILE *fp = fopen(argv[1], "r");
  while (fgets(line, sizeof (line), fp) != NULL) 
  {
    //printf("%s", line);
    len = strlen(line) - 1;
    if (line[len] == '\n')
      line[len] = 0;
      
    if (strcmp(line, "(version)") == 0)
      fgets(line, sizeof line, fp);
      
    if (strcmp(line, "(box)") == 0) 
    {
      fgets(line, sizeof line, fp);
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
      fgets(line, sizeof (line), fp);
      len = strlen(line) - 1;
      if (line[len] == '\n')
	line[len] = 0;
      count = 0;
      variable_count = 0;
      
      while (strcmp(line, "(logic_to_physic)") != 0 && strcmp(line, "(version)") != 0 && strcmp(line, "(box)") != 0 && strcmp(line, "(bits)") && strcmp(line, "(bitsperblock)") != 0 && strcmp(line, "(blocksperfile)") != 0 && strcmp(line, "(filename_template)") != 0 && strcmp(line, "(time)") != 0)
      {
        pch1 = strtok(line, " *+");
        while (pch1 != NULL)
        {
          //printf("");
          if (count == 0)
            strcpy(variable_name[variable_count], strdup(pch1));

          if (count == 1)
            values_per_sample[variable_count] = atoi(pch1);

          if (count == 2)
          {
            len = strlen(pch1) - 1;
            if (pch1[len] == '\n')
              pch1[len] = 0;
            if (strcmp(pch1, "float64") != 0)
            {
              fprintf(stderr, "Currently supporting only float64 types\n");
              return 0;
            }
          }
          count++;
          pch1 = strtok(NULL, " *+");
        }
        count = 0;

        fgets(line, sizeof (line), fp);
        len = strlen(line) - 1;
        if (line[len] == '\n')
          line[len] = 0;
        variable_count++;
      }
    }
    if (strcmp(line, "(bits)") == 0)
      fgets(line, sizeof line, fp);

    if (strcmp(line, "(bitsperblock)") == 0) 
    {
      fgets(line, sizeof line, fp);
      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      bits_per_block = atoi(line);
      samples_per_block = pow(2, bits_per_block);
    }
    if (strcmp(line, "(blocksperfile)") == 0) 
    {
      fgets(line, sizeof line, fp);
      len = strlen(line) - 1;
      if (line[len] == '\n')
	line[len] = 0;
      blocks_per_file= atoi(line);
    }
    if (strcmp(line, "(filename_template)") == 0) 
    {
      fgets(line, sizeof line, fp);
      len = strlen(line) - 1;
      if (line[len] == '\n')
	line[len] = 0;
    }
    if (strcmp(line, "(time)") == 0)
    {
      fgets(line, sizeof line, fp);
      len = strlen(line) - 1;
      if (line[len] == '\n')
        line[len] = 0;
      
      pch1 = strtok(line, " ");
      count = 0;
      while (pch1 != NULL)
      {
        //printf("");
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
  
  printf("Finished Parsing %s\n", argv[1]);
  printf("Starting time step %d and Ending time step %d\n", start_time_step, end_time_step);
  printf("(box)\n0 %d 0 %d 0 %d 0 %d 0 %d\n", global_bounds[0], global_bounds[1], global_bounds[2], global_bounds[3], global_bounds[4]);
  printf("(fields)\n");
  for(var = 0 ; var < variable_count ; var++)
  {
      printf("%s %d*", variable_name[var], values_per_sample[var]);
      printf("float64 ");
      if(var != variable_count - 1)
	printf(" + \n");
  }
  printf("(bitsperblock)\n%d\n(blocksperfile)\n%d\n", bits_per_block, blocks_per_file);
  printf("(filename_template)\n%s\n", filename_template);
    
  
  int bounding_box[2][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
  };

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++) 
  {
    bounding_box[0][i] = 0;
    bounding_box[1][i] = global_bounds[i];
  }
    
  PointND extents;
  extents.x = global_bounds[0];
  extents.y = global_bounds[1];
  extents.z = global_bounds[2];
  extents.u = global_bounds[3];
  extents.v = global_bounds[4];
  GuessBitmaskPattern(bitSequence, extents);
  maxh = strlen(bitSequence);

  for (i = 0; i <= maxh; i++)
    bitPattern[i] = RegExBitmaskBit(bitSequence, i);

  block_layout* global_block_layout =  malloc(sizeof (global_block_layout));
  createBlockBitmap(bounding_box, blocks_per_file, bits_per_block, maxh, bitPattern, global_block_layout);
  
  k = 1;
  for (i = 1; i < (global_block_layout->levels); i++)
  {
    counter = 0;
    for(j = 0 ; j < k ; j++)
    {
      if(global_block_layout->hz_block_number_array[i][j] != 0)
      {
        global_block_layout->hz_block_number_array[i][counter] = global_block_layout->hz_block_number_array[i][j];
        counter++;
      }
    }
    k = k * 2;
  }
  
  /// maximum number of files possible
  max_files = (getPowerOf2(global_bounds[0]) * getPowerOf2(global_bounds[1]) * getPowerOf2(global_bounds[2]) * getPowerOf2(global_bounds[3]) * getPowerOf2(global_bounds[4])) / ((long long) pow(2, bits_per_block) * (long long) blocks_per_file);
  if ((getPowerOf2(global_bounds[0]) * getPowerOf2(global_bounds[1]) * getPowerOf2(global_bounds[2]) * getPowerOf2(global_bounds[3]) * getPowerOf2(global_bounds[4])) % ((long long) pow(2, bits_per_block) * (long long) blocks_per_file))
    max_files++;
  assert(max_files != 0);

  int file_no;
  int *existing_file_index;
  existing_file_index = (int*) malloc(max_files * sizeof (int));
  memset(existing_file_index, 0, max_files * sizeof (int));
  existing_file_index[0] = 1;
  
  for (i = 1; i < global_block_layout->levels; i++) 
  {
    for (j = 0; j < global_block_layout->hz_block_count_array[i]; j++) 
    {
      file_no = global_block_layout->hz_block_number_array[i][j] / blocks_per_file;
      existing_file_index[file_no] = 1;
    }
  }
    
  uint32_t* binheader;
  int binheader_count;
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader = (uint32_t*) malloc(sizeof (*binheader)*(binheader_count));
    
  int fd;
  long long ZYX[PIDX_MAX_DIMENSIONS];
  
  for (t = start_time_step; t <= end_time_step; t++)
  {
    long long lost_element_count = 0, element_count = 0;
    for (i = 0; i < max_files; i++) 
    {
      if (existing_file_index[i] == 1) 
      {
        generate_file_name_template(maxh, bits_per_block, argv[1], t, filename_template);
       
        char bin_file[PATH_MAX];
        ret = generate_file_name(blocks_per_file, filename_template, i, bin_file, PATH_MAX);
        if (ret == -1)
        {
          fprintf(stderr, "[File : %s] [Line : %d] generate_file_name\n", __FILE__, __LINE__);
          return 0;
        }
        
        printf("Parsing File : %s\n", bin_file);
        
        fd = open(bin_file, O_RDONLY);
        if (fd < 0) 
        {
          fprintf(stderr, "[File : %s] [Line : %d] open\n", __FILE__, __LINE__);
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
        unsigned long long* data_buffer = NULL;
        int check_bit = 1, s = 0;
        long long hz_index, hz_val;
        
        size_t data_size;
        off_t data_offset;
        for (var = 0; var < variable_count; var++) 
        {
          for (bpf = 0; bpf < blocks_per_file; bpf++) 
          {
            if (is_block_present((bpf + (blocks_per_file * i)), global_block_layout)) 
            {
              data_offset = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 12]);
              data_size = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 14]);

              //if(var == 2 || var == 1)
              //printf("[%d] Offset %ld Count %ld\n", bpf, data_offset, data_size);
              data_buffer = malloc(data_size);
              memset(data_buffer, 0, data_size);

              ret = pread(fd, data_buffer, data_size, data_offset);
              //printf("[%d] %ld and %ld\n", bpf, data_size, ret);
              assert(ret == data_size);
              
              for (hz_val = 0; hz_val < data_size/(sizeof(unsigned long long) * values_per_sample[var]); hz_val++) 
              {
                hz_index = (blocks_per_file * i * samples_per_block) + (bpf * samples_per_block) + hz_val;
                Hz_to_xyz(bitPattern, maxh - 1, hz_index, ZYX);
                
                //if (data_buffer[hz_val * values_per_sample[var] + 0] == 0 && (ZYX[2] < global_bounds[2]))
                if (ZYX[2] >= global_bounds[2] || ZYX[1] >= global_bounds[1] || ZYX[0] >= global_bounds[0])
                  continue;

                check_bit = 1, s = 0;
                for (s = 0; s < values_per_sample[var]; s++)
                  check_bit = check_bit && (data_buffer[hz_val * values_per_sample[var] + s] == 100 + (global_bounds[0] * global_bounds[1] * ZYX[2])+(global_bounds[0]*(ZYX[1])) + ZYX[0]);

                if (check_bit == 0)
                {
                  lost_element_count++;
                  printf("[%d] [%lld %lld %lld] [%lld : %lld %lld %lld] Actual: %lld Should Be %lld\n", 
                         var,
                         lost_element_count, hz_index/samples_per_block, hz_val, 
                         hz_index, ZYX[0], ZYX[1], ZYX[2], 
                         (unsigned long long)data_buffer[hz_val * values_per_sample[var] + 0], (100 + (global_bounds[0] * global_bounds[1] * ZYX[2])+(global_bounds[0]*(ZYX[1])) + ZYX[0]));
                }
                else 
                {
                  element_count++;
                }
              }
              
              free(data_buffer);
              data_buffer = 0;
            }
          }
        }
      }
    }
    
    int total_samples = 0;
    for (var = 0; var < variable_count; var++)
      total_samples = total_samples + values_per_sample[var];

    printf("%lld + %lld (%lld) : %lld\n", (long long) (element_count), (long long)lost_element_count, element_count + lost_element_count, (long long) global_bounds[0] * global_bounds[1] * global_bounds[2] * global_bounds[3] * global_bounds[4] * variable_count);
    assert(element_count == (long long) global_bounds[0] * global_bounds[1] * global_bounds[2] * global_bounds[3] * global_bounds[4] * variable_count);
    
  }
  
  destroyBlockBitmap(global_block_layout);
  free(global_block_layout);
  global_block_layout = 0;
    
  
  
  return 0;
}

static int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen) 
{
  long long address = 0;
  unsigned int segs[MAX_TEMPLATE_DEPTH] = {0};
  int seg_count = 0;
  char* pos;
  int ret;

  //printf("[generate_file_name]: %d %s %d :: %s\n", file_number, filename, maxlen, filename_template);
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
  sprintf(data_set_path, "%s/time%06d.idx", directory_path, current_time_step);
  
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

static void revstr(char* str)
{
  long long i;
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

static int is_block_present(int block_number, block_layout* layout)
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


static void destroyBlockBitmap(block_layout* layout)
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

static int createBlockBitmap(int bounding_box[2][5], int blocks_per_file, int bits_per_block, int maxH, const char* bitPattern, block_layout* layout)
{
  long long hz_from = 0, hz_to = 0, block_number = 1;
  int i, j, m, n_blocks = 1, ctr = 1;
  long long *ZYX_from, *ZYX_to;
  
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
  
  hz_from = (long long)(block_number - 1) * pow(2, bits_per_block);
  hz_to = (long long)(block_number * pow(2, bits_per_block)) - 1;
  
  for(m = 1 ; m < (maxH - bits_per_block); m++)
  {
    n_blocks = pow(2, (m - 1));
    int t = 0;
    for(t = 0 ; t < n_blocks ; t++)
    {
      block_number = block_number + 1;
      
      hz_from = (long long)(block_number - 1) * pow(2, bits_per_block);
      hz_to = (long long)(block_number * pow(2, bits_per_block)) - 1;
      
      ZYX_to = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
      ZYX_from = malloc(sizeof(long long) * PIDX_MAX_DIMENSIONS);
      
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


static void Hz_to_xyz(const char* bitmask,  int maxh, long long hzaddress, long long* xyz)
{
  long long lastbitmask=((long long)1)<<maxh;
  
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