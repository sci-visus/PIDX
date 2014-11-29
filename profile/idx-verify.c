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

#include "PIDX_inc.h"

#define MAX_VARIABLE_COUNT 1024
#define MAX_TEMPLATE_DEPTH 6

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

  block_layout* block_layout =  malloc(sizeof (block_layout));
  createBlockBitmap(bounding_box, blocks_per_file, bits_per_block, maxh, bitPattern, block_layout);
  
  k = 1;
  for (i = 1; i < (block_layout->levels); i++)
  {
    counter = 0;
    for(j = 0 ; j < k ; j++)
    {
      if(block_layout->hz_block_number_array[i][j] != 0)
      {
        block_layout->hz_block_number_array[i][counter] = block_layout->hz_block_number_array[i][j];
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
  
  for (i = 1; i < block_layout->levels; i++) 
  {
    for (j = 0; j < block_layout->hz_block_count_array[i]; j++) 
    {
      file_no = block_layout->hz_block_number_array[i][j] / blocks_per_file;
      existing_file_index[file_no] = 1;
    }
  }
    
  uint32_t* binheader;
  int binheader_count;
  binheader_count = 10 + 10 * blocks_per_file * variable_count;
  binheader = (uint32_t*) malloc(sizeof (*binheader)*(binheader_count));
    
  int fd;
  long long ZYX[PIDX_MAX_DIMENSIONS], lost_element_count = 0, element_count = 0;
  
  for (t = start_time_step; t <= end_time_step; t++)
  {
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
        double* data_buffer = NULL;
        int check_bit = 1, s = 0;
        long long hz_index, hz_val;
        
        size_t data_size;
        off_t data_offset;
        for (var = 0; var < variable_count; var++) 
        {
          for (bpf = 0; bpf < blocks_per_file; bpf++) 
          {
            if (is_block_present((bpf + (blocks_per_file * i)), block_layout)) 
            {
              data_offset = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 12]);
              data_size = ntohl(binheader[(bpf + var * blocks_per_file)*10 + 14]);

              data_buffer = malloc(data_size);
              memset(data_buffer, 0, data_size);

              ret = pread(fd, data_buffer, data_size, data_offset);
              assert(ret == data_size);
              
              for (hz_val = 0; hz_val < data_size/(sizeof(double) * values_per_sample[var]); hz_val++) 
              {
                hz_index = (blocks_per_file * i * samples_per_block) + (bpf * samples_per_block) + hz_val;
                Hz_to_xyz(bitPattern, maxh - 1, hz_index, ZYX);
                
                if (data_buffer[hz_val * values_per_sample[var] + 0] == 0)
                  continue;

                check_bit = 1, s = 0;
                for (s = 0; s < values_per_sample[var]; s++)
                  check_bit = check_bit && (data_buffer[hz_val * values_per_sample[var] + s] == s + var + 100 + (global_bounds[0] * global_bounds[1] * ZYX[2])+(global_bounds[0]*(ZYX[1])) + ZYX[0]);

                if (check_bit == 0)
                {
                  lost_element_count++;
                  printf("[%d %d] [%lld : %lld %lld %lld] Actual: %d Should Be %d\n", lost_element_count, hz_index/samples_per_block, hz_index, ZYX[0], ZYX[1], ZYX[2], (int)data_buffer[hz_val * values_per_sample[var] + 0], (var + 100 + (global_bounds[0] * global_bounds[1] * ZYX[2])+(global_bounds[0]*(ZYX[1])) + ZYX[0]));
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
  }
    
  int total_samples = 0;
  for (var = 0; var < variable_count; var++)
    total_samples = total_samples + values_per_sample[var];

  printf("%lld + %lld (%lld) : %lld\n", (long long) (element_count), (long long)lost_element_count, element_count + lost_element_count, (long long) global_bounds[0] * global_bounds[1] * global_bounds[2] * global_bounds[3] * global_bounds[4] * variable_count);
  assert(element_count == (long long) global_bounds[0] * global_bounds[1] * global_bounds[2] * global_bounds[3] * global_bounds[4] * variable_count);
  
}