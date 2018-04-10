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

#include <byteswap.h>
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
#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>

#include "hdf5.h"

static int process_count = 1, rank = 0;
static int partition_count[3];
static int global_box_size[3];
static int local_box_size[3];
static int time_step_count = 1;
static int variable_count = 1;
static char output_file_template[512];
static char output_file_name[512];

static void init_mpi(int argc, char **argv);
static void parse_args(int argc, char **argv);
static void terminate_with_error_msg(const char *format, ...);
static void terminate();


int main(int argc, char **argv)
{

  init_mpi(argc, argv);
  parse_args(argc, argv);

  int color = 0;
  int i = 0, j = 0, k = 0;
  int ts, var;
  int slice;
  int sub_div[3], local_offset[3];

  /// IDX File Name
  const char *output_file; 
  
  double **double_data;
  
  MPI_Comm global_comm, local_comm;
  global_comm = MPI_COMM_WORLD;
  
  /// Calculating every process's offset and count  
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * local_box_size[1];
  local_offset[0] = (slice % sub_div[0]) * local_box_size[0];

  double_data = (double**)malloc(sizeof(*double_data) * variable_count);
  memset(double_data, 0, sizeof(*double_data) * variable_count);
  
  for(var = 0; var < variable_count; var++)
  {
    double_data[var] = (double*)malloc(sizeof (uint64_t) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
    
    //printf("Local Extents %d %d %d\n", local_box_size[0], local_box_size[1], local_box_size[2]);
    for (k = 0; k < local_box_size[2]; k++)
      for (j = 0; j < local_box_size[1]; j++)
        for (i = 0; i < local_box_size[0]; i++)
        {
          int64_t index = (int64_t) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
          double_data[var][index] = ((global_box_size[0] * global_box_size[1]*(local_offset[2] + k))+(global_box_size[0]*(local_offset[1] + j)) + (local_offset[0] + i));
        }
  }
  
  double sim_start = 0, sim_end = 0, total_time = 0, max_time = 0, global_max_time = 0;

  printf ("%d %d %d\n", partition_count[0], partition_count[1], partition_count[2]);
  if (partition_count[0] != 1 || partition_count[1] != 1 || partition_count[2] != 1 )
  {
    global_box_size[0] = global_box_size[0] / partition_count[0];
    global_box_size[1] = global_box_size[1] / partition_count[1];
    global_box_size[2] = global_box_size[2] / partition_count[2];
  
    for (j = 0; j < global_box_size[0] * partition_count[0]; j = j + (global_box_size[0]))
      if (local_offset[0] >= j && local_offset[0] < (j + global_box_size[0]))
      {
        local_offset[0] = local_offset[0] - j;
        break;
      }
        
    for (j = 0; j <= global_box_size[1] * partition_count[1]; j = j + (global_box_size[1]))
      if (local_offset[1] >= j && local_offset[1] < j + (global_box_size[1] ))
      {
        local_offset[1] = local_offset[1] - j;
        break;
      }
        
    for (j = 0; j < global_box_size[2] * partition_count[2]; j = j + (global_box_size[2]))
      if (local_offset[2] >= j && local_offset[2] < j + (global_box_size[2] /* partition_count */))
      {
        local_offset[2] = local_offset[2] - j;
        break;
      }
        
    int rank_x, rank_y, rank_z, rank_slice;
    int i = 0, j = 0, k = 0;
    rank_z = rank / (sub_div[0] * sub_div[1]);
    rank_slice = rank % (sub_div[0] * sub_div[1]);
    rank_y = (rank_slice / sub_div[0]);
    rank_x = (rank_slice % sub_div[0]);
    

    int* colors = (int*)malloc(sizeof(*colors) * partition_count[0] * partition_count[1] * partition_count[2]);
    memset(colors, 0, sizeof(*colors) * partition_count[0] * partition_count[1] * partition_count[2]);
    
    for (k = 0; k < partition_count[2]; k++)
      for (j = 0; j < partition_count[1]; j++)
        for (i = 0; i < partition_count[0]; i++)
          colors[(partition_count[0] * partition_count[1] * k) + (partition_count[0] * j) + i] = (partition_count[0] * partition_count[1] * k) + (partition_count[0] * j) + i;
          
    int index_x = 0, index_y = 0, index_z = 0;
    for (i = 0; i < sub_div[0]; i = i + (sub_div[0] / partition_count[0]))
    {
      if (rank_x >= i && rank_x < i + (sub_div[0] / partition_count[0]))
      {
        index_x = i;
        break;
      }
    }
    for (i = 0; i < sub_div[1]; i = i + (sub_div[1] / partition_count[1]))
    {
      if (rank_y >= i && rank_y < i + (sub_div[1] / partition_count[1]))
      {
        index_y = i;
        break;
      }
    }
    for (i = 0; i < sub_div[2]; i = i + (sub_div[2] / partition_count[2]))
    {
      if (rank_z >= i && rank_z < i + (sub_div[2] / partition_count[2]))
      {
        index_z = i;
        break;
      }
    }
    
    //(partition_count[0] * partition_count[1] * (index_z/(sub_div[2] / partition_count[2]))) + (partition_count[0] * (index_y/ (sub_div[1] / partition_count[1]))) + (index_x / (sub_div[0] / partition_count[0]));
    color = colors[(partition_count[0] * partition_count[1] * (index_z/(sub_div[2] / partition_count[2]))) + (partition_count[0] * (index_y/ (sub_div[1] / partition_count[1]))) + (index_x / (sub_div[0] / partition_count[0]))];
    
    free(colors);
    
    //printf("%d = %d %d %d :: %d\n", rank, rank_x, rank_y, rank_z, args.idx_derived_ptr->color);
    MPI_Comm_split(global_comm, color, rank, &local_comm);
  }
  else
    MPI_Comm_dup(global_comm, &local_comm);


  
  hid_t file_id;
  hid_t plist_id;     
  hid_t dataset_id;
  
  hid_t file_dataspace;
  hid_t mem_dataspace;
  hid_t filespace;
  
  hsize_t dimsf[3];
  dimsf[0] = global_box_size[2];
  dimsf[1] = global_box_size[1];
  dimsf[2] = global_box_size[0];
  
  hsize_t count[3];
  count[0] = local_box_size[2];
  count[1] = local_box_size[1];
  count[2] = local_box_size[0];
  
  hsize_t offset[3];
  offset[0] = local_offset[2];
  offset[1] = local_offset[1];
  offset[2] = local_offset[0];
  
  MPI_Info info  = MPI_INFO_NULL;
#if 1
  for (ts = 0; ts < time_step_count; ts++)
  {
    /// Creating the filename 
    sprintf(output_file_name, "%s_%d_%d%s", output_file_template, color, ts, ".h5");
    output_file = output_file_name;
    
    sim_start = MPI_Wtime();
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, local_comm, info);
    file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
  
    char variable_name[512];
    for(var = 0; var < variable_count; var++)
    {
      sprintf(variable_name, "variable_%d", var);
      
      filespace = H5Screate_simple(3, dimsf, NULL); 
      dataset_id = H5Dcreate(file_id, variable_name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(filespace);
            
      mem_dataspace = H5Screate_simple (3, count, NULL);
      
      file_dataspace = H5Dget_space (dataset_id);
      H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, plist_id, double_data[var]);
      
      H5Dclose(dataset_id);
      H5Sclose(file_dataspace);
      H5Sclose(mem_dataspace);
      H5Pclose(plist_id);
    }
    H5Fclose(file_id);
    
    sim_end = MPI_Wtime();
    total_time = sim_end - sim_start;
    MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, local_comm);
    
    if (total_time == max_time)
    {
      printf("[Local] [%d] [%d] %s Variables %d IDXs %d %d %d Global Dimension %d %d %d Time %f seconds Throughput %f MB/sec\n", rank, process_count, output_file, variable_count, partition_count[0], partition_count[1], partition_count[2], global_box_size[0], global_box_size[1], global_box_size[2], max_time, (double)(global_box_size[0] * global_box_size[1] * global_box_size[2] * variable_count * sizeof(double) * partition_count[0] * partition_count[1] * partition_count[2])/(max_time*1000*1000));
    }
    
    //if (partition_count[0] != 1 || partition_count[1] != 1 || partition_count[2] != 1 )
    //{
      int global_rank, global_process_count;
      MPI_Allreduce(&max_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, global_comm);
      MPI_Comm_rank(global_comm, &global_rank);
      MPI_Comm_size(global_comm, &global_process_count);
      
      if (total_time == global_max_time)
      {
        printf("[Global][%d] [%d] %s Variables %d IDXs %d %d %d Global Dimension %d %d %d Time %f seconds Throughput %f MB/sec\n", rank, process_count, output_file, variable_count, partition_count[0], partition_count[1], partition_count[2], global_box_size[0], global_box_size[1], global_box_size[2], global_max_time, (double)(global_box_size[0] * global_box_size[1] * global_box_size[2] * variable_count * sizeof(double) * partition_count[0] * partition_count[1] * partition_count[2])/(global_max_time*1000*1000));
        printf("=============================================================================================================================\n\n");
      }
    //if (rank == 0)
  }

#endif
  for(var = 0; var < variable_count; var++)
  {
    free(double_data[var]);
    double_data[var] = 0;
  }
  free(double_data);
  double_data = 0;


  
  MPI_Finalize();

  return 0;
}


//----------------------------------------------------------------
static void init_mpi(int argc, char **argv)
{
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
}

//----------------------------------------------------------------
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:r:f:t:v:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
    case('g'): // global dimension
      if ((sscanf(optarg, "%dx%dx%d", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) || (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
        terminate_with_error_msg("Invalid global dimensions\n%s", "wrong usage");
      break;

    case('l'): // local dimension
      if ((sscanf(optarg, "%dx%dx%d", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||(local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
        terminate_with_error_msg("Invalid local dimension\n%s", "wrong usage");
      break;

    case('r'): // local dimension
      if ((sscanf(optarg, "%dx%dx%d", &partition_count[0], &partition_count[1], &partition_count[2]) == EOF) ||(partition_count[0] < 1 || partition_count[1] < 1 || partition_count[2] < 1))
        terminate_with_error_msg("Invalid partition dimension\n%s", "wrong usage");
      break;

    case('f'): // output file name
      if (sprintf(output_file_template, "%s", optarg) < 0)
        terminate_with_error_msg("Invalid output file name template\n%s", "wrong usage");
      sprintf(output_file_name, "%s%s", output_file_template, ".idx");
      break;

    case('t'): // number of timesteps
      if (sscanf(optarg, "%d", &time_step_count) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", "wrong usage");
      break;

    case('v'): // number of variables
        if (sscanf(optarg, "%d", &variable_count) < 0)
          terminate_with_error_msg("Invalid variable file\n%s", "wrong usage");
        break;

    default:
      terminate_with_error_msg("Wrong arguments\n%s", "wrong usage");
    }
  }
}


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
