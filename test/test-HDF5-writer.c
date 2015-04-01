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


#include "pidxtest.h"

#include <PIDX.h>

#if PIDX_OPTION_HDF5
#include "hdf5.h"
#endif

int test_hdf5_writer(struct Args args, int rank, int nprocs)
{
#if PIDX_OPTION_HDF5
#if PIDX_HAVE_MPI
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
  
  /// The command line arguments are shared by all processes
  MPI_Bcast(args.extents, 5, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(args.count_local, 5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.time_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.variable_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.idx_count, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&args.output_file_template, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  /// Calculating every process's offset and count  
  sub_div[0] = (args.extents[0] / args.count_local[0]);
  sub_div[1] = (args.extents[1] / args.count_local[1]);
  sub_div[2] = (args.extents[2] / args.count_local[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * args.count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * args.count_local[1];
  local_offset[0] = (slice % sub_div[0]) * args.count_local[0];

  double_data = (double**)malloc(sizeof(*double_data) * args.variable_count);
  memset(double_data, 0, sizeof(*double_data) * args.variable_count);
  
  for(var = 0; var < args.variable_count; var++)
  {
    double_data[var] = (double*)malloc(sizeof (uint64_t) * args.count_local[0] * args.count_local[1] * args.count_local[2]);
    
    //printf("Local Extents %d %d %d\n", args.count_local[0], args.count_local[1], args.count_local[2]);
    for (k = 0; k < args.count_local[2]; k++)
      for (j = 0; j < args.count_local[1]; j++)
        for (i = 0; i < args.count_local[0]; i++)
        {
          int64_t index = (int64_t) (args.count_local[0] * args.count_local[1] * k) + (args.count_local[0] * j) + i;
          double_data[var][index] = ((args.extents[0] * args.extents[1]*(local_offset[2] + k))+(args.extents[0]*(local_offset[1] + j)) + (local_offset[0] + i));
          //if (rank == 0)
          //printf("[%d] Samples are %f \n", rank, double_data[var][index]);  
        }
  }
  
  double sim_start = 0, sim_end = 0, total_time = 0, max_time = 0, global_max_time = 0;
  
  //printf ("%d %d %d\n", args.idx_count[0], args.idx_count[1], args.idx_count[2]);
  if (args.idx_count[0] != 1 || args.idx_count[1] != 1 || args.idx_count[2] != 1 )
  {
    args.extents[0] = args.extents[0] / args.idx_count[0];
    args.extents[1] = args.extents[1] / args.idx_count[1];
    args.extents[2] = args.extents[2] / args.idx_count[2];
  
    for (j = 0; j < args.extents[0] * args.idx_count[0]; j = j + (args.extents[0]))
      if (local_offset[0] >= j && local_offset[0] < (j + args.extents[0]))
      {
        local_offset[0] = local_offset[0] - j;
        break;
      }
        
    for (j = 0; j <= args.extents[1] * args.idx_count[1]; j = j + (args.extents[1]))
      if (local_offset[1] >= j && local_offset[1] < j + (args.extents[1] ))
      {
        local_offset[1] = local_offset[1] - j;
        break;
      }
        
    for (j = 0; j < args.extents[2] * args.idx_count[2]; j = j + (args.extents[2]))
      if (local_offset[2] >= j && local_offset[2] < j + (args.extents[2] /* args.idx_count */))
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
    
    int* colors = malloc(sizeof(*colors) * args.idx_count[0] * args.idx_count[1] * args.idx_count[2]);
    memset(colors, 0, sizeof(*colors) * args.idx_count[0] * args.idx_count[1] * args.idx_count[2]);
    
    for (k = 0; k < args.idx_count[2]; k++)
      for (j = 0; j < args.idx_count[1]; j++)
        for (i = 0; i < args.idx_count[0]; i++)
          colors[(args.idx_count[0] * args.idx_count[1] * k) + (args.idx_count[0] * j) + i] = (args.idx_count[0] * args.idx_count[1] * k) + (args.idx_count[0] * j) + i;
          
    int index_x = 0, index_y = 0, index_z = 0;
    for (i = 0; i < sub_div[0]; i = i + (sub_div[0] / args.idx_count[0]))
    {
      if (rank_x >= i && rank_x < i + (sub_div[0] / args.idx_count[0]))
      {
        index_x = i;
        break;
      }
    }
    for (i = 0; i < sub_div[1]; i = i + (sub_div[1] / args.idx_count[1]))
    {
      if (rank_y >= i && rank_y < i + (sub_div[1] / args.idx_count[1]))
      {
        index_y = i;
        break;
      }
    }
    for (i = 0; i < sub_div[2]; i = i + (sub_div[2] / args.idx_count[2]))
    {
      if (rank_z >= i && rank_z < i + (sub_div[2] / args.idx_count[2]))
      {
        index_z = i;
        break;
      }
    }
    
    //(args.idx_count[0] * args.idx_count[1] * (index_z/(sub_div[2] / args.idx_count[2]))) + (args.idx_count[0] * (index_y/ (sub_div[1] / args.idx_count[1]))) + (index_x / (sub_div[0] / args.idx_count[0]));
    color = colors[(args.idx_count[0] * args.idx_count[1] * (index_z/(sub_div[2] / args.idx_count[2]))) + (args.idx_count[0] * (index_y/ (sub_div[1] / args.idx_count[1]))) + (index_x / (sub_div[0] / args.idx_count[0]))];
    
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
  dimsf[0] = args.extents[2];
  dimsf[1] = args.extents[1];
  dimsf[2] = args.extents[0];
  
  hsize_t count[3];
  count[0] = args.count_local[2];
  count[1] = args.count_local[1];
  count[2] = args.count_local[0];
  
  hsize_t offset[3];
  offset[0] = local_offset[2];
  offset[1] = local_offset[1];
  offset[2] = local_offset[0];
  
  MPI_Info info  = MPI_INFO_NULL;
  
  for (ts = 0; ts < args.time_step; ts++) 
  {
    /// Creating the filename 
    args.output_file_name = (char*) malloc(sizeof (char) * 512);
    sprintf(args.output_file_name, "%s_%d_%d%s", args.output_file_template, color, ts, ".h5");
    output_file = args.output_file_name;
    
    sim_start = MPI_Wtime();
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, local_comm, info);
    file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
  
    char variable_name[512];
    for(var = 0; var < args.variable_count; var++)
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
      printf("[Local] [%d] [%d] %s Variables %d IDXs %lld %lld %lld Global Dimension %lld %lld %lld Time %f seconds Throughput %f MB/sec\n", rank, nprocs, output_file, args.variable_count, (long long)args.idx_count[0], (long long)args.idx_count[1], (long long)args.idx_count[2], (long long)args.extents[0], (long long)args.extents[1], (long long)args.extents[2], max_time, (double)(args.extents[0] * args.extents[1] * args.extents[2] * args.variable_count * sizeof(double))/(max_time*1000*1000));
    }
    
    //if (args.idx_count[0] != 1 || args.idx_count[1] != 1 || args.idx_count[2] != 1 )
    //{
      int global_rank, global_nprocs;
      MPI_Allreduce(&max_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, global_comm);
      MPI_Comm_rank(global_comm, &global_rank);
      MPI_Comm_size(global_comm, &global_nprocs);
      
      if (total_time == global_max_time)
      {
        printf("[Global][%d] [%d] %s Variables %d IDXs %lld %lld %lld Global Dimension %lld %lld %lld Time %f seconds Throughput %f MB/sec\n", rank, nprocs, output_file, args.variable_count, (long long)args.idx_count[0], (long long)args.idx_count[1], (long long)args.idx_count[2], (long long)args.extents[0], (long long)args.extents[1], (long long)args.extents[2], global_max_time, (double)(args.extents[0] * args.extents[1] * args.extents[2] * args.variable_count * sizeof(double))/(global_max_time*1000*1000));
        printf("=============================================================================================================================\n\n");
      }
    
    //if (rank == 0)
      
  }
  
  for(var = 0; var < args.variable_count; var++)
  {
    free(double_data[var]);
    double_data[var] = 0;
  }
  free(double_data);
  double_data = 0;
  
  free(args.output_file_name);
  
#endif
#endif
  return 0;
}