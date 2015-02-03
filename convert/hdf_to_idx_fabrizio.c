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

#include "hdf_to_idx.h"

#if PIDX_OPTION_HDF5

#define HDF_IO 1
#include <PIDX.h>
#include "hdf5.h"


static char *output_file_name;
static double* buffer;

int main(int argc, char **argv) 
{
  int sub_div[3], local_offset[3], count_local[3];
  int t = 0, time_count = 0;
  
#if HDF_IO
  hid_t file_id;
  hid_t plist_id;     
  hid_t group_id;
  hid_t dataset_id;
  hid_t file_dataspace;
  hid_t mem_dataspace;
#endif
  
  PIDX_file file;
  PIDX_access access;
  const int bits_per_block = 15;
  const int blocks_per_file = 256;
  PIDX_variable variable;
      
  int nprocs=1, rank=0, slice;
  
  MPI_Comm comm  = MPI_COMM_WORLD;
  
#if HDF_IO
  MPI_Info info  = MPI_INFO_NULL;
#endif
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank); 
  
  output_file_name = (char*) malloc(sizeof (char) * 1024);
  sprintf(output_file_name, "%s%s", "/scratch/project/visus/datasets/flame", ".idx");
  
  if (nprocs == 320)
  {
    count_local[0] = 256;
    count_local[1] = 188;
    count_local[2] = 192;
  }
  else if (nprocs == 640)
  {
    count_local[0] = 256;
    count_local[1] = 188;
    count_local[2] = 96;
  }
  else if (nprocs == 1280)
  {
    count_local[0] = 256;
    count_local[1] = 94;
    count_local[2] = 96;
  }
  else if (nprocs == 2560)
  {
    count_local[0] = 128;
    count_local[1] = 94;
    count_local[2] = 96;
  }

  sub_div[0] = (1024 / count_local[0]);
  sub_div[1] = (940 / count_local[1]);
  sub_div[2] = (3072 / count_local[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * count_local[1];
  local_offset[0] = (slice % sub_div[0]) * count_local[0];
  
  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  PIDX_set_point_5D(1024, 940, 3072, 1, 1, global_bounding_box);
  PIDX_set_point_5D(local_offset[0], local_offset[1], local_offset[2], 0, 0, local_offset_point);
  PIDX_set_point_5D(count_local[0], count_local[1], count_local[2], 1, 1, local_box_count_point);
  
#if HDF_IO
  hsize_t count[3];
  count[0] = count_local[0];
  count[1] = count_local[1];
  count[2] = count_local[2];
  
  hsize_t offset[3];
  offset[0] = local_offset[0];
  offset[1] = local_offset[1];
  offset[2] = local_offset[2];
#endif
  
  // P
  // U
  // V
  // W
  // ZMIX
  
#if HDF_IO
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
#endif
 
  FILE *fp;
  fp = fopen("list", "r");
  char file_name[1][1024];
  while (!feof(fp)) 
  {
    if (fscanf(fp, "%s", file_name[time_count]) != 1)
      break;
    if(rank == 0)
      printf("%s\n", file_name[time_count]);
    time_count++;
  }
  fclose(fp);
  
  PIDX_create_access(&access);
  PIDX_set_mpi_access(access, 1, 1, 1 comm);
  PIDX_set_process_extent(access, sub_div[0], sub_div[1], sub_div[2]);
  
  if (rank == 0)
    printf("Number of timesteps = %d\n", time_count);
  for (t = 0; t < time_count; t++)
  {
#if HDF_IO
    file_id = H5Fopen(file_name[t], H5F_ACC_RDONLY, plist_id);
#endif
     
    PIDX_file_create(output_file_name, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, t);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_aggregation_factor(file, 1);
    PIDX_set_block_count(file, blocks_per_file);
    PIDX_set_variable_count(file, 1);
    
    PIDX_debug_rst(file, 0);
    PIDX_debug_hz(file, 0);
    PIDX_dump_agg_info(file, 0);
    
    PIDX_enable_hz(file, 1;
    PIDX_enable_agg(file, 1);
    PIDX_enable_io(file, 1);
    
    
#if HDF_IO
    group_id = H5Gopen(file_id, "/data", H5P_DEFAULT);
    dataset_id = H5Dopen2(group_id, "P", H5P_DEFAULT);
    
    mem_dataspace = H5Screate_simple (3, count, NULL);
    file_dataspace = H5Dget_space (dataset_id);
    H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
#endif
    
    buffer = malloc(sizeof(double) * 1024 * (940 * 3072/nprocs));
    memset(buffer, 0, sizeof(double) * 1024 * (940 * 3072/nprocs));
    
#if HDF_IO
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer);
    H5Sclose(mem_dataspace);
    H5Sclose(file_dataspace);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
#else
    
#endif
    
    PIDX_variable_create(file, "P", sizeof(double) * 8, "1*float64", &variable);
    PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
    
    if (rank == 0)
      printf("Writing for time %d\n", t);
    
    PIDX_close(file);
    
    free(buffer);
    buffer = 0;
  }
  
  //////////
#if HDF_IO
  H5Pclose(plist_id);
#endif
  
  PIDX_close_access(access);
  
  free(output_file_name);
  output_file_name = 0;
  
  MPI_Finalize();
  
  return 0;
}
#endif //PIDX_OPTION_HDF5