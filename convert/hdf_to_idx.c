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

#include <PIDX.h>
#include <Generic_data_structs.h>

#include "hdf5.h"
#include "hdf_to_idx.h"

static char *output_file_name;
static double* buffer;

int main(int argc, char **argv) 
{
  int sub_div[3], local_offset[3], count_local[3];
  hid_t file_id;
  hid_t plist_id;     
  hid_t group_id;
  hid_t dataset_id;
  
  hid_t file_dataspace;
  hid_t mem_dataspace;
  
  PIDX_file file;
  PIDX_access access;
  const char *output_file;
  const int bits_per_block = 15;
  const int blocks_per_file = 256;
  PIDX_variable variable;
  PIDX_variable variable2;
  PIDX_variable variable3;
  
  int time = 0;
    
  int nprocs, rank, slice;
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank); 
  
  output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(output_file_name, "%s%s", "ALL", ".idx");
  output_file = output_file_name;
  
  if (nprocs == 4)
  {
    count_local[0] = 128;
    count_local[1] = 256;
    count_local[2] = 256;
  }
  else if (nprocs == 2)
  {
    count_local[0] = 256;
    count_local[1] = 256;
    count_local[2] = 256;
  }
  else
  {
    count_local[0] = 512;
    count_local[1] = 256;
    count_local[2] = 256;
  }
  
  sub_div[0] = (512 / count_local[0]);
  sub_div[1] = (256 / count_local[1]);
  sub_div[2] = (256 / count_local[2]);
  local_offset[2] = (rank / (sub_div[0] * sub_div[1])) * count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_offset[1] = (slice / sub_div[0]) * count_local[1];
  local_offset[0] = (slice % sub_div[0]) * count_local[0];
  
  PIDX_point global_bounding_box, local_offset_point, local_box_count_point;
  PIDX_set_point_5D(512, 256, 256, 1, 1, global_bounding_box);
  PIDX_set_point_5D(local_offset[0], local_offset[1], local_offset[2], 0, 0, local_offset_point);
  PIDX_set_point_5D(count_local[0], count_local[1], count_local[2], 1, 1, local_box_count_point);
  
  hsize_t count[3];
  count[0] = count_local[0];
  count[1] = count_local[1];
  count[2] = count_local[2];
  
  hsize_t offset[3];
  offset[0] = local_offset[0];
  offset[1] = local_offset[1];
  offset[2] = local_offset[2];
  
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  
  PIDX_create_access(&access);
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
 
  file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, plist_id);
  PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
  PIDX_set_dims(file, global_bounding_box);
  PIDX_set_current_time_step(file, 0);
  PIDX_set_block_size(file, bits_per_block);
  PIDX_set_block_count(file, blocks_per_file);
  
  
  //////////////////////////////////
  buffer = (double*)malloc(sizeof(double) * 512 * 256 * 256/nprocs);
  memset(buffer, 0, sizeof(double) * 512 * 256 * 256/nprocs);
  
  group_id = H5Gopen(file_id, "/Species", H5P_DEFAULT);
  dataset_id = H5Dopen2(group_id, "O2", H5P_DEFAULT);
  
  mem_dataspace = H5Screate_simple (3, count, NULL);
  file_dataspace = H5Dget_space (dataset_id);
  H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace/*H5S_ALL*/, file_dataspace/*H5S_ALL*/, H5P_DEFAULT, buffer);
  
  PIDX_variable_create(file, "Species-O2", sizeof(double) * 8, "1*float64", &variable);
  PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
  PIDX_flush(file);
  
  H5Sclose(mem_dataspace);
  H5Sclose(file_dataspace);
  H5Dclose(dataset_id);
  H5Gclose(group_id);
  /////////////////////////////////
  
  
  ////////////////////////////////
  group_id = H5Gopen(file_id, "/Flow properties", H5P_DEFAULT);
  dataset_id = H5Dopen2(group_id, "Pressure (atm)", H5P_DEFAULT);
  
  mem_dataspace = H5Screate_simple (3, count, NULL);
  file_dataspace = H5Dget_space (dataset_id);
  H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer);
  
  PIDX_variable_create(file, "Flow-properties-pressure", sizeof(double) * 8, "1*float64", &variable2);
  PIDX_append_and_write_variable(variable2, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
  PIDX_flush(file);
  
  H5Sclose(mem_dataspace);
  H5Sclose(file_dataspace);
  H5Dclose(dataset_id);
  H5Gclose(group_id);
  /////////////////////////////////
  
  ////////////////////////////////
  group_id = H5Gopen(file_id, "/Flow properties", H5P_DEFAULT);
  dataset_id = H5Dopen2(group_id, "Vorticity (s^-1)", H5P_DEFAULT);
  
  mem_dataspace = H5Screate_simple (3, count, NULL);
  file_dataspace = H5Dget_space (dataset_id);
  H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer);
  
  PIDX_variable_create(file, "Flow-properties-vorticity", sizeof(double) * 8, "1*float64", &variable3);
  PIDX_append_and_write_variable(variable3, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
  PIDX_flush(file);
  
  H5Sclose(mem_dataspace);
  H5Sclose(file_dataspace);
  H5Dclose(dataset_id);
  H5Gclose(group_id);
  /////////////////////////////////
  
  
  //////////
  H5Pclose(plist_id);
  H5Fclose(file_id);
  
  PIDX_close(&file);
  PIDX_close_access(&access);
  //////////
  
  free(buffer);
  buffer = 0;
  
  MPI_Finalize();
  return 0;
}