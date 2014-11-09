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

#include <PIDX.h>
#include "hdf5.h"


static char *output_file_name;
static double** buffer;

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
  const int bits_per_block = 15;
  const int blocks_per_file = 256;
  PIDX_variable variable;
      
  int nprocs=1, rank=0, slice;
  
#if PIDX_HAVE_MPI
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank); 
#endif
  
  output_file_name = (char*) malloc(sizeof (char) * 1024);
  sprintf(output_file_name, "%s%s", "/scratch/project/visus/datasets/flame", ".idx");
  
  if (nprocs == 2)
  {
    count_local[0] = 256;
    count_local[1] = 256;
    count_local[2] = 256;
  }
  else if (nprocs == 4)
  {
    count_local[0] = 128;
    count_local[1] = 256;
    count_local[2] = 256;
  }
  else if (nprocs == 8)
  {
    count_local[0] = 128;
    count_local[1] = 128;
    count_local[2] = 256;
  }
  else if (nprocs == 16)
  {
    count_local[0] = 128;
    count_local[1] = 128;
    count_local[2] = 128;
  }
  else if (nprocs == 32)
  {
    count_local[0] = 64;
    count_local[1] = 128;
    count_local[2] = 128;
  }
  else if (nprocs == 64)
  {
    count_local[0] = 64;
    count_local[1] = 64;
    count_local[2] = 128;
  }
  else if (nprocs == 128)
  {
    count_local[0] = 64;
    count_local[1] = 64;
    count_local[2] = 64;
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
  
  int group_count = 1;
  char *group_name;
  
  char dataset_count = 5;
  char *dataset_name[5];
  char *pidx_dataset_name[5];
  
  int variable_count = 4;
  
  group_name = strdup("/data");
  
  dataset_name[0] = strdup("P");
  dataset_name[1] = strdup("U");
  dataset_name[2] = strdup("V");
  dataset_name[3] = strdup("W");
  dataset_name[4] = strdup("ZMIX");
  

  pidx_dataset_name[0] = strdup("P");
  pidx_dataset_name[1] = strdup("U");
  pidx_dataset_name[2] = strdup("V");
  pidx_dataset_name[3] = strdup("W");
  pidx_dataset_name[4] = strdup("ZMIX");
  
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  PIDX_create_access(&access);

#if PIDX_HAVE_MPI  
  H5Pset_fapl_mpio(plist_id, comm, info);
  PIDX_set_mpi_access(access, comm);
#endif
 
  
  FILE *fp;
  fp = fopen("list", "r");
  char file_name[1][1024];
  int time_count = 0;
  while (!feof(fp)) 
  {
    if (fscanf(fp, "%s", file_name[time_count]) != 1)
      break;
    if(rank == 0)
      printf("%s\n", file_name[time_count]);
    time_count++;
  }
  fclose(fp);
  
  int g = 0, d = 0, t = 0, time_step, var_count, var = 0;
  time_step = time_count;
  if (rank == 0)
    printf("Number of timesteps = %d\n", time_step);
  for (t = 0; t < time_step; t++)
  {
    var_count = 0;
    buffer = malloc(sizeof(double*) * variable_count);
    memset(buffer, 0, sizeof(double*) * variable_count);
    
    file_id = H5Fopen(file_name[t], H5F_ACC_RDONLY, plist_id);
     
    PIDX_file_create(output_file_name, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, t);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    PIDX_set_variable_count(file, 34);
    
    group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
    
    for(d = 0; d < dataset_count; d++)
    {
      dataset_id = H5Dopen2(group_id, dataset_name[d], H5P_DEFAULT);
    
      mem_dataspace = H5Screate_simple (3, count, NULL);
      file_dataspace = H5Dget_space (dataset_id);
      H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      
      if ( var_count == 0)
      {
	for (var = 0; var < variable_count; var++)
	{
	  buffer[var] = malloc(sizeof(double) * 512 * 256 * 256/nprocs);
	  memset(buffer[var], 0, sizeof(double) * 512 * 256 * 256/nprocs);
	}
      }
      //printf("[%d] Writing %d variable from %d group: %s\n", t, g, d, dataset_name[g][d] );
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer[var_count]);
      PIDX_variable_create(file, pidx_dataset_name[d], sizeof(double) * 8, "1*float64", &variable);
      PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, buffer[var_count], PIDX_column_major);
      
      var_count++;
      
      if (rank == 0)
	printf("[%d] Writing %d variable from 1 group: %s\n", t, d, dataset_name[d]);
      
      if (var_count == variable_count)
      {
	PIDX_flush(file);
	if(rank == 0)
	  printf("Flushed %d variables\n", variable_count);
	for (var = 0; var < variable_count; var++)
	{
	  free(buffer[var]);
	  buffer[var] = 0;
	}
	var_count = 0;
      }
      
      H5Sclose(mem_dataspace);
      H5Sclose(file_dataspace);
      H5Dclose(dataset_id);
    }
    
    H5Gclose(group_id);
    
    PIDX_close(file);
    H5Fclose(file_id);
    
    free(buffer);
    buffer = 0;
  }
  
  //////////
  H5Pclose(plist_id);
  
  
  
  PIDX_close_access(access);
  
  //////////
  
  free(output_file_name);
  output_file_name = 0;
  
#if PIDX_HAVE_MPI    
  MPI_Finalize();
#endif
  
  return 0;
}
#endif //PIDX_OPTION_HDF5
