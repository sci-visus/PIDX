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
  
  int time = 0;
    
  int nprocs, rank, slice;
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank); 
  
  output_file_name = (char*) malloc(sizeof (char) * 512);
  sprintf(output_file_name, "%s%s", "/media/My Passport/Kaust_data/F", ".idx");
  output_file = output_file_name;
  
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
  
  int group_count = 6;
  char *group_name[6];
  
  char dataset_count[6] = {13, 3, 3, 3, 9, 3};
  char *dataset_name[6][13];
  char *pidx_dataset_name[6][13];
  
  group_name[0] = strdup("/Flow properties");
  group_name[1] = strdup("/Intermediate strain eigenvector");
  group_name[2] = strdup("/Most compressive eigenvector");
  group_name[3] = strdup("/Most extensive eigenvector");
  group_name[4] = strdup("/Species");
  group_name[5] = strdup("/Velocity");
  
  dataset_name[0][0] = strdup("C - Scalar Dissipation Rate (s^-1)");
  dataset_name[0][1] = strdup("Density (g.cm^-3)");
  dataset_name[0][2] = strdup("Displacement speed (cm.s^-1)");
  dataset_name[0][3] = strdup("Divergence (s^-1)");
  dataset_name[0][4] = strdup("Flame thickness - thermal (mm)");
  dataset_name[0][5] = strdup("Heat Release Rate (W.cm^-3)");
  dataset_name[0][6] = strdup("Pressure (atm)");
  dataset_name[0][7] = strdup("Progress Variable");
  dataset_name[0][8] = strdup("Q-criterion (s^-2)");
  dataset_name[0][9] = strdup("Temperature (K)");
  dataset_name[0][10] = strdup("Thermal diffusivity (m^2.s^-1)");
  dataset_name[0][11] = strdup("Viscosity (g.cm^-1.s^-1)");
  dataset_name[0][12] = strdup("Vorticity (s^-1)");
  
  dataset_name[1][0] = strdup("Eig-12 (s^-1)");
  dataset_name[1][1] = strdup("Eig-22 (s^-1)");
  dataset_name[1][2] = strdup("Eig-32 (s^-1)");
  
  dataset_name[2][0] = strdup("Eig-11 (s^-1)");
  dataset_name[2][1] = strdup("Eig-21 (s^-1)");
  dataset_name[2][2] = strdup("Eig-31 (s^-1)");
  
  dataset_name[3][0] = strdup("Eig-13 (s^-1)");
  dataset_name[3][1] = strdup("Eig-13 (s^-1)");
  dataset_name[3][2] = strdup("Eig-13 (s^-1)");

  dataset_name[4][0] = strdup("H");
  dataset_name[4][1] = strdup("H2");
  dataset_name[4][2] = strdup("H2)");
  dataset_name[4][3] = strdup("H2O2");
  dataset_name[4][4] = strdup("HO2");
  dataset_name[4][5] = strdup("N2");
  dataset_name[4][6] = strdup("O");
  dataset_name[4][7] = strdup("O2");
  dataset_name[4][8] = strdup("OH");

  dataset_name[5][0] = strdup("U (cm.s^-1)");
  dataset_name[5][1] = strdup("V (cm.s^-1)");
  dataset_name[5][2] = strdup("W (cm.s^-1)");
  
  pidx_dataset_name[0][0] = strdup("Flow-properties-C-Scalar-Dissipation-Rate");
  pidx_dataset_name[0][1] = strdup("Flow-properties-Density");
  pidx_dataset_name[0][2] = strdup("Flow-properties-Displacement-speed");
  pidx_dataset_name[0][3] = strdup("Flow-properties-Divergence");
  pidx_dataset_name[0][4] = strdup("Flow-properties-Flame-thickness-thermal");
  pidx_dataset_name[0][5] = strdup("Flow-properties-Heat-Release-Rate");
  pidx_dataset_name[0][6] = strdup("Flow-properties-Pressure");
  pidx_dataset_name[0][7] = strdup("Flow-properties-Progress-Variable");
  pidx_dataset_name[0][8] = strdup("Flow-properties-Q-criterion");
  pidx_dataset_name[0][9] = strdup("Flow-properties-Temperature");
  pidx_dataset_name[0][10] = strdup("Flow-properties-Thermal-diffusivity");
  pidx_dataset_name[0][11] = strdup("Flow-properties-Viscosity");
  pidx_dataset_name[0][12] = strdup("Flow-properties-Vorticity");
  
  pidx_dataset_name[1][0] = strdup("Intermediate-strain-eigenvector/Eig-12");
  pidx_dataset_name[1][1] = strdup("Intermediate-strain-eigenvector/Eig-22");
  pidx_dataset_name[1][2] = strdup("Intermediate-strain-eigenvector/Eig-32");
  
  pidx_dataset_name[2][0] = strdup("Most-compressive-eigenvector/Eig-11");
  pidx_dataset_name[2][1] = strdup("Most-compressive-eigenvector/Eig-21");
  pidx_dataset_name[2][2] = strdup("Most-compressive-eigenvector/Eig-31");
  
  pidx_dataset_name[3][0] = strdup("Most-extensive-eigenvector/Eig-13");
  pidx_dataset_name[3][1] = strdup("Most-extensive-eigenvector/Eig-13");
  pidx_dataset_name[3][2] = strdup("Most-extensive-eigenvector/Eig-13");

  pidx_dataset_name[4][0] = strdup("Species/H");
  pidx_dataset_name[4][1] = strdup("Species/H2");
  pidx_dataset_name[4][2] = strdup("Species/H2)");
  pidx_dataset_name[4][3] = strdup("Species/H2O2");
  pidx_dataset_name[4][4] = strdup("Species/HO2");
  pidx_dataset_name[4][5] = strdup("Species/N2");
  pidx_dataset_name[4][6] = strdup("Species/O");
  pidx_dataset_name[4][7] = strdup("Species/O2");
  pidx_dataset_name[4][8] = strdup("Species/OH");

  pidx_dataset_name[5][0] = strdup("Velocity/U");
  pidx_dataset_name[5][1] = strdup("Velocity/V");
  pidx_dataset_name[5][2] = strdup("Velocity/W");

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  
  PIDX_create_access(&access);
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
 
  buffer = (double*)malloc(sizeof(double) * 512 * 256 * 256/nprocs);
  memset(buffer, 0, sizeof(double) * 512 * 256 * 256/nprocs);
  /*
  FILE *fp;
  fp = fopen("list", "r");
  char file_name[1024];
  while (!feof(fp)) 
  {
    if (fscanf(fp, "%s", file_name) != 1)
      break;
    printf("%s\n", file_name);
  }
  fclose(fp);
  */
  int g = 0, d = 0, t = 0, time_step = 1;
  for (t = 0; t < time_step; t++)
  {
    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, plist_id);
    PIDX_file_create(output_file, PIDX_file_trunc, access, &file);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, t);
    PIDX_set_block_size(file, bits_per_block);
    PIDX_set_block_count(file, blocks_per_file);
    
    for(g = 0; g < group_count; g++)
    {
      group_id = H5Gopen(file_id, group_name[g], H5P_DEFAULT);
      
      for(d = 0; d < dataset_count[g]; d++)
      {
	//printf("Opening [%d] Group %s with [%d] Dataset %s\n", g, group_name[g], d, dataset_name[g][d]);
	dataset_id = H5Dopen2(group_id, dataset_name[g][d], H5P_DEFAULT);
      
	mem_dataspace = H5Screate_simple (3, count, NULL);
	file_dataspace = H5Dget_space (dataset_id);
	H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer);
	
	PIDX_variable_create(file, pidx_dataset_name[g][d], sizeof(double) * 8, "1*float64", &variable);
	PIDX_append_and_write_variable(variable, local_offset_point, local_box_count_point, buffer, PIDX_column_major);
	PIDX_flush(file);
	
	H5Sclose(mem_dataspace);
	H5Sclose(file_dataspace);
	H5Dclose(dataset_id);
	if(rank == 0)
	  printf("Done writing [%d] Group %s [%d] Dataset %s\n", g, group_name[g], d, dataset_name[g][d]);
      }
      
      H5Gclose(group_id);
    }
  }  
  
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