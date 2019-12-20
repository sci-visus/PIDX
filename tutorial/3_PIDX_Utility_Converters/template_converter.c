/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

#define PIDX_IO    0
#define HDF5_IO    0 && PIDX_HAVE_HDF5
#define PNETCDF_IO 0 && PIDX_HAVE_PNETCDF
#define NETCDF_IO  1 && PIDX_HAVE_NETCDF
#define NVISUS_IO  1 && PIDX_HAVE_NVISUSIO

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif

#if PIDX_HAVE_HDF5
  #include <hdf5.h>
#endif
  
#if PIDX_HAVE_PNETCDF
  #include <pnetcdf.h>
#endif
  
#if PIDX_HAVE_NETCDF
  #include <netcdf.h>
  #include <netcdf_par.h>
#endif

static int parse_args(int argc, char **argv);
static void usage(void);
static void delete_buffers();

static size_t global_box_size[4] = {0, 0, 0, 0}; ///< global dimensions of 3D volume
static size_t local_box_size[4] = {0, 0, 0, 0};  ///< local dimensions of the per-process block
static int time_step_count = 1;                           ///< Number of time-steps
static int variable_count = 1;                            ///< Number of fields
static char **var_name;
static char output_file_template[512] = "test_idx";       ///< output IDX file Name Template
static char var_file[512];
static char input_file[512];
static char **file_name;
static float **buffer;
static int *values_per_sample;                            ///< Example: 1 for scalar 3 for vector

///< Print error and exit program
static void handle_error(int status, char *message, char *filename, int line)
{
#if PIDX_HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(stderr,"(%d) ",rank);
#endif

#if PNETCDF_IO
  fprintf(stderr, "Error at %s:%d: %s (%s)\n", filename, line, ncmpi_strerror(status), message);
#elif HDF5_IO
  fprintf(stderr, "Error at %s:%d: %s (%s)\n", filename, line, hdf5_strerror(status), message); //or whatever is hdf5 error function
#elif NETCDF_IO
  fprintf(stderr, "Error at %s:%d: %s (%s)\n", filename, line, nc_strerror(status), message);
#else
  fprintf(stderr, "Error at %s:%d: status = %d, message = %s\n", filename, line, status, message);
#endif

#if PIDX_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#endif
  exit(-1);
}


int main(int argc, char **argv)
{
  int ret,i;
  int t;
  int var;
  int slice = 0;
  int nprocs = 1, rank = 0;
  char output_file_name[1024];
  size_t local_box_offset[4];

  // MPI initialization
#if PIDX_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#if NETCDF_IO
  const char *str=nc_inq_libvers();
  printf("NetCDF version: %s\n", str);
#endif

  ret = parse_args(argc, argv);
  if (ret < 0)
  {
    usage();
    handle_error(-1,"syntax error",__FILE__,__LINE__);
  }

#if 0
  // check if the num procs is appropriate
  int num_bricks = (global_box_size[0] / local_box_size[0]) * (global_box_size[1] / local_box_size[1]) * (global_box_size[2] / local_box_size[2]);

  if(num_bricks != nprocs)
  {
    fprintf(stderr, "Error: number of sub-blocks (%d) doesn't match number of procs (%d)\n", num_bricks, nprocs);
    fprintf(stderr, "Incorrect distribution of data across processes i.e.\n(global_x / local_x) X (global_x / local_x) X (global_x / local_x) != nprocs\n(%d/%d) X (%d/%d) X (%d/%d) != %d\n", (int)global_box_size[0], (int)local_box_size[0], (int)global_box_size[1], (int)local_box_size[1], (int)global_box_size[2], (int)local_box_size[2], nprocs);

#if PIDX_HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(-1);
#endif
  }
#endif

#if NETCDF_IO
  int ncid, varid;
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
  //int open_mode=NC_MPIIO;
  int open_mode=NC_MPIPOSIX;   //<ctc> very slightly faster on my mac laptop
  if ((ret = nc_open_par(file_name[0], open_mode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid)))
    handle_error(ret,"",__FILE__,__LINE__);
  MPI_Barrier(MPI_COMM_WORLD);

  /* Get the varid of the data variable, based on its name. */
  if ((ret = nc_inq_varid(ncid, var_name[0], &varid)))
    handle_error(ret,"",__FILE__,__LINE__);

  /* Get global_box_size */
  //todo: use nc_inq_varndims, etc...
  global_box_size[0]=5760;
  global_box_size[1]=2881;
  global_box_size[2]=72;
  global_box_size[3]=1;
  local_box_size[3]=1;

  // set collective io
  if (ret = ret=nc_var_par_access(ncid,varid,NC_COLLECTIVE))
    handle_error(ret,"",__FILE__,__LINE__);

  /* Check the counts of dims, vars, and atts. */
  int ndims,nvars,natts,unlimdimid;
  if ((ret = nc_inq(ncid, &ndims, &nvars, &natts, &unlimdimid)))
    handle_error(ret,"",__FILE__,__LINE__);

  /* Check dims. */
  int d;
  for (d = 0; d < ndims; d++)
  {
    char name[NC_MAX_NAME + 1];
    size_t len;
    if ((ret=nc_inq_dim(ncid, d, name, &len)))
      handle_error(ret,"",__FILE__,__LINE__);
    printf("(%d) d %d: %s (len=%zu)\n",rank,d,name,len);
  }
#endif

  // Creating the filename
  //sprintf(output_file_name, "%s%s", output_file_template,".idx");
  sprintf(output_file_name, "%s", output_file_template);

  // Calculating every process data's offset and size
  int sub_div[4];
  sub_div[0] = (global_box_size[0] / local_box_size[0]);
  sub_div[1] = (global_box_size[1] / local_box_size[1]);
  sub_div[2] = (global_box_size[2] / local_box_size[2]);
  sub_div[3] = (global_box_size[3] / local_box_size[3]);
  local_box_offset[3] = 0;
  local_box_offset[2] = (rank / (sub_div[0] * sub_div[1])) * local_box_size[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  local_box_offset[1] = (slice / sub_div[0]) * local_box_size[1];
  local_box_offset[0] = (slice % sub_div[0]) * local_box_size[0];
#if HDF5_IO
  hid_t file_id, plist_id, dataset_id;
  hid_t file_dataspace, mem_dataspace;
#endif


#if NETCDF_IO
  /* Read the data. */  
  printf("(%d) offset: %d %d %d %d, size: %d %d %d %d\n",rank,local_box_offset[0],local_box_offset[1],local_box_offset[2],local_box_offset[3],local_box_size[0],local_box_size[1],local_box_size[2],local_box_size[3]);

  //reverse order
  size_t reversed_local_box_offset[4];
  size_t reversed_local_box_size[4];
  for (i=0;i<=3;i++)
  {
    reversed_local_box_size[i]=local_box_size[3-i];
    reversed_local_box_offset[i]=local_box_offset[3-i];
  }

  if ((ret = nc_get_vara(ncid, varid, reversed_local_box_offset, reversed_local_box_size, *buffer)))
    handle_error(ret,"",__FILE__,__LINE__);
  printf("(%d) read data successfully!\n",rank);

  /* Close the file, freeing all resources. */
  if ((ret = nc_close(ncid)))
    handle_error(ret,"",__FILE__,__LINE__);
#endif

#if PNETCDF_IO
  //char *filename="/Users/cam/data/uvcdat/c1440_NR.inst30mn_3d_CLOUD_Nv.20060616_1100z.nc4";
  char *filename="/Users/cam/data/uvcdat/tas_Amon_CESM1-CAM5-1-FV2_historical_r1i1p1_185001-200512.nc";
  char *var_name="CLOUD";
  nc_type xtypep;
  int varidp;

  int ncfile,ndims,nvars,ngatts,unlimited;
  MPI_Offset *dim_sizes, var_size;

  ret = ncmpi_open(MPI_COMM_WORLD,filename, NC_NOWRITE, MPI_INFO_NULL, &ncfile);
  if (ret != NC_NOERR) handle_error(ret, "", __FILE__, __LINE__);

  /* no commnunication needed after ncmpi_open: all processors have a cached
   * veiw of the metadata once ncmpi_open returns */

  /* reader knows nothing about dataset, but we can interrogate with query
   * routines: ncmpi_inq tells us how many of each kind of "thing"
   * (dimension, variable, attribute) we will find in the file  */

  ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
  if (ret != NC_NOERR) handle_error(ret,"", __FILE__, __LINE__);

  printf("ndims: %d, nvars: %d, ngatts: %d, unlimited: %d\n",ndims,nvars,ngatts,unlimited);

  /* we do not really need the name of the dimension or the variable for
   * reading in this example.  we could, in a different example, take the
   * name of a variable on the command line and read just that one */

  dim_sizes = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
  /* netcdf dimension identifiers are allocated sequentially starting
   * at zero; same for variable identifiers */
  for(i=0; i<ndims; i++)  
  {
    ret = ncmpi_inq_dimlen(ncfile, i, &(dim_sizes[i]) );
    if (ret != NC_NOERR) handle_error(ret,"", __FILE__, __LINE__);
  }

#if 0 //finish me!
  for(i=0; i<nvars; i++) { 
    /* much less coordination in this case compared to rank 0 doing all
     * the i/o: everyone already has the necessary information */
    ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                        &var_natts);
    if (ret != NC_NOERR) handle_error(ret,"", __FILE__, __LINE__);

    start = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));
    count = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));

    /* we will simply decompose along one dimension.  Generally the
     * application has some algorithim for domain decomposistion.  Note
     * that data decomposistion can have an impact on i/o performance.
     * Often it's best just to do what is natural for the application,
     * but something to consider if performance is not what was
     * expected/desired */

    start[0] = (dim_sizes[dimids[0]]/nprocs)*rank;
    count[0] = (dim_sizes[dimids[0]]/nprocs);
    var_size = count[0];

    for (j=1; j<var_ndims; j++) {
      start[j] = 0;
      count[j] = dim_sizes[dimids[j]];
      var_size *= count[j];
    }

    switch(type) {
      case NC_INT:
        data = (int*) calloc(var_size, sizeof(int));
        ret = ncmpi_get_vara_int_all(ncfile, i, start, count, data);
        if (ret != NC_NOERR) handle_error(ret,"", __FILE__, __LINE__);
        break;
      default:
        /* we can do this for all the known netcdf types but this
         * example is already getting too long  */
        fprintf(stderr, "unsupported NetCDF type \n");
    }

    free(start);
    free(count);
    if (data != NULL) free(data);
  }
#endif
#endif

#if HDF5_IO
  plist_id = H5Pcreate(H5P_FILE_ACCESS);

#if PIDX_HAVE_MPI
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#endif

#if PIDX_IO
  PIDX_point global_bounding_box, local_offset, local_size;
  PIDX_set_point_5D(global_bounding_box, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);
  PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);
  PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);

  PIDX_file file;
  PIDX_access access;
  PIDX_variable *variable;

  PIDX_create_access(&access);

#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  PIDX_time_step_caching_ON();
#endif
  for (t = 0; t < time_step_count; t++)
  {
#if HDF5_IO
    file_id = H5Fopen(file_name[t], H5F_ACC_RDONLY, plist_id);
    for(var = 0; var < variable_count; var++)
    {
      memset(buffer[var], 0, sizeof(double) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
      dataset_id = H5Dopen2(file_id, var_name[var], H5P_DEFAULT);

      mem_dataspace = H5Screate_simple (3, local_box_size, NULL);
      file_dataspace = H5Dget_space (dataset_id);
      H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, local_box_offset, NULL, local_box_size, NULL);

      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, buffer[var]);

      H5Sclose(mem_dataspace);
      H5Sclose(file_dataspace);
      H5Dclose(dataset_id);
    }
    H5Fclose(file_id);
#endif

#if PIDX_IO
    variable = malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);
    PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, &file);
    PIDX_set_variable_count(file, variable_count);
    PIDX_set_dims(file, global_bounding_box);
    PIDX_set_current_time_step(file, t);
    int64_t restructured_box_size[5] = {32, 32, 32, 1, 1};
    ret = PIDX_set_restructuring_box(file, restructured_box_size);
    if (ret != PIDX_success)  
      handle_error(-1,"PIDX_set_restructuring_box",__FILE__,__LINE__);

    for(var = 0; var < variable_count; var++)
    {
      printf("PIDX create variable %s\n",var_name[var]);
      ret = PIDX_variable_create(var_name[var], sizeof(float) * 8, "1*float32", &variable[var]);
      if (ret != PIDX_success)  handle_error(-1, "PIDX_variable_data_layout", __FILE__, __LINE__);

      printf("PIDX variable write data layout %s\n",var_name[var]);
      ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, buffer[var], PIDX_row_major);
      if (ret != PIDX_success)  handle_error(-1, "PIDX_variable_data_layout", __FILE__, __LINE__);

      printf("PIDX append and write variable %s\n",var_name[var]);
      ret = PIDX_append_and_write_variable(file, variable[var]);
      if (ret != PIDX_success)  handle_error(-1, "PIDX_append_and_write_variable", __FILE__, __LINE__);
    }
    printf("closing idx...\n");
    PIDX_close(file);
    printf("done!\n");
    free(variable);
#endif

  }
#if PIDX_IO
  printf("closing pidx access...\n");
  PIDX_time_step_caching_OFF();
  PIDX_close_access(access);
  printf("closed!\n");
#endif

#if HDF5_IO
  H5Pclose(plist_id);
#endif

  delete_buffers();

#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

static void delete_buffers()
{
  int i;

  for (i = 0; i < variable_count; i++)
  {
    free(buffer[i]);
    buffer[i] = 0;
  }
  free(buffer);
  buffer = 0;

  free(values_per_sample);
  values_per_sample = 0;

  for (i = 0; i < variable_count; i++)
  {
    free(var_name[i]);
    var_name[i] = 0;
  }
  free(var_name);
  var_name = 0;

  for (i = 0; i < time_step_count; i++)
  {
    free(file_name[i]);
    file_name[i] = 0;
  }
  free(file_name);
  file_name = 0;

  return;
}

///< Parse the input arguments
static int parse_args(int argc, char **argv)
{
  int i;
  for (i=0;i<argc;i++)
    printf("%s\n",argv[i]);

  int ret;
  //char flags[] = "v:l:f:i:t:";
  char flags[] = "l:i:v:f:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
      // case('g'):
      //     sscanf(optarg, "%zux%zux%zu", &global_box_size[0], &global_box_size[1], &global_box_size[2]);
      //     break;
      case('l'):
          sscanf(optarg, "%zux%zux%zu", &local_box_size[0], &local_box_size[1], &local_box_size[2]);
          break;
      case('f'):
          sprintf(output_file_template, "%s", optarg);
          break;
      case('i'):
          sprintf(input_file, "%s", optarg);
          break;
      case('v'):
          sprintf(var_file, "%s", optarg);
          break;
      case('?'):
      default:
          return (-1);
    }
  }
#if 0
    /* need positive dimensions */
  if (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1 || local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1)
  {
    fprintf(stderr, "Error: bad dimension specification.\n");
    return (-1);
  }

  /* need global dimension to be larger than the local */
  if (global_box_size[0] < local_box_size[0] || global_box_size[1] < local_box_size[1] || global_box_size[2] < local_box_size[2])
  {
    fprintf(stderr, "Error: Per-process local box size cannot be greater than the global box\n");
    return (-1);
  }

  if (global_box_size[0] == 0 || global_box_size[1] == 0 || global_box_size[2] == 0)
  {
    fprintf(stderr, "Global Dimension cannot be 0!!!!!!!!!\n");
    return (-1);
  }
#endif

  if (local_box_size[0] == 0 || local_box_size[1] == 0 || local_box_size[2] == 0)
  {
    fprintf(stderr, "Local Dimension cannot be 0!!!!!!!!!\n");
    return (-1);
  }

  FILE *fp = fopen(var_file, "r");
  ret = fscanf(fp, "%d", &variable_count);
  printf("var_file: %s, var_count: %d\n",var_file,variable_count);
  if (ret != EOF && ret != 1)
    return (-1);

  var_name = malloc(sizeof(char*) * variable_count);
  memset(var_name, 0, sizeof(char*) * variable_count);

  buffer = malloc(sizeof(float*) * variable_count);
  memset(buffer, 0, sizeof(float*) * variable_count);
  for (i = 0; i < variable_count; i++)
  {
    buffer[i] = malloc(sizeof(float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
    //memset(buffer[i], 0, sizeof(float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
  }

  values_per_sample = malloc(sizeof(*values_per_sample) * variable_count);
  memset(values_per_sample, 0, sizeof(*values_per_sample) * variable_count);

  for (i = 0; i < variable_count; i++)
  {
    char temp_var_name[1024];
    ret = fscanf(fp, "%s %d", temp_var_name, &values_per_sample[i]);
    if (ret != 2 || ret == EOF)
      return (-1);
    var_name[i] = strdup(temp_var_name);
  }
  fclose(fp);


  fp = fopen(input_file, "r");
  ret = fscanf(fp, "%d", &time_step_count);
  if (ret != EOF && ret != 1)
    return (-1);
  file_name = malloc(sizeof(char*) * time_step_count);
  memset(file_name, 0, sizeof(char*) * time_step_count);

  for (i = 0; i < time_step_count; i++)
  {
    char temp_file_name[1024];
    ret = fscanf(fp, "%s", temp_file_name);
    if (ret != 1 || ret == EOF)
      return (-1);
    file_name[i] = strdup(temp_file_name);
  }
  fclose(fp);

  return (0);
}


///< How to use this progam
static void usage(void)
{
  printf("Serial Usage: ./template_converter -l 4x4x4 -v var_list -i input_file_list -f output_idx_file_name\n");
  printf("Parallel Usage: mpirun -n 8 ./template_converter -v var_list -i input_file_list -f output_idx_file_name -l 2x2x2 -f filename.idx -t 1 -v 1\n");
  printf("  -i: list of input files (space-separated)\n");
  printf("  -v: list of input fields (space-separated, must all have same size domain)\n");
  printf("  -l: local (per-process) dimensions\n");
  printf("  -t: timesteps (space-separated range + step, e.g. 0 100 10)\n");
  printf("  -f: IDX Filename\n");
  printf("\n");

  return;
}

