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

// TODO: allow the user to specify the number of variables read before flushing
// TODO: check the dimension for each dataset (variable) in the HDF5 file

#include <PIDX.h>
#include <pnetcdf.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#if PIDX_HAVE_MPI
#include <mpi.h>
#endif

#define PIDX_ROW_OR_COLUMN_MAJOR PIDX_row_major

enum { X, Y, Z, NUM_DIMS };
enum AtomicType { CHAR = 1, UCHAR = 1, INT = 4, UINT = 4, FLOAT = 4, DOUBLE = 8, INVALID };
struct Type
{
  enum AtomicType atomic_type; ///< e.g. int, float, etc
  unsigned long long num_values; ///< e.g. 1, 2, 3
};
static unsigned long long global_box_size[NUM_DIMS] = { 0 }; ///< global dimensions of 3D volume
static unsigned long long local_box_size[NUM_DIMS] = { 0 }; ///< local dimensions of the per-process block
static unsigned long long local_box_offset[NUM_DIMS] = { 0 };
static PIDX_point pidx_global_box_size = { 0 };
static PIDX_point pidx_local_box_offset = { 0 };
static PIDX_point pidx_local_box_size = { 0 };
static int process_count = 1; ///< Number of processes
static int rank = 0;
static int time_step_count = 1; ///< Number of time-steps
static int var_count = 1; ///< Number of fields
static char output_file_template[512] = "test"; ///< output IDX file Name Template
static char output_file_name[512] = "test.idx";
static char var_file[512];
static char netcdf_file_list[512];
static char **netcdf_var_names = 0;
static char **pidx_var_names = 0;
static char **netcdf_file_names = 0;
static struct Type *var_types = 0;
static PIDX_variable *pidx_vars = 0;
static void *var_data = 0;

static char *usage = "Serial Usage: ./pnetcdf-to-idx -g 4x4x4 -l 4x4x4 -v var_list -i pnetcdf_file_names_list -f output_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./netcdf-to-idx -g 4x4x4 -l 2x2x2 -f Filename_ -v var_list -i netcdf_file_names_list\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -i: file containing list of input netcdf files\n"
                     "  -v: file containing list of input fields\n";

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

//----------------------------------------------------------------
static void rank_0_print(const char *format, ...)
{
  if (rank != 0) return;
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
}

//----------------------------------------------------------------
static void init_mpi(int argc, char **argv)
{
#if PIDX_HAVE_MPI
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
#endif
}

//----------------------------------------------------------------
static void shutdown_mpi()
{
#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
}

//----------------------------------------------------------------
// Read all the lines in a file into a list.
// Return the number of lines read.
// A line starting with // is treated as a comment and does not increase the count.
static int read_list_in_file(const char *file_name, char ***list_ptr)
{
  assert(file_name != 0);
  assert(list_ptr != 0);

  rank_0_print("Opening file %s\n", file_name);
  FILE *fp = 0;
  if ((fp = fopen(file_name, "r")) == 0) {
    terminate_with_error_msg("ERROR: Cannot open file here %s\n", file_name);
}

  // first pass, count the number of non-comment lines
  int line_count = 0;
  char line[512]; // a line cannot be longer than 512 characters
  //char *line = malloc(sizeof(char) * 512);
  while (fgets(line, 512, fp) != 0)
  {
    if (line[0] != '/' || line[1] != '/')
      ++line_count;
  }
  // second pass, actually read the data
  *list_ptr = (char **)calloc(line_count, sizeof(*list_ptr));
  char **list = *list_ptr;
  if (list == 0)
    terminate_with_error_msg("ERROR: Failed to allocate memory for a list of names. Bytes requested = %d (items) * %u (bytes)\n", line_count, sizeof(*list_ptr));
  rewind(fp);
  int i = 0;
  while (fgets(line, sizeof(line), fp) != 0)
  {
    if (line[0] != '/' || line[1] != '/')
    {
      line[strcspn(line, "\r\n")] = 0; // trim the newline character at the end if any
      list[i] = strdup(line);
      ++i;
    }
  }

  fclose(fp);
  return line_count;
}

//----------------------------------------------------------------
static void free_memory()
{
//  free(var_data);
  var_data = 0;

  int i = 0;
  for (i = 0; i < var_count; i++)
  {
    free(netcdf_var_names[i]);
    netcdf_var_names[i] = 0;
    free(pidx_var_names[i]);
    pidx_var_names[i] = 0;
  }
  free(netcdf_var_names);
  netcdf_var_names = 0;
  free(pidx_var_names);
  pidx_var_names = 0;
  free(var_types);
  free(pidx_vars);

  for (i = 0; i < time_step_count; i++)
  {
    free(netcdf_file_names[i]);
    netcdf_file_names[i] = 0;
  }
  free(netcdf_file_names);
  netcdf_file_names = 0;
}

//----------------------------------------------------------------
///< Parse the input arguments
static void parse_args(int argc, char **argv)
{
  if (argc != 11)
    terminate_with_error_msg("ERROR: Wrong number of arguments.\n%s", usage);

  char flags[] = "g:l:f:i:v:";
  int opt = 0;
  while ((opt = getopt(argc, argv, flags)) != -1)
  {
    switch (opt)
    {
      case('g'): // global dimension
        if ((sscanf(optarg, "%lldx%lldx%lld", &global_box_size[0], &global_box_size[1], &global_box_size[2]) == EOF) ||
            (global_box_size[0] < 1 || global_box_size[1] < 1 || global_box_size[2] < 1))
          terminate_with_error_msg("Invalid global dimensions\n%s", usage);
        break;
      case('l'): // local dimension
        if ((sscanf(optarg, "%lldx%lldx%lld", &local_box_size[0], &local_box_size[1], &local_box_size[2]) == EOF) ||
            (local_box_size[0] < 1 || local_box_size[1] < 1 || local_box_size[2] < 1))
          terminate_with_error_msg("Invalid local dimension\n%s", usage);
        break;
      case('f'): // output file name
        if (sprintf(output_file_template, "%s", optarg) < 0)
          terminate_with_error_msg("Invalid output file name template\n%s", usage);
        sprintf(output_file_name, "%s%s", output_file_template, ".idx");
        break;
      case('i'): // a file with a list of NetCDF files
        if (sprintf(netcdf_file_list, "%s", optarg) < 0)
          terminate_with_error_msg("Invalid input file\n%s", usage);
        break;
      case('v'): // a file with a list of variables
        if (sprintf(var_file, "%s", optarg) < 0)
          terminate_with_error_msg("Invalid variable file\n%s", usage);
        break;
      default:
        terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

//----------------------------------------------------------------
static void check_args()
{
  if (global_box_size[X] < local_box_size[X] || global_box_size[Y] < local_box_size[Y] || global_box_size[Z] < local_box_size[Z])
    terminate_with_error_msg("ERROR: Global box is smaller than local box in one of the dimensions\n");

  // check if the number of processes given by the user is consistent with the actual number of processes needed
  int brick_count = (int)((global_box_size[X] + local_box_size[X] - 1) / local_box_size[X]) *
                    (int)((global_box_size[Y] + local_box_size[Y] - 1) / local_box_size[Y]) *
                    (int)((global_box_size[Z] + local_box_size[Z] - 1) / local_box_size[Z]);
  if(brick_count != process_count)
    terminate_with_error_msg("ERROR: Number of sub-blocks (%d) doesn't match number of processes (%d)\n", brick_count, process_count);
}

//----------------------------------------------------------------
static void calculate_per_process_offsets()
{
  int sub_div[NUM_DIMS];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
}

//----------------------------------------------------------------
// Convert an atomic NETCDF datatype to native C datatype
static enum AtomicType from_netcdf_atomic_type(nc_type input)
{
  switch(input) {
   case NC_CHAR:
   return CHAR;
   case NC_INT:
    return INT;
   case NC_FLOAT:
   return FLOAT;
   case NC_DOUBLE:
    return DOUBLE;
   default:
  return INVALID;
 }
}

//----------------------------------------------------------------
// Convert a NETCDF datatype to C native or array datatype
/*
static struct Type from_netcdf_type(nc_type xtypep)
{
  struct Type type;
  int *ndimsp;
//  H5T_class_t type_class = H5Tget_class(type_id);
  if (type_class == H5T_ARRAY)
  {

    int num_dims = ncmpi_inq_varndims(int ncid,varid,ndimsp);

    int num_dims =  ncmpi_inq_vardimid  (int ncid, int varid, int dimids[]);
    if (num_dims != 1)
    {
      type.atomic_type = INVALID; // we don't support arrays of more than 1 dimension
      return type;
    }
    if (H5Tget_array_dims2(type_id, &type.num_values) < 0)
    {
      type.atomic_type = INVALID;
      return type;
    }
  }
  else if (xtypep == NC_FLOAT || xtypep == NC_INTEGER)
  {
    type.num_values = 1;
  }
  else // we don't support HD5_COMPOUND datatype for example
  {
    type.atomic_type = INVALID;
  }

  hid_t atomic_type_id = H5Tget_native_type(type_id, H5T_DIR_DESCEND);
  type.atomic_type = from_netcdf_atomic_type(atomic_type_id);
  H5Tclose(atomic_type_id);
  return type;
}
*/
//----------------------------------------------------------------
// Open the first NETCDF file and query all the types for all the variables
static void determine_var_types()
{
  assert(netcdf_file_names != 0);
  assert(netcdf_var_names != 0);
  struct Type type;
  int plist_id,i = 0,varidp,ndimsp;

  int file_id = ncmpi_open(MPI_COMM_WORLD,netcdf_file_names[0], NC_NOWRITE, MPI_INFO_NULL, &plist_id);

 if (file_id != NC_NOERR)
   terminate_with_error_msg("ERROR: Cannot open file %s\n", netcdf_file_names[0]);

  var_types = (struct Type *)calloc(var_count, sizeof(*var_types));

  nc_type xtypep;
  for (i = 0; i < var_count; ++i)
  {
   int dataset_id = ncmpi_inq_varid (plist_id, netcdf_var_names[i], &varidp);
   dataset_id = ncmpi_inq_vartype(plist_id,varidp, &xtypep);
   
   if (dataset_id != NC_NOERR) 
        terminate_with_error_msg("ERROR: Cannot read the datatype of the variable %s\n", netcdf_var_names[i]);


  int num_dims = ncmpi_inq_varndims(plist_id,varidp,&ndimsp);

    if (ndimsp > 3 )
    {
//TODO:probably we can make it more clever to handle more dimensions as they are not related to the variable itself on netcdf
      type.atomic_type = INVALID; // we don't support arrays of more than 3 dimension
      return;
    }
   else if (xtypep == NC_FLOAT || xtypep == NC_INT)
   {
    type.num_values = 1;
   }
  else // we don't support HD5_COMPOUND datatype for example
  {
    type.atomic_type = INVALID;
  }

  type.atomic_type = from_netcdf_atomic_type(xtypep);

    var_types[i] = type;

    if (var_types[i].atomic_type == INVALID)
      terminate_with_error_msg("ERROR: The datatype of the %s variable is not supported\n", netcdf_var_names[i]);
  }
  ncmpi_close(plist_id);
}

//----------------------------------------------------------------
// Return a negative value when failed, otherwise return 0
int read_var_from_netcdf(int file_id, const char *var_name, struct Type type)
{
  int TIMES=1;
  int LATS=local_box_size[1];
  int LONS=local_box_size[0];
  assert(var_name != 0);
  assert(var_data != 0);
  int varidp,ndims,nvars,ngatts,unlimited;
  nc_type xtypep;
//int dataset_id = H5Dopen2(file_id, var_name, H5P_DEFAULT);
 
  int dataset_id = ncmpi_inq_varid(file_id, var_name, &varidp);
  MPI_Offset start[]={local_box_offset[2],local_box_offset[1],local_box_offset[0]};
  MPI_Offset count[]={TIMES,LATS,LONS};

  dataset_id = ncmpi_inq_vartype(file_id,varidp, &xtypep);


//  if (dataset_id !=0)
//    terminate_with_error_msg("ERROR: Failed to open NetCDF dataset for variable %s\n", var_name);

  int read_error = 0;

  if (type.atomic_type == DOUBLE)
    
    read_error = ncmpi_get_vara_double(file_id, varidp, start, count, var_data);
  else if (type.atomic_type == FLOAT){

    ncmpi_begin_indep_data(file_id);
   read_error = ncmpi_get_vara_float(file_id, varidp, start, count, (float *)var_data);
   ncmpi_end_indep_data(file_id);
 if (read_error != NC_NOERR) 
      terminate_with_error_msg("ERROR: Can not read the data for the variable %s \n", var_name);

}

  else if (type.atomic_type == INT)
    read_error = ncmpi_get_vara_int(file_id, varidp, start, count, var_data);
//  else if (type.atomic_type == UINT)
//    read_error = ncmpi_get_vara_uint(file_id, varidp, start, count, var_data);

//  else if (type.atomic_type == CHAR)
//    read_error = ncmpi_get_vara_char(file_id, varidp, start, count, var_data);

  else if (type.atomic_type == UCHAR)
    read_error = ncmpi_get_vara_uchar(file_id, varidp, start, count, var_data);
  else
    terminate_with_error_msg("ERROR: Unsupported type. Type = %d\n", type.atomic_type);

  if (read_error < 0)
    return -1;
  return 0;
}

//----------------------------------------------------------------
static void to_idx_type_string(struct Type type, char *type_string)
{
  assert(type_string != 0);

  if (type.atomic_type == DOUBLE)
    sprintf(type_string, "%lld*float64", type.num_values);
  else if (type.atomic_type == FLOAT)
      sprintf(type_string, "%lld*float32", type.num_values);
  else if (type.atomic_type == INT)
    sprintf(type_string, "%lld*int32", type.num_values);
  else if (type.atomic_type == UINT)
    sprintf(type_string, "%lld*uint32", type.num_values);
  else if (type.atomic_type == CHAR)
    sprintf(type_string, "%lld*int8", type.num_values);
  else if (type.atomic_type == UCHAR)
    sprintf(type_string, "%lld*uint8", type.num_values);
  else
    terminate_with_error_msg("ERROR: Unsupported type. Type = %d\n", type.atomic_type);
}

//----------------------------------------------------------------
// Convert an NetCDF variable name to a PIDX variable name.
// Remove the leading "/" character.
// Replace any "/" and " " character with "-".
static void netcdf_var_name_to_pidx_var_name(const char *netcdf_name, char *pidx_name)
{
  assert(netcdf_name != 0);
  assert(pidx_name != 0);

  int i = 0;
  while (netcdf_name[i] == '/')
    ++i;
  int j = 0;
  while (netcdf_name[i] != '\0')
  {
    pidx_name[j] = (netcdf_name[i] == '/' || netcdf_name[i] == ' ') ? '-' : netcdf_name[i];
    ++j;
    ++i;
  }
  pidx_name[j] = '\0';
}

//----------------------------------------------------------------
static void create_pidx_var_names()
{
  assert(netcdf_var_names != 0);

  pidx_var_names = (char **)calloc(var_count, sizeof(*pidx_var_names));
  if (pidx_var_names == 0)
    terminate_with_error_msg("ERROR: Failed to allocate memory to store PIDX var names\n. Bytes requested = %d (values) * %u (bytes)\n", var_count, sizeof(*pidx_var_names));

  int i = 0;
  for (i = 0; i < var_count; ++i)
  {
    pidx_var_names[i] = (char *)calloc(strlen(netcdf_var_names[i]) + 1, sizeof(*pidx_var_names[i]));
    if (pidx_var_names == 0)
      terminate_with_error_msg("ERROR: Failed to allocate memory to store the PIDX var name for %s\n. Bytes requested = %d (values) * %u (bytes)\n", netcdf_var_names[i], strlen(netcdf_var_names[i]) + 1, sizeof(*pidx_var_names[i]));
    netcdf_var_name_to_pidx_var_name(netcdf_var_names[i], pidx_var_names[i]);
    rank_0_print("NetCDF variable %s becomes\n PIDX variable %s\n", netcdf_var_names[i], pidx_var_names[i]);
  }
}

//----------------------------------------------------------------
static void create_pidx_vars()
{
  assert(netcdf_var_names != 0);
  assert(var_types != 0);

  pidx_vars = (PIDX_variable *)calloc(var_count, sizeof(*pidx_vars));
  int i = 0;
  for (i = 0; i < var_count; ++i)
  {

    char type_string[32] = { 0 };
    to_idx_type_string(var_types[i], type_string);
    int ret = PIDX_variable_create(pidx_var_names[i], var_types[i].atomic_type * 8, type_string, &pidx_vars[i]);
    if (ret != PIDX_success)
      terminate_with_error_msg("ERROR: PIDX failed to create PIDX variable %s\n", pidx_var_names[i]);
  }
}

//----------------------------------------------------------------
static void write_var_to_idx(PIDX_file pidx_file, const char *var_name, PIDX_variable pidx_var)
{
  assert(var_name != 0);
  assert(var_data != 0);
  PIDX_set_point_5D(pidx_local_box_offset, local_box_offset[X], local_box_offset[Y], local_box_offset[Z], 0, 0);
  PIDX_set_point_5D(pidx_local_box_size, local_box_size[X], local_box_size[Y], local_box_size[Z], 1, 1);
  int ret = PIDX_success,i;

  ret = PIDX_variable_write_data_layout(pidx_var, pidx_local_box_offset, pidx_local_box_size, var_data, PIDX_ROW_OR_COLUMN_MAJOR);
  if (ret != PIDX_success)
    terminate_with_error_msg("ERROR: PIDX failed to specify variable data layout for %s\n", var_name);
  ret = PIDX_append_and_write_variable(pidx_file, pidx_var);
  if (ret != PIDX_success)
    terminate_with_error_msg("ERROR: PIDX failed to append and write variable %s\n", var_name);
}

//----------------------------------------------------------------
static void create_pidx_access(PIDX_access *pidx_access)
{
  int ret = PIDX_create_access(pidx_access);
  if (ret != PIDX_success)
    terminate_with_error_msg("ERROR: Failed to create PIDX access\n");
#if PIDX_HAVE_MPI
  ret = PIDX_set_mpi_access(*pidx_access, MPI_COMM_WORLD);
  if (ret != PIDX_success)
    terminate_with_error_msg("ERROR: Failed to create PIDX MPI access\n");
#endif
}

//----------------------------------------------------------------
static void set_pidx_params(PIDX_file pidx_file)
{
  int ret = PIDX_set_point_5D(pidx_global_box_size, global_box_size[X], global_box_size[Y], global_box_size[Z], 1, 1);
  ret = PIDX_set_dims(pidx_file, pidx_global_box_size);
  if (ret != PIDX_success)
    terminate_with_error_msg("ERROR: Failed to set PIDX global dimensions\n");
  ret = PIDX_set_variable_count(pidx_file, var_count);
  if (ret != 0)
    terminate_with_error_msg("ERROR: Failed to set PIDX variable count\n");
//      PIDX_set_block_size(pidx_file,16);
//      PIDX_set_block_count(pidx_file,16);

}

//----------------------------------------------------------------
int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  printf("start %d %d\n",rank,process_count);
   int i;
  if(rank==0) {
    parse_args(argc, argv);
//    var_count = read_list_in_file(var_file, &netcdf_var_names);
//    rank_0_print("Number of variables = %d\n", var_count);
//    time_step_count = read_list_in_file(netcdf_file_list, &netcdf_file_names);
//    rank_0_print("Number of timesteps = %d\n", time_step_count);
    check_args();
  }

   
  //  The command line arguments are shared by all processes
#if PIDX_HAVE_MPI
  MPI_Bcast(global_box_size, 3, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_box_size, 3, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
//  MPI_Bcast(&time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&var_file, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netcdf_file_list, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

//  MPI_Bcast(&var_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&output_file_name, 512, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

//TODO: Only the rank 0 should read the input files and broadcast the data
    var_count = read_list_in_file(var_file, &netcdf_var_names);
    rank_0_print("Number of variables = %d\n", var_count);


    time_step_count = read_list_in_file(netcdf_file_list, &netcdf_file_names);
    rank_0_print("Number of timesteps = %d\n", time_step_count);

  calculate_per_process_offsets();
  
  create_pidx_var_names();

  PIDX_access pidx_access;
  create_pidx_access(&pidx_access);
  PIDX_time_step_caching_ON();

  determine_var_types();
  create_pidx_vars();

  int t = 0,plist_id;
  for (t = 0; t < time_step_count; ++t)
  {
    rank_0_print("Processing time step %d (file %s)\n", t, netcdf_file_names[t]);

    PIDX_file pidx_file;
    int ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, pidx_access, &pidx_file);
  

    if (ret != PIDX_success)
      terminate_with_error_msg("ERROR: Failed to create PIDX file\n");
    set_pidx_params(pidx_file);
    PIDX_set_current_time_step(pidx_file, t);

    int file_id = ncmpi_open(MPI_COMM_WORLD,netcdf_file_names[t], NC_NOWRITE, MPI_INFO_NULL, &plist_id);

    if (file_id !=0)
      terminate_with_error_msg("ERROR: Failed to open file %s\n", netcdf_file_names[t]);

    int v = 0;
    for(v = 0; v < var_count; ++v)
    {
      rank_0_print("Processing variable %s\n", netcdf_var_names[v]);
      var_data = malloc(var_types[v].atomic_type * var_types[v].num_values * local_box_size[0] * local_box_size[1] * local_box_size[2]);

      read_var_from_netcdf(plist_id, netcdf_var_names[v], var_types[v]);
     
      write_var_to_idx(pidx_file, pidx_var_names[v], pidx_vars[v]);
      if (PIDX_flush(pidx_file) != PIDX_success)
        terminate_with_error_msg("ERROR: Failed to flush variable %s, time step %d\n", pidx_var_names[v], t);
      free(var_data);
    }

    ncmpi_close(plist_id);
    PIDX_close(pidx_file);
  }

  PIDX_time_step_caching_OFF();
  PIDX_close_access(pidx_access);

  free_memory();
  shutdown_mpi();

  return 0;
}
