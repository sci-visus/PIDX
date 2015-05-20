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
#include "testdefs.h"

#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#if PIDX_HAVE_MPI
  #include <mpi.h>
#endif


/// usage
static void usage(enum Kind kind)
{
  printf("Usage for pidxtest -k %s:\n\n",kindToStr(kind));

  switch (kind)
  {
  //case SERIAL_WRITER:                      usage_serial();              break;
  case PARALLEL_WRITER:                    usage_multi_idx_writer();    break;
  case SIMPLE_WRITER:                      usage_simple_idx_writer();    break;
  //case PARALLEL_MULTI_PATCH_WRITER:        usage_multi_var_writer();    break;
  //case SERIAL_READER:                      /*usage_serial_reader();*/              break;
  //case PARALLEL_READER:                    usage_reader();              break;
  case DEFAULT:
  default:                                 usage_multi_idx_writer();
  }
}

/// main
int main(int argc, char **argv)
{
  int ret, nprocs = 1, rank = 0;
  struct Args args;

#if PIDX_HAVE_MPI
  /// MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /// Rank 0 parses the command line arguments
  if (rank == 0)
  {
    ret = parse_args(&args, argc, argv);
    if (ret < 0)
    {
      usage(args.kind);
      fprintf(stderr, "File [%s] Line [%d] Error [usage]\n", __FILE__, __LINE__);
#if PIDX_HAVE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    }
    if (args.count_local[0] == 0 || args.count_local[1] == 0 || args.count_local[2] == 0)
    {
      usage(args.kind);
      fprintf(stderr, "File [%s] Line [%d] Error [Local Dimension cannot be 0]\n", __FILE__, __LINE__);
#if PIDX_HAVE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    }
    else
    {
      if ((args.extents[0] / args.count_local[0]) * (args.extents[1] / args.count_local[1]) * (args.extents[2] / args.count_local[2]) != nprocs)
      {
	usage(args.kind);
        fprintf(stderr, "File [%s] Line [%d] Error [Wrong Number of Processes]\n", __FILE__, __LINE__);
#if PIDX_HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
      }
    }
  }

#if PIDX_HAVE_MPI
  MPI_Bcast(&args.kind, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  /// run the specified test
  switch (args.kind)
  {

    case PARALLEL_READER:
      if(rank == 0)
        printf("Performing Parallel Read....\n");
      test_reader(args, rank, nprocs);
      break;
    /*
    case SERIAL_READER:
      if(rank == 0)
	printf("Performing Serial Read....\n");
      //serial_reader(args);
      break;
    */

    case PARALLEL_WRITER:
      if(rank == 0)
        printf("Performing Parallel Write....\n");
      test_multi_idx_writer(args, rank, nprocs);
      break;

    case SIMPLE_WRITER:
      if(rank == 0)
        printf("Performing Simple Parallel Write....\n");
      test_simple_idx_writer(args, rank, nprocs);
      break;

    case PARALLEL_CONVERTER:
      if (rank == 0)
        printf("Performing Parallel Convert....\n");
      test_converter(args, rank);
      break;
    case PARALLEL_CONVERTER2:
      if (rank == 0)
        printf("Performing Parallel Convert 2....\n");
      test_converter2(args, rank);
      break;

    case HDF5_WRITER:
      if(rank == 0)
        printf("Performing Parallel Write....\n");
      test_hdf5_writer(args, rank, nprocs);
      break;

    case HDF5_READER:
      if(rank == 0)
        printf("Performing Parallel Write....\n");
      test_hdf5_reader(args, rank, nprocs);
      break;

    /*
    case PARALLEL_MULTI_PATCH_WRITER:
      if(rank == 0)
	printf("Performing Parallel Write....\n");
      test_multi_patch_writer(args, rank, nprocs);
      break;

    case SERIAL_WRITER:
      if(rank == 0)
	printf("Performing Serial Write....\n");
      serial_writer(args);
      break;

    case SERIAL_READER:
      if(rank == 0)
        printf("Performing Serial Read....\n");
      serial_reader(args);
      break;
    */
    default:
      test_multi_idx_writer(args, rank, nprocs);
  }

#if PIDX_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

/*parse_args*/
int parse_args(struct Args *args, int argc, char **argv)
{
  char strkind[128];
  FILE* config_file;

  if (!args)
    return -1;

  memset(args, 0, sizeof(struct Args));

  config_file = fopen(argv[1], "r");
  if (!config_file)
  {
    fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
    return 1;
  }

  fscanf(config_file, "(mode)\n");
  fscanf(config_file, "%s\n", strkind);
  args->kind = strToKind(strkind);

  fscanf(config_file, "(global box)\n");
  fscanf(config_file, "%lld %lld %lld\n", (long long*)&args->extents[0], (long long*)&args->extents[1], (long long*)&args->extents[2]);

  fscanf(config_file, "(local box)\n");
  fscanf(config_file, "%d %d %d\n", &args->count_local[0], &args->count_local[1], &args->count_local[2]);

  fscanf(config_file, "(restructured box size)\n");
  fscanf(config_file, "%d %d %d\n", &args->restructured_box_size[0], &args->restructured_box_size[1], &args->restructured_box_size[2]);

  fscanf(config_file, "(file name)\n");
  fscanf(config_file, "%s\n", args->output_file_template);

  fscanf(config_file, "(time steps)\n");
  fscanf(config_file, "%d", &args->time_step_count);
  args->time_steps = (int *)malloc(args->time_step_count * sizeof(int));
  int ts = 0;
  for (ts = 0; ts < args->time_step_count; ++ts)
  {
    args->time_steps[ts] = ts; // default
  }
  char c;
  if ((c = fgetc(config_file)) != '\n') // there is a list of time steps following
  {
    for (ts = 0; ts < args->time_step_count; ++ts)
    {
      fscanf(config_file, "%d", &args->time_steps[ts]);
    }
    fscanf(config_file, "\n");
  }

  fscanf(config_file, "(fields)\n");
  fscanf(config_file, "%d", &args->variable_count);
  args->variables = (int *)malloc(args->variable_count * sizeof(int));
  int v = 0;
  for (v = 0; v < args->variable_count; ++v)
  {
    args->variables[v] = v; // default
  }
  if ((c = fgetc(config_file)) != '\n') // there is a list of variables following
  {
    for (v = 0; v < args->variable_count; ++v)
    {
      fscanf(config_file, "%d", &args->variables[v]);
    }
    fscanf(config_file, "\n");
  }

  fscanf(config_file, "(idx count x:y:z)\n");
  fscanf(config_file, "%d %d %d\n", &args->idx_count[0], &args->idx_count[1], &args->idx_count[2]);

  if (strcmp(strkind, "hdf5-writer") != 0 && strcmp(strkind, "hdf5-reader") != 0)
  {
    //if (strcmp(strkind, "parallel-writer") == 0 || strcmp(strkind, "serial-writer") == 0)
    {
      fscanf(config_file, "(debug rst:hz:agg)\n");
      fscanf(config_file, "%d %d %d\n", &args->debug_rst, &args->debug_hz, &args->dump_agg);

      fscanf(config_file, "(perform brst:hz:comp:agg:io)\n");
      fscanf(config_file, "%d %d %d %d %d\n", &args->perform_brst, &args->perform_hz, &args->perform_compression, &args->perform_agg, &args->perform_io);

      fscanf(config_file, "(compression block size)\n");
      fscanf(config_file, "%lld %lld %lld\n", (long long*)&args->chunk_size[0], (long long*)&args->chunk_size[1], (long long*)&args->chunk_size[2]);

      fscanf(config_file, "(compression type:bits)\n");
      fscanf(config_file, "%d %d\n", &args->compression_type, &args->compression_bit_rate);

      fscanf(config_file, "(blocks per file)\n");
      fscanf(config_file, "%d\n", &args->blocks_per_file);

      fscanf(config_file, "(samples per block)\n");
      fscanf(config_file, "%d\n", &args->bits_per_block);

      fscanf(config_file, "(aggregation factor)\n");
      fscanf(config_file, "%d\n", &args->aggregation_factor);

      fscanf(config_file, "(rank z ordering)\n");
      fscanf(config_file, "%d\n", &args->is_rank_z_ordering);

      fscanf(config_file, "(global indexing order)\n");
      fscanf(config_file, "%d\n", &args->is_global_indexing);

      fscanf(config_file, "(hz from:to)\n");
      fscanf(config_file, "%d %d\n", &args->hz_from, &args->hz_to);
    }
  }

  fclose(config_file);

  /* need positive dimensions */
  if (args->extents[0] < 1 || args->extents[1] < 1 || args->extents[2] < 1 || args->count_local[0] < 1 || args->count_local[1] < 1 || args->count_local[2] < 1)
  {
    printf("Error: bad dimension specification.\n");
    return (-1);
  }

  /* need global dimension to be larger than the local */
  if (args->extents[0] < args->count_local[0] || args->extents[1] < args->count_local[1] || args->extents[2] < args->count_local[2])
  {
    printf("Error: global dimensions and local dimensions aren't evenly divisible\n");
    return (-1);
  }

  args->extents[3] = 1;
  args->extents[4] = 1;

  return (0);
}


/// kindToStr
char* kindToStr(enum Kind k)
{
  switch (k)
  {
    case PARALLEL_READER:                return "parallel-reader";
    //case SERIAL_READER:                  return "serial-reader";
    case PARALLEL_WRITER:                return "parallel-writer";
    case SIMPLE_WRITER:                  return "simple-writer";
    case PARALLEL_CONVERTER:             return "parallel-converter";
    case PARALLEL_CONVERTER2:            return "parallel-converter2";
    case HDF5_WRITER:                    return "hdf5-writer";
    case HDF5_READER:                    return "hdf5-reader";
    //case PARALLEL_MULTI_PATCH_WRITER:    return "parallel-multi-patch-writer";
    //case SERIAL_WRITER:                  return "serial-writer";
    case DEFAULT:
    default:                             return "default";
  }
}

/// strToKind
enum Kind strToKind(const char *str)
{
  if (strcmp(str,"parallel-reader")   == 0)             return PARALLEL_READER;
  //if (strcmp(str,"serial-reader")     == 0)             return SERIAL_READER;
  if (strcmp(str,"parallel-writer")   == 0)             return PARALLEL_WRITER;
  if (strcmp(str,"simple-writer")   == 0)             return SIMPLE_WRITER;
  if (strcmp(str, "parallel-converter") == 0)           return PARALLEL_CONVERTER;
  if (strcmp(str, "parallel-converter2") == 0)          return PARALLEL_CONVERTER2;
  if (strcmp(str,"hdf5-writer")   == 0)                 return HDF5_WRITER;
  if (strcmp(str,"hdf5-reader")   == 0)                 return HDF5_READER;
  //if (strcmp(str,"parallel-multi-patch-writer")   == 0) return PARALLEL_MULTI_PATCH_WRITER;
  //if (strcmp(str,"serial-writer")     == 0)             return SERIAL_WRITER;
  else                                                  return DEFAULT;
}
