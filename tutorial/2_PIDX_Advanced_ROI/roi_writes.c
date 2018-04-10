/*
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

#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <PIDX.h>

enum { X, Y, Z, NUM_DIMS };
static int rank = 0;
static unsigned long long global_box_size[3] = {0, 0, 0};
static unsigned long long local_box_offset[3];
static unsigned long long local_box_size[3] = {0, 0, 0};
static int time_step_count = 1;
static char output_file_template[512] = "test";
static float *data;
static char output_file_name[512] = "test.idx";
static PIDX_point global_size, local_offset, local_size;
static int roi_type = 0;
static int current_time_step = 0;
static char *usage = "Serial Usage: ./roi_writes -g 128x128x128 -l 128x128x128 -f out_file_name -r 1\n"
                     "Parallel Usage: mpirun -n 64 ./roi_writes -g 128x128x128 -l 32x32x32 -f out_file_name -r 1\n"
                     "  -g: global dimensions\n"
                     "  -l: local (per-process) dimensions\n"
                     "  -f: IDX filename\n"
                     "  -r: ROI type (1, 2, 4)\n";


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

static void calculate_per_process_offsecurrent_time_step()
{
#if 0
  int sub_div[NUM_DIMS];
  sub_div[X] = (global_box_size[X] / local_box_size[X]);
  sub_div[Y] = (global_box_size[Y] / local_box_size[Y]);
  sub_div[Z] = (global_box_size[Z] / local_box_size[Z]);
  local_box_offset[Z] = (rank / (sub_div[X] * sub_div[Y])) * local_box_size[Z];
  int slice = rank % (sub_div[X] * sub_div[Y]);
  local_box_offset[Y] = (slice / sub_div[X]) * local_box_size[Y];
  local_box_offset[X] = (slice % sub_div[X]) * local_box_size[X];
#endif
#if 1
  if (rank == 0)
  {
    local_box_offset[0] = 0;
    local_box_offset[1] = 5;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 11;
    local_box_size[1] = 17;
    local_box_size[2] = 32;
  }
  
  if (rank == 1)
  {
    local_box_offset[0] = 11;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 18;
    local_box_size[1] = 32;
    local_box_size[2] = 32;
  }
  
  if (rank == 2)
  {
    local_box_offset[0] = 29;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 18;
    local_box_size[1] = 32;
    local_box_size[2] = 32;
  }

  if (rank == 3)
  {
    local_box_offset[0] = 47;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 18;
    local_box_size[1] = 32;
    local_box_size[2] = 32;
  }
  
  if (rank == 4)
  {
    local_box_offset[0] = 65;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 18;
    local_box_size[1] = 32;
    local_box_size[2] = 32;
  }
  
  if (rank == 5)
  {
    local_box_offset[0] = 83;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 18;
    local_box_size[1] = 32;
    local_box_size[2] = 32;
  }
  
  if (rank == 6)
  {
    local_box_offset[0] = 101;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 10;
    local_box_size[1] = 42;
    local_box_size[2] = 32;
  }
  
  if (rank == 7)
  {
    local_box_offset[0] = 111;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;
    
    local_box_size[0] = 21;
    local_box_size[1] = 52;
    local_box_size[2] = 32;
  }
#endif

#if 0
  [PIDX 0] [0] Offset and Count -1 9 -1 : 11 11 32
  [PIDX 1] [0] Offset and Count -1 20 -1 : 11 11 32
  [PIDX 2] [0] Offset and Count 9 -1 -1 : 11 11 32
  [PIDX 3] [0] Offset and Count 10 10 -1 : 10 10 32
  [PIDX 4] [0] Offset and Count 10 20 -1 : 10 11 32
  [PIDX 5] [0] Offset and Count 20 -1 -1 : 10 11 32
  [PIDX 6] [0] Offset and Count 20 10 -1 : 10 10 32
  [PIDX 7] [0] Offset and Count 20 20 -1 : 10 11 32
  [PIDX 8] [0] Offset and Count 30 -1 -1 : 10 11 32
  [PIDX 9] [0] Offset and Count 30 10 -1 : 10 10 32
  [PIDX 10] [0] Offset and Count 30 20 -1 : 10 11 32
  [PIDX 11] [0] Offset and Count 40 -1 -1 : 10 11 32
  [PIDX 12] [0] Offset and Count 40 10 -1 : 10 10 32
  [PIDX 13] [0] Offset and Count 40 20 -1 : 10 11 32
  [PIDX 14] [0] Offset and Count 50 -1 -1 : 10 11 32
  [PIDX 15] [0] Offset and Count 50 10 -1 : 10 10 32
  [PIDX 16] [0] Offset and Count 50 20 -1 : 10 10 32
  [PIDX 17] [0] Offset and Count 49 30 -1 : 11 11 32
  [PIDX 18] [0] Offset and Count 60 -1 -1 : 11 11 32
  [PIDX 19] [0] Offset and Count 60 10 -1 : 11 10 32
  [PIDX 20] [0] Offset and Count 60 20 -1 : 11 10 32
  [PIDX 21] [0] Offset and Count 60 30 -1 : 11 10 32
  [PIDX 22] [0] Offset and Count 59 40 -1 : 12 11 32

  
  //[PIDX 0] [0] Offset and Count -1 9 -1 : 11 11 32
  if (rank == 0)
  {
    local_box_offset[0] = 0;
    local_box_offset[1] = 10;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 1] [0] Offset and Count -1 20 -1 : 11 11 32
  if (rank == 1)
  {
    local_box_offset[0] = 0;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 2] [0] Offset and Count 9 -1 -1 : 11 11 32
  if (rank == 2)
  {
    local_box_offset[0] = 10;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 3] [0] Offset and Count 10 10 -1 : 10 10 32
  if (rank == 3)
  {
    local_box_offset[0] = 11;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 4] [0] Offset and Count 10 20 -1 : 10 11 32
  if (rank == 4)
  {
    local_box_offset[0] = 11;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 5] [0] Offset and Count 20 -1 -1 : 10 11 32
  if (rank == 5)
  {
    local_box_offset[0] = 21;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 6] [0] Offset and Count 20 10 -1 : 10 10 32
  if (rank == 6)
  {
    local_box_offset[0] = 21;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 7] [0] Offset and Count 20 20 -1 : 10 11 32
  if (rank == 7)
  {
    local_box_offset[0] = 21;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 8] [0] Offset and Count 30 -1 -1 : 10 11 32
  if (rank == 8)
  {
    local_box_offset[0] = 31;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 9] [0] Offset and Count 30 10 -1 : 10 10 32
  if (rank == 9)
  {
    local_box_offset[0] = 31;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 10] [0] Offset and Count 30 20 -1 : 10 11 32
  if (rank == 10)
  {
    local_box_offset[0] = 31;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 11] [0] Offset and Count 40 -1 -1 : 10 11 32
  if (rank == 11)
  {
    local_box_offset[0] = 41;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 12] [0] Offset and Count 40 10 -1 : 10 10 32
  if (rank == 12)
  {
    local_box_offset[0] = 41;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 13] [0] Offset and Count 40 20 -1 : 10 11 32
  if (rank == 13)
  {
    local_box_offset[0] = 41;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 14] [0] Offset and Count 50 -1 -1 : 10 11 32
  if (rank == 14)
  {
    local_box_offset[0] = 51;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 15] [0] Offset and Count 50 10 -1 : 10 10 32
  if (rank == 15)
  {
    local_box_offset[0] = 51;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 16] [0] Offset and Count 50 20 -1 : 10 10 32
  if (rank == 16)
  {
    local_box_offset[0] = 51;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 10;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 17] [0] Offset and Count 49 30 -1 : 11 11 32
  if (rank == 17)
  {
    local_box_offset[0] = 50;
    local_box_offset[1] = 31;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 18] [0] Offset and Count 60 -1 -1 : 11 11 32
  if (rank == 18)
  {
    local_box_offset[0] = 61;
    local_box_offset[1] = 0;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }

  //[PIDX 19] [0] Offset and Count 60 10 -1 : 11 10 32
  if (rank == 19)
  {
    local_box_offset[0] = 61;
    local_box_offset[1] = 11;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 20] [0] Offset and Count 60 20 -1 : 11 10 32
  if (rank == 20)
  {
    local_box_offset[0] = 61;
    local_box_offset[1] = 21;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 21] [0] Offset and Count 60 30 -1 : 11 10 32
  if (rank == 21)
  {
    local_box_offset[0] = 61;
    local_box_offset[1] = 31;
    local_box_offset[2] = 0;

    local_box_size[0] = 11;
    local_box_size[1] = 10;
    local_box_size[2] = 32;
  }

  //[PIDX 22] [0] Offset and Count 59 40 -1 : 12 11 32
  if (rank == 22)
  {
    local_box_offset[0] = 60;
    local_box_offset[1] = 41;
    local_box_offset[2] = 0;

    local_box_size[0] = 12;
    local_box_size[1] = 11;
    local_box_size[2] = 32;
  }
#endif

}

//(11*17*32)+(18*32*32)+(18*32*32)+(18*32*32)+(18*32*32)+(18*32*32)+(10*42*32)+(21*52*32)
static void destroy_synthetic_simulation_data()
{
  switch (roi_type)
  {
    case 0:
      free(data);
      break;

    case 1:
      if (rank % 2 == 0)
        free(data);
      break;

    case 2:
      if ((rank / (global_box_size[1]/local_box_size[1])) % 2 == 0)
      {
        if (rank % 2 == 0)
          free(data);
      }
      else
      {
        if (rank % 2 == 1)
          free(data);
      }
      break;
  case 4:
    if (rank == current_time_step * (global_box_size[X] / local_box_size[X]) * (global_box_size[Y] / local_box_size[Y]))
      free (data);

    default:
      break;
  }
}

///< Parse the input argumencurrent_time_step
static void parse_args(int argc, char **argv)
{
  char flags[] = "g:l:f:r:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
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

    case('r'): // ROI type
      if (sscanf(optarg, "%d", &roi_type) < 0)
        terminate_with_error_msg("Invalid variable file\n%s", usage);
      break;

    default:
      terminate_with_error_msg("Wrong argumencurrent_time_step\n%s", usage);
    }
  }
}


static void create_synthetic_simulation_data_and_ROI_extents()
{
  switch (roi_type)
  {
    case 0:
      data = malloc(sizeof (float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
      //for (i = 0; i < local_box_size[0] * local_box_size[1] * local_box_size[2]; i++)
      //  data[i] = rank;

      int k,j,i;
      for (k = 0; k < local_box_size[2]; k++)
        for (j = 0; j < local_box_size[1]; j++)
          for (i = 0; i < local_box_size[0]; i++)
          {
            unsigned long long index = (unsigned long long) (local_box_size[0] * local_box_size[1] * k) + (local_box_size[0] * j) + i;
            data[index] = 100 + ((global_box_size[0] * global_box_size[1]*(local_box_offset[2] + k))+(global_box_size[0]*(local_box_offset[1] + j)) + (local_box_offset[0] + i));
          }

      PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
      break;

    case 1:
      if (rank % 2 == 0)
      {
        data = malloc(sizeof (float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
        for (i = 0; i < local_box_size[0] * local_box_size[1] * local_box_size[2]; i++)
          data[i] = 100;

        PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
      }
      else
      {
        data = NULL;
        PIDX_set_point_5D(local_size, 0, 0, 0, 1, 1);
      }
      break;

    case 2:
      if ((rank / (global_box_size[1]/local_box_size[1])) % 2 == 0)
      {
        if (rank % 2 == 0)
        {
          data = malloc(sizeof (float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
          for (i = 0; i < local_box_size[0] * local_box_size[1] * local_box_size[2]; i++)
            data[i] = 100;
          PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
        }
        else
        {
          data = NULL;
            PIDX_set_point_5D(local_size, 0, 0, 0, 1, 1);
        }
      }
      else
      {
        if (rank % 2 == 1)
        {
          data = malloc(sizeof (float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
          for (i = 0; i < local_box_size[0] * local_box_size[1] * local_box_size[2]; i++)
            data[i] = 100;

          PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
        }
        else
        {
          data = NULL;
          PIDX_set_point_5D(local_size, 0, 0, 0, 1, 1);
        }
      }
      break;

    case 4:
      if (rank == current_time_step * (global_box_size[X] / local_box_size[X]) * (global_box_size[Y] / local_box_size[Y]))
      {
        data = malloc(sizeof (float) * local_box_size[0] * local_box_size[1] * local_box_size[2]);
        for (i = 0; i < local_box_size[0] * local_box_size[1] * local_box_size[2]; i++)
          data[i] = 100;

        PIDX_set_point_5D(local_size, local_box_size[0], local_box_size[1], local_box_size[2], 1, 1);
      }
      else
      {
        data = NULL;
        PIDX_set_point_5D(local_size, 0, 0, 0, 1, 1);
      }

      break;

    default:
      break;
  }
}

int main(int argc, char **argv)
{
  init_mpi(argc, argv);
  parse_args(argc, argv);
  calculate_per_process_offsecurrent_time_step();

  rank_0_print("Simulation Data Created\n");

  int ret;
  PIDX_file file;            // IDX file descriptor
  PIDX_variable variable;   // variable descriptor

  PIDX_set_point_5D(global_size, global_box_size[0], global_box_size[1], global_box_size[2], 1, 1);
  PIDX_set_point_5D(local_offset, local_box_offset[0], local_box_offset[1], local_box_offset[2], 0, 0);

  //  Creating access
  PIDX_access access;
  PIDX_create_access(&access);
#if PIDX_HAVE_MPI
  PIDX_set_mpi_access(access, MPI_COMM_WORLD);
#endif

  if (roi_type == 4)
    time_step_count = global_box_size[Z]/local_box_size[Z];

  for (current_time_step = 0; current_time_step < time_step_count; current_time_step++)
  {
    //  PIDX mandatory calls
    create_synthetic_simulation_data_and_ROI_extents();

    ret = PIDX_file_create(output_file_name, PIDX_MODE_CREATE, access, global_size, &file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

    PIDX_set_block_size(file, 16);

    //PIDX_debug_rst(file, 1);
    //PIDX_set_compression_type(file, PIDX_CHUNKING_ONLY);
    PIDX_set_compression_type(file, PIDX_CHUNKING_ZFP);
    PIDX_set_lossy_compression_bit_rate(file, 1);

    ret = PIDX_set_current_time_step(file, current_time_step);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");

    ret = PIDX_set_variable_count(file, 1);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

    //ret = PIDX_set_resolution(file, 0, 2);
    //if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_resolution");

    //ret = PIDX_set_ROI_writes(file);
    //if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_ROI_writes");

    //unsigned long long restructured_box_size[5] = {128, 128, 128, 1, 1};
    //ret = PIDX_set_restructuring_box(file, restructured_box_size);
    //if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_restructuring_box");

    ret = PIDX_variable_create("ROI_Var", sizeof(float) * 8, FLOAT32, &variable);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

    ret = PIDX_variable_write_data_layout(variable, local_offset, local_size, data, PIDX_row_major);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_data_layout");

    ret = PIDX_append_and_write_variable(file, variable);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

    ret = PIDX_close(file);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");

    destroy_synthetic_simulation_data();
  }

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

  shutdown_mpi();

  return 0;
}
