#include "PIDX_io.h"


struct PIDX_io_descriptor
{

#if PIDX_HAVE_MPI
  MPI_Comm comm;                               ///< MPI sub-communicator (including all processes per IDX file)
#endif

  PIDX_partition_merge_idx_io partition_merge_idx_io;
  PIDX_partitioned_idx_io partitioned_idx_io;
  PIDX_idx_io idx_io;
  PIDX_raw_io raw_io;


  idx_dataset idx;                             ///< Contains all relevant IDX file info
                                               ///< Blocks per file, samples per block, bitmask, box, file name template

  idx_dataset_derived_metadata idx_d;          ///< Contains all derieved IDX file info
                                               ///< number of files, files that are ging to be populated

  idx_debug idx_dbg;

};

/// Returns elapsed time
double PIDX_get_time()
{
#if PIDX_HAVE_MPI
  return MPI_Wtime();
#else
  struct timeval temp;
  gettimeofday(&temp, NULL);
  return (double)(temp.tv_sec) + (double)(temp.tv_usec)/1000000.0;
#endif
}

void PIDX_init_timming_buffers1(PIDX_time time, int variable_count)
{
  variable_count = variable_count * 2;
  time->init_start = malloc (sizeof(double) * variable_count);            memset(time->init_start, 0, sizeof(double) * variable_count);
  time->init_end = malloc (sizeof(double) * variable_count);              memset(time->init_end, 0, sizeof(double) * variable_count);
  time->write_init_start = malloc (sizeof(double) * variable_count);      memset(time->write_init_start, 0, sizeof(double) * variable_count);
  time->write_init_end = malloc (sizeof(double) * variable_count);        memset(time->write_init_end, 0, sizeof(double) * variable_count);
}



void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count)
{
  variable_count = variable_count * 2;
  time->header_counter = 0;
  time->variable_counter = 0;
  time->startup_start = malloc (sizeof(double) * variable_count);                 memset(time->startup_start, 0, sizeof(double) * variable_count);
  time->startup_end = malloc (sizeof(double) * variable_count);                   memset(time->startup_end, 0, sizeof(double) * variable_count);

  time->rst_start = malloc (sizeof(double) * variable_count);                     memset(time->rst_start, 0, sizeof(double) * variable_count);
  time->rst_end = malloc (sizeof(double) * variable_count);                       memset(time->rst_end, 0, sizeof(double) * variable_count);

  time->rst_io_start = malloc (sizeof(double) * variable_count);                  memset(time->rst_io_start, 0, sizeof(double) * variable_count);
  time->rst_io_end = malloc (sizeof(double) * variable_count);                    memset(time->rst_io_end, 0, sizeof(double) * variable_count);

  time->hz_start = malloc (sizeof(double) * variable_count);                      memset(time->hz_start, 0, sizeof(double) * variable_count);
  time->hz_end = malloc (sizeof(double) * variable_count);                        memset(time->hz_end, 0, sizeof(double) * variable_count);

  time->cleanup_start = malloc (sizeof(double) * variable_count);                 memset(time->cleanup_start, 0, sizeof(double) * variable_count);
  time->cleanup_end = malloc (sizeof(double) * variable_count);                   memset(time->cleanup_end, 0, sizeof(double) * variable_count);
  time->finalize_start = malloc (sizeof(double) * variable_count);                memset(time->finalize_start, 0, sizeof(double) * variable_count);
  time->finalize_end = malloc (sizeof(double) * variable_count);                  memset(time->finalize_end, 0, sizeof(double) * variable_count);

  time->buffer_start = malloc (sizeof(double) * variable_count);                  memset(time->buffer_start, 0, sizeof(double) * variable_count);
  time->buffer_end = malloc (sizeof(double) * variable_count);                    memset(time->buffer_end, 0, sizeof(double) * variable_count);

  time->chunk_start =  malloc (sizeof(double) * variable_count);                  memset(time->chunk_start, 0, sizeof(double) * variable_count);
  time->chunk_end =  malloc (sizeof(double) * variable_count);                    memset(time->chunk_end, 0, sizeof(double) * variable_count);
  time->compression_start =  malloc (sizeof(double) * variable_count);            memset(time->compression_start, 0, sizeof(double) * variable_count);
  time->compression_end =  malloc (sizeof(double) * variable_count);              memset(time->compression_end, 0, sizeof(double) * variable_count);

  time->agg_start = malloc (sizeof(double*) * variable_count);       memset(time->agg_start, 0, sizeof(double*) * variable_count);
  time->agg_end = malloc (sizeof(double*) * variable_count);         memset(time->agg_end, 0, sizeof(double*) * variable_count);
  time->agg_buf_start = malloc (sizeof(double*) * variable_count);   memset(time->agg_buf_start, 0, sizeof(double*) * variable_count);
  time->agg_buf_end = malloc (sizeof(double*) * variable_count);     memset(time->agg_buf_end, 0, sizeof(double*) * variable_count);
  time->io_start = malloc (sizeof(double*) * variable_count);        memset(time->io_start, 0, sizeof(double*) * variable_count);
  time->io_end = malloc (sizeof(double*) * variable_count);          memset(time->io_end, 0, sizeof(double*) * variable_count);
  time->io_per_process_start = malloc (sizeof(double*) * variable_count);          memset(time->io_per_process_start, 0, sizeof(double*) * variable_count);
  time->io_per_process_end = malloc (sizeof(double*) * variable_count);          memset(time->io_per_process_end, 0, sizeof(double*) * variable_count);

  int i = 0;
  for (i = 0; i < variable_count; i++)
  {
    time->agg_start[i] = malloc (sizeof(double) * layout_count);       memset(time->agg_start[i], 0, sizeof(double) * layout_count);
    time->agg_end[i] = malloc (sizeof(double) * layout_count);         memset(time->agg_end[i], 0, sizeof(double) * layout_count);
    time->agg_buf_start[i] = malloc (sizeof(double) * layout_count);   memset(time->agg_buf_start[i], 0, sizeof(double) * layout_count);
    time->agg_buf_end[i] = malloc (sizeof(double) * layout_count);     memset(time->agg_buf_end[i], 0, sizeof(double) * layout_count);
    time->io_start[i] = malloc (sizeof(double) * layout_count);        memset(time->io_start[i], 0, sizeof(double) * layout_count);
    time->io_end[i] = malloc (sizeof(double) * layout_count);          memset(time->io_end[i], 0, sizeof(double) * layout_count);

    time->io_per_process_start[i] = malloc (sizeof(double) * layout_count);        memset(time->io_per_process_start[i], 0, sizeof(double) * layout_count);
    time->io_per_process_end[i] = malloc (sizeof(double) * layout_count);          memset(time->io_per_process_end[i], 0, sizeof(double) * layout_count);
  }
}


void PIDX_delete_timming_buffers1(PIDX_time time)
{
  free(time->init_start);
  free(time->init_end);
  free(time->write_init_start);
  free(time->write_init_end);
}



void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count)
{

  free(time->startup_start);
  free(time->startup_end);

  free(time->rst_start);
  free(time->rst_end);

  free(time->rst_io_end);
  free(time->rst_io_start);

  free(time->hz_start);
  free(time->hz_end);

  free(time->cleanup_start);
  free(time->cleanup_end);
  free(time->finalize_start);
  free(time->finalize_end);

  free(time->buffer_start);
  free(time->buffer_end);

  free(time->chunk_start);
  free(time->chunk_end);
  free(time->compression_start);
  free(time->compression_end);

  int i = 0;
  variable_count = variable_count * 2;
  for (i = 0; i < variable_count; i++)
  {
    free(time->agg_start[i]);
    free(time->agg_end[i]);
    free(time->agg_buf_start[i]);
    free(time->agg_buf_end[i]);
    free(time->io_start[i]);
    free(time->io_end[i]);

    free(time->io_per_process_start[i]);
    free(time->io_per_process_end[i]);
  }
  free(time->agg_start);
  free(time->agg_end);
  free(time->agg_buf_start);
  free(time->agg_buf_end);
  free(time->io_start);
  free(time->io_end);
  free(time->io_per_process_start);
  free(time->io_per_process_end);

}



void PIDX_print_raw_io_timing(MPI_Comm comm, PIDX_time time, int var_count, int layout_count)
{

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  int var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
#else
  total_time = max_time;
#endif

  if (max_time == total_time)
  {
    fprintf(stdout, "Time Taken: %f Seconds\n", max_time);
    fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
    printf("Block layout creation time %f\n", time->populate_idx_end_time - time->populate_idx_start_time);
    printf("Meta data IO time %f\n", time->meta_data_end_io - time->meta_data_start_io);

    double header_io_time = 0;
    for (var = 0; var < time->header_counter; var++)
    {
      header_io_time = header_io_time + (time->write_init_end[var] - time->write_init_start[var]);
      fprintf(stdout, "File Create time (+ header IO) %f\n", (time->write_init_end[var] - time->write_init_start[var]));
    }

    double total_time_rch = 0;
    for (var = 0; var < var_count; var++)
    {
      fprintf(stdout, "[%d] RST + RST IO + FINALIZE = %f + %f + %f = %f\n", var, (time->rst_end[var] - time->rst_start[var]), (time->rst_io_end[var] - time->rst_io_start[var]), (time->finalize_end[var] - time->finalize_start[var]), ((time->rst_end[var] - time->rst_start[var]) + (time->rst_io_end[var] - time->rst_io_start[var]) + (time->finalize_end[var] - time->finalize_start[var])) );
      total_time_rch = total_time_rch + ((time->rst_end[var] - time->rst_start[var]) + (time->rst_io_end[var] - time->rst_io_start[var]) + (time->finalize_end[var] - time->finalize_start[var]));
    }

    fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f] [%f]\n", (time->populate_idx_end_time - time->populate_idx_start_time) + (time->meta_data_end_io - time->meta_data_start_io) + total_time_rch + header_io_time, (time->populate_idx_end_time - time->populate_idx_start_time), (time->meta_data_end_io - time->meta_data_start_io), header_io_time, total_time_rch, max_time);

    fprintf(stdout, "==========================================================================================================\n");
  }
}

void PIDX_print_idx_io_timing(MPI_Comm comm, PIDX_time time, int var_count, int layout_count)
{
  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  int var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
#else
  total_time = max_time;
#endif
  if (max_time == total_time)
  {
    fprintf(stdout, "Time Taken: %f Seconds\n", max_time);
    fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
    printf("Block layout creation time %f\n", time->populate_idx_end_time - time->populate_idx_start_time);
    fprintf(stdout, "File Create Time: %f Seconds\n", (time->file_create_time - time->sim_start));

    double header_io_time = 0;
    for (var = 0; var < time->header_counter; var++)
    {
      header_io_time = header_io_time + (time->write_init_end[var] - time->write_init_start[var]);
      fprintf(stdout, "File Create time (+ header IO) %f\n", (time->write_init_end[var] - time->write_init_start[var]));
    }
    double total_time_ai = 0, total_time_bc = 0, total_time_a = 0, total_time_i = 0, total_time_pi = 0;
    int p = 0;
    for (var = 0; var < var_count; var++)
    {
      for (p = 0; p < layout_count; p++)
      {
        fprintf(stdout, "[%d %d] Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time = %f + %f + %f + %f = %f\n", var, p, (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]), (time->agg_end[var][p] - time->agg_start[var][p]), (time->io_end[var][p] - time->io_start[var][p]), (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]), (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]) + (time->agg_end[var][p] - time->agg_start[var][p]) + (time->io_end[var][p] - time->io_start[var][p]) + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]));

        total_time_bc = total_time_bc + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]);
        total_time_a = total_time_a + (time->agg_end[var][p] - time->agg_start[var][p]);
        total_time_i = total_time_i + (time->io_end[var][p] - time->io_start[var][p]);
        total_time_pi = total_time_pi + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]);
      }
    }
    total_time_ai = total_time_bc + total_time_a + total_time_i + total_time_pi;
    fprintf(stdout, "Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time : %f + %f + %f + %f = %f\n", total_time_bc, total_time_a, total_time_i, total_time_pi, total_time_ai);



    //printf("static_var_counter = %d\n", static_var_counter);
    double total_time_rch = 0;
    for (var = 0; var < var_count; var++)
    {
      fprintf(stdout, "[%d] STARTUP + RST + BRST + HZ = %f + %f + %f + %f = %f\n", var, (time->startup_end[var] - time->startup_start[var]), (time->rst_end[var] - time->rst_start[var]), (time->chunk_end[var] - time->chunk_start[var]), (time->hz_end[var] - time->hz_start[var]), (time->startup_end[var] - time->startup_start[var]) + (time->rst_end[var] - time->rst_start[var]) + (time->chunk_end[var] - time->chunk_start[var]) + (time->hz_end[var] - time->hz_start[var]));
      total_time_rch = total_time_rch + (time->startup_end[var] - time->startup_start[var]) + (time->rst_end[var] - time->rst_start[var]) + (time->chunk_end[var] - time->chunk_start[var]) + (time->hz_end[var] - time->hz_start[var]);
    }

    fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f + %f] [%f]\n", total_time_ai + total_time_rch + (time->file_create_time - time->sim_start) + (time->populate_idx_end_time - time->populate_idx_start_time) + header_io_time, (time->populate_idx_end_time - time->populate_idx_start_time), (time->file_create_time - time->sim_start), header_io_time, total_time_rch, total_time_ai, max_time);

    fprintf(stdout, "==========================================================================================================\n");
  }
}


void PIDX_print_partition_timing(MPI_Comm comm, PIDX_time time, int var_count, int layout_count)
{

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  int var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
#else
  total_time = max_time;
#endif

  if (max_time == total_time)
  {
    fprintf(stdout, "Time Taken: %f Seconds\n", max_time);
    fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
    printf("Partition time %f\n", time->populate_idx_end_time - time->populate_idx_start_time);
    fprintf(stdout, "File Create Time: %f Seconds\n", (time->file_create_time - time->sim_start));

    double header_io_time = 0;
    for (var = 0; var < time->header_counter; var++)
    {
      header_io_time = header_io_time + (time->write_init_end[var] - time->write_init_start[var]);
      fprintf(stdout, "File Create time (+ header IO) %f\n", (time->write_init_end[var] - time->write_init_start[var]));
    }
    double total_time_ai = 0, total_time_bc = 0, total_time_a = 0, total_time_i = 0, total_time_pi = 0;
    int p = 0;
    for (var = 0; var <  var_count * 2; var++)
    {
      for (p = 0; p < layout_count; p++)
      {
        fprintf(stdout, "[%d %d] Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time = %f + %f + %f + %f = %f\n", var, p, (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]), (time->agg_end[var][p] - time->agg_start[var][p]), (time->io_end[var][p] - time->io_start[var][p]), (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]), (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]) + (time->agg_end[var][p] - time->agg_start[var][p]) + (time->io_end[var][p] - time->io_start[var][p]) + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]));

        total_time_bc = total_time_bc + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]);
        total_time_a = total_time_a + (time->agg_end[var][p] - time->agg_start[var][p]);
        total_time_i = total_time_i + (time->io_end[var][p] - time->io_start[var][p]);
        total_time_pi = total_time_pi + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]);
      }
    }
    total_time_ai = total_time_bc + total_time_a + total_time_i + total_time_pi;
    fprintf(stdout, "Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time : %f + %f + %f + %f = %f\n", total_time_bc, total_time_a, total_time_i, total_time_pi, total_time_ai);


    double total_time_rch = 0;
    for (var = 0; var < var_count * 2; var++)
    {
      fprintf(stdout, "[%d] STARTUP + HZ + FINALIZE = %f + %f + %f = %f\n", var, (time->startup_end[var] - time->startup_start[var]), (time->hz_end[var] - time->hz_start[var]), (time->finalize_end[var] - time->finalize_start[var]), (time->startup_end[var] - time->startup_start[var]) + (time->hz_end[var] - time->hz_start[var]));
      total_time_rch = total_time_rch + (time->startup_end[var] - time->startup_start[var]) + (time->finalize_end[var] - time->finalize_start[var]) + (time->hz_end[var] - time->hz_start[var]);
    }

    fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f + %f] [%f]\n", total_time_ai + total_time_rch + (time->file_create_time - time->sim_start) + (time->populate_idx_end_time - time->populate_idx_start_time) + header_io_time, (time->populate_idx_end_time - time->populate_idx_start_time), (time->file_create_time - time->sim_start), header_io_time, total_time_rch, total_time_ai, max_time);

    fprintf(stdout, "==========================================================================================================\n");
  }
}


void PIDX_print_partition_merge_timing(MPI_Comm comm, PIDX_time time, int var_count, int layout_count)
{

  double total_time = time->sim_end - time->sim_start;
  double max_time = total_time;
  int var = 0, rank = 0, nprocs = 1;

#if PIDX_HAVE_MPI
  MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
#else
  total_time = max_time;
#endif

  if (max_time == total_time)
  {
    fprintf(stdout, "Time Taken: %f Seconds\n", max_time);
    fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");
    printf("Block layout creation time %f\n", time->populate_idx_end_time - time->populate_idx_start_time);
    printf("Partiton time %f\n", time->partition_end_time - time->partition_start_time);
    fprintf(stdout, "File Create Time: %f Seconds\n", (time->file_create_time - time->sim_start));

    double header_io_time = 0;
    for (var = 0; var < time->header_counter; var++)
    {
      header_io_time = header_io_time + (time->write_init_end[var] - time->write_init_start[var]);
      fprintf(stdout, "File Create time (+ header IO) %f\n", (time->write_init_end[var] - time->write_init_start[var]));
    }
    double total_time_ai = 0, total_time_bc = 0, total_time_a = 0, total_time_i = 0, total_time_pi = 0;
    int p = 0;
    for (var = 0; var <  var_count * 2; var++)
    {
      for (p = 0; p < layout_count; p++)
      {
        fprintf(stdout, "[%d %d] Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time = %f + %f + %f + %f = %f\n", var, p, (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]), (time->agg_end[var][p] - time->agg_start[var][p]), (time->io_end[var][p] - time->io_start[var][p]), (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]), (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]) + (time->agg_end[var][p] - time->agg_start[var][p]) + (time->io_end[var][p] - time->io_start[var][p]) + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]));

        total_time_bc = total_time_bc + (time->agg_buf_end[var][p] - time->agg_buf_start[var][p]);
        total_time_a = total_time_a + (time->agg_end[var][p] - time->agg_start[var][p]);
        total_time_i = total_time_i + (time->io_end[var][p] - time->io_start[var][p]);
        total_time_pi = total_time_pi + (time->io_per_process_end[var][p] - time->io_per_process_start[var][p]);
      }
    }
    total_time_ai = total_time_bc + total_time_a + total_time_i + total_time_pi;
    fprintf(stdout, "Agg Buf Time + Agg time + AGG I/O time + Per-Process I/O time : %f + %f + %f + %f = %f\n", total_time_bc, total_time_a, total_time_i, total_time_pi, total_time_ai);


    double total_time_rch = 0;
    for (var = 0; var < var_count * 2; var++)
    {
      fprintf(stdout, "[%d] STARTUP + HZ + FINALIZE = %f + %f + %f = %f\n", var, (time->startup_end[var] - time->startup_start[var]), (time->hz_end[var] - time->hz_start[var]), (time->finalize_end[var] - time->finalize_start[var]), (time->startup_end[var] - time->startup_start[var]) + (time->hz_end[var] - time->hz_start[var]));
      total_time_rch = total_time_rch + (time->startup_end[var] - time->startup_start[var]) + (time->finalize_end[var] - time->finalize_start[var]) + (time->hz_end[var] - time->hz_start[var]);
    }

    fprintf(stdout, "PIDX Total Time = %f [%f + %f + %f + %f + %f + %f] [%f]\n", total_time_ai + total_time_rch + (time->file_create_time - time->sim_start) + (time->populate_idx_end_time - time->populate_idx_start_time) + (time->partition_end_time - time->partition_start_time) + header_io_time, (time->populate_idx_end_time - time->populate_idx_start_time), (time->partition_end_time - time->partition_start_time), (time->file_create_time - time->sim_start), header_io_time, total_time_rch, total_time_ai, max_time);

    fprintf(stdout, "==========================================================================================================\n");
  }
}

PIDX_io PIDX_io_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, idx_debug idx_dbg)
{
  //Creating the restructuring ID
  PIDX_io io;
  io = malloc(sizeof (*io));
  memset(io, 0, sizeof (*io));

  io->idx = idx_meta_data;
  io->idx_d = idx_derived_ptr;
  io->idx_dbg = idx_dbg;

  return (io);
}

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_io_set_communicator(PIDX_io file, MPI_Comm comm)
{
  if (file == NULL)
    return PIDX_err_id;

  file->comm = comm;

  return PIDX_success;
}
#endif

PIDX_return_code PIDX_io_io(PIDX_io file, int mode, int io_type, int start_var_index, int end_var_index)
{
  int ret;
  if (mode == PIDX_MODE_CREATE)
  {
    if (io_type == PIDX_PARTITIONED_IDX_IO)
    {
      file->partitioned_idx_io = PIDX_partitioned_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->partitioned_idx_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_partitioned_idx_io_set_communicator(file->partitioned_idx_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_partitioned_idx_write(file->partitioned_idx_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_partitioned_idx_io_finalize(file->partitioned_idx_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }

    else if (io_type == PIDX_IDX_IO)
    {
      if (start_var_index == file->idx->variable_count)
        return PIDX_success;

      file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->idx_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_idx_write(file->idx_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_idx_io_finalize(file->idx_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }

    else if (io_type == PIDX_RAW_IO)
    {
      if (start_var_index == file->idx->variable_count)
        return PIDX_success;

      file->raw_io = PIDX_raw_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->raw_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_raw_io_set_communicator(file->raw_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_raw_write(file->raw_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_raw_io_finalize(file->raw_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
    else
    {
      file->partition_merge_idx_io = PIDX_partition_merge_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->partition_merge_idx_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_partition_merge_idx_io_set_communicator(file->partition_merge_idx_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_partition_merge_idx_write(file->partition_merge_idx_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_partition_merge_idx_io_finalize(file->partition_merge_idx_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
  }

  else if (mode == PIDX_MODE_RDONLY)
  {
    if (io_type == PIDX_IDX_IO)
    {
      if (start_var_index == file->idx->variable_count)
        return PIDX_success;

      file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->idx_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_idx_read(file->idx_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_idx_io_finalize(file->idx_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
    else
    {
      if (start_var_index == file->idx->variable_count)
        return PIDX_success;

      file->raw_io = PIDX_raw_io_init(file->idx, file->idx_d, file->idx_dbg);
      if (file->raw_io == NULL)
        return PIDX_err_flush;

      ret = PIDX_raw_io_set_communicator(file->raw_io, file->comm);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      int ncores = 1;
#if PIDX_HAVE_MPI
      MPI_Comm_size(file->comm, &ncores);
#endif
      if (file->idx_d->data_core_count == ncores)
        ret = PIDX_raw_read(file->raw_io, start_var_index, end_var_index);
      else
        ret = PIDX_forced_raw_read(file->raw_io, start_var_index, end_var_index);
      if (ret != PIDX_success)
        return PIDX_err_flush;

      ret = PIDX_raw_io_finalize(file->raw_io);
      if (ret != PIDX_success)
        return PIDX_err_flush;
    }
  }

  /*
  else if (mode == PIDX_MODE_RDWR)
  {
    int state = file->idx->variable[file->local_variable_index]->io_state;
    int state_index = file->local_variable_index;
    int new_state, same_state_count = 0;

    for (i = file->local_variable_index; i < file->local_variable_index + file->local_variable_count; i++)
    {
      new_state = file->idx->variable[i]->io_state;
      if (state == new_state)
      {
        same_state_count++;
        if (i == file->local_variable_index + file->local_variable_count - 1)
        {
          if (state == 1)
          {
            if (file->local_variable_index == file->idx->variable_count)
              return PIDX_success;

            file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
            if (file->idx_io == NULL)
              return PIDX_err_flush;

            ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
            if (ret != PIDX_success)
              return PIDX_err_flush;

            PIDX_idx_write(file->idx_io, state_index, state_index + same_state_count);
          }
          else if (state == 0)
          {
            file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
            if (file->idx_io == NULL)
              return PIDX_err_flush;

            ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
            if (ret != PIDX_success)
              return PIDX_err_flush;

            PIDX_idx_read(file->idx_io, state_index, state_index + same_state_count);
          }
        }
      }
      else
      {
        if (state == 1)
        {
          if (file->local_variable_index == file->idx->variable_count)
            return PIDX_success;

          PIDX_idx_write(file->idx_io, state_index, state_index + same_state_count);

          file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
          if (file->idx_io == NULL)
            return PIDX_err_flush;

          ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
          if (ret != PIDX_success)
            return PIDX_err_flush;
        }
        else if (state == 0)
        {
          PIDX_idx_read(file->idx_io, state_index, state_index + same_state_count);
          file->idx_io = PIDX_idx_io_init(file->idx, file->idx_d, file->idx_dbg);
          if (file->idx_io == NULL)
            return PIDX_err_flush;

          ret = PIDX_idx_io_set_communicator(file->idx_io, file->comm);
          if (ret != PIDX_success)
            return PIDX_err_flush;
        }

        state = new_state;
        state_index = i;
        same_state_count = 1;
      }
    }
  }
  */
  return PIDX_success;
}

PIDX_return_code PIDX_io_finalize(PIDX_io io)
{
  free(io);
  io = 0;

  return PIDX_success;
}
