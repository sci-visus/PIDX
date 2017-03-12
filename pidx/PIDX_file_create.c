#include "PIDX_file_handler.h"


/// Function to create IDX file descriptor (based on flags and access)
PIDX_return_code PIDX_file_create(const char* filename, PIDX_flags flags, PIDX_access access_type, PIDX_point dims, PIDX_file* file)
{
  if (flags != PIDX_MODE_CREATE && flags != PIDX_MODE_EXCL)
    return PIDX_err_unsupported_flags;

  if (flags == PIDX_MODE_EXCL)
  {
    struct stat buffer;
    if (stat(filename, &buffer) != 0)
      return PIDX_err_file_exists;
  }

  unsigned long long i = 0;
  int ret;
  char file_name_skeleton[1024];

  if (strncmp(".idx", &filename[strlen(filename) - 4], 4) != 0 && !filename)
    return PIDX_err_name;

  *file = malloc(sizeof (*(*file)));
  memset(*file, 0, sizeof (*(*file)));

  (*file)->idx = malloc(sizeof (*((*file)->idx)));
  memset((*file)->idx, 0, sizeof (*((*file)->idx)));

  (*file)->idx_c = malloc(sizeof (*((*file)->idx_c)));
  memset((*file)->idx_c, 0, sizeof (*((*file)->idx_c)));

  // IDX derived pointer (everything derived from idx)
  (*file)->idx_d = malloc(sizeof (*((*file)->idx_d)));
  memset((*file)->idx_d, 0, sizeof (*((*file)->idx_d)));

  (*file)->idx_d->time = malloc(sizeof (*((*file)->idx_d->time)));
  memset((*file)->idx_d->time, 0, sizeof (*((*file)->idx_d->time)));

  (*file)->idx_d->time->sim_start = PIDX_get_time();

  (*file)->idx_dbg = malloc(sizeof (*((*file)->idx_dbg)));
  memset((*file)->idx_dbg, 0, sizeof (*((*file)->idx_dbg)));

  (*file)->flags = flags;

  memcpy((*file)->idx->bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
  memcpy((*file)->idx->box_bounds, dims, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  for (i = 0; i < PIDX_MAX_DIMENSIONS; i++)
  {
    (*file)->idx_d->partition_count[i] = 1;
    (*file)->idx_d->partition_size[i] = getPowerOf2(dims[i]);
    (*file)->idx_d->partition_offset[i] = 0;
  }

  (*file)->idx_dbg->debug_do_rst = 1;
  (*file)->idx_dbg->debug_do_chunk = 1;
  (*file)->idx_dbg->debug_do_compress = 1;
  (*file)->idx_dbg->debug_do_hz = 1;
  (*file)->idx_dbg->debug_do_agg = 1;
  (*file)->idx_dbg->debug_do_io = 1;
  (*file)->idx_dbg->debug_rst = 0;
  (*file)->idx_dbg->debug_hz = 0;

  (*file)->local_group_index = 0;
  (*file)->local_group_count = 1;
  (*file)->flush_used = 0;
  (*file)->write_on_close = 0;

  (*file)->idx_d->color = 0;
  (*file)->idx->io_type = PIDX_IDX_IO;
  (*file)->idx_d->data_core_count = -1;

  (*file)->idx_d->reduced_res_from = 0;
  (*file)->idx_d->reduced_res_to = 0;

  (*file)->idx_d->parallel_mode = access_type->parallel;
  (*file)->idx_d->raw_io_pipe_length = 0;

  (*file)->idx_c->global_comm = access_type->comm;
  (*file)->idx_c->local_comm = access_type->comm;
  MPI_Comm_rank((*file)->idx_c->global_comm, &((*file)->idx_c->grank));
  MPI_Comm_size((*file)->idx_c->global_comm, &((*file)->idx_c->gnprocs));
  MPI_Comm_rank((*file)->idx_c->local_comm, &((*file)->idx_c->lrank));
  MPI_Comm_size((*file)->idx_c->local_comm, &((*file)->idx_c->lnprocs));

  (*file)->idx_c->gnproc_x = -1;
  (*file)->idx_c->gnproc_y = -1;
  (*file)->idx_c->gnproc_z = -1;


  (*file)->idx->enable_rst = 1;
  (*file)->idx->enable_agg = 1;

  (*file)->idx->current_time_step = 0;
  (*file)->idx->variable_count = -1;
  (*file)->idx->group_index_tracker = 0;

  (*file)->idx->agg_counter = 0;

  (*file)->idx->compression_type = PIDX_NO_COMPRESSION;

  strncpy(file_name_skeleton, filename, strlen(filename) - 4);
  file_name_skeleton[strlen(filename) - 4] = '\0';
  sprintf((*file)->idx->filename, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_global, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_partition, "%s.idx", file_name_skeleton);
  sprintf((*file)->idx->filename_file_zero, "%s.idx", file_name_skeleton);

  (*file)->idx->bits_per_block = PIDX_default_bits_per_block;
  (*file)->idx->blocks_per_file = PIDX_default_blocks_per_file;

  //initialize logic_to_physic transform to identity
  (*file)->idx->transform[0]  = 1.0;
  (*file)->idx->transform[5]  = 1.0;
  (*file)->idx->transform[10] = 1.0;
  (*file)->idx->transform[15] = 1.0;

  (*file)->idx->variable_group_count = 1;

  memset((*file)->idx->bitPattern, 0, 512);
  memset((*file)->idx->bitSequence, 0, 512);
  memset((*file)->idx->reg_patch_size, 0, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  (*file)->idx->reg_box_set = PIDX_BOX_PER_PROCESS;
  (*file)->idx->compression_factor = 1;
  (*file)->idx->compression_bit_rate = 64;
  for (i=0;i<PIDX_MAX_DIMENSIONS;i++)
    (*file)->idx->chunk_size[i] = 1;

  (*file)->idx_d->dimension = 0;
  (*file)->idx_d->samples_per_block = (int)pow(2, PIDX_default_bits_per_block);
  (*file)->idx_d->maxh = 0;
  (*file)->idx_d->max_file_count = 0;
  (*file)->idx_d->fs_block_size = 0;
  (*file)->idx_d->start_fs_block = 0;

  (*file)->idx_dbg->state_dump = PIDX_NO_META_DATA_DUMP;

  (*file)->idx_d->layout_count = 0;
  (*file)->idx_d->reduced_res_from = 0;
  (*file)->idx_d->reduced_res_to = 0;

  (*file)->idx->shared_face = 0;

  for (i = 0; i < 16; i++)
  {
    (*file)->idx->variable_grp[i] = malloc(sizeof(*((*file)->idx->variable_grp[i])));
    memset((*file)->idx->variable_grp[i], 0, sizeof(*((*file)->idx->variable_grp[i])));
  }

  if ((*file)->idx_c->grank == 0)
  {
    //TODO: close and delete the file (there is a way to do this automatically by fopen...)
    struct stat stat_buf;
    FILE *dummy = fopen(".dummy.txt", "w");
    fclose(dummy);
    ret = stat(".dummy.txt", &stat_buf);
    if (ret != 0)
    {
      fprintf(stderr, "[%s] [%d] Unable to identify File-System block size\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    (*file)->idx_d->fs_block_size = stat_buf.st_blksize;
  }

#if PIDX_HAVE_MPI
  if ((*file)->idx_d->parallel_mode == 1)
    MPI_Bcast(&((*file)->idx_d->fs_block_size), 1, MPI_INT, 0, (*file)->idx_c->global_comm);
#endif

  (*file)->idx->flip_endian = 0;

  unsigned int endian = 1;
  char *c = (char*)&endian;
  if (*c)
    (*file)->idx->endian = 1;
  else
    (*file)->idx->endian = 0;

  return PIDX_success;
}
