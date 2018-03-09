#include "../../../PIDX_inc.h"

static int maximum_neighbor_count = 256;
static int intersectNDChunk(PIDX_patch A, PIDX_patch B);
static int pointInChunk(PIDX_patch p, const double *pos);
static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi);
static PIDX_return_code PIDX_meta_data_write(PIDX_io file, int gi, int svi);
static PIDX_return_code PIDX_particle_raw_read(PIDX_io file, int gi, int svi, int evi);


PIDX_return_code PIDX_particle_write(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_time time = file->idx_d->time;

  if (group_meta_data_init(file, gi, svi, evi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  time->particle_meta_data_io_start = MPI_Wtime();
  if (PIDX_meta_data_write(file, gi, svi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->particle_meta_data_io_end = MPI_Wtime();


  time->particle_data_io_start = MPI_Wtime();
  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  PIDX_variable var0 = var_grp->variable[svi];
  for (int p = 0; p < var0->sim_patch_count; p++)
  {
    char file_name[PATH_MAX];
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));
    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, file->idx_c->grank, p);

    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    off_t data_offset = 0;
    for (int si = svi; si < evi; si++)
    {
      PIDX_variable var = var_grp->variable[si];

      int sample_count;
      int bits_per_sample;
      PIDX_get_datatype_details(var->type_name, &sample_count, &bits_per_sample);

      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

#if 0
      fprintf(stderr, "[Inside PIDX] Rank %d has %d particles [%lld %lld %lld - %lld %lld %lld]\n", file->idx_c->grank, var->sim_patch[p]->particle_count, var->sim_patch[p]->offset[0], var->sim_patch[p]->offset[1], var->sim_patch[p]->offset[2], var->sim_patch[p]->size[0], var->sim_patch[p]->size[1], var->sim_patch[p]->size[2]);
      if (si == 0)
      {
        double particle_x = 0, particle_y, particle_z = 0;
        for (int pt = 0; pt < var->sim_patch[p]->particle_count; pt++)
        {
          memcpy(&particle_x, var->sim_patch[p]->buffer + (pt * sample_count + 0) * sizeof(double), sizeof (double));
          memcpy(&particle_y, var->sim_patch[p]->buffer + (pt * sample_count + 1) * sizeof(double), sizeof (double));
          memcpy(&particle_z, var->sim_patch[p]->buffer + (pt * sample_count + 2) * sizeof(double), sizeof (double));
          printf("%f\t%f\t%f\n", particle_x, particle_y, particle_z);
        }
      }
      printf("[%s] [%d] = %d %d %d\n", var->type_name, var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT), var->sim_patch[p]->particle_count, sample_count, (bits_per_sample/CHAR_BIT));
#endif

      ssize_t buffer_size = var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT);
      ssize_t write_count = pwrite(fp, var->sim_patch[p]->buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
      data_offset = data_offset + buffer_size;
    }

    close(fp);
  }

  free (directory_path);
  time->particle_data_io_end = MPI_Wtime();

  return PIDX_success;
}



PIDX_return_code PIDX_particle_read(PIDX_io file, int gi, int svi, int evi)
{

  PIDX_particle_raw_read(file, gi, svi, evi);

  return PIDX_success;
}



static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi)
{
  PIDX_time time = file->idx_d->time;

  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  int ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();

  return PIDX_success;
}


/* Will:
 * The metadata file structure is as follows:
 * Header:
 *   global number of procs: f64
 *   max patch count on a rank: f64
 *
 * List of patch data, with nprocs * max_count entries
 * Each entry consists of:
 *   local patch count: f64
 *   box.lower: 3*f64
 *   box.upper: 3*f64
 *   local particle count: f64
 *
 * TODO: Wouldn't this be a lot easier to construct and write if we had
 * some structs for the types? Then we also wouldn't need to store the
 * counts as f64 (they should be u64).
 */
static PIDX_return_code PIDX_meta_data_write(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  double *global_patch;
  PIDX_variable var0 = var_grp->variable[svi];
  int max_patch_count;
  int patch_count = var0->sim_patch_count;
  // TODO WILL: Why the max patch count across all ranks? Ranks could have differing numbers of patches
  // right? So would we leave some unused space in the file here?
  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);

  double *local_patch = malloc(sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));
  memset(local_patch, 0, sizeof(double) * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1));

  int pcounter = 0;
  local_patch[0] = (double)patch_count;
  for (int i = 0; i < patch_count; i++)
  {
    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * PIDX_MAX_DIMENSIONS + d + 1] = var0->sim_patch[i]->physical_offset[d];

    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
      local_patch[i * PIDX_MAX_DIMENSIONS + PIDX_MAX_DIMENSIONS + d + 1] = var0->sim_patch[i]->physical_size[d];

    local_patch[i * PIDX_MAX_DIMENSIONS + 2*PIDX_MAX_DIMENSIONS + 1] = var0->sim_patch[i]->particle_count;

    pcounter++;
  }

  global_patch = malloc((file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));
  memset(global_patch, 0,(file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double));

  MPI_Allgather(local_patch, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, global_patch + 2, (2 * PIDX_MAX_DIMENSIONS + 1) * max_patch_count + 1, MPI_DOUBLE, file->idx_c->global_comm);

  global_patch[0] = (double)file->idx_c->gnprocs;
  global_patch[1] = (double)max_patch_count;

  char file_path[PATH_MAX];

  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(file_path, "%s_OFFSET_SIZE", directory_path);
  free(directory_path);
  if (file->idx_c->grank == 1 || file->idx_c->gnprocs == 1)
  {
    int fp = open(file_path, O_CREAT | O_WRONLY, 0664);

    //for (int i = 0; i < (file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2); i++)
    //  printf("[%d] [np %d] ----> %f\n", i, file->idx_c->gnprocs, global_patch[i]);

    ssize_t write_count = pwrite(fp, global_patch, (file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double), 0);
    if (write_count != (file->idx_c->gnprocs * (max_patch_count * (2 * PIDX_MAX_DIMENSIONS + 1) + 1) + 2) * sizeof(double))
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
      return PIDX_err_io;
    }
    close(fp);
  }

  free(local_patch);
  free(global_patch);

  return PIDX_success;
}



static PIDX_return_code PIDX_particle_raw_read(PIDX_io file, int gi, int svi, int evi)
{
  int max_dim = 3;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  char size_path[PATH_MAX];

  char *idx_directory_path = malloc(sizeof(*idx_directory_path) * PATH_MAX);
  memset(idx_directory_path, 0, sizeof(*idx_directory_path) * PATH_MAX);
  strncpy(idx_directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  sprintf(size_path, "%s_OFFSET_SIZE", idx_directory_path);
  free(idx_directory_path);

  double number_cores = 0;
  int fp = open(size_path, O_RDONLY);
  // TODO WILL: This would be a lot easier to follow if we pread into a structure of some kind
  ssize_t read_count = pread(fp, &number_cores, sizeof(double), 0);
  if (read_count != sizeof(double))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  double max_patch_count = 0;
  // TODO WILL: This would be a lot easier to follow if we pread into a structure of some kind
  read_count = pread(fp, &max_patch_count, sizeof(double), sizeof(double));
  if (read_count != sizeof(double))
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  int buffer_read_size = (number_cores * ((int)max_patch_count * (2 * max_dim + 1) + 1)) * sizeof(double);

  double *size_buffer = malloc(buffer_read_size);
  memset(size_buffer, 0, buffer_read_size);

  // TODO WILL: This would be a lot easier to follow if we pread into a structure of some kind
  read_count = pread(fp, size_buffer, buffer_read_size, 2 * sizeof(double));
  if (read_count != buffer_read_size)
  {
    fprintf(stderr, "[%s] [%d] pread() failed.\n", __FILE__, __LINE__);
    return PIDX_err_io;
  }

  close(fp);

  char *file_name = malloc(PATH_MAX * sizeof(*file_name));
  memset(file_name, 0, PATH_MAX * sizeof(*file_name));

  char *directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  char *data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, file->idx->current_time_step);

  PIDX_buffer tmp_patch_read_buf = PIDX_buffer_create_empty();

  for (int pc1 = 0; pc1 < var_grp->variable[svi]->sim_patch_count; pc1++)
  {
    PIDX_patch local_proc_patch = (PIDX_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->physical_offset[d] = var_grp->variable[svi]->sim_patch[pc1]->physical_offset[d];
      local_proc_patch->physical_size[d] = var_grp->variable[svi]->sim_patch[pc1]->physical_size[d];
    }

    // PC - - - - - - - - - -   PC - - - - - - - - - -      PC - - - - -
    // 0  1 2 3 4 5 6 7 8 9 10   11 12 13 14 15 16 17 18 19 20 21  22
    PIDX_patch n_proc_patch = (PIDX_patch)malloc(sizeof (*n_proc_patch));
    memset(n_proc_patch, 0, sizeof (*n_proc_patch));
    int p_counter = 1;

    int patch_count = 0;
    for (int n = 0; n < number_cores; n++)
    {
      int pc = (int)size_buffer[n * ((int)max_patch_count * (2 * max_dim + 1) + 1)];
      int pc_index = n * ((int)max_patch_count * (2 * max_dim + 1) + 1);
      //if (file->idx_c->grank == 0)
      //  printf("Index %d PC %d\n", n * ((int)max_patch_count * (2 * max_dim + 1) + 1), pc);
      for (int m = 0; m < pc; m++)
      {
        for (int d = 0; d < PIDX_MAX_DIMENSIONS; d++)
        {
          n_proc_patch->physical_offset[d] = size_buffer[pc_index + m * (2 * max_dim + 1) + d + 1];
          n_proc_patch->physical_size[d] = size_buffer[pc_index + m * (2 * max_dim + 1) + max_dim + d + 1];
          n_proc_patch->particle_count = (int)size_buffer[pc_index + m * (2 * max_dim + 1) + max_dim + d + 1 + 1];
        }

        if (file->idx_c->grank == 0)
          printf("[PC %lu] OC %f %f %f - %f %f %f\n", n_proc_patch->particle_count, n_proc_patch->physical_offset[0],
              n_proc_patch->physical_offset[1], n_proc_patch->physical_offset[2], n_proc_patch->physical_size[0],
              n_proc_patch->physical_size[1], n_proc_patch->physical_size[2]);

        if (intersectNDChunk(local_proc_patch, n_proc_patch))
        {
          printf("Reading from n proc patch %d\n", n);
          sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, n, m);
          int fpx = open(file_name, O_RDONLY);

          // TODO WILL: For particles we need to rethink how we do this loop, we'll need to
          // check that the particle position is inside our query box, and only then get
          // the variable from the file. So there's a few options: we could store some mask
          // of which particles in the patch we're reading and compute this in a pre-pass
          // where we read (and optionally keep) the positions, and then read all the other
          // vars based on this mask. Or, we can go through for each particle and do the
          // test and read its data on a particle-by-particle basis (that is copying from the patch,
          // though this will then mean we load all vars the user wants for the patch up front,
          // and then do this loop through and copy over)
          // The first option may be easier to extend on to doing some acceleration structures
          // where we can then now that a whole subtree of particles are inside/outside the
          // query. We can also compute once up front the indices to write the particles
          // within the region to by doing a scan, then we can read the whole patch for
          // that attrib, and copy over the data.
          // The latter will be easier to implement, and since the patches are not large
          // it should be fine to do.
          for (int start_index = svi; start_index < evi; start_index = start_index + 1)
          {
            int other_offset = 0;
            for (int v1 = 0; v1 < start_index; v1++)
            {
              PIDX_variable var1 = var_grp->variable[v1];
              other_offset = other_offset + ((var1->bpv/8) * var1->vps * n_proc_patch->particle_count);
            }
            PIDX_variable var = var_grp->variable[start_index];
            const size_t bytes_per_sample = var->vps * var->bpv/8;

            const size_t proc_particle_read_size = n_proc_patch->particle_count * bytes_per_sample;
            PIDX_buffer_resize(&tmp_patch_read_buf, proc_particle_read_size);

            const size_t preadc = pread(fpx, tmp_patch_read_buf.buffer, proc_particle_read_size, other_offset);
            if (preadc != proc_particle_read_size)
            {
              fprintf(stderr, "[%s] [%d] Error in pread [%d %d]\n", __FILE__, __LINE__, (int)preadc,
                  (int)proc_particle_read_size);
              return PIDX_err_rst;
            }

            // TODO WILL: How do we know which variable the user has written as the "position" of
            // particles, so that we can filter them by the box query?
            // QUESTION for SID: the setup for the sim patch count info doesn't make much sense
            // it seems. I'm not sure how the patches are added? 
            // TODO FURTHERMORE: When we're reading and filtering the non-positional vars in
            // a box query we still need to know the positions of the particles to do the filtering
            // correctly. So regardless of what the user requests to read, we always must read the
            // positions
            size_t patch_particle_offset = *var->sim_patch[pc1]->read_particle_count * bytes_per_sample;

            // Always alloc space for at least half the particles in the patch to be stored
            // TODO: With better acceleration structures built on the particles within a patch,
            // we can make better estimates (or know exactly) how many particles this query will take
            // from this patch.
            if (var->sim_patch[pc1]->read_particle_buffer_capacity == 0
                || var->sim_patch[pc1]->read_particle_buffer_capacity - patch_particle_offset
                    < n_proc_patch->particle_count * bytes_per_sample / 2)
            {
              // Round up to a particle when allocating
              var->sim_patch[pc1]->read_particle_buffer_capacity += n_proc_patch->particle_count * bytes_per_sample / 2
                + bytes_per_sample;
              *var->sim_patch[pc1]->read_particle_buffer = realloc(*var->sim_patch[pc1]->read_particle_buffer,
                  var->sim_patch[pc1]->read_particle_buffer_capacity);
            }

            // TODO WILL: This assumes there's only on attribute and that it's the position and that
            // the position is a vec3d
            // TODO: Will: what we need is a std::vector style struct.
            size_t patch_particles_read = 0;
            for (size_t i = 0; i < n_proc_patch->particle_count; ++i) {
              // TODO WILL: This assumes var->vps == PIDX_MAX_DIMENSIONS
              if (pointInChunk(local_proc_patch, (double*)(tmp_patch_read_buf.buffer + i * bytes_per_sample))) {

                memcpy(*var->sim_patch[pc1]->read_particle_buffer + patch_particle_offset,
                    tmp_patch_read_buf.buffer + i * bytes_per_sample, bytes_per_sample);

                ++patch_particles_read;
                patch_particle_offset += bytes_per_sample;

                // Have we got enough room to store the particles we may need to store if we take
                // half of all the rest?
                if (var->sim_patch[pc1]->read_particle_buffer_capacity - patch_particle_offset
                    < (n_proc_patch->particle_count - i) * bytes_per_sample / 2)
                {
                  var->sim_patch[pc1]->read_particle_buffer_capacity +=
                    (n_proc_patch->particle_count - i) * bytes_per_sample / 2 + bytes_per_sample;

                  *var->sim_patch[pc1]->read_particle_buffer = realloc(*var->sim_patch[pc1]->read_particle_buffer,
                      var->sim_patch[pc1]->read_particle_buffer_capacity);
                }
              }
            }

            *var->sim_patch[pc1]->read_particle_count += patch_particles_read;
            var->sim_patch[pc1]->particle_count = *var->sim_patch[pc1]->read_particle_count;

            /*
            printf("Capacity %lu particles, currently have read %d\n",
                var->sim_patch[pc1]->read_particle_buffer_capacity / bytes_per_sample,
                *var->sim_patch[pc1]->read_particle_count);
                */
          }
          close(fpx);

          patch_count++;
        }
        p_counter++;
      }
    }

    free(local_proc_patch);
    local_proc_patch = 0;
    free(n_proc_patch);
    n_proc_patch = 0;

  }

  PIDX_buffer_free(&tmp_patch_read_buf);
  free(file_name);
  free(data_set_path);
  free(directory_path);

  free(size_buffer);

  return PIDX_success;
}


// TODO WILL: Correct this function for intersecting the chunks
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int check_bit = 0;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    check_bit = check_bit
      || (A->physical_offset[d] + A->physical_size[d]) < B->physical_offset[d]
      || (B->physical_offset[d] + B->physical_size[d]) < A->physical_offset[d]; 

  return !check_bit;
}
static int pointInChunk(PIDX_patch p, const double *pos)
{
  int contains_point = 1;
  for (int d = 0; d < PIDX_MAX_DIMENSIONS; ++d)
    contains_point = contains_point
      && pos[d] >= p->physical_offset[d]
      && pos[d] <= p->physical_offset[d] + p->physical_size[d];

  return contains_point;
}

