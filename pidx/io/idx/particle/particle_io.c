#include "../../../PIDX_inc.h"

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi);

PIDX_return_code PIDX_particle_write(PIDX_io file, int gi, int svi, int evi)
{

  int si = 0, p = 0, pt = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  if (group_meta_data_init(file, gi, svi, evi) != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, file->idx->filename, strlen(file->idx->filename) - 4);

  PIDX_variable var0 = var_grp->variable[svi];
  for (p = 0; p < var0->sim_patch_count; p++)
  {
    char file_name[PATH_MAX];
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));
    sprintf(file_name, "%s/time%09d/%d_%d", directory_path, file->idx->current_time_step, file->idx_c->grank, p);

    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    off_t data_offset = 0;
    for (si = svi; si < evi; si++)
    {
      PIDX_variable var = var_grp->variable[si];

      int sample_count;
      int bits_per_sample;
      PIDX_get_datatype_details(var->type_name, &sample_count, &bits_per_sample);

      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      double particle_x = 0, particle_y, particle_z = 0;

      //fprintf(stderr, "[Inside PIDX] Rank %d has %d particles [%lld %lld %lld - %lld %lld %lld]\n", file->idx_c->grank, var->sim_patch[p]->particle_count, var->sim_patch[p]->offset[0], var->sim_patch[p]->offset[1], var->sim_patch[p]->offset[2], var->sim_patch[p]->size[0], var->sim_patch[p]->size[1], var->sim_patch[p]->size[2]);
      if (si == 0)
      {
        for (pt = 0; pt < var->sim_patch[p]->particle_count; pt++)
        {
          memcpy(&particle_x, var->sim_patch[p]->buffer + (pt * sample_count + 0) * sizeof(double), sizeof (double));
          memcpy(&particle_y, var->sim_patch[p]->buffer + (pt * sample_count + 1) * sizeof(double), sizeof (double));
          memcpy(&particle_z, var->sim_patch[p]->buffer + (pt * sample_count + 2) * sizeof(double), sizeof (double));
          printf("%f\t%f\t%f\n", particle_x, particle_y, particle_z);
        }
      }

      ssize_t buffer_size = var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT);
      //printf("[%d] = %d %d %d\n", var->sim_patch[p]->particle_count * sample_count * (bits_per_sample/CHAR_BIT), var->sim_patch[p]->particle_count, sample_count, (bits_per_sample/CHAR_BIT));

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

  return PIDX_success;
}


PIDX_return_code PIDX_particle_read(PIDX_io file, int gi, int svi, int evi)
{

  return PIDX_success;
}


static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->header_io_start = PIDX_get_time();
  // Creates the file heirarchy and writes the header info for all binary files
  ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
  if (ret != PIDX_success)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->header_io_end = PIDX_get_time();

  return PIDX_success;
}
