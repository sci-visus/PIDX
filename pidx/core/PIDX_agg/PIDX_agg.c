#include "../../PIDX_inc.h"

#define PIDX_ACTIVE_TARGET 1

static PIDX_return_code report_error(PIDX_return_code ret, char* file, int line);
static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab);
static PIDX_return_code create_open_log_file (PIDX_agg_id id);
static PIDX_return_code close_log_file (PIDX_agg_id id);
static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode);
static PIDX_return_code aggregate(PIDX_agg_id id, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer ab, PIDX_block_layout lbl, int MODE, int layout_id);

struct PIDX_agg_struct
{

  //MPI_Comm comm;
  //MPI_Comm global_comm;

  MPI_Win win;


  idx_comm idx_c;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int gi;
  int fi;
  int li;

  int ***agg_r;
};



PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int fi, int li)
{
  PIDX_agg_id id;

  id = malloc(sizeof (*id));
  memset(id, 0, sizeof (*id));

  id->idx = idx_meta_data;
  id->idx_d = idx_d;
  id->idx_c = idx_c;

  id->gi = 0;
  id->fi = fi;
  id->li = li;

  return id;
}


PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl)
{
  int i = 0, j = 0, j1 = 0, d = 0;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];

  ab->aggregator_interval = id->idx_c->nprocs / ((id->li - id->fi + 1) * lbl->efc * ab->agg_f);
  assert(ab->aggregator_interval != 0);

  ab->buffer_size = 0;
  ab->sample_number = -1;
  ab->var_number = -1;
  ab->file_number = -1;

  id->agg_r = malloc(lbl->efc * sizeof (int**));
  memset(id->agg_r, 0, lbl->efc * sizeof (int**));
  for (i = 0; i < lbl->efc; i++)
  {
    id->agg_r[i] = malloc((id->li - id->fi + 1) * sizeof (int*));
    memset(id->agg_r[i], 0, (id->li - id->fi + 1) * sizeof (int*));
    for (j = id->fi; j <= id->li; j++)
    {
      j1 = j - id->fi;
      PIDX_variable var = var_grp->variable[j];
      id->agg_r[i][j1] = malloc(var->vps * sizeof (int) * ab->agg_f);
      memset(id->agg_r[i][j1], 0, var->vps * sizeof (int) * ab->agg_f);
      for (d = 0 ; d < var->vps * ab->agg_f; d++)
        id->agg_r[i][j1][d] = -1;
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_agg_buf_create_multiple_level(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset, int var_offset, int file_status)
{
  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  int i = 0, j = 0, k = 0;
  for (k = 0; k < lbl->efc; k++)
  {
    for (i = id->fi; i <= id->li; i++)
    {
      for (j = 0; j < var_grp->variable[i]->vps * ab->agg_f; j++)
      {
        unsigned long long file_index = 0;

        int *first;
        first = malloc(sizeof(*first) * PIDX_MAX_DIMENSIONS);
        memset(first, 0, sizeof(*first) * PIDX_MAX_DIMENSIONS);

        int file_count = 1;

        int negative_file_index = 0;
        if (agg_offset == 0 || agg_offset == 1)
        {
          file_count = 1;
          file_index = 0 * ab->agg_f + j;
        }
        else
        {
          file_count = (int)pow(2, agg_offset - 1);
          negative_file_index = (int)pow(2, agg_offset - 1);
          file_index = (lbl->existing_file_index[k] - negative_file_index) * ab->agg_f + j;
        }

        int bits = (id->idx_d->maxh - 1 - (int)log2(file_count * ab->agg_f));
        file_index = file_index << bits;

        Deinterleave(id->idx->bitPattern, (id->idx_d->maxh - 1), file_index, first);

        int calculated_rank = 0;
#if 0
        int rank_x = first[0] / (var_grp->variable[id->fi]->sim_patch[0]->size[0]);
        int rank_y = first[1] / (var_grp->variable[id->fi]->sim_patch[0]->size[1]);
        int rank_z = first[2] / (var_grp->variable[id->fi]->sim_patch[0]->size[2]);
        int nrank_x = (id->idx->bounds[0] / var_grp->variable[id->fi]->sim_patch[0]->size[0]);
        int nrank_y = (id->idx->bounds[1] / var_grp->variable[id->fi]->sim_patch[0]->size[1]);
#endif
        int rank_x = first[0] / (id->idx->reg_patch_size[0]);
        int rank_y = first[1] / (id->idx->reg_patch_size[1]);
        int rank_z = first[2] / (id->idx->reg_patch_size[2]);
        int nrank_x = ((id->idx_d->partition_size[0] * id->idx_d->partition_count[0]) / id->idx->reg_patch_size[0]);
        int nrank_y = ((id->idx_d->partition_size[1] * id->idx_d->partition_count[1]) / id->idx->reg_patch_size[1]);

#if 0
        int nrank_x = (id->idx_d->partition_size[0] / id->idx->reg_patch_size[0]);
        int nrank_y = (id->idx->bounds[1] / id->idx->reg_patch_size[1]);
#endif
        calculated_rank = rank_x + (rank_y * nrank_x) + (rank_z * nrank_x * nrank_y);

        int trank = 0;
        int interval = (id->idx_c->nprocs/ (lbl->efc * ab->agg_f * id->idx->variable_count));

        if (file_status == 1)
          trank = var_grp->rank_buffer[calculated_rank + var_offset * interval + (interval/2)];
        else if (file_status == 0)
          trank = var_grp->rank_buffer[calculated_rank + var_offset * interval];
        else if (file_status == 2)
          trank = id->idx_c->gnprocs - 1;

        //if (rank == 0)
        //  printf("%d FC %d Bits %d FI %d :: %d %d [%d(%d) %d %d] CR %d [%d %d %d : %d %d %d : %d %d] AI %d trank %d\n", interval, file_count, bits, file_index, file_status, interval, k, lbl->efc, i, j, calculated_rank, first[0], first[1], first[2], rank_x, rank_y, rank_z, nrank_x, nrank_y, calculated_rank + var_offset * interval + (interval/2), trank);

        id->agg_r[k][i - id->fi][j] = trank;

        if(id->idx_c->rank == id->agg_r[k][i - id->fi][j])
        {
          ab->file_number = lbl->existing_file_index[k];
          ab->var_number = i;
          ab->sample_number = j;

          unsigned long long sample_count = lbl->bcpf[ab->file_number] * id->idx_d->samples_per_block / ab->agg_f;

          int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

          int bpdt = 0;
          bpdt = (chunk_size * var_grp->variable[ab->var_number]->bpv/8) / (id->idx->compression_factor);

          ab->buffer_size = sample_count * bpdt;

          if (i == 0)// || i == id->idx->variable_count - 1)
            printf("[G %d] [L %d] [Lid %d] [V %d] [LFi %d] [GFi %d] [Si %d] [F/S/N %d] -> [AGG %d [CR %d (%d %d %d) Rank (%d %d %d)]] [Buffer %lld (%d x %d x %d)]\n", id->idx_c->grank, id->idx_c->rank, agg_offset, i, k, lbl->existing_file_index[k], j, file_status, trank, calculated_rank, first[0], first[1], first[2], rank_x, rank_y, rank_z, ab->buffer_size, lbl->bcpf[ab->file_number], id->idx_d->samples_per_block, bpdt);//, first[0], first[1], first[2], rank_x, rank_y, rank_z);

          ab->buffer = malloc(ab->buffer_size);
          memset(ab->buffer, 0, ab->buffer_size);
          if (ab->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) ab->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
        free(first);
      }
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_global_and_local(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl,  int MODE)
{
  int ret;
  ret = create_window(id, ab);
  if (ret != PIDX_success) report_error(PIDX_err_agg, __FILE__, __LINE__);

  ret = create_open_log_file(id);
  if (ret != PIDX_success) report_error(PIDX_err_agg, __FILE__, __LINE__);

#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, id->win);
  if (ret != MPI_SUCCESS) report_error(PIDX_err_agg, __FILE__, __LINE__);
#endif

  ret = one_sided_data_com(id, ab, layout_id, lbl, MODE);
  if (ret != PIDX_success) report_error(PIDX_err_agg, __FILE__, __LINE__);

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  ret = MPI_Win_fence(0, id->win);
  if (ret != MPI_SUCCESS) report_error(PIDX_err_agg, __FILE__, __LINE__);
#endif

  ret = MPI_Win_free(&(id->win));
  if (ret != MPI_SUCCESS) report_error(PIDX_err_agg, __FILE__, __LINE__);
#endif

  ret = close_log_file(id);
  if (ret != PIDX_success) report_error(PIDX_err_agg, __FILE__, __LINE__);

  return PIDX_success;
}


PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id id, PIDX_block_layout lbl)
{
  int i = 0, j = 0;
  for (i = 0; i < lbl->efc; i++)
  {
    for (j = id->fi; j <= id->li; j++)
    {
      free(id->agg_r[i][j - id->fi]);
      id->agg_r[i][j - id->fi] = 0;
    }
    free(id->agg_r[i]);
  }
  free(id->agg_r);
  id->agg_r = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_destroy(Agg_buffer ab)
{
  if (ab->buffer_size != 0)
  {
    free(ab->buffer);
    ab->buffer = 0;
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_finalize(PIDX_agg_id id)
{
  free(id);
  id = 0;

  return PIDX_success;
}


static PIDX_return_code create_window(PIDX_agg_id id, Agg_buffer ab)
{
  int ret = 0;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var = var_grp->variable[ab->var_number];

  if (ab->buffer_size != 0)
  {
    int tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
    int bpdt = tcs * (var->bpv/8) / (id->idx->compression_factor);

    ret = MPI_Win_create(ab->buffer, ab->buffer_size, bpdt, MPI_INFO_NULL, id->idx_c->comm, &(id->win));
    if (ret != MPI_SUCCESS) report_error(PIDX_err_agg, __FILE__, __LINE__);
  }
  else
  {
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, id->idx_c->comm, &(id->win));
    if (ret != MPI_SUCCESS) report_error(PIDX_err_agg, __FILE__, __LINE__);
  }

  return PIDX_success;
}


static PIDX_return_code create_open_log_file (PIDX_agg_id id)
{
#ifdef PIDX_DUMP_AGG
  if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
  {
    char agg_file_name[1024];
    ret = mkdir(id->idx_d->agg_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, id->idx_d->agg_dump_dir_name);
      return PIDX_err_agg;
    }

#if PIDX_HAVE_MPI
    MPI_Barrier(id->comm);
#endif

    sprintf(agg_file_name, "%s/rank_%d", id->idx_d->agg_dump_dir_name, id->idx_c->rank);
    agg_dump_fp = fopen(agg_file_name, "a+");
    if (!agg_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] agg_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, agg_file_name);
      return PIDX_err_agg;
    }
  }
#endif

  return PIDX_success;
}

static PIDX_return_code close_log_file (PIDX_agg_id id)
{
#ifdef PIDX_DUMP_AGG
  if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
  {
    fprintf(agg_dump_fp, "\n");
    fclose(agg_dump_fp);
  }
#endif

  return PIDX_success;
}


static PIDX_return_code one_sided_data_com(PIDX_agg_id id, Agg_buffer ab, int layout_id, PIDX_block_layout lbl, int mode)
{
  int i, p, v, ret = 0;
  unsigned long long index = 0, count = 0;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var0 = var_grp->variable[id->fi];

  for(v = id->fi; v <= id->li; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for (p = 0; p < var->patch_group_count; p++)
    {
      index = 0, count = 0;
      HZ_buffer hz_buf = var->hz_buffer[p];

      if (hz_buf->type == 1)
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
          fflush(agg_dump_fp);
        }
#endif

        for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
        {
          if (hz_buf->nsamples_per_level[i][0] * hz_buf->nsamples_per_level[i][1] * hz_buf->nsamples_per_level[i][2] != 0)
          {
            index = 0;
            count =  hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1;

#ifdef PIDX_DUMP_AGG
            if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif
            //printf("A [Level %d] : Offset %d Count %d\n", i, hz_buf->start_hz_index[i], count);
            ret = aggregate(id, v, hz_buf->start_hz_index[i], count, hz_buf->buffer[i], 0, ab, lbl, mode, layout_id);
            if (ret != PIDX_success)
            {
              fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
              return PIDX_err_agg;
            }
          }
        }
      }
      else if (hz_buf->type == 2)
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "Type %d Variable %d Patch %d\n",hz_buf->type, v, p);
          fflush(agg_dump_fp);
        }
#endif

        for (i = lbl->resolution_from; i < lbl->resolution_to; i++)
        {
          if (var0->hz_buffer[p]->nsamples_per_level[i][0] * var0->hz_buffer[p]->nsamples_per_level[i][1] * var0->hz_buffer[p]->nsamples_per_level[i][2] != 0)
          {
#ifdef PIDX_DUMP_AGG
            if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
            {
              fprintf(agg_dump_fp, "[%d]: ", i);
              fflush(agg_dump_fp);
            }
#endif
            int start_block_index = hz_buf->start_hz_index[i] / id->idx_d->samples_per_block;
            int end_block_index = hz_buf->end_hz_index[i] / id->idx_d->samples_per_block;
            assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);

            if (end_block_index == start_block_index)
            {
              count = (hz_buf->end_hz_index[i] - hz_buf->start_hz_index[i] + 1);
              //printf("B [Level %d] : Offset %d Count %d\n", i, hz_buf->start_hz_index[i], count);
              ret = aggregate(id, v, var0->hz_buffer[p]->start_hz_index[i], count, hz_buf->buffer[i], 0, ab, lbl, mode, layout_id);
              if (ret != PIDX_success)
              {
                fprintf(stderr, " Error in aggregate Line %d File %s\n", __LINE__, __FILE__);
                return PIDX_err_agg;
              }
            }
            //
            else
            {
              int send_index = 0;
              int bl;
              for (bl = start_block_index; bl <= end_block_index; bl++)
              {
                if (PIDX_blocks_is_block_present(bl, lbl))
                {
                  if (bl == start_block_index)
                  {
                    index = 0;
                    count = ((start_block_index + 1) * id->idx_d->samples_per_block) - hz_buf->start_hz_index[i];
                  }
                  else if (bl == end_block_index)
                  {
                    index = (end_block_index * id->idx_d->samples_per_block - hz_buf->start_hz_index[i]);
                    count = hz_buf->end_hz_index[i] - ((end_block_index) * id->idx_d->samples_per_block) + 1;
                  }
                  else
                  {
                    index = (bl * id->idx_d->samples_per_block - hz_buf->start_hz_index[i]);
                    count = id->idx_d->samples_per_block;
                  }

                  //printf("C [Level %d] : Offset %d Count %d\n", i, hz_buf->start_hz_index[i], count);
                  ret = aggregate(id, v, index + hz_buf->start_hz_index[i], count, hz_buf->buffer[i], send_index, ab, lbl, mode, layout_id);
                  if (ret != PIDX_success)
                  {
                    fprintf(stderr, "[%s] [%d] write_read_samples() failed.\n", __FILE__, __LINE__);
                    return PIDX_err_agg;
                  }
                  send_index = send_index + count;
                }
                else
                  send_index = send_index + id->idx_d->samples_per_block;
              }
            }
            //
          }
        }
      }
    }
  }

  return PIDX_success;
}

static PIDX_return_code aggregate(PIDX_agg_id id, int variable_index, unsigned long long hz_start, unsigned long long hz_count, unsigned char* hz_buffer, int buffer_offset, Agg_buffer ab, PIDX_block_layout lbl, int MODE, int layout_id)
{
  int ret;
  int itr;
  int bpdt;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, vps;
  int target_rank = 0;
  unsigned long long start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, samples_in_file = 0;
  unsigned long long samples_per_file = (unsigned long long) id->idx_d->samples_per_block * id->idx->blocks_per_file;

  unsigned long long tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];


  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var = var_grp->variable[variable_index];

  vps = var->vps; //number of samples for variable j

  // hz_start is the starting HZ index for the data buffer at level "level" and for regular patch number "patch"
  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / id->idx_d->samples_per_block;


  //number of empty blocks befor block "block_no" in the file "file_no"
  //negative_block_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, block_no, id->idx->variable[id->ini]->global_block_layout);
  negative_block_offset = PIDX_blocks_find_negative_offset(id->idx->blocks_per_file, block_no, lbl);
  if (negative_block_offset < 0)
    return PIDX_err_agg;

  //number of samples in file "file_no"
  //samples_in_file = id->idx->variable[id->ini]->bcpf[file_no] * id->idx_d->samples_per_block;
  samples_in_file = lbl->bcpf[file_no] * id->idx_d->samples_per_block;
  if (samples_in_file > samples_per_file)
    return PIDX_err_agg;

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * id->idx_d->samples_per_block))) * vps)
      %
      (samples_in_file * vps);

  sample_index = target_disp / (samples_in_file / ab->agg_f);
  if (sample_index >= var->vps * ab->agg_f)
    return PIDX_err_agg;

  target_disp = target_disp % (samples_in_file / ab->agg_f);

  //if (rank == 60)
  //  printf("%d ----> %d TD %d NBO %d\n", hz_start, block_no, target_disp, negative_block_offset);

  target_rank = id->agg_r[lbl->inverse_existing_file_index[file_no]][variable_index - id->fi][sample_index];
  //printf("[%d] File no %d IFI %d TR %d\n", rank, file_no, lbl->inverse_existing_file_index[file_no], target_rank);

  /*
  if (layout_id != 0 && id->idx->current_time_step == 0)
  {
  MPI_Comm agg_comm;
  int max_rank = 0;
  int min_rank = 0;

  MPI_Comm_split(id->comm, target_rank, rank, &agg_comm);
  int nrank = 0;
  MPI_Comm_rank(agg_comm, &nrank);

  MPI_Allreduce(&rank, &max_rank, 1, MPI_INT, MPI_MAX, agg_comm);
  MPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, agg_comm);

  MPI_Comm_free(&agg_comm);

  if (target_rank < min_rank || target_rank > max_rank || rank < min_rank || rank > max_rank)
  {
    printf("[TR %d] [%d] V %d P %d A %d FN %d [%d %d]\n", target_rank, rank, variable_index, layout_id, lbl->inverse_existing_file_index[file_no], file_no, min_rank, max_rank);
  }

  assert(target_rank >= min_rank);
  assert(target_rank <= max_rank);
  assert(rank >= min_rank);
  assert(rank <= max_rank);
  }
  */

  target_count = hz_count * vps;
  bpdt = ((var->bpv / 8) * tcs) / (id->idx->compression_factor);
  hz_buffer = hz_buffer + buffer_offset * bpdt * vps;

  start_agg_index = target_disp / (unsigned long long) (samples_in_file / ab->agg_f);
  end_agg_index = ((target_disp + target_count - 1) / (unsigned long long) (samples_in_file / ab->agg_f));

  //printf("SAI : EAI TR %d FN %d IFN %d :: %d : %d\n", target_rank, file_no, lbl->inverse_existing_file_index[file_no], start_agg_index, end_agg_index);
  if (start_agg_index != end_agg_index)
  {
    if (target_rank != id->idx_c->rank)
    {
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , id->win);
#endif
      //target_disp_address = target_disp;
      if (MODE == PIDX_WRITE)
      {
        ret = MPI_Put(hz_buffer, ((samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        printf("[1] TR %d\n", target_rank);
        ret = MPI_Get(hz_buffer, ((samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, target_rank, target_disp, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }

#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, id->win);
#endif
    }
    else
    {
      if (MODE == PIDX_WRITE)
      {
        memcpy( ab->buffer + target_disp * bpdt, hz_buffer, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt);
      }
      else
      {
        memcpy( hz_buffer, ab->buffer + target_disp * bpdt, ( (samples_in_file / ab->agg_f) - target_disp) * bpdt);
      }
    }

    for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++)
    {
      if (target_rank != id->idx_c->rank)
      {
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank + ab->aggregator_interval, 0, id->win);
#endif
        if (MODE == PIDX_WRITE)
        {
          ret = MPI_Put(hz_buffer + (( (samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
        else
        {
          printf("[2] TR %d\n", target_rank + ab->aggregator_interval);
          ret = MPI_Get(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (samples_in_file / ab->agg_f) * bpdt, MPI_BYTE, id->win);
          if (ret != MPI_SUCCESS)
          {
            fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
#ifndef PIDX_ACTIVE_TARGET
        MPI_Win_unlock(target_rank + ab->aggregator_interval, id->win);
#endif
      }
      else
      {
        if (MODE == PIDX_WRITE)
        {
          memcpy( ab->buffer, hz_buffer + (( (samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, (samples_in_file / ab->agg_f) * bpdt);
        }
        else
        {
          memcpy( hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + (itr * (samples_in_file / ab->agg_f))) * bpdt, ab->buffer, (samples_in_file / ab->agg_f) * bpdt);
        }
      }
    }

    if (target_rank + ab->aggregator_interval != id->idx_c->rank)
    {
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank + ab->aggregator_interval, 0, id->win);
#endif
      if (MODE == PIDX_WRITE)
      {
        ret = MPI_Put(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / ab->agg_f))) + (((samples_in_file / ab->agg_f)) - target_disp))) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
        printf("[3] TR %d (%d %d)\n", target_rank + ab->aggregator_interval, target_rank, ab->aggregator_interval);
        ret = MPI_Get(hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - (((end_agg_index - start_agg_index - 1) * ((samples_in_file / ab->agg_f))) + (((samples_in_file / ab->agg_f)) - target_disp))) * bpdt, MPI_BYTE, target_rank + ab->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Get Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank + ab->aggregator_interval, id->win);
#endif
    }
    else
    {
      if(MODE == PIDX_WRITE)
      {
        memcpy( ab->buffer, hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt);
      }
      else
      {
        memcpy( hz_buffer + (((samples_in_file / ab->agg_f) - target_disp) + ((end_agg_index - start_agg_index - 1) * (samples_in_file / ab->agg_f))) * bpdt, ab->buffer, (target_count - ((end_agg_index - start_agg_index) * (samples_in_file / ab->agg_f) - target_disp)) * bpdt);
      }
    }
  }
  //{
  //  printf("read vs write conflict!\n");
  //}
  else
  {
    if(target_rank != id->idx_c->rank)
    {
#if PIDX_HAVE_MPI
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , id->win);
#endif
      //target_disp_address = target_disp;
      if(MODE == PIDX_WRITE)
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Put(hz_buffer, hz_count * vps * bpdt, MPI_BYTE, target_rank, target_disp, hz_count * vps * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[D] Target Rank %d Count %lld Local Disp %d Target Disp %lld\n", target_rank,  (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif

        ret = MPI_Get(hz_buffer, hz_count * vps * bpdt, MPI_BYTE, target_rank, target_disp, hz_count * vps * bpdt, MPI_BYTE, id->win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
#ifndef PIDX_ACTIVE_TARGET
      MPI_Win_unlock(target_rank, id->win);
#endif
#endif
    }
    else
    {
      if(MODE == PIDX_WRITE)
      {

#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( ab->buffer + target_disp * bpdt, hz_buffer, hz_count * vps * bpdt);
      }
      else
      {
#ifdef PIDX_DUMP_AGG
        if (id->idx_d->dump_agg_info == 1 && id->idx->current_time_step == 0)
        {
          fprintf(agg_dump_fp, "[MD] Count %lld Local Disp %d Target Disp %lld\n", (long long)hz_count, 0, (long long)target_disp);
          fflush(agg_dump_fp);
        }
#endif
        memcpy( hz_buffer, ab->buffer + target_disp * bpdt, hz_count * vps * bpdt);
      }
    }
  }

  return PIDX_success;
}

static PIDX_return_code report_error(PIDX_return_code ret, char* file, int line)
{
  fprintf(stdout,"File %s Line %d\n", file, line);
  return ret;
}
