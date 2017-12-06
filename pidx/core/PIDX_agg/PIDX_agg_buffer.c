#include "../../PIDX_inc.h"

static int intersectNDChunk(PIDX_patch A, PIDX_patch B);

PIDX_return_code PIDX_agg_create_randomized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset, int var_offset, int file_status)
{
  return PIDX_success;
}


PIDX_return_code PIDX_agg_buf_create_local_uniform_dist(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl)
{
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];

  int rank_counter = 0;
  int aggregator_interval = id->idx_c->lnprocs / ((id->fi - id->li + 1) * lbl->efc);

  for (k = 0; k < lbl->efc; k++)
  {
    for (i = id->fi; i <= id->li; i++)
    {
      for (j = 0; j < var_grp->variable[i]->vps * ab->agg_f; j++)
      {
        id->agg_r[k][i - id->fi][j] = rank_counter;
        rank_counter = rank_counter + aggregator_interval;

        if(id->idx_c->lrank == id->agg_r[k][i - id->fi][j])
        {
          ab->file_number = lbl->existing_file_index[k];
          ab->var_number = i;
          ab->sample_number = j;

          unsigned long long sample_count = lbl->bcpf[ab->file_number] * id->idx_d->samples_per_block / ab->agg_f;

          int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];

          int bpdt = 0;
          bpdt = (chunk_size * var_grp->variable[ab->var_number]->bpv/8) / (id->idx->compression_factor);

          ab->buffer_size = sample_count * bpdt;

          ab->buffer = malloc(ab->buffer_size);
          memset(ab->buffer, 0, ab->buffer_size);
          if (ab->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) ab->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_create_global_partition_localized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset)
{
  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var0 = var_grp->variable[id->fi];
  int i = 0, j = 0, k = 0;
  unsigned long long local_patch_offset[PIDX_MAX_DIMENSIONS];
  unsigned long long local_patch_size[PIDX_MAX_DIMENSIONS];

  int d = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_patch_offset[d] = var0->restructured_super_patch->restructured_patch->offset[d];
    local_patch_size[d] = var0->restructured_super_patch->restructured_patch->size[d];
  }

  int wc = id->idx_c->lnprocs * (PIDX_MAX_DIMENSIONS);

  unsigned long long* global_patch_offset = malloc(wc * sizeof(*global_patch_offset));
  memset(global_patch_offset, 0, wc * sizeof(*global_patch_offset));

  unsigned long long* global_patch_size = malloc(wc * sizeof(*global_patch_size));
  memset(global_patch_size, 0, wc * sizeof(*global_patch_size));

  MPI_Allgather(&local_patch_offset[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, id->idx_c->local_comm);

  MPI_Allgather(&local_patch_size[0], PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, id->idx_c->local_comm);

#if 0
  if (id->idx_c->grank == 0 && agg_offset == 0)
  {
    for (i = 0; i < id->idx_c->lnprocs; i++)
      fprintf(stderr, "[%d] (%d) -----> %d %d %d - %d %d %d\n", i, id->idx_c->lnprocs, (int)global_patch_offset[PIDX_MAX_DIMENSIONS * i + 0], (int)global_patch_offset[PIDX_MAX_DIMENSIONS * i + 1], (int)global_patch_offset[PIDX_MAX_DIMENSIONS * i + 2], (int)global_patch_size[PIDX_MAX_DIMENSIONS * i + 0], (int)global_patch_size[PIDX_MAX_DIMENSIONS * i + 1], (int)global_patch_size[PIDX_MAX_DIMENSIONS * i + 2]);
  }
#endif

  for (k = 0; k < lbl->efc; k++)
  {
    for (i = id->fi; i <= id->li; i++)
    {
      for (j = 0; j < var_grp->variable[i]->vps * ab->agg_f; j++)
      {
        int start_rank = -1, end_rank = -1;
        unsigned long long global_file_index = lbl->existing_file_index[k];

        int first_block = -1, last_block = -1;
        int b = 0;
        for (b = 0; b < id->idx->blocks_per_file; b++)
        {
          if (id->idx_d->block_bitmap[global_file_index][b] == 1)
          {
            first_block = b;
            break;
          }
        }
        for (b = id->idx->blocks_per_file - 1; b >= 0; b--)
        {
          if (id->idx_d->block_bitmap[global_file_index][b] == 1)
          {
            last_block = b;
            break;
          }
        }
        assert(last_block == lbl->lbi[global_file_index]);

        int s = 0;
        int last_index = -1;
        int first_index = -1;
        unsigned long long ZYX[PIDX_MAX_DIMENSIONS];

#if MULTI_BOX
          int break_counter1 = 0;
          for (s = 0; s < id->idx_d->samples_per_block; s++)
          {
            // HZ index of last sample of the block.
            last_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (lbl->lbi[global_file_index] + 1) * id->idx_d->samples_per_block - 1 - s;

            // xyz index of the last sample of the block.
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, last_index, ZYX);

            // check to see if the sample is within bounds.
            for (i1 = 0; i1 < id->idx_c->lnprocs * max_patch_count; i1++)
            {
              if (ZYX[0] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 0] && ZYX[0] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 0] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 0]) &&
                  ZYX[1] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 1] && ZYX[1] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 1] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 1]) &&
                  ZYX[2] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 2] && ZYX[2] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 2] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 2]))
              {
                break_counter1 = 1;
                break;
              }
            }
            if (break_counter1 == 1)
              break;
          }
#else
          for (s = 0; s < id->idx_d->samples_per_block; s++)
          {
            // HZ index of last sample of the block.
            last_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (lbl->lbi[global_file_index] + 1) * id->idx_d->samples_per_block - 1 - s;

            // xyz index of the last sample of the block.
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, last_index, ZYX);

            // check to see if the sample is within bounds.
            if (ZYX[0] < PIDX_MIN(id->idx->box_bounds[0], id->idx_d->partition_offset[0] + id->idx_d->partition_size[0]) &&
                ZYX[1] < PIDX_MIN(id->idx->box_bounds[1], id->idx_d->partition_offset[1] + id->idx_d->partition_size[1]) &&
                ZYX[2] < PIDX_MIN(id->idx->box_bounds[2], id->idx_d->partition_offset[2] + id->idx_d->partition_size[2]))
              break;
          }
#endif

#if MULTI_BOX
          int break_counter2 = 0;
          for (s = 0; s < id->idx_d->samples_per_block; s++)
          {
             // HZ index of first sample of the block.
            first_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (first_block) * id->idx_d->samples_per_block + s;

            // xyz index of the first sample of the block.
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, first_index, ZYX);

            // check to see if the sample is within bounds.
            for (i1 = 0; i1 < id->idx_c->lnprocs * max_patch_count; i1++)
            {
              if (ZYX[0] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 0] && ZYX[0] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 0] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 0]) &&
                  ZYX[1] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 1] && ZYX[1] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 1] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 1]) &&
                  ZYX[2] >= id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 2] && ZYX[2] < (id->idx->all_offset[PIDX_MAX_DIMENSIONS * i1 + 2] + id->idx->all_size[PIDX_MAX_DIMENSIONS * i1 + 2]))
              {
                break_counter2 = 1;
                break;
              }
            }
            if (break_counter2 == 1)
              break;
          }
#else
          for (s = 0; s < id->idx_d->samples_per_block; s++)
          {
             // HZ index of first sample of the block.
            first_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (first_block) * id->idx_d->samples_per_block + s;

            // xyz index of the first sample of the block.
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, first_index, ZYX);

            // check to see if the sample is within bounds.
            if (ZYX[0] >= id->idx_d->partition_offset[0] && ZYX[0] < PIDX_MIN(id->idx->box_bounds[0], id->idx_d->partition_offset[0] + id->idx_d->partition_size[0]) &&
                ZYX[1] >= id->idx_d->partition_offset[1] && ZYX[1] < PIDX_MIN(id->idx->box_bounds[1], id->idx_d->partition_offset[1] + id->idx_d->partition_size[1]) &&
                ZYX[2] >= id->idx_d->partition_offset[2] && ZYX[2] < PIDX_MIN(id->idx->box_bounds[2], id->idx_d->partition_offset[2] + id->idx_d->partition_size[2]))
              break;
          }
#endif

        unsigned long long global_start_hz = first_index;
        unsigned long long global_end_hz = last_index;
        unsigned long long global_start_ZYX[PIDX_MAX_DIMENSIONS], global_end_ZYX[PIDX_MAX_DIMENSIONS];

        Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_start_hz, global_start_ZYX);
        Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_end_hz, global_end_ZYX);

        PIDX_patch global_start_point = (PIDX_patch)malloc(sizeof (*global_start_point));
        memset(global_start_point, 0, sizeof (*global_start_point));
        global_start_point->offset[0] = global_start_ZYX[0];
        global_start_point->offset[1] = global_start_ZYX[1];
        global_start_point->offset[2] = global_start_ZYX[2];
        global_start_point->size[0] = 1;
        global_start_point->size[1] = 1;
        global_start_point->size[2] = 1;

        //Extent of process with rank r
        PIDX_patch rank_r_patch = malloc(sizeof (*rank_r_patch));
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));

        int r = 0, d = 0, m = 0;
        int break_counter = 0;
        for (r = 0; r < id->idx_c->lnprocs; r++)
        {
          if (global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 0] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 1] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 2] == 0)
          {
            continue;
          }

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            rank_r_patch->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
            rank_r_patch->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
          }

          if (intersectNDChunk(global_start_point, rank_r_patch))
          {
            start_rank = r;
            break_counter = 1;
            break;
          }

          if (break_counter == 1)
            break;
        }
        free(global_start_point);

        PIDX_patch global_end_point = (PIDX_patch)malloc(sizeof (*global_end_point));
        memset(global_end_point, 0, sizeof (*global_end_point));
        global_end_point->offset[0] = global_end_ZYX[0];
        global_end_point->offset[1] = global_end_ZYX[1];
        global_end_point->offset[2] = global_end_ZYX[2];
        global_end_point->size[0] = 1;
        global_end_point->size[1] = 1;
        global_end_point->size[2] = 1;

        break_counter = 0;
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));
        for (r = 0; r < id->idx_c->lnprocs; r++)
        {
          if (global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 0] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 1] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 2] == 0)
            continue;

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            rank_r_patch->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
            rank_r_patch->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
          }

          if (intersectNDChunk(global_end_point, rank_r_patch))
          {
            end_rank = r;
            break_counter = 1;
            break;
          }

          if (break_counter == 1)
            break;
        }
        free(rank_r_patch);
        free(global_end_point);

        float range = (float)(end_rank - start_rank + 1) / (id->idx_d->variable_pipe_length + 1);

#if 0
        if (agg_offset < var_grp->shared_end_layout_index)
          id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->lvi) * range);
        else
          id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->lvi) * range) + (range/2);
#endif

        if (agg_offset < var_grp->shared_end_layout_index)
          id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->fi) * range);
        else
          id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->fi) * range) + (range/2);


        if (id->idx_c->lrank == id->agg_r[k][i - id->fi][j])
        {
          ab->file_number = lbl->existing_file_index[k];
          ab->var_number = i;
          ab->sample_number = j;

          unsigned long long sample_count = lbl->bcpf[ab->file_number] * id->idx_d->samples_per_block / ab->agg_f;
          int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
          int bpdt = (chunk_size * var_grp->variable[ab->var_number]->bpv/8) / (id->idx->compression_factor);

          ab->buffer_size = sample_count * bpdt;

#if DETAIL_OUTPUT
          fprintf(stderr, "[Lid %d] [TS %d] [%d] [C %d] [G %d %d] [L %d %d] [S E R %d (%lld : %lld %lld %lld) - %d (%lld : %lld %lld %lld) R %f] [V %d P %d] [LFi %d] [GFi %d] [Si %d] [Buffer %lld (%d x %d x %d)]\n",
               agg_offset, id->idx->current_time_step,
               id->agg_r[k][i - id->fi][j],
               id->idx_d->color,
               id->idx_c->grank, id->idx_c->gnprocs,
               id->idx_c->lrank, id->idx_c->lnprocs,
               start_rank, global_start_hz, global_start_ZYX[0], global_start_ZYX[1], global_start_ZYX[2],
               end_rank, global_end_hz, global_end_ZYX[0], global_end_ZYX[1], global_end_ZYX[2],
               range,
               i, (id->idx_d->variable_pipe_length + 1),
               k,
               lbl->existing_file_index[k],
               j,
               ab->buffer_size, lbl->bcpf[ab->file_number], id->idx_d->samples_per_block, bpdt);
#endif

          ab->buffer = malloc(ab->buffer_size);
          memset(ab->buffer, 0, ab->buffer_size);
          if (ab->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) ab->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
  }

  free(global_patch_size);
  free(global_patch_offset);

  return PIDX_success;
}



PIDX_return_code PIDX_agg_create_local_partition_localized_aggregation_buffer(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int agg_offset)
{
  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var0 = var_grp->variable[id->fi];
  int i = 0, j = 0, k = 0;
  unsigned long long local_patch_offset[PIDX_MAX_DIMENSIONS];
  unsigned long long local_patch_size[PIDX_MAX_DIMENSIONS];

  int d = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_patch_offset[d] = var0->restructured_super_patch->restructured_patch->offset[d];
    assert(local_patch_offset[d] % id->idx->chunk_size[d] == 0);
    local_patch_offset[d] = local_patch_offset[d] / id->idx->chunk_size[d];

    local_patch_size[d] = var0->restructured_super_patch->restructured_patch->size[d];
    if (local_patch_size[d] % id->idx->chunk_size[d] == 0)
      local_patch_size[d] = local_patch_size[d] / id->idx->chunk_size[d];
    else
      local_patch_size[d] = (local_patch_size[d] / id->idx->chunk_size[d]) + 1;
  }

  int wc = id->idx_c->lnprocs * PIDX_MAX_DIMENSIONS;

  //if (id->idx_d->color == 1)
  //fprintf(stderr, "[%d] [%d] O : S :: %d %d %d - %d %d %d\n", id->idx_d->color, id->idx_c->lrank, local_patch_offset[0], local_patch_offset[1], local_patch_offset[2], local_patch_size[0], local_patch_size[1], local_patch_size[2]);

  unsigned long long* global_patch_offset = malloc(wc * sizeof(*global_patch_offset));
  memset(global_patch_offset, 0, wc * sizeof(*global_patch_offset));

  unsigned long long* global_patch_size = malloc(wc * sizeof(*global_patch_size));
  memset(global_patch_size, 0, wc * sizeof(*global_patch_size));

  MPI_Allgather(local_patch_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, id->idx_c->local_comm);

  MPI_Allgather(local_patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, global_patch_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, id->idx_c->local_comm);

  //fprintf(stderr, "[%d : %d] - %d %d %d - %d %d %d\n", id->idx_c->lrank, id->idx_c->grank, (int)local_patch_offset[0], (int)local_patch_offset[1], (int)local_patch_offset[2], (int)local_patch_size[0], (int)local_patch_size[1], (int)local_patch_size[2]);

  //fprintf(stderr, "Fi Li %d %d\n", id->fi, id->li);
  for (k = 0; k < lbl->efc; k++)
  {
    for (i = id->fi; i <= id->li; i++)
    {
      for (j = 0; j < var_grp->variable[i]->vps * ab->agg_f; j++)
      {
        int start_rank = -1, end_rank = -1;
        unsigned long long global_file_index = lbl->existing_file_index[k];

        int first_block = -1, last_block = -1;
        int b = 0;
        for (b = 0; b < id->idx->blocks_per_file; b++)
        {
          if (id->idx_d->block_bitmap[global_file_index][b] == 1)
          {
            first_block = b;
            break;
          }
        }
        for (b = id->idx->blocks_per_file - 1; b >= 0; b--)
        {
          if (id->idx_d->block_bitmap[global_file_index][b] == 1)
          {
            last_block = b;
            break;
          }
        }
        assert(last_block == lbl->lbi[global_file_index]);

        int s = 0;
        int last_index = -1;
        int first_index = -1;
        unsigned long long ZYX[PIDX_MAX_DIMENSIONS];

#if 0
        int x = 0;
        for (x = 0; x <= last_block; x++)
        {
          for (s = 0; s < id->idx_d->samples_per_block; s++)
          {
            // HZ index of last sample of the block.
            last_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (/*lbl->lbi[global_file_index]*/x + 1) * id->idx_d->samples_per_block - 1 - s;

            // xyz index of the last sample of the block.
            Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, last_index, ZYX);

            if (id->idx_c->lrank == 0)
              fprintf(stderr, "[%d] ZYX: %d %d %d\n", x, ZYX[0], ZYX[1], ZYX[2]);

            // check to see if the sample is within bounds.
            if (ZYX[0] < id->idx->box_bounds[0] / id->idx->chunk_size[0] &&
                ZYX[1] < id->idx->box_bounds[1] / id->idx->chunk_size[1] &&
                ZYX[2] < id->idx->box_bounds[2] / id->idx->chunk_size[2])
              break;
          }
        },ik
#endif
        last_index = -1;
        for (s = 0; s < id->idx_d->samples_per_block; s++)
        {
          // HZ index of last sample of the block.
          last_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (lbl->lbi[global_file_index] + 1) * id->idx_d->samples_per_block - 1 - s;

          // xyz index of the last sample of the block.
          Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, last_index, ZYX);

          //if (id->idx_c->lrank == 0)
          //  fprintf(stderr, "ZYX: %d %d %d\n", ZYX[0], ZYX[1], ZYX[2]);

          // check to see if the sample is within bounds.
          if (ZYX[0] < id->idx->box_bounds[0] / id->idx->chunk_size[0] &&
              ZYX[1] < id->idx->box_bounds[1] / id->idx->chunk_size[1] &&
              ZYX[2] < id->idx->box_bounds[2] / id->idx->chunk_size[2])
            break;
        }

        for (s = 0; s < id->idx_d->samples_per_block; s++)
        {
          // HZ index of first sample of the block.
          first_index = global_file_index * id->idx->blocks_per_file * id->idx_d->samples_per_block + (first_block) * id->idx_d->samples_per_block + s;

          // xyz index of the first sample of the block.
          Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, first_index, ZYX);

          // check to see if the sample is within bounds.
          if (ZYX[0] >= 0 && ZYX[0] < id->idx->box_bounds[0] / id->idx->chunk_size[0] &&
              ZYX[1] >= 0 && ZYX[1] < id->idx->box_bounds[1] / id->idx->chunk_size[1] &&
              ZYX[2] >= 0 && ZYX[2] < id->idx->box_bounds[2] / id->idx->chunk_size[2])
            break;
        }


        unsigned long long global_start_hz = first_index;
        unsigned long long global_end_hz = last_index;

        unsigned long long global_start_ZYX[PIDX_MAX_DIMENSIONS], global_end_ZYX[PIDX_MAX_DIMENSIONS];

        Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_start_hz, global_start_ZYX);
        Hz_to_xyz(id->idx->bitPattern, id->idx_d->maxh - 1, global_end_hz, global_end_ZYX);

        PIDX_patch global_start_point = (PIDX_patch)malloc(sizeof (*global_start_point));
        memset(global_start_point, 0, sizeof (*global_start_point));
        global_start_point->offset[0] = global_start_ZYX[0];
        global_start_point->offset[1] = global_start_ZYX[1];
        global_start_point->offset[2] = global_start_ZYX[2];
        global_start_point->size[0] = 1;
        global_start_point->size[1] = 1;
        global_start_point->size[2] = 1;

        //Extent of process with rank r
        PIDX_patch rank_r_patch = malloc(sizeof (*rank_r_patch));
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));

        int r = 0, d = 0, m = 0;
        int break_counter = 0;
        for (r = 0; r < id->idx_c->lnprocs; r++)
        {
          if (global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 0] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 1] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 2] == 0)
          {
            continue;
          }

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            rank_r_patch->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
            rank_r_patch->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
          }

          if (intersectNDChunk(global_start_point, rank_r_patch))
          {
            start_rank = r;
            break_counter = 1;
            break;
          }

          if (break_counter == 1)
            break;
        }
        free(global_start_point);

        PIDX_patch global_end_point = (PIDX_patch)malloc(sizeof (*global_end_point));
        memset(global_end_point, 0, sizeof (*global_end_point));
        global_end_point->offset[0] = global_end_ZYX[0];
        global_end_point->offset[1] = global_end_ZYX[1];
        global_end_point->offset[2] = global_end_ZYX[2];
        global_end_point->size[0] = 1;
        global_end_point->size[1] = 1;
        global_end_point->size[2] = 1;

        break_counter = 0;
        memset(rank_r_patch, 0, sizeof (*rank_r_patch));
        for (r = 0; r < id->idx_c->lnprocs; r++)
        {
          if (global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 0] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 1] == 0 &&
              global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + 2] == 0)
            continue;

          for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
          {
            rank_r_patch->offset[d] = global_patch_offset[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
            rank_r_patch->size[d] = global_patch_size[PIDX_MAX_DIMENSIONS * r + m * PIDX_MAX_DIMENSIONS + d];
          }

          if (intersectNDChunk(global_end_point, rank_r_patch))
          {
            end_rank = r;
            break_counter = 1;
            break;
          }

          if (break_counter == 1)
            break;
        }
        free(rank_r_patch);
        free(global_end_point);

        float range = (float)(end_rank - start_rank + 1) / (id->idx_d->variable_pipe_length + 1);

#if 0
        if (agg_offset < var_grp->shared_end_layout_index)
            id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->lvi) * range);
        else
            id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->lvi) * range) + (range/2);
#endif
        if (agg_offset < var_grp->shared_end_layout_index)
            id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->fi) * range);
        else
            id->agg_r[k][i - id->fi][j] = start_rank + (int)((float)(i - id->fi) * range) + (range/2);

        if (id->idx_c->lrank == id->agg_r[k][i - id->fi][j])
        {
          ab->file_number = lbl->existing_file_index[k];
          ab->var_number = i;
          ab->sample_number = j;

          unsigned long long sample_count = lbl->bcpf[ab->file_number] * id->idx_d->samples_per_block / ab->agg_f;
          int chunk_size = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
          int bpdt = (chunk_size * var_grp->variable[ab->var_number]->bpv/8) / (id->idx->compression_factor);

          ab->buffer_size = sample_count * bpdt;

#if DETAIL_OUTPUT
          fprintf(stderr, "[Lid %d] [TS %d] [%d] [C %d] [G %d %d] [L %d %d] [S E R %d (%lld : %lld %lld %lld) - %d (%lld : %lld %lld %lld) R %f] [V %d LVI %d P %d] [LFi %d] [GFi %d %d] [Si %d] [Buffer %lld (%d x %d x %d)]\n",
               agg_offset, id->idx->current_time_step,
               id->agg_r[k][i - id->fi][j],
               id->idx_d->color,
               id->idx_c->grank, id->idx_c->gnprocs,
               id->idx_c->lrank, id->idx_c->lnprocs,
               start_rank, global_start_hz, global_start_ZYX[0], global_start_ZYX[1], global_start_ZYX[2],
               end_rank, global_end_hz, global_end_ZYX[0], global_end_ZYX[1], global_end_ZYX[2],
               range,
               i, (i - id->fi), (id->idx_d->variable_pipe_length + 1),
               k,
               lbl->existing_file_index[k], global_file_index,
               j,
               ab->buffer_size, lbl->bcpf[ab->file_number], id->idx_d->samples_per_block, bpdt);
#endif

          //fprintf(stderr, "Agg buffer size %d (%d x %d (%d x %d / %d))\n", ab->buffer_size, sample_count, bpdt, chunk_size, var_grp->variable[ab->var_number]->bpv/8, id->idx->compression_factor);
          ab->buffer = malloc(ab->buffer_size);
          memset(ab->buffer, 0, ab->buffer_size);
          if (ab->buffer == NULL)
          {
            fprintf(stderr, " Error in malloc %lld: Line %d File %s\n", (long long) ab->buffer_size, __LINE__, __FILE__);
            return PIDX_err_agg;
          }
        }
      }
    }
  }

  free(global_patch_size);
  free(global_patch_offset);

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


/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}
