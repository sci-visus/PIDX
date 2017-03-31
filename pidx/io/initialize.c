#include "../PIDX_inc.h"

static int maximum_neighbor_count = 256;
static int getPowerOftwo(int x);
static int calculate_patch_group_count_for_patch_per_process(PIDX_io file, int gi, int svi);
static int calculate_patch_group_count_for_multi_patch_per_process(PIDX_io file, int gi, int svi, unsigned long long sim_max_patch_group_count);
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);
static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file);
static int contains_patch(Ndim_patch reg_patch, Ndim_patch* patches, int count);

PIDX_return_code set_rst_box_size(PIDX_io file, int gi, int svi)
{
  int factor = 1;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  if (file->idx->reg_box_set == PIDX_CLOSEST_POWER_TWO)
  {
    file->idx->reg_patch_size[0] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[0]);
    file->idx->reg_patch_size[1] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[1]);
    file->idx->reg_patch_size[2] = getPowerOf2(file->idx->variable_grp[gi]->variable[svi]->sim_patch[0]->size[2]);
  }
  else if (file->idx->reg_box_set == PIDX_USER_RST_BOX)
  {
    assert(file->idx->reg_patch_size[0] != 0);
    assert(file->idx->reg_patch_size[1] != 0);
    assert(file->idx->reg_patch_size[2] != 0);
  }
  else if (file->idx->reg_box_set == PIDX_BOX_PER_PROCESS)
  {
    int rst_case_type = 0;
    int patch_count = var0->sim_patch_count;
    int max_patch_count = 0;
    MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);
    if (max_patch_count > 1)
      rst_case_type = 1;

    restructure_tag:

    file->idx->reg_patch_size[0] = factor * getPowerOftwo(var0->sim_patch[0]->size[0]);
    file->idx->reg_patch_size[1] = factor * getPowerOftwo(var0->sim_patch[0]->size[1]);
    file->idx->reg_patch_size[2] = factor * getPowerOftwo(var0->sim_patch[0]->size[2]);

    int grp_count = 0;
    if (rst_case_type == 0)
      grp_count = calculate_patch_group_count_for_patch_per_process(file, gi, svi);
    else
      grp_count = calculate_patch_group_count_for_multi_patch_per_process(file, gi, svi, max_patch_count);
    int max_grp_count = 0;

    MPI_Allreduce(&grp_count, &max_grp_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm );
    if (max_grp_count > 1)
    {
      factor = factor * 2;
      goto restructure_tag;
    }
  }
  else if (file->idx->reg_box_set == PIDX_BOX_FROM_BITSTRING)
    set_reg_patch_size_from_bit_string(file);

  else if (file->idx->reg_box_set == PIDX_UNIFORMLY_DISTRIBUTED_BOX)
  {
    file->idx->number_processes[0] = ceil(file->idx->box_bounds[0] / file->idx->reg_patch_size[0]);
    file->idx->number_processes[1] = ceil(file->idx->box_bounds[1] / file->idx->reg_patch_size[1]);
    file->idx->number_processes[2] = ceil(file->idx->box_bounds[2] / file->idx->reg_patch_size[2]);

    file->idx->new_box_set = malloc(file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2] * sizeof(*file->idx->new_box_set));
    memset(file->idx->new_box_set, 0, (file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2] * sizeof(*file->idx->new_box_set)));

    int i = 0, j = 0, k = 0;
    for (i = 0; i < file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]; i++)
    {
      file->idx->new_box_set[i] = malloc(sizeof(*(file->idx->new_box_set[i])));
      memset(file->idx->new_box_set[i], 0, sizeof(*(file->idx->new_box_set[i])));

      file->idx->new_box_set[i]->rank = -1;
    }

    int rank_count = 0;
    int index = 0;
    //int color = 0;
    int nx = 0, ny = 0, nz = 0;
    int int_x = file->idx_c->gnproc_x / file->idx->number_processes[0];
    int int_y = file->idx_c->gnproc_y / file->idx->number_processes[1];
    int int_z = file->idx_c->gnproc_z / file->idx->number_processes[2];

    for (k = 0; k < file->idx->box_bounds[2]; k = k + file->idx->reg_patch_size[2])
      for (j = 0; j < file->idx->box_bounds[1]; j = j + file->idx->reg_patch_size[1])
        for (i = 0; i < file->idx->box_bounds[0]; i = i + file->idx->reg_patch_size[0])
        {
          //Interior regular patches
          index = ((k / file->idx->reg_patch_size[2]) * file->idx->number_processes[0] * file->idx->number_processes[1]) + ((j / file->idx->reg_patch_size[1]) * file->idx->number_processes[0]) + (i / file->idx->reg_patch_size[0]);

          file->idx->new_box_set[index]->offset[0] = i;
          file->idx->new_box_set[index]->offset[1] = j;
          file->idx->new_box_set[index]->offset[2] = k;

          file->idx->new_box_set[index]->size[0] = file->idx->reg_patch_size[0];
          file->idx->new_box_set[index]->size[1] = file->idx->reg_patch_size[1];
          file->idx->new_box_set[index]->size[2] = file->idx->reg_patch_size[2];

          file->idx->new_box_set[index]->edge = 1;

          //Edge regular patches
          if ((i + file->idx->reg_patch_size[0] /*+ 1*/) > file->idx->box_bounds[0])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[0] = file->idx->box_bounds[0] - i;
          }

          if ((j + file->idx->reg_patch_size[1] /*+ 1*/) > file->idx->box_bounds[1])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[1] = file->idx->box_bounds[1] - j;
          }

          if ((k + file->idx->reg_patch_size[2] /*+ 1*/) > file->idx->box_bounds[2])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[2] = file->idx->box_bounds[2] - k;
          }

          file->idx->new_box_set[index]->rank = rank_count * (file->idx_c->gnprocs / (file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]));
          rank_count++;

          //if (file->idx_c->grank == file->idx->new_box_set[index]->rank)
          //  color = 1;
        }


    if (file->idx_c->gnproc_x != -1 && file->idx_c->gnproc_y != -1 && file->idx_c->gnproc_z != -1)
    {
      for (nz = 0; nz < file->idx_c->gnproc_z; nz = nz + int_z)
        for (ny = 0; ny < file->idx_c->gnproc_y; ny = ny + int_y)
          for (nx = 0; nx < file->idx_c->gnproc_x; nx = nx + int_x)
          {
            index = ((nz / int_z) * file->idx->number_processes[0] * file->idx->number_processes[1]) + ((ny / int_y) * file->idx->number_processes[0]) + (nx / int_x);

            //file->idx->new_box_set[index]->rank = (nz * file->idx->number_processes[0] * file->idx->number_processes[1]) + (ny * file->idx->number_processes[0]) + nx;
            file->idx->new_box_set[index]->rank = (nz * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (ny * file->idx_c->gnproc_x) + nx;
            if (file->idx_c->grank == file->idx->new_box_set[index]->rank)
            {
              file->idx_c->grank_x = nx;
              file->idx_c->grank_y = ny;
              file->idx_c->grank_z = nz;
            }
          }
    }

    /*
    if (file->idx_c->grank == 0)
    {
      for (k = 0; k < file->idx->number_processes[2]; k++)
        for (j = 0; j < file->idx->number_processes[1]; j++)
          for (i = 0; i < file->idx->number_processes[0]; i++)
          {
            index = (k * file->idx->number_processes[0] * file->idx->number_processes[1]) + (j * file->idx->number_processes[0]) + i;

            printf("[%d %d %d] [%d %d %d] Rank %d\n", i, j, k, file->idx_c->grank_x, file->idx_c->grank_y, file->idx_c->grank_z, file->idx->new_box_set[index]->rank);
          }
    }
    */

    //assert(rank_count == file->idx_c->gnprocs);
    assert(rank_count == file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]);

  }

  else if (file->idx->reg_box_set == PIDX_WAVELET_BOX)
  {
    file->idx->number_processes[0] = ceil(file->idx->box_bounds[0] / file->idx->reg_patch_size[0]);
    file->idx->number_processes[1] = ceil(file->idx->box_bounds[1] / file->idx->reg_patch_size[1]);
    file->idx->number_processes[2] = ceil(file->idx->box_bounds[2] / file->idx->reg_patch_size[2]);

    file->idx->new_box_set = malloc(file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2] * sizeof(*file->idx->new_box_set));
    memset(file->idx->new_box_set, 0, (file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2] * sizeof(*file->idx->new_box_set)));

    int i = 0, j = 0, k = 0;
    for (i = 0; i < file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]; i++)
    {
      file->idx->new_box_set[i] = malloc(sizeof(*(file->idx->new_box_set[i])));
      memset(file->idx->new_box_set[i], 0, sizeof(*(file->idx->new_box_set[i])));

      file->idx->new_box_set[i]->rank = -1;
    }

    int l = 0, left = 0, right = 0;
    for (l = 0; l < file->idx_d->wavelet_levels; l++)
    {
      left = left + (int)pow(2, l + 1);
      right = right + (int)pow(2, l);
    }


    int rank_count = 0;
    int index = 0;
    //int color = 0;
    int nx = 0, ny = 0, nz = 0;
    int int_x = file->idx_c->gnproc_x / file->idx->number_processes[0];
    int int_y = file->idx_c->gnproc_y / file->idx->number_processes[1];
    int int_z = file->idx_c->gnproc_z / file->idx->number_processes[2];

    for (k = 0; k < file->idx->box_bounds[2]; k = k + file->idx->reg_patch_size[2])
      for (j = 0; j < file->idx->box_bounds[1]; j = j + file->idx->reg_patch_size[1])
        for (i = 0; i < file->idx->box_bounds[0]; i = i + file->idx->reg_patch_size[0])
        {
          //Interior regular patches
          index = ((k / file->idx->reg_patch_size[2]) * file->idx->number_processes[0] * file->idx->number_processes[1]) + ((j / file->idx->reg_patch_size[1]) * file->idx->number_processes[0]) + (i / file->idx->reg_patch_size[0]);

          file->idx->new_box_set[index]->offset[0] = i;
          file->idx->new_box_set[index]->offset[1] = j;
          file->idx->new_box_set[index]->offset[2] = k;

          file->idx->new_box_set[index]->size[0] = file->idx->reg_patch_size[0];
          file->idx->new_box_set[index]->size[1] = file->idx->reg_patch_size[1];
          file->idx->new_box_set[index]->size[2] = file->idx->reg_patch_size[2];

          file->idx->new_box_set[index]->edge = 1;

          //Edge regular patches
          if ((i + file->idx->reg_patch_size[0] /*+ 1*/) > file->idx->box_bounds[0])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[0] = file->idx->box_bounds[0] - i;
          }

          if ((j + file->idx->reg_patch_size[1] /*+ 1*/) > file->idx->box_bounds[1])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[1] = file->idx->box_bounds[1] - j;
          }

          if ((k + file->idx->reg_patch_size[2] /*+ 1*/) > file->idx->box_bounds[2])
          {
            file->idx->new_box_set[index]->edge = 2;
            file->idx->new_box_set[index]->size[2] = file->idx->box_bounds[2] - k;
          }

          file->idx->new_box_set[index]->rank = rank_count * (file->idx_c->gnprocs / (file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]));
          rank_count++;

          //if (file->idx_c->grank == file->idx->new_box_set[index]->rank)
          //  color = 1;
        }

    if (file->idx_c->gnproc_x != -1 && file->idx_c->gnproc_y != -1 && file->idx_c->gnproc_z != -1)
    {
      for (nz = 0; nz < file->idx_c->gnproc_z; nz = nz + int_z)
        for (ny = 0; ny < file->idx_c->gnproc_y; ny = ny + int_y)
          for (nx = 0; nx < file->idx_c->gnproc_x; nx = nx + int_x)
          {
            index = ((nz / int_z) * file->idx->number_processes[0] * file->idx->number_processes[1]) + ((ny / int_y) * file->idx->number_processes[0]) + (nx / int_x);

            file->idx->new_box_set[index]->rank = (nz * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (ny * file->idx_c->gnproc_x) + nx;

            // TODO: Check this again
            if (nx != 0)
              file->idx->new_box_set[index]->size_nx = left;
            if (nx != file->idx_c->gnproc_x - int_x)
              file->idx->new_box_set[index]->size_px = right;

            if (ny != 0)
              file->idx->new_box_set[index]->size_ny = left;
            if (ny != file->idx_c->gnproc_y - int_y)
              file->idx->new_box_set[index]->size_py = right;

            if (nz != 0)
              file->idx->new_box_set[index]->size_nz = left;
            if (nz != file->idx_c->gnproc_z - int_z)
              file->idx->new_box_set[index]->size_pz = right;

            //if (file->idx_c->grank == 5)
            //  printf("[INIT R %d] : %d %d -- %d %d -- %d %d\n", file->idx->new_box_set[index]->rank, file->idx->new_box_set[index]->size_nx, file->idx->new_box_set[index]->size_px, file->idx->new_box_set[index]->size_ny, file->idx->new_box_set[index]->size_py, file->idx->new_box_set[index]->size_nz, file->idx->new_box_set[index]->size_pz );

            if (file->idx_c->grank == file->idx->new_box_set[index]->rank)
            {
              file->idx_c->grank_x = nx;
              file->idx_c->grank_y = ny;
              file->idx_c->grank_z = nz;

              file->idx_d->w_nx = file->idx->new_box_set[index]->size_nx;
              file->idx_d->w_px = file->idx->new_box_set[index]->size_px;

              file->idx_d->w_ny = file->idx->new_box_set[index]->size_ny;
              file->idx_d->w_py = file->idx->new_box_set[index]->size_py;

              file->idx_d->w_nz = file->idx->new_box_set[index]->size_nz;
              file->idx_d->w_pz = file->idx->new_box_set[index]->size_pz;
            }
          }
    }

    /*
    if (file->idx_c->grank == 0)
    {
      for (k = 0; k < file->idx->number_processes[2]; k++)
        for (j = 0; j < file->idx->number_processes[1]; j++)
          for (i = 0; i < file->idx->number_processes[0]; i++)
          {
            index = (k * file->idx->number_processes[0] * file->idx->number_processes[1]) + (j * file->idx->number_processes[0]) + i;

            printf("[%d %d %d] [%d %d %d] Rank %d\n", i, j, k, file->idx_c->grank_x, file->idx_c->grank_y, file->idx_c->grank_z, file->idx->new_box_set[index]->rank);
          }
    }
    */

    //assert(rank_count == file->idx_c->gnprocs);
    assert(rank_count == file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]);

  }

  return PIDX_success;
}

PIDX_return_code free_rst_box_size(PIDX_io file)
{
  if (file->idx->reg_box_set == PIDX_WAVELET_BOX || file->idx->reg_box_set == PIDX_UNIFORMLY_DISTRIBUTED_BOX)
  {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < file->idx->number_processes[0] * file->idx->number_processes[1] * file->idx->number_processes[2]; i++)
      free(file->idx->new_box_set[i]);
    free(file->idx->new_box_set);
  }
  return PIDX_success;
}




static int calculate_patch_group_count_for_multi_patch_per_process(PIDX_io file, int gi, int svi, unsigned long long sim_max_patch_group_count)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  int j = 0;
  int patch_group_count = 0;
  unsigned long long* sim_multi_patch_r_count;
  unsigned long long* sim_multi_patch_r_offset;

  int r, d, c;
  unsigned long long i, k, max_vol, patch_count, pc;
  int reg_patch_count, edge_case = 0;
  int reg_multi_patch_grp_count;
  Ndim_multi_patch_group* reg_multi_patch_grp;

  sim_multi_patch_r_count = malloc(sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count);
  memset(sim_multi_patch_r_count, 0, (sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count));
  sim_multi_patch_r_offset = malloc(sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count);
  memset(sim_multi_patch_r_offset, 0, (sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count));

  for(pc=0; pc < var0->sim_patch_count; pc++)
  {
    unsigned long long* tempoff = var0->sim_patch[pc]->offset;
    unsigned long long* tempsize = var0->sim_patch[pc]->size;

    unsigned long long index = file->idx_c->grank * (PIDX_MAX_DIMENSIONS * sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
    unsigned long long* curr_patch_offset = &sim_multi_patch_r_offset[index];
    unsigned long long* curr_patch_size = &sim_multi_patch_r_count[index];

    memcpy(curr_patch_offset, tempoff,sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
    memcpy(curr_patch_size, tempsize,sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  }

  unsigned long long* count_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*count_buffer_copy));
  memset(count_buffer_copy, 0, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*count_buffer_copy));

  memcpy(count_buffer_copy, &sim_multi_patch_r_count[file->idx_c->grank * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*count_buffer_copy));

  MPI_Allgather(count_buffer_copy, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count, MPI_LONG_LONG, sim_multi_patch_r_count, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count, MPI_LONG_LONG, file->idx_c->global_comm);
  free(count_buffer_copy);

  unsigned long long* offset_buffer_copy = malloc(PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*offset_buffer_copy));
  memset(offset_buffer_copy, 0, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*offset_buffer_copy));

  memcpy(offset_buffer_copy, &sim_multi_patch_r_offset[file->idx_c->grank * PIDX_MAX_DIMENSIONS * sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*sim_max_patch_group_count * sizeof(*offset_buffer_copy));

  MPI_Allgather(offset_buffer_copy, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count, MPI_LONG_LONG, sim_multi_patch_r_offset, PIDX_MAX_DIMENSIONS*sim_max_patch_group_count, MPI_LONG_LONG, file->idx_c->global_comm);
  free(offset_buffer_copy);

  /// extents for the local process(file->idx_c->grank)
  unsigned long long adjusted_bounds[PIDX_MAX_DIMENSIONS];
  memcpy(adjusted_bounds, file->idx->bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  int max_found_reg_patches = 1;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
#if 0
    adjusted_bounds[d] = file->idx->bounds[d];
    if (file->idx->bounds[d] % file->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((file->idx->bounds[d] / file->idx->chunk_size[d]) + 1) * file->idx->chunk_size[d];
    max_found_reg_patches *= ceil((float)file->idx->bounds[d]/(float)file->idx->reg_patch_size[d]);
#endif
    adjusted_bounds[d] = file->idx->box_bounds[d];
    if (file->idx->box_bounds[d] % file->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((file->idx->box_bounds[d] / file->idx->chunk_size[d]) + 1) * file->idx->chunk_size[d];
    max_found_reg_patches *= ceil((float)file->idx->box_bounds[d]/(float)file->idx->reg_patch_size[d]);
  }

  reg_multi_patch_grp_count = 0;

  unsigned long long pc0 = 0, d0 = 0;
  Ndim_patch* found_reg_patches = malloc(sizeof(Ndim_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(Ndim_patch*)*max_found_reg_patches);

  int found_reg_patches_count = 0;
  for(pc0 = 0; pc0 < var0->sim_patch_count; pc0++)
  {
    Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
    {
      local_proc_patch->offset[d0] = var0->sim_patch[pc0]->offset[d0];
      local_proc_patch->size[d0] = var0->sim_patch[pc0]->size[d0];
    }

    for (i = 0; i < adjusted_bounds[0]; i = i + file->idx->reg_patch_size[0])
      for (j = 0; j < adjusted_bounds[1]; j = j + file->idx->reg_patch_size[1])
        for (k = 0; k < adjusted_bounds[2]; k = k + file->idx->reg_patch_size[2])
        {
          Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
          memset(reg_patch, 0, sizeof (*reg_patch));

          //Interior regular patches
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = file->idx->reg_patch_size[0];
          reg_patch->size[1] = file->idx->reg_patch_size[1];
          reg_patch->size[2] = file->idx->reg_patch_size[2];

          //Edge regular patches
          if ((i + file->idx->reg_patch_size[0]) > adjusted_bounds[0])
            reg_patch->size[0] = adjusted_bounds[0] - i;
          if ((j + file->idx->reg_patch_size[1]) > adjusted_bounds[1])
            reg_patch->size[1] = adjusted_bounds[1] - j;
          if ((k + file->idx->reg_patch_size[2]) > adjusted_bounds[2])
            reg_patch->size[2] = adjusted_bounds[2] - k;

          if (intersectNDChunk(reg_patch, local_proc_patch))
          {
            if(!contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
            {
              found_reg_patches[found_reg_patches_count] = (Ndim_patch)malloc(sizeof (*reg_patch));
              memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));

              found_reg_patches_count++;
              reg_multi_patch_grp_count++;
            }
          }
          free(reg_patch);
        }
    free(local_proc_patch);
  }

  for(i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  found_reg_patches = malloc(sizeof(Ndim_patch*)*max_found_reg_patches);
  memset(found_reg_patches, 0, sizeof(Ndim_patch*)*max_found_reg_patches);

  found_reg_patches_count = 0;

  reg_multi_patch_grp = (Ndim_multi_patch_group*)malloc(sizeof(*reg_multi_patch_grp) * reg_multi_patch_grp_count);
  memset(reg_multi_patch_grp, 0, sizeof(*reg_multi_patch_grp) * reg_multi_patch_grp_count);

  reg_patch_count = 0;


  /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with

  for (i = 0; i < adjusted_bounds[0]; i = i + file->idx->reg_patch_size[0])
    for (j = 0; j < adjusted_bounds[1]; j = j + file->idx->reg_patch_size[1])
      for (k = 0; k < adjusted_bounds[2]; k = k + file->idx->reg_patch_size[2])
      {
        Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));

        //Interior regular patches
        reg_patch->offset[0] = i;
        reg_patch->offset[1] = j;
        reg_patch->offset[2] = k;
        reg_patch->size[0] = file->idx->reg_patch_size[0];
        reg_patch->size[1] = file->idx->reg_patch_size[1];
        reg_patch->size[2] = file->idx->reg_patch_size[2];

        //Edge regular patches
        edge_case = 0;
        if ((i + file->idx->reg_patch_size[0]) > adjusted_bounds[0])
        {
          reg_patch->size[0] = adjusted_bounds[0] - i;
          edge_case = 1;
        }
        if ((j + file->idx->reg_patch_size[1]) > adjusted_bounds[1])
        {
          reg_patch->size[1] = adjusted_bounds[1] - j;
          edge_case = 1;
        }
        if ((k + file->idx->reg_patch_size[2]) > adjusted_bounds[2])
        {
          reg_patch->size[2] = adjusted_bounds[2] - k;
          edge_case = 1;
        }

        for(pc0 = 0; pc0 < var0->sim_patch_count; pc0++)
        {
          Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
          memset(local_proc_patch, 0, sizeof (*local_proc_patch));
          for (d0 = 0; d0 < PIDX_MAX_DIMENSIONS; d0++)
          {
            local_proc_patch->offset[d0] = var0->sim_patch[pc0]->offset[d0];
            local_proc_patch->size[d0] = var0->sim_patch[pc0]->size[d0];
          }

          /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
          if (intersectNDChunk(reg_patch, local_proc_patch) && !contains_patch(reg_patch, found_reg_patches, found_reg_patches_count))
          {
            found_reg_patches[found_reg_patches_count] = (Ndim_patch)malloc(sizeof (*reg_patch));
            memcpy(found_reg_patches[found_reg_patches_count], reg_patch, sizeof (*reg_patch));
            found_reg_patches_count++;

            reg_multi_patch_grp[reg_patch_count] = malloc(sizeof(*(reg_multi_patch_grp[reg_patch_count])));
            memset(reg_multi_patch_grp[reg_patch_count], 0, sizeof(*(reg_multi_patch_grp[reg_patch_count])));

            Ndim_multi_patch_group patch_grp = reg_multi_patch_grp[reg_patch_count];

            patch_grp->source_patch = (PIDX_source_patch_index*)malloc(sizeof(PIDX_source_patch_index) * maximum_neighbor_count);
            patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
            patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
            memset(patch_grp->source_patch, 0, sizeof(PIDX_source_patch_index) * maximum_neighbor_count);
            memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
            memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

            patch_count = 0;
            patch_grp->count = 0;
            if(edge_case == 0)
              patch_grp->type = 1;
            else
              patch_grp->type = 2;

            //Iterate through all processes
            for (r = 0; r < file->idx_c->gnprocs; r++)
            {
              for(pc = 0; pc < sim_max_patch_group_count; pc++)
              {
                //Extent of process with rank r
                Ndim_patch curr_patch = malloc(sizeof (*curr_patch));
                memset(curr_patch, 0, sizeof (*curr_patch));

                for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                {
                  unsigned long long index = r * (PIDX_MAX_DIMENSIONS * sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
                  curr_patch->offset[d] = sim_multi_patch_r_offset[index+d];
                  curr_patch->size[d] = sim_multi_patch_r_count[index+d];
                }

                if(curr_patch->size[0] == 0)
                {
                  free(curr_patch);
                  continue;
                }

                if (intersectNDChunk(reg_patch, curr_patch))
                {
                  patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
                  memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

                  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                  {
                    //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
                    if (curr_patch->offset[d] <= reg_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                    {
                      patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = (curr_patch->offset[d] + curr_patch->size[d] - 1) - reg_patch->offset[d] + 1;
                    }
                    else if (reg_patch->offset[d] <= curr_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                    {
                      patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - curr_patch->offset[d] + 1;
                    }
                    else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (curr_patch->offset[d] + curr_patch->size[d] - 1) && reg_patch->offset[d] >= curr_patch->offset[d])
                    {
                      patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
                    }
                    else if (( curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && curr_patch->offset[d] >= reg_patch->offset[d])
                    {
                      patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                      patch_grp->patch[patch_count]->size[d] = curr_patch->size[d];
                    }

                    //offset and count of intersecting regular patch
                    patch_grp->reg_patch->offset[d] = reg_patch->offset[d];
                    patch_grp->reg_patch->size[d] = reg_patch->size[d];
                  }

                  patch_grp->source_patch[patch_count].rank = r;
                  patch_grp->source_patch[patch_count].index = pc;
                  patch_count++;

                  if (patch_count >= maximum_neighbor_count)
                  {
                    maximum_neighbor_count = maximum_neighbor_count * 2;

                    PIDX_source_patch_index *temp_buffer2 = realloc(patch_grp->source_patch, maximum_neighbor_count * sizeof(PIDX_source_patch_index));
                    if (temp_buffer2 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->source_patch = temp_buffer2;

                    Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
                    if (temp_buffer3 == NULL)
                    {
                      fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                      return PIDX_err_rst;
                    }
                    else
                      patch_grp->patch = temp_buffer3;

                    if (file->idx_c->grank == 0)
                      printf("[ERROR] maximum_neighbor_count needs to be increased\n");
                  }

                  patch_grp->count = patch_count;
                }
                free(curr_patch);
              }
            }


            patch_grp->max_patch_rank = patch_grp->source_patch[0].rank;
            max_vol = 1;
            for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              max_vol = max_vol * patch_grp->patch[0]->size[d];

            unsigned long long c_vol = 1;
            for(c = 1; c < patch_grp->count ; c++)
            {
              c_vol = 1;
              for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                c_vol = c_vol * patch_grp->patch[c]->size[d];
              if(c_vol > max_vol)
              {
                max_vol = c_vol;
                patch_grp->max_patch_rank = patch_grp->source_patch[c].rank;
              }
            }

            if(file->idx_c->grank == patch_grp->max_patch_rank)
              patch_group_count = patch_group_count + 1;

            reg_patch_count++;
          }

          free(local_proc_patch);
        }
        free(reg_patch);
      }

  for(i=0; i<found_reg_patches_count; i++)
  {
    free(found_reg_patches[i]);
    found_reg_patches[i] = 0;
  }
  free(found_reg_patches);

  free(sim_multi_patch_r_offset);
  free(sim_multi_patch_r_count);


  for (i = 0; i < reg_multi_patch_grp_count; i++)
  {
    for (j = 0; j < reg_multi_patch_grp[i]->count ; j++ )
    {
      free(reg_multi_patch_grp[i]->patch[j]);
      reg_multi_patch_grp[i]->patch[j] = 0;
    }

    free(reg_multi_patch_grp[i]->source_patch);
    reg_multi_patch_grp[i]->source_patch = 0;

    free(reg_multi_patch_grp[i]->patch);
    reg_multi_patch_grp[i]->patch = 0;

    free(reg_multi_patch_grp[i]->reg_patch);
    reg_multi_patch_grp[i]->reg_patch = 0;

    free(reg_multi_patch_grp[i]);
    reg_multi_patch_grp[i] = 0;
  }

  free(reg_multi_patch_grp);
  reg_multi_patch_grp = 0;

  return patch_group_count;
}



static int calculate_patch_group_count_for_patch_per_process(PIDX_io file, int gi, int svi)
{
  unsigned long long *rank_r_offset;
  unsigned long long *rank_r_count;
  int patch_group_count = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  Ndim_patch_group* reg_patch_grp;
  int j = 0;
  int r, d, c;
  unsigned long long i, k, max_vol, patch_count;
  int reg_patch_count, edge_case = 0;


  patch_group_count = 0;

  rank_r_offset = malloc(sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_offset, 0, (sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS));

  rank_r_count =  malloc(sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS);
  memset(rank_r_count, 0, (sizeof (unsigned long long) * file->idx_c->gnprocs * PIDX_MAX_DIMENSIONS));

  MPI_Allgather(var0->sim_patch[0]->offset , PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rank_r_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->global_comm);

  MPI_Allgather(var0->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, rank_r_count, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->global_comm);

  /// extents for the local process(rank)
  Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
  memset(local_proc_patch, 0, sizeof (*local_proc_patch));
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    local_proc_patch->offset[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * file->idx_c->grank + d];
    local_proc_patch->size[d] = rank_r_count[PIDX_MAX_DIMENSIONS * file->idx_c->grank + d];
  }

  unsigned long long adjusted_bounds[PIDX_MAX_DIMENSIONS];
  memcpy(adjusted_bounds, file->idx->bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));

  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
#if 0
    adjusted_bounds[d] = file->idx->bounds[d];
    if (file->idx->bounds[d] % file->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((file->idx->bounds[d] / file->idx->chunk_size[d]) + 1) * file->idx->chunk_size[d];
#endif
    adjusted_bounds[d] = file->idx->box_bounds[d];
    if (file->idx->box_bounds[d] % file->idx->chunk_size[d] != 0)
      adjusted_bounds[d] = ((file->idx->box_bounds[d] / file->idx->chunk_size[d]) + 1) * file->idx->chunk_size[d];
  }

  int reg_patch_grp_count = 0;
  for (i = 0; i < adjusted_bounds[0]; i = i + file->idx->reg_patch_size[0])
    for (j = 0; j < adjusted_bounds[1]; j = j + file->idx->reg_patch_size[1])
      for (k = 0; k < adjusted_bounds[2]; k = k + file->idx->reg_patch_size[2])
      {
        Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));

        //Interior regular patches
        reg_patch->offset[0] = i;
        reg_patch->offset[1] = j;
        reg_patch->offset[2] = k;
        reg_patch->size[0] = file->idx->reg_patch_size[0];
        reg_patch->size[1] = file->idx->reg_patch_size[1];
        reg_patch->size[2] = file->idx->reg_patch_size[2];

        //Edge regular patches
        if ((i + file->idx->reg_patch_size[0]) > adjusted_bounds[0])
          reg_patch->size[0] = adjusted_bounds[0] - i;
        if ((j + file->idx->reg_patch_size[1]) > adjusted_bounds[1])
          reg_patch->size[1] = adjusted_bounds[1] - j;
        if ((k + file->idx->reg_patch_size[2]) > adjusted_bounds[2])
          reg_patch->size[2] = adjusted_bounds[2] - k;

        if (intersectNDChunk(reg_patch, local_proc_patch))
          reg_patch_grp_count++;

        free(reg_patch);
      }

  reg_patch_grp = (Ndim_patch_group*)malloc(sizeof(*reg_patch_grp) * reg_patch_grp_count);
  memset(reg_patch_grp, 0, sizeof(*reg_patch_grp) * reg_patch_grp_count);

  reg_patch_count = 0;
  /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with

  for (i = 0; i < adjusted_bounds[0]; i = i + file->idx->reg_patch_size[0])
    for (j = 0; j < adjusted_bounds[1]; j = j + file->idx->reg_patch_size[1])
      for (k = 0; k < adjusted_bounds[2]; k = k + file->idx->reg_patch_size[2])
      {
        Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
        memset(reg_patch, 0, sizeof (*reg_patch));

        //Interior regular patches
        reg_patch->offset[0] = i;
        reg_patch->offset[1] = j;
        reg_patch->offset[2] = k;
        reg_patch->size[0] = file->idx->reg_patch_size[0];
        reg_patch->size[1] = file->idx->reg_patch_size[1];
        reg_patch->size[2] = file->idx->reg_patch_size[2];

        //Edge regular patches
        edge_case = 0;
        if ((i + file->idx->reg_patch_size[0]) > adjusted_bounds[0])
        {
          reg_patch->size[0] = adjusted_bounds[0] - i;
          edge_case = 1;
        }
        if ((j + file->idx->reg_patch_size[1]) > adjusted_bounds[1])
        {
          reg_patch->size[1] = adjusted_bounds[1] - j;
          edge_case = 1;
        }
        if ((k + file->idx->reg_patch_size[2]) > adjusted_bounds[2])
        {
          reg_patch->size[2] = adjusted_bounds[2] - k;
          edge_case = 1;
        }

        /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
        if (intersectNDChunk(reg_patch, local_proc_patch))
        {
          reg_patch_grp[reg_patch_count] = malloc(sizeof(*(reg_patch_grp[reg_patch_count])));
          memset(reg_patch_grp[reg_patch_count], 0, sizeof(*(reg_patch_grp[reg_patch_count])));

          Ndim_patch_group patch_grp = reg_patch_grp[reg_patch_count];

          patch_grp->source_patch_rank = (int*)malloc(sizeof(int) * maximum_neighbor_count);
          patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
          patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
          memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
          memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
          memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));

          patch_count = 0;
          patch_grp->count = 0;
          if(edge_case == 0)
            patch_grp->type = 1;
          else
            patch_grp->type = 2;

          //Iterate through all processes
          for (r = 0; r < file->idx_c->gnprocs; r++)
          {
            //Extent of process with rank r
            Ndim_patch rank_r_patch = malloc(sizeof (*rank_r_patch));
            memset(rank_r_patch, 0, sizeof (*rank_r_patch));

            for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            {
              rank_r_patch->offset[d] = rank_r_offset[PIDX_MAX_DIMENSIONS * r + d];
              rank_r_patch->size[d] = rank_r_count[PIDX_MAX_DIMENSIONS * r + d];
            }

            //If process with rank r intersects with the regular patch, then calculate the offset, count and volume of the intersecting volume
            if (intersectNDChunk(reg_patch, rank_r_patch))
            {
              patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
              memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));

              for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              {
                //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
                if (rank_r_patch->offset[d] <= reg_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                {
                  patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                  patch_grp->patch[patch_count]->size[d] = (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) - reg_patch->offset[d] + 1;
                }
                else if (reg_patch->offset[d] <= rank_r_patch->offset[d] && (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                {
                  patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
                  patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - rank_r_patch->offset[d] + 1;
                }
                else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) && reg_patch->offset[d] >= rank_r_patch->offset[d])
                {
                  patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                  patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
                }
                else if (( rank_r_patch->offset[d] + rank_r_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && rank_r_patch->offset[d] >= reg_patch->offset[d])
                {
                  patch_grp->patch[patch_count]->offset[d] = rank_r_patch->offset[d];
                  patch_grp->patch[patch_count]->size[d] = rank_r_patch->size[d];
                }

                //offset and count of intersecting regular patch
                patch_grp->reg_patch->offset[d] = reg_patch->offset[d];
                patch_grp->reg_patch->size[d] = reg_patch->size[d];
              }

              patch_grp->source_patch_rank[patch_count] = r;
              patch_count++;

              if (patch_count >= maximum_neighbor_count)
              {
                maximum_neighbor_count = maximum_neighbor_count * 2;

                int *temp_buffer2 = realloc(patch_grp->source_patch_rank, maximum_neighbor_count * sizeof(int));
                if (temp_buffer2 == NULL)
                {
                  fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_rst;
                }
                else
                  patch_grp->source_patch_rank = temp_buffer2;

                Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
                if (temp_buffer3 == NULL)
                {
                  fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                  return PIDX_err_rst;
                }
                else
                  patch_grp->patch = temp_buffer3;

                if (file->idx_c->grank == 0)
                  printf("maximum_neighbor_count needs to be increased to %d\n", maximum_neighbor_count);
              }

              patch_grp->count = patch_count;
            }
            free(rank_r_patch);
          }

          patch_grp->max_patch_rank = patch_grp->source_patch_rank[0];
          max_vol = 1;
          for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
            max_vol = max_vol * patch_grp->patch[0]->size[d];
          unsigned long long c_vol = 1;
          for(c = 1; c < patch_grp->count ; c++)
          {
            c_vol = 1;
            for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
              c_vol = c_vol * patch_grp->patch[c]->size[d];
            if(c_vol > max_vol)
            {
              max_vol = c_vol;
              patch_grp->max_patch_rank = patch_grp->source_patch_rank[c];
            }
          }

          if(file->idx_c->grank == patch_grp->max_patch_rank)
            patch_group_count = patch_group_count + 1;
          reg_patch_count++;
        }
        free(reg_patch);
      }

  free(local_proc_patch);
  free(rank_r_offset);
  free(rank_r_count);

  for (i = 0; i < reg_patch_grp_count; i++)
  {
    for (j = 0; j < reg_patch_grp[i]->count ; j++ )
    {
      free(reg_patch_grp[i]->patch[j]);
      reg_patch_grp[i]->patch[j] = 0;
    }

    free(reg_patch_grp[i]->source_patch_rank);
    reg_patch_grp[i]->source_patch_rank = 0;

    free(reg_patch_grp[i]->patch);
    reg_patch_grp[i]->patch = 0;

    free(reg_patch_grp[i]->reg_patch);
    reg_patch_grp[i]->reg_patch = 0;

    free(reg_patch_grp[i]);
    reg_patch_grp[i] = 0;
  }

  free(reg_patch_grp);
  reg_patch_grp = 0;

  return patch_group_count;
}


PIDX_return_code idx_init(PIDX_io file, int gi, int svi)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  var_grp->rank_buffer = malloc(file->idx_c->gnprocs * sizeof(*var_grp->rank_buffer));
  memset(var_grp->rank_buffer, 0, file->idx_c->gnprocs * sizeof(*var_grp->rank_buffer));
  MPI_Allgather(&(file->idx_c->grank), 1, MPI_INT, var_grp->rank_buffer, 1, MPI_INT, file->idx_c->global_comm);
#if 0
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];

  MPI_Allgather(var0->sim_patch[0]->offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx->all_offset, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->local_comm);
  MPI_Allgather(var0->sim_patch[0]->size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx->all_size, PIDX_MAX_DIMENSIONS, MPI_UNSIGNED_LONG_LONG, file->idx_c->local_comm);
#endif
  return PIDX_success;
}



/// Function to find the power of 2 of an integer value (example 5->8)
static int getPowerOftwo(int x)
{
  int n = 1;
  while (n < x)
    n <<= 1;
  return n;
}


/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}


static PIDX_return_code set_reg_patch_size_from_bit_string(PIDX_io file)
{
  int bits = log2(getPowerOf2(file->idx_c->gnprocs));
  int counter = 1;
  unsigned long long power_two_bound[PIDX_MAX_DIMENSIONS];
  power_two_bound[0] = file->idx_d->partition_count[0] * file->idx_d->partition_size[0];
  power_two_bound[1] = file->idx_d->partition_count[1] * file->idx_d->partition_size[1];
  power_two_bound[2] = file->idx_d->partition_count[2] * file->idx_d->partition_size[2];

  memcpy(file->idx->reg_patch_size, power_two_bound, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  while (bits != 0)
  {
    if (file->idx->bitSequence[counter] == '0')
      file->idx->reg_patch_size[0] = file->idx->reg_patch_size[0] / 2;

    else if (file->idx->bitSequence[counter] == '1')
      file->idx->reg_patch_size[1] = file->idx->reg_patch_size[1] / 2;

    else if (file->idx->bitSequence[counter] == '2')
      file->idx->reg_patch_size[2] = file->idx->reg_patch_size[2] / 2;

    counter++;
    bits--;
  }

  return PIDX_success;
}

static int contains_patch(Ndim_patch reg_patch, Ndim_patch* patches, int count)
{
  int i=0;

  for(i=0; i<count; i++)
  {
    int d=0;
    int matches = 0;
    for(d=0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      if(reg_patch->offset[d] == patches[i]->offset[d] && reg_patch->size[d] == patches[i]->size[d])
      matches++;
    }

    if(matches == PIDX_MAX_DIMENSIONS)
      return 1;
  }

  return 0;
}
