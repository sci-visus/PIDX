#include "../../PIDX_inc.h"

static int intersectNDChunk(PIDX_patch A, PIDX_patch B);


PIDX_return_code partition_setup(PIDX_io file, int gi, int svi)
{
  int *colors;
  int index_i = 0, index_j = 0, index_k = 0;
  int i = 0, j = 0, k = 0, d = 0;

  colors = malloc(sizeof(*colors) * file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]);
  memset(colors, 0, sizeof(*colors) * file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]);
  file->idx_d->color = (file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * file->idx_d->partition_count[2]) + 1;

  for (k = 0; k < file->idx_d->partition_count[2]; k++)
    for (j = 0; j < file->idx_d->partition_count[1]; j++)
      for (i = 0; i < file->idx_d->partition_count[0]; i++)
      {
        colors[(file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * k) + (file->idx_d->partition_count[0] * j) + i] = (file->idx_d->partition_count[0] * file->idx_d->partition_count[1] * k) + (file->idx_d->partition_count[0] * j) + i;
      }

  PIDX_patch local_p = (PIDX_patch)malloc(sizeof (*local_p));
  memset(local_p, 0, sizeof (*local_p));

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[svi];

  if (var->restructured_super_patch_count == 1)
  {
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_p->offset[d] = var->restructured_super_patch->restructured_patch->offset[d];
      local_p->size[d] = var->restructured_super_patch->restructured_patch->size[d];
    }
    PIDX_patch reg_patch = (PIDX_patch)malloc(sizeof (*reg_patch));
    memset(reg_patch, 0, sizeof (*reg_patch));

    Point3D bounds_point;
    int maxH = 0;
    bounds_point.x = (int) file->idx_d->partition_count[0];
    bounds_point.y = (int) file->idx_d->partition_count[1];
    bounds_point.z = (int) file->idx_d->partition_count[2];
    char bitSequence[512];
    char bitPattern[512];
    GuessBitmaskPattern(bitSequence, bounds_point);
    maxH = strlen(bitSequence);

    for (i = 0; i <= maxH; i++)
      bitPattern[i] = RegExBitmaskBit(bitSequence, i);

    int z_order = 0;
    int number_levels = maxH - 1;

    for (i = 0, index_i = 0; i < file->idx->bounds[0]; i = i + file->idx_d->partition_size[0], index_i++)
    {
      for (j = 0, index_j = 0; j < file->idx->bounds[1]; j = j + file->idx_d->partition_size[1], index_j++)
      {
        for (k = 0, index_k = 0; k < file->idx->bounds[2]; k = k + file->idx_d->partition_size[2], index_k++)
        {
          reg_patch->offset[0] = i;
          reg_patch->offset[1] = j;
          reg_patch->offset[2] = k;
          reg_patch->size[0] = file->idx_d->partition_size[0];
          reg_patch->size[1] = file->idx_d->partition_size[1];
          reg_patch->size[2] = file->idx_d->partition_size[2];

         if (intersectNDChunk(reg_patch, local_p))
         {
            Point3D xyzuv_Index;
            xyzuv_Index.x = index_i;
            xyzuv_Index.y = index_j;
            xyzuv_Index.z = index_k;

            z_order = 0;
            Point3D zero;
            zero.x = 0;
            zero.y = 0;
            zero.z = 0;
            memset(&zero, 0, sizeof (Point3D));

            int cnt = 0;
            for (cnt = 0; memcmp(&xyzuv_Index, &zero, sizeof (Point3D)); cnt++, number_levels--)
            {
              int bit = bitPattern[number_levels];
              z_order |= ((unsigned long long) PGET(xyzuv_Index, bit) & 1) << cnt;
              PGET(xyzuv_Index, bit) >>= 1;
            }

            file->idx_d->color = colors[z_order];

            file->idx_d->partition_offset[0] = i;
            file->idx_d->partition_offset[1] = j;
            file->idx_d->partition_offset[2] = k;

            break;
          }
        }
      }
    }
    free(reg_patch);
  }
  else if (var->restructured_super_patch_count > 1)
    fprintf(stderr, "RST artifact %d\n", var->restructured_super_patch_count);

  free(colors);

  char file_name_skeleton[1024];
  strncpy(file_name_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  file_name_skeleton[strlen(file->idx->filename) - 4] = '\0';

  //if (file->idx_d->partition_count[0] != 1 || file->idx_d->partition_count[1] != 1 || file->idx_d->partition_count[2] != 1)
  sprintf(file->idx->filename_partition, "%s_%d.idx", file_name_skeleton, file->idx_d->color);
  //else
  //  strcpy(file->idx->filename_partition, file->idx->filename);

  free(local_p);

  return PIDX_success;
}


PIDX_return_code create_local_comm(PIDX_io file)
{
  int ret;
  ret = MPI_Comm_split(file->idx_c->rst_comm, file->idx_d->color, file->idx_c->rrank, &(file->idx_c->local_comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  MPI_Comm_size(file->idx_c->local_comm, &(file->idx_c->lnprocs));
  MPI_Comm_rank(file->idx_c->local_comm, &(file->idx_c->lrank));

  //printf("[%d %d] [%d %d]\n", file->idx_c->gnprocs, file->idx_c->grank, file->idx_c->lnprocs, file->idx_c->lrank);

  return PIDX_success;
}


PIDX_return_code destroy_local_comm(PIDX_io file)
{
  int ret;
  ret = MPI_Comm_free(&(file->idx_c->local_comm));
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  return PIDX_success;
}


PIDX_return_code find_partition_count(PIDX_io file)
{
  int d = 0;
  // calculate number of partitions
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
  {
    file->idx_d->partition_size[d] = file->idx->bounds[d] / file->idx_d->partition_count[d];
    if (file->idx->bounds[d] % file->idx_d->partition_count[d] != 0)
      file->idx_d->partition_size[d]++;

    file->idx_d->partition_size[d] = pow(2, (int)ceil(log2(file->idx_d->partition_size[d])));
  }

  unsigned long long dims;
  unsigned long long adjusted_bounds[PIDX_MAX_DIMENSIONS];
  adjusted_bounds[0] = file->idx_d->partition_size[0] / file->idx->chunk_size[0];
  adjusted_bounds[1] = file->idx_d->partition_size[1] / file->idx->chunk_size[1];
  adjusted_bounds[2] = file->idx_d->partition_size[2] / file->idx->chunk_size[2];

  if (PIDX_inner_product(&dims, adjusted_bounds))
    return PIDX_err_size;

  //printf("dims %d spb %d\n", dims, file->idx_d->samples_per_block);

  if (dims < file->idx_d->samples_per_block)
  {
    // ensure blocksize is a subset of the total volume.
    file->idx_d->samples_per_block = getPowerOf2(dims) >> 1;
    file->idx->bits_per_block = getNumBits(file->idx_d->samples_per_block) - 1;
  }

  if (file->idx->bits_per_block == 0)
  {
    file->idx->bits_per_block = 0;
    file->idx_d->samples_per_block = 1;
  }

  //printf("Rank %d Partition size %d %d %d\n", file->idx_c->grank, file->idx_d->partition_size[0], file->idx_d->partition_size[1], file->idx_d->partition_size[2]);

  return PIDX_success;
}


static int intersectNDChunk(PIDX_patch A, PIDX_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];

  return !(check_bit);
}
