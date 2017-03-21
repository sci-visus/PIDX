#include "../PIDX_inc.h"

#define FINAL_RESULT 0
#define INTERMEDIATE_RESULT 0

static MPI_Comm verifyComm;


static PIDX_return_code print_global_data (PIDX_io file, int gi, int v);
static PIDX_return_code create_debug_comm (PIDX_io file, int gi, int v);

static PIDX_return_code wavelet_rst_x (PIDX_io file, int gi, int v, int l, int nx);
static PIDX_return_code wavelet_rst_y (PIDX_io file, int gi, int v, int l, int ny);
static PIDX_return_code wavelet_rst_z (PIDX_io file, int gi, int v, int l, int nz);


PIDX_return_code rst_wavelet(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int v = 0;
  int l = 0;
  int i1 = 0, j1 = 0, k1 = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_time time = file->idx_d->time;

  for (v = svi; v <= evi; v++)
  {
    PIDX_variable var = var_grp->variable[v];

#if FINAL_RESULT || INTERMEDIATE_RESULT
    if (create_debug_comm(file, gi, v) != PIDX_success)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return PIDX_err_wavelet;
    }
#endif

    if (var->patch_group_count == 0)
      goto comm_cleanup;

    Ndim_patch_group patch_group = var->rst_patch_group[0];
    int nx = file->idx_d->w_nx;//patch_group->size_nx;
    int px = file->idx_d->w_px;//patch_group->size_px;
    int ny = file->idx_d->w_ny;//patch_group->size_ny;
    int py = file->idx_d->w_py;//patch_group->size_py;
    int nz = file->idx_d->w_nz;//patch_group->size_nz;
    int pz = file->idx_d->w_pz;//patch_group->size_pz;

    int s_x = (int) patch_group->reg_patch->size[0];
    int s_y = (int) patch_group->reg_patch->size[1];
    int s_z = (int) patch_group->reg_patch->size[2];

#if INTERMEDIATE_RESULT
    if (file->idx_c->grank == 1)
      printf("ORIGINAL [%d] SIZE %d %d %d\n", file->idx_c->grank, s_x, s_y, s_z);

    for (k1 = 0; k1 < s_z - 0; k1++)
    {
      for (j1 = 0; j1 < s_y - 0; j1++)
      {
        for (i1 = 0; i1 < s_x - 0; i1++)
        {
          if (var->bpv/8 == 4)
          {
            float x;
            int index = (s_x * s_y * k1) + (s_x * j1) + i1;
            memcpy(&x, patch_group->reg_patch->buffer + index * var->bpv/8, var->bpv/8);
            if (file->idx_c->grank == 1)
              printf("%3.3f\t", x);
          }
        }
        if (file->idx_c->grank == 1)
          printf("\n");
      }
    }
#endif

    for (l = 0; l < file->idx_d->wavelet_levels; l++)
    {
      //if (file->idx_c->grank == 0)
      //{
      //printf("nx px ny py - %d %d %d %d\n", nx, px, ny, py);
      if ((int) patch_group->reg_patch->size[0] >= (int)pow(2, l + 1))
      {
        time->w_rst_comp_x_start[gi][v][l] = MPI_Wtime();
        if (wavelet_rst_x (file, gi, v, l, nx) != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
        time->w_rst_comp_x_end[gi][v][l] = MPI_Wtime();
      }


      if ((int) patch_group->reg_patch->size[1] >= (int)pow(2, l + 1))
      {
        time->w_rst_comp_y_start[gi][v][l] = MPI_Wtime();
        if (wavelet_rst_y (file, gi, v, l, ny) != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
        time->w_rst_comp_y_end[gi][v][l] = MPI_Wtime();
      }


      if ((int) patch_group->reg_patch->size[2] >= (int)pow(2, l + 1))
      {
        time->w_rst_comp_z_start[gi][v][l] = MPI_Wtime();
        if (wavelet_rst_z (file, gi, v, l, nz) != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
        time->w_rst_comp_z_end[gi][v][l] = MPI_Wtime();
      }
      //}
    }


    patch_group->reg_patch->size[0] = s_x - px - nx;
    patch_group->reg_patch->size[1] = s_y - py - ny;
    patch_group->reg_patch->size[2] = s_z - pz - nz;


    patch_group->reg_patch->offset[0] = patch_group->reg_patch->offset[0] + nx;
    patch_group->reg_patch->offset[1] = patch_group->reg_patch->offset[1] + ny;
    patch_group->reg_patch->offset[2] = patch_group->reg_patch->offset[2] + nz;


  unsigned char* temp_buffer = malloc (patch_group->reg_patch->size[0] * patch_group->reg_patch->size[1] * patch_group->reg_patch->size[2] * var->bpv/8);


  int write_count = 0;
  for (k1 = nz; k1 < s_z - pz; k1++)
  {
    for (j1 = ny; j1 < s_y - py; j1++)
    {
      i1 = nx;
      int index = (s_x * s_y * k1) + (s_x * j1) + i1;
      memcpy(temp_buffer + write_count * var->bpv/8, patch_group->reg_patch->buffer + index * var->bpv/8, patch_group->reg_patch->size[0] * var->bpv/8);
      write_count = write_count + patch_group->reg_patch->size[0];
    }
  }

  /*
  if (file->idx_c->grank == 1)
    printf("SIZE %d (%d - %d - %d) %d (%d - %d - %d) %d (%d - %d - %d)\n", patch_group->reg_patch->size[0], s_x, px, nx ,patch_group->reg_patch->size[1], s_y, py, ny, patch_group->reg_patch->size[2], s_z, pz, nz);

  for (k1 = 0; k1 < patch_group->reg_patch->size[2]; k1++)
  {
    for (j1 = 0; j1 < patch_group->reg_patch->size[1]; j1++)
    {
      for (i1 = 0; i1 < patch_group->reg_patch->size[0]; i1++)
      {
        int index = (patch_group->reg_patch->size[0] * patch_group->reg_patch->size[1] * k1) + (patch_group->reg_patch->size[0] * j1) + i1;
        float x;
        memcpy (&x, temp_buffer + index * var->bpv/8, var->bpv/8);
        if (file->idx_c->grank == 1)
        printf("%3.3f\t", x);
      }
      if (file->idx_c->grank == 1)
      printf("\n");
    }
  }
  */

  unsigned char *temp_buffer2 = realloc(patch_group->reg_patch->buffer, patch_group->reg_patch->size[0] * patch_group->reg_patch->size[1] * patch_group->reg_patch->size[2] * var->bpv/8);
  if (temp_buffer2 == NULL)
  {
    fprintf(stderr, "[%s] [%d] realloc() failed %d %d %d.\n", __FILE__, __LINE__, (int)patch_group->reg_patch->size[0], (int)patch_group->reg_patch->size[1], (int)patch_group->reg_patch->size[2]);
    return PIDX_err_rst;
  }
  else
    patch_group->reg_patch->buffer = temp_buffer2;


  memcpy(patch_group->reg_patch->buffer, temp_buffer, patch_group->reg_patch->size[0] * patch_group->reg_patch->size[1] * patch_group->reg_patch->size[2] * var->bpv/8);
  free(temp_buffer);

#if FINAL_RESULT
  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif



comm_cleanup:
    ;
#if FINAL_RESULT || INTERMEDIATE_RESULT
    MPI_Comm_free(&verifyComm);
#endif
  }

    return PIDX_success;
}




static PIDX_return_code create_debug_comm (PIDX_io file, int gi, int v)
{
  int ret;
  int color = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];

  if (var->patch_group_count == 0)
    color = 0;
  else
    color = 1;

  ret = MPI_Comm_split(file->idx_c->global_comm, color, file->idx_c->grank, &verifyComm);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }

  return PIDX_success;
}




static PIDX_return_code wavelet_rst_x (PIDX_io file, int gi, int v, int l, int nx)
{
  int stride = (int)pow(2, l);
  int odd_start_offset = (int)pow(2, l);
  int even_start_offset = 0;
  int interval = (int)pow(2,l);
  int stride_x = (int)pow(2, l + 1);
  int index = 0;
  int i = 0, j = 0, k = 0;
#if 0
  if (nx != 0)
  {
  //if (l % 2 == 0)
  if (l == 0)
  {
    odd_start_offset = (int)pow(2,l);
    even_start_offset = 0;
  }
  else
  {
    even_start_offset = (int)pow(2,l);
    odd_start_offset = 0;
  }
  }
#endif

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  for (i = odd_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[0] + (int) var->rst_patch_group[0]->reg_patch->size[0] - 1; i = i + stride_x)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[0])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[0]);
      odd_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[0];
      break;
    }
  }

  for (i = even_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[0] + (int) var->rst_patch_group[0]->reg_patch->size[0] - 1; i = i + stride_x)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[0])
    {
      //if (file->idx_c->grank == 1)
      //  printf("EVEN start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[0]);
      even_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[0];
      break;
    }
  }

  int y_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[1] + (int) var->rst_patch_group[0]->reg_patch->size[1] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[1])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      y_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[1];
      break;
    }
  }
  int z_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[2] + (int) var->rst_patch_group[0]->reg_patch->size[2] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[2])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[2]);
      z_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[2];
      break;
    }
  }

  //Ndim_patch_group patch_group = var->rst_patch_group[0];
  /*
  int nx = patch_group->size_nx;
  int px = patch_group->size_px;
  int ny = patch_group->size_ny;
  int py = patch_group->size_py;
  int nz = patch_group->size_nz;
  int pz = patch_group->size_pz;
  */

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  for (k = z_start_offset; k < s_z - 0; k = k + stride)
  {
    for (j = y_start_offset; j < s_y - 0; j =  j + stride)
    {
      // every sample but the last
      for (i = odd_start_offset ; i < s_x; i = i + stride_x)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           if (i - interval < 0)
             continue;

           float left, right, new_val;
           memcpy (&left, wb + (index - interval) * bytes_for_datatype, bytes_for_datatype);

           if (i + interval >= s_x)
             memcpy (&right, &left, bytes_for_datatype);
           else
             memcpy (&right, wb + (index + interval) * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           //if (file->idx_c->grank == 1)
           //  printf("[O] %f - %f - %f\n", left, new_val, right);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

#if INTERMEDIATE_RESULT
  int i1 = 0, j1 = 0, k1 = 0;
  if (file->idx_c->grank == 1)
    printf("X ODD [%d] SIZE %d %d %d\n", file->idx_c->grank, s_x, s_y, s_z);

  Ndim_patch_group patch_group = var->rst_patch_group[0];
  for (k1 = 0; k1 < s_z - 0; k1++)
  {
    for (j1 = 0; j1 < s_y - 0; j1++)
    {
      for (i1 = 0; i1 < s_x - 0; i1++)
      {
        if (var->bpv/8 == 4)
        {
          float x;
          int index = (s_x * s_y * k1) + (s_x * j1) + i1;
          memcpy(&x, patch_group->reg_patch->buffer + index * var->bpv/8, var->bpv/8);
          if (file->idx_c->grank == 1)
            printf("%3.3f\t", x);
        }
      }
      if (file->idx_c->grank == 1)
        printf("\n");
    }
  }
#endif
  //


  //
  //printf("Size = %d %d %d\n", s_x, s_y, s_z);
  for (k = z_start_offset; k < s_z - 0; k = k + stride)
  {
    for (j = y_start_offset; j < s_y - 0; j =  j + stride)
    {
      // every sample but the last
      for (i = even_start_offset ; i < s_x; i = i + stride_x)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;

         if (i + interval >= s_x)
           continue;

         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;

           //printf("[%d] [%d %d %d] index %d interval %d\n", file->idx_c->grank, i, j, k, index, interval);
           memcpy (&right, wb + (index + interval) * bytes_for_datatype, bytes_for_datatype);

           if (i - interval < 0)
             memcpy (&left, &right, bytes_for_datatype);
           else
             memcpy (&left, wb + (index - interval) * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           //if (file->idx_c->grank == 1 && l == 1)
           //  printf("[E index %d interval %d i %d] %f - %f - %f\n", index, interval, i, left, new_val, right);
           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

#if INTERMEDIATE_RESULT
  if (file->idx_c->grank == 1)
    printf("X EVEN ORIGINAL [%d] SIZE %d %d %d\n", file->idx_c->grank, s_x, s_y, s_z);

  for (k1 = 0; k1 < s_z - 0; k1++)
  {
    for (j1 = 0; j1 < s_y - 0; j1++)
    {
      for (i1 = 0; i1 < s_x - 0; i1++)
      {
        if (var->bpv/8 == 4)
        {
          float x;
          int index = (s_x * s_y * k1) + (s_x * j1) + i1;
          memcpy(&x, patch_group->reg_patch->buffer + index * var->bpv/8, var->bpv/8);
          if (file->idx_c->grank == 1)
            printf("%3.3f\t", x);
        }
      }
      if (file->idx_c->grank == 1)
        printf("\n");
    }
  }
#endif

  return PIDX_success;
}



static PIDX_return_code wavelet_rst_y (PIDX_io file, int gi, int v, int l, int ny)
{
  int stride = (int)pow(2, l);
  int odd_start_offset = (int)pow(2, l);
  int even_start_offset = 0;
  int interval = (int)pow(2,l);
  int stride_y = (int)pow(2, l + 1);
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  for (i = odd_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[1] + (int) var->rst_patch_group[0]->reg_patch->size[1] - 1; i = i + stride_y)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[1])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      odd_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[1];
      break;
    }
  }

  for (i = even_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[1] + (int) var->rst_patch_group[0]->reg_patch->size[1] - 1; i = i + stride_y)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[1])
    {
      //if (file->idx_c->grank == 1)
      //  printf("EVEN start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      even_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[1];
      break;
    }
  }

  int x_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[0] + (int) var->rst_patch_group[0]->reg_patch->size[0] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[0])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      x_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[0];
      break;
    }
  }

  int z_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[2] + (int) var->rst_patch_group[0]->reg_patch->size[2] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[2])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[2]);
      z_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[2];
      break;
    }
  }

#if 0
  if (ny != 0)
  {
  //if (l % 2 == 0)
  if (l == 0)
  {
    odd_start_offset = (int)pow(2,l);
    even_start_offset = 0;
  }
  else
  {
    even_start_offset = (int)pow(2,l);
    odd_start_offset = 0;
  }
  }
#endif

  //Ndim_patch_group patch_group = var->rst_patch_group[0];
  /*
  int nx = patch_group->size_nx;
  int px = patch_group->size_px;
  int ny = patch_group->size_ny;
  int py = patch_group->size_py;
  int nz = patch_group->size_nz;
  int pz = patch_group->size_pz;
  */

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  for (k = z_start_offset; k < s_z - 0; k = k + stride)
  {
    for (i = x_start_offset ; i < s_x - 0; i = i + stride)
    {
      // every sample but the last
      for (j = odd_start_offset; j < s_y; j =  j + stride_y)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;

         if (j - interval < 0)
           continue;

         int up_index = (s_x * s_y * k) + (s_x * (j + interval)) + i;
         int down_index = (s_x * s_y * k) + (s_x * (j - interval)) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + down_index * bytes_for_datatype, bytes_for_datatype);

           if (j + interval >= s_y)
             memcpy (&right, &left, bytes_for_datatype);
           else
             memcpy (&right, wb + up_index * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           //if (file->idx_c->grank == 1)
           //  printf("%f -- %f -- %f\n", left, new_val, right);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

#if INTERMEDIATE_RESULT
  int i1 = 0, j1 = 0, k1 = 0;
  if (file->idx_c->grank == 1)
    printf("Y ODD [%d] SIZE %d %d %d\n", file->idx_c->grank, s_x, s_y, s_z);

  Ndim_patch_group patch_group = var->rst_patch_group[0];
  for (k1 = 0; k1 < s_z - 0; k1++)
  {
    for (j1 = 0; j1 < s_y - 0; j1++)
    {
      for (i1 = 0; i1 < s_x - 0; i1++)
      {
        if (var->bpv/8 == 4)
        {
          float x;
          int index = (s_x * s_y * k1) + (s_x * j1) + i1;
          memcpy(&x, patch_group->reg_patch->buffer + index * var->bpv/8, var->bpv/8);
          if (file->idx_c->grank == 1)
            printf("%3.3f\t", x);
        }
      }
      if (file->idx_c->grank == 1)
        printf("\n");
    }
  }
#endif


  for (k = z_start_offset; k < s_z - 0; k = k + stride)
  {
    for (i = x_start_offset ; i < s_x - 0; i = i + stride)
    {
      // every sample but the last
      for (j = even_start_offset; j < s_y; j =  j + stride_y)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;

         if (j + interval >= s_y)
           continue;

         int up_index = (s_x * s_y * k) + (s_x * (j + interval)) + i;
         int down_index = (s_x * s_y * k) + (s_x * (j - interval)) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&right, wb + up_index * bytes_for_datatype, bytes_for_datatype);

           if (j - interval < 0)
             memcpy (&left, &right, bytes_for_datatype);
           else
             memcpy (&left, wb + down_index * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           //if (file->idx_c->grank == 1 && l == 2)
           //  printf("[E index %d interval %d i %d] %f - %f - %f\n", index, interval, j, left, new_val, right);

           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

#if INTERMEDIATE_RESULT
  if (file->idx_c->grank == 1)
    printf("Y EVEN [%d] SIZE %d %d %d\n", file->idx_c->grank, s_x, s_y, s_z);

  if (file->idx_c->grank == 1)
  {
  for (k1 = 0; k1 < s_z - 0; k1++)
  {
    for (j1 = 0; j1 < s_y - 0; j1++)
    {
      for (i1 = 0; i1 < s_x - 0; i1++)
      {
        if (var->bpv/8 == 4)
        {
          float x;
          int index = (s_x * s_y * k1) + (s_x * j1) + i1;
          memcpy(&x, patch_group->reg_patch->buffer + index * var->bpv/8, var->bpv/8);

            printf("%3.3f\t", x);
        }
      }
      if (file->idx_c->grank == 1)
        printf("\n");
    }
  }
  }
#endif
//
  return PIDX_success;
}


static PIDX_return_code wavelet_rst_z (PIDX_io file, int gi, int v, int l, int nz)
{
  int stride = (int)pow(2, l);
  int odd_start_offset = (int)pow(2, l);
  int even_start_offset = 0;
  int interval = (int)pow(2,l);
  int stride_z = (int)pow(2, l + 1);
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

#if 0
  if (nz != 0)
  {
  //if (l % 2 == 0)
  if (l == 0)
  {
    odd_start_offset = (int)pow(2,l);
    even_start_offset = 0;
  }
  else
  {
    even_start_offset = (int)pow(2,l);
    odd_start_offset = 0;
  }
  }
#endif

  for (i = odd_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[2] + (int) var->rst_patch_group[0]->reg_patch->size[2] - 1; i = i + stride_z)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[2])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[0]);
      odd_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[2];
      break;
    }
  }

  for (i = even_start_offset; i < (int) var->rst_patch_group[0]->reg_patch->offset[2] + (int) var->rst_patch_group[0]->reg_patch->size[2] - 1; i = i + stride_z)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[2])
    {
      //if (file->idx_c->grank == 1)
      //  printf("EVEN start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[0]);
      even_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[2];
      break;
    }
  }

  int x_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[0] + (int) var->rst_patch_group[0]->reg_patch->size[0] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[0])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      x_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[0];
      break;
    }
  }

  int y_start_offset = 0;
  for (i = 0; i < (int) var->rst_patch_group[0]->reg_patch->offset[1] + (int) var->rst_patch_group[0]->reg_patch->size[1] - 1; i = i + stride)
  {
    if (i >= (int) var->rst_patch_group[0]->reg_patch->offset[1])
    {
      //if (file->idx_c->grank == 1)
      //  printf("ODD start index at level %d = %d - %d\n", l, i, (int)var->rst_patch_group[0]->reg_patch->offset[1]);
      y_start_offset = i - (int)var->rst_patch_group[0]->reg_patch->offset[1];
      break;
    }
  }

  //Ndim_patch_group patch_group = var->rst_patch_group[0];
  /*
  int nx = patch_group->size_nx;
  int px = patch_group->size_px;
  int ny = patch_group->size_ny;
  int py = patch_group->size_py;
  int nz = patch_group->size_nz;
  int pz = patch_group->size_pz;
  */


  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  for (j = y_start_offset; j < s_y - 0; j = j + stride)
  {
    for (i = x_start_offset ; i < s_x - 0; i = i + stride)
    {
      // every sample but the last
      for (k = odd_start_offset ; k < s_z; k =  k + stride_z)
      {
         if (k - interval < 0)
           continue;
         index = (s_x * s_y * k) + (s_x * j) + i;
         int up_index = (s_x * s_y * (k + interval)) + (s_x * j) + i;
         int down_index = (s_x * s_y * (k - interval)) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           //if (file->idx_c->grank == 1)
           //  printf ("[%d %d] -> %d %d %d\n", i, j, index, up_index, down_index);
           float left, right, new_val;
           memcpy (&left, wb + down_index * bytes_for_datatype, bytes_for_datatype);

           if (k + interval >= s_z)
             memcpy (&right, &left, bytes_for_datatype);
           else
             memcpy (&right, wb + up_index * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

  for (j = y_start_offset; j < s_y - 0; j = j + stride)
  {
    for (i = x_start_offset ; i < s_x - 0; i = i + stride)
    {
      // every sample but the last
      for (k = even_start_offset ; k < s_z; k =  k + stride_z)
      {
         if (k + interval >= s_z)
           continue;
         index = (s_x * s_y * k) + (s_x * j) + i;
         int up_index = (s_x * s_y * (k + interval)) + (s_x * j) + i;
         int down_index = (s_x * s_y * (k - interval)) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           //if (file->idx_c->grank == 1)
           //  printf ("[%d %d] -> %d %d %d\n", i, j, index, up_index, down_index);
           float left, right, new_val;
           memcpy (&right, wb + up_index * bytes_for_datatype, bytes_for_datatype);

           if (k - interval < 0)
             memcpy (&left, &right, bytes_for_datatype);
           else
             memcpy (&left, wb + down_index * bytes_for_datatype, bytes_for_datatype);

           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }
    }
  }

  return PIDX_success;
}



static PIDX_return_code print_global_data (PIDX_io file, int gi, int v)
{
  int ret;
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];

  if (var->patch_group_count == 0)
    goto exit;

  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];

  int o_x = (int) var->rst_patch_group[0]->reg_patch->offset[0];
  int o_y = (int) var->rst_patch_group[0]->reg_patch->offset[1];
  int o_z = (int) var->rst_patch_group[0]->reg_patch->offset[2];

  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  MPI_Win win;
  unsigned char* global_buffer;

  if (file->idx_c->grank == 0)
  {
    int global_buffer_size = file->idx->bounds[0] * file->idx->bounds[1] * file->idx->bounds[2] * bytes_for_datatype;
    global_buffer = malloc(global_buffer_size);
    memset(global_buffer, 0, global_buffer_size);

    ret = MPI_Win_create(global_buffer, global_buffer_size, bytes_for_datatype, MPI_INFO_NULL, verifyComm, &win);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " Error in Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_wavelet;
    }
  }
  else
  {
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, verifyComm, &win);
    if (ret != MPI_SUCCESS)
    {
      fprintf(stderr, " Error in Line %d File %s\n", __LINE__, __FILE__);
      return PIDX_err_wavelet;
    }
  }

  ret = MPI_Win_fence(0, win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " Error in Line %d File %s\n", __LINE__, __FILE__);
    return PIDX_err_wavelet;
  }

  int i = 0, j = 0, k = 0;
  for (k = 0; k < s_z; k++)
  {
    for (j = 0; j < s_y; j++)
    {
      for (i = 0; i < s_x; i++)
      {
        int index = (s_x * s_y * k) + (s_x * j) + i;
        float val;
        int target_disp = ((file->idx->bounds[0] * file->idx->bounds[1]*(o_z + k))+(file->idx->bounds[0]*(o_y + j)) + (o_x + i));

        memcpy(&val, wb + index * bytes_for_datatype, bytes_for_datatype);
        ret = MPI_Put(&val, bytes_for_datatype, MPI_BYTE, 0, target_disp, bytes_for_datatype, MPI_BYTE, win);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, " Error in MPI_Put Line %d File %s\n", __LINE__, __FILE__);
          return PIDX_err_agg;
        }
      }
    }
  }

  ret = MPI_Win_fence(0, win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " Error in Line %d File %s\n", __LINE__, __FILE__);
    return PIDX_err_wavelet;
  }

  ret = MPI_Win_free(&win);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, " Error in Line %d File %s\n", __LINE__, __FILE__);
    return PIDX_err_wavelet;
  }

#if 1
  if (file->idx_c->grank == 0)
  {
    printf("\n");
    for (k = 0; k < file->idx->bounds[2]; k++)
    {
      for (j = 0; j < file->idx->bounds[1]; j++)
      //for (j = file->idx->bounds[1] -1; j >= 0 ; j--)
      {
        for (i = 0; i < file->idx->bounds[0]; i++)
        {
          int index = (file->idx->bounds[0] * file->idx->bounds[1] * k) + (file->idx->bounds[0] * j) + i;
          float v;
          memcpy(&v, global_buffer + index * bytes_for_datatype, bytes_for_datatype);
          printf("%3.3f\t", v);
        }
        printf("\n");
      }
    }
    printf("\n");
  }
#endif

  if (file->idx_c->grank == 0)
    free(global_buffer);

  exit:
    ;

  return PIDX_success;
}
