#include "../PIDX_inc.h"

#define SORIGINAL 0
#define SINTERMEDIATE_RESULT 0
#define FINAL_RESULT 0


static MPI_Comm verifyComm;

static int sent_x = 0, receive_x = 0;
static int sent_y = 0, receive_y = 0;
static int sent_z = 0, receive_z = 0;

static int positive_rank_x = -1, negative_rank_x = -1;
static unsigned char *positive_buffer_x, *negative_buffer_x;

static int positive_rank_y = -1, negative_rank_y = -1;
static unsigned char *positive_buffer_y, *negative_buffer_y;

static int positive_rank_z = -1, negative_rank_z = -1;
static unsigned char *positive_buffer_z, *negative_buffer_z;

static PIDX_return_code print_global_data (PIDX_io file, int gi, int v);
static PIDX_return_code create_debug_comm (PIDX_io file, int gi, int v);
//static PIDX_return_code create_wavelet_buffers (PIDX_io file, int gi, int v);

static PIDX_return_code calculate_neighbor_ranks (PIDX_io file);
static PIDX_return_code calculate_neighbor_ranks_x (PIDX_io file);
static PIDX_return_code calculate_neighbor_ranks_y (PIDX_io file);
static PIDX_return_code calculate_neighbor_ranks_z (PIDX_io file);

static PIDX_return_code create_stencil_buffers (PIDX_io file, int gi, int v, int level);
static PIDX_return_code destroy_stencil_buffers();

static PIDX_return_code wavelet_odd_x (PIDX_io file, int gi, int v, int level, int x_stride_x, int x_stride_y, int x_stride_z);
static PIDX_return_code wavelet_comm_p2p_odd_x (PIDX_io file, int gi, int v, int x_stride_y, int x_stride_z);
static PIDX_return_code wavelet_comp_odd_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z);

static PIDX_return_code wavelet_even_x (PIDX_io file, int gi, int v, int level, int x_stride_x, int x_stride_y, int x_stride_z);
static PIDX_return_code wavelet_comm_p2p_even_x (PIDX_io file, int gi, int v, int l, int x_stride_y, int x_stride_z);
static PIDX_return_code wavelet_comp_even_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z);

static PIDX_return_code wavelet_odd_y (PIDX_io file, int gi, int v, int level, int y_stride_x, int y_stride_y, int y_stride_z);
static PIDX_return_code wavelet_comm_p2p_odd_y (PIDX_io file, int gi, int v, int y_stride_y, int y_stride_z);
static PIDX_return_code wavelet_comp_odd_y (PIDX_io file, int gi, int v, int l, int y_stride_x, int y_stride_y, int y_stride_z);

static PIDX_return_code wavelet_even_y (PIDX_io file, int gi, int v, int level, int y_stride_x, int y_stride_y, int y_stride_z);
static PIDX_return_code wavelet_comm_p2p_even_y (PIDX_io file, int gi, int v, int level, int y_stride_y, int y_stride_z);
static PIDX_return_code wavelet_comp_even_y (PIDX_io file, int gi, int v, int l, int y_stride_x, int y_stride_y, int y_stride_z);

static PIDX_return_code wavelet_odd_z (PIDX_io file, int gi, int v, int level, int stride_x, int stride_y, int stride_z);
static PIDX_return_code wavelet_comm_p2p_odd_z (PIDX_io file, int gi, int v, int stride_y, int stride_z);
static PIDX_return_code wavelet_comp_odd_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z);

static PIDX_return_code wavelet_even_z (PIDX_io file, int gi, int v, int level, int stride_x, int stride_y, int stride_z);
static PIDX_return_code wavelet_comm_p2p_even_z (PIDX_io file, int gi, int v, int level, int stride_y, int stride_z);
static PIDX_return_code wavelet_comp_even_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z);


PIDX_return_code compute_average(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int v = 0, p = 0, b = 0, i = 0, j = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (v = svi; v <= evi; v++)
  //for (v = 1; v < 2; v++)
  {
    PIDX_variable var = var_grp->variable[v];

    for (p = 0; p < var->patch_group_count; p++)
    {
      for (b = 0; b < var->chunk_patch_group[p]->count; b++)
      {
        Ndim_patch patch = var->chunk_patch_group[p]->patch[b];
        unsigned char* buffer = patch->buffer;
        int element_count = patch->size[0] * patch->size[1] * patch->size[2] * var->vps;
        unsigned char* temp = malloc(element_count / 64 * sizeof (float));
        int count = 0;
        float sum = 0;
        float average = 0;

        for (i = 0; i < element_count; i = i + 64)
        {
          sum = 0;
          for (j = 0; j < 64; j++)
          {
            float value;
            //value = ((float*) buffer)[i + j];
            memcpy(&value, buffer + (i + j) * sizeof (float), sizeof(float));
            sum = sum + value;
            //sum = sum + ((float*) buffer)[i + j];
          }
          average = sum / 64;

          //if (file->idx_c->grank == 1)
          //  printf("Average = %f\n", average);

          memcpy(temp + count * sizeof (float), &average, sizeof (float));
          count++;
        }
        memcpy(patch->buffer, temp, element_count / 64 * sizeof (float));

        unsigned char* temp_buffer = realloc(patch->buffer, count* sizeof (float));
        if (temp_buffer == NULL)
          return PIDX_err_compress;
        else
          patch->buffer = temp_buffer;
      }
    }
  }
  return PIDX_success;
}


PIDX_return_code idx_stencil_wavelet(PIDX_io file, int gi, int svi, int evi, int mode)
{
#if 1
  int ret = 0;
  int v = 0;
  int l = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (v = svi; v <= evi; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    assert(var->patch_group_count <= 1);

    if (file->idx->current_time_step == 0)
    {
      ret = create_debug_comm(file, gi, v);
      if (ret != PIDX_success)
      {
        fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
        return PIDX_err_wavelet;
      }
    }

    if (var->patch_group_count == 0)
      goto comm_cleanup;

#if 0
    ret = create_wavelet_buffers(file, gi, v);
    if (ret != PIDX_success)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return PIDX_err_wavelet;
    }
#endif

#if SORIGINAL
    if (file->idx_c->grank == 0)
      printf("SORIGINAL");

    ret = print_global_data(file, gi, v);
    if (ret != PIDX_success)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return PIDX_err_wavelet;
    }
#endif
#if 1
    for (l = 0; l < file->idx_d->wavelet_levels; l++)
    {
      sent_x = sent_y = sent_z = 0;
      receive_x = receive_y = receive_z = 0;

      ret = calculate_neighbor_ranks(file);
      if (ret != PIDX_success)
      {
        fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
        return PIDX_err_wavelet;
      }

      ret = create_stencil_buffers(file, gi, v, l);
      if (ret != PIDX_success)
      {
        fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
        return PIDX_err_wavelet;
      }

      int x_stride_x, x_stride_y, x_stride_z;
      int y_stride_x, y_stride_y, y_stride_z;
      int z_stride_x, z_stride_y, z_stride_z;
      if ((int) var->rst_patch_group[0]->reg_patch->size[0] >= (int)pow(2, l + 1))
      {
        x_stride_x = (int)pow(2, l + 1);
        x_stride_y = (int)pow(2, l);
        x_stride_z = (int)pow(2, l);
        ret = wavelet_odd_x (file, gi, v, l, x_stride_x, x_stride_y, x_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }

        ret = wavelet_even_x (file, gi, v, l, x_stride_x, x_stride_y, x_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
      }

      y_stride_x = x_stride_x;
      y_stride_y = (int)pow(2, l + 1);
      y_stride_z = (int)pow(2, l);
      if ((int) var->rst_patch_group[0]->reg_patch->size[1] >= (int)pow(2, l + 1))
      {
        ret = wavelet_odd_y (file, gi, v, l, y_stride_x, y_stride_y, y_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }

        ret = wavelet_even_y (file, gi, v, l, y_stride_x, y_stride_y, y_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
      }

      z_stride_x = y_stride_x;
      z_stride_y = y_stride_y;
      z_stride_z = (int)pow(2, l + 1);
      if ((int) var->rst_patch_group[0]->reg_patch->size[2] >= (int)pow(2, l + 1))
      {
        ret = wavelet_odd_z (file, gi, v, l, z_stride_x, z_stride_y, z_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }

        ret = wavelet_even_z (file, gi, v, l, z_stride_x, z_stride_y, z_stride_z);
        if (ret != PIDX_success)
        {
          fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
          return PIDX_err_wavelet;
        }
      }
#if 1
      if (file->idx->current_time_step == 0)
      {
        int vrank;
        int total_receive = 0, total_sent = 0;
        int local_receive = receive_x + receive_y + receive_z;
        int local_sent = sent_x + sent_y + sent_z;
        MPI_Comm_rank(verifyComm, &vrank);

        MPI_Allreduce(&local_receive, &total_receive, 1, MPI_INT, MPI_SUM, verifyComm);
        MPI_Allreduce(&local_sent, &total_sent, 1, MPI_INT, MPI_SUM, verifyComm);

        if (vrank == 0)
          printf("[%d] Total Data sent %d Total Data received %d\n", l, total_sent, total_receive);
      }
#endif
      //printf("[%d] Rank %d Sent %d [%d %d %d] Receive %d [%d %d %d]\n", l, file->idx_c->grank, sent_x + sent_y + sent_z, sent_x, sent_y, sent_z, receive_x + receive_y + receive_z, receive_x, receive_y, receive_z);
      destroy_stencil_buffers();
    }
#endif
#if FINAL_RESULT
    ret = print_global_data(file, gi, v);
    if (ret != PIDX_success)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return PIDX_err_wavelet;
    }
#endif


    comm_cleanup:
      ;
      if (file->idx->current_time_step == 0)
        MPI_Comm_free(&verifyComm);



  }


#endif
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


#if 0
static PIDX_return_code create_wavelet_buffers (PIDX_io file, int gi, int v)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];

  var->rst_patch_group[0]->wavelet_reg_patch = malloc(sizeof(*(var->rst_patch_group[0]->wavelet_reg_patch)));

  memcpy(var->rst_patch_group[0]->wavelet_reg_patch->offset, var->rst_patch_group[0]->reg_patch->offset, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);
  memcpy(var->rst_patch_group[0]->wavelet_reg_patch->size, var->rst_patch_group[0]->reg_patch->size, sizeof(unsigned long long) * PIDX_MAX_DIMENSIONS);

  int s_x = (int) var->rst_patch_group[0]->wavelet_reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->wavelet_reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->wavelet_reg_patch->size[2];

  var->rst_patch_group[0]->wavelet_reg_patch->buffer = malloc(s_x * s_y * s_z * (var->bpv/8) * var->vps);
  if (var->rst_patch_group[0]->wavelet_reg_patch->buffer == NULL)
  {
    fprintf(stderr, "[%s] [%d] malloc() failed for buffer size %d (%d x %d x %d x %d x %d).\n", __FILE__, __LINE__, (int) (s_x * s_y * s_z * (var->bpv/8) * var->vps), (int)s_x, (int)s_y, (int)s_z, (int)(var->bpv/8), (int)var->vps);
    return PIDX_err_chunk;
  }

  memcpy(var->rst_patch_group[0]->wavelet_reg_patch->buffer, var->rst_patch_group[0]->reg_patch->buffer, s_x * s_y * s_z * (var->bpv/8) * var->vps);


  return PIDX_success;
}
#endif


static PIDX_return_code calculate_neighbor_ranks (PIDX_io file)
{
  int ret = 0;

  ret = calculate_neighbor_ranks_x(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }

  ret = calculate_neighbor_ranks_y(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }

   ret = calculate_neighbor_ranks_z(file);
  if (ret != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }

  return PIDX_success;
}



static PIDX_return_code calculate_neighbor_ranks_x (PIDX_io file)
{
  int p_rank = 0, n_rank = 0;

  p_rank = file->idx_c->grank_x + file->idx_c->gnproc_x / file->idx->number_processes[0];
  if (p_rank >= file->idx_c->gnproc_x)
    p_rank = -1;

  n_rank = file->idx_c->grank_x - file->idx_c->gnproc_x / file->idx->number_processes[0];
  if (n_rank < 0)
    n_rank = -1;

  if (p_rank != -1)
    positive_rank_x = (file->idx_c->grank_z * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (file->idx_c->grank_y * file->idx_c->gnproc_x) + p_rank;
  else
    positive_rank_x = -1;

  if (n_rank != -1)
    negative_rank_x = (file->idx_c->grank_z * file->idx_c->gnproc_x * file->idx_c->gnproc_y) + (file->idx_c->grank_y * file->idx_c->gnproc_x) + n_rank;
  else
    negative_rank_x = -1;

  return PIDX_success;
}




static PIDX_return_code calculate_neighbor_ranks_y (PIDX_io file)
{
  int p_rank = 0, n_rank = 0;

  p_rank = file->idx_c->grank_y + file->idx_c->gnproc_y / file->idx->number_processes[1];
  if (p_rank >= file->idx_c->gnproc_y)
    p_rank = -1;

  n_rank = file->idx_c->grank_y - file->idx_c->gnproc_y / file->idx->number_processes[1];
  if (n_rank < 0)
    n_rank = -1;

  if (p_rank != -1)
    positive_rank_y = (file->idx_c->grank_z * file->idx_c->gnproc_x * file->idx_c->gnproc_y) +
                      (p_rank * file->idx_c->gnproc_x) +
                       file->idx_c->grank_x;
  else
    positive_rank_y = -1;

  if (n_rank != -1)
    negative_rank_y = (file->idx_c->grank_z * file->idx_c->gnproc_x * file->idx_c->gnproc_y) +
                      (n_rank * file->idx_c->gnproc_x) +
                       file->idx_c->grank_x;
  else
    negative_rank_y = -1;

  //if (file->idx_c->grank == 0)
  //printf("My Rank %d My Up rank = %d My Down rank %d\n", file->idx_c->grank, positive_rank_y, negative_rank_y);

  return PIDX_success;
}



static PIDX_return_code calculate_neighbor_ranks_z (PIDX_io file)
{
  int p_rank = 0, n_rank = 0;

  p_rank = file->idx_c->grank_z + file->idx_c->gnproc_z / file->idx->number_processes[2];
  if (p_rank >= file->idx_c->gnproc_z)
    p_rank = -1;

  n_rank = file->idx_c->grank_z - file->idx_c->gnproc_z / file->idx->number_processes[2];
  if (n_rank < 0)
    n_rank = -1;

  if (p_rank != -1)
    positive_rank_z = (p_rank * file->idx_c->gnproc_x * file->idx_c->gnproc_y) +
                      (file->idx_c->grank_y * file->idx_c->gnproc_x) +
                       file->idx_c->grank_x;
  else
    positive_rank_z = -1;

  if (n_rank != -1)
    negative_rank_z = (n_rank * file->idx_c->gnproc_x * file->idx_c->gnproc_y) +
                      (file->idx_c->grank_y * file->idx_c->gnproc_x) +
                       file->idx_c->grank_x;
  else
    negative_rank_z = -1;

  return PIDX_success;
}


static PIDX_return_code create_stencil_buffers (PIDX_io file, int gi, int v, int l)
{
  int stride = (int)pow(2, l);
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];

  positive_buffer_x = malloc((int)ceil((double)s_z/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype);
  memset(positive_buffer_x, 0, (int)ceil((double)s_z/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype);

  negative_buffer_x = malloc((int)ceil((double)s_z/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype);
  memset(negative_buffer_x, 0, (int)ceil((double)s_z/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype);

  positive_buffer_y = malloc(((int)ceil((double)s_z/stride) * (int)ceil((double)s_x/stride) * bytes_for_datatype));
  memset(positive_buffer_y, 0, ((int)ceil((double)s_z/stride) * (int)ceil((double)s_x/stride) * bytes_for_datatype));

  negative_buffer_y = malloc(((int)ceil((double)s_z/stride) * (int)ceil((double)s_x/stride) * bytes_for_datatype));
  memset(negative_buffer_y, 0, ((int)ceil((double)s_z/stride) * (int)ceil((double)s_x/stride) * bytes_for_datatype));

  positive_buffer_z = malloc(((int)ceil((double)s_x/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype));
  memset(positive_buffer_z, 0, ((int)ceil((double)s_x/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype));

  negative_buffer_z = malloc(((int)ceil((double)s_x/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype));
  memset(negative_buffer_z, 0, ((int)ceil((double)s_x/stride) * (int)ceil((double)s_y/stride) * bytes_for_datatype));

  return PIDX_success;
}


static PIDX_return_code destroy_stencil_buffers()
{
  free(positive_buffer_x);
  free(negative_buffer_x);
  free(positive_buffer_y);
  free(negative_buffer_y);
  free(positive_buffer_z);
  free(negative_buffer_z);

  return PIDX_success;
}



static PIDX_return_code wavelet_odd_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_x_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_odd_x (file, gi, v, x_stride_y, x_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_x_odd_end[gi][v][l] = MPI_Wtime();


  time->w_stencil_comp_x_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_odd_x (file, gi, v, l, x_stride_x, x_stride_y, x_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_x_odd_end[gi][v][l] = MPI_Wtime();


#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet XXXX ---- 1", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif


  return PIDX_success;
}


static PIDX_return_code wavelet_even_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_x_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_even_x (file, gi, v, l, x_stride_y, x_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_x_even_end[gi][v][l] = MPI_Wtime();

  time->w_stencil_comp_x_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_even_x (file, gi, v, l, x_stride_x, x_stride_y, x_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_x_even_end[gi][v][l] = MPI_Wtime();

#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet XXXX ---- 2", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code wavelet_odd_y (PIDX_io file, int gi, int v, int l, int y_stride_x, int y_stride_y, int y_stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_y_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_odd_y (file, gi, v, y_stride_x, y_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_y_odd_end[gi][v][l] = MPI_Wtime();

  time->w_stencil_comp_y_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_odd_y (file, gi, v, l, y_stride_x, y_stride_y, y_stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_y_odd_end[gi][v][l] = MPI_Wtime();

#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet YYYY ---- 1", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code wavelet_even_y (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_y_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_even_y (file, gi, v, l, stride_x, stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_y_even_end[gi][v][l] = MPI_Wtime();

  time->w_stencil_comp_y_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_even_y (file, gi, v, l, stride_x, stride_y, stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_y_even_end[gi][v][l] = MPI_Wtime();

#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet YYYY ---- 2", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code wavelet_odd_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_z_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_odd_z (file, gi, v, stride_x, stride_y) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_z_odd_end[gi][v][l] = MPI_Wtime();

  time->w_stencil_comp_z_odd_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_odd_z (file, gi, v, l, stride_x, stride_y, stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_z_odd_end[gi][v][l] = MPI_Wtime();

#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet ZZZZ ---- 1", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code wavelet_even_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  PIDX_time time = file->idx_d->time;

  time->w_stencil_comm_z_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comm_p2p_even_z (file, gi, v, l, stride_x, stride_y) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comm_z_even_end[gi][v][l] = MPI_Wtime();

  time->w_stencil_comp_z_even_start[gi][v][l] = MPI_Wtime();
  if (wavelet_comp_even_z (file, gi, v, l, stride_x, stride_y, stride_z) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
  time->w_stencil_comp_z_even_end[gi][v][l] = MPI_Wtime();

#if SINTERMEDIATE_RESULT
  if (file->idx_c->grank == 0)
    printf("[%d] Wavelet ZZZZ ---- 2", l);

  if (print_global_data(file, gi, v) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif

  return PIDX_success;
}


static PIDX_return_code wavelet_comm_p2p_odd_x (PIDX_io file, int gi, int v, int x_stride_y, int x_stride_z)
{
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  //int stride = (int)pow(2, l);
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  //printf("My rank %d [%d + (%d/%d)] Rank I get data from %d Rank I send data to %d\n", file->idx_c->grank, file->idx_c->grank_x, file->idx_c->gnproc_x, file->idx->number_processes[0], positive_rank_x, negative_rank_x);

  int req_count = 0;
  if (positive_rank_x != -1)
  {
    receive_x = receive_x + ((int)ceil((double)s_z/x_stride_z) * (int)ceil((double)s_y/x_stride_y) * bytes_for_datatype);
    MPI_Irecv(positive_buffer_x, ((int)ceil((double)s_z/x_stride_z) * (int)ceil((double)s_y/x_stride_y) * bytes_for_datatype), MPI_BYTE, positive_rank_x, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  if (negative_rank_x != -1)
  {
    int count = 0;
    int lindex = 0;
    for (k = 0; k < s_z; k = k + x_stride_z)
    {
      for (j = 0; j < s_y; j =  j + x_stride_y)
      {
        i = 0;
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(negative_buffer_x + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_x = sent_x + ((int)ceil((double)s_z/x_stride_z) * (int)ceil((double)s_y/x_stride_y) * bytes_for_datatype);
    MPI_Isend(negative_buffer_x, ((int)ceil((double)s_z/x_stride_z) * (int)ceil((double)s_y/x_stride_y) * bytes_for_datatype), MPI_BYTE, negative_rank_x, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  int ret;
  MPI_Status status[2];
  ret = MPI_Waitall(req_count, request, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }

  return PIDX_success;
}



static PIDX_return_code wavelet_comp_odd_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z)
{
  int odd_start_offset = (int)pow(2, l);
  int interval = (int)pow(2,l);
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  int sample_count = 0;
  for (k = 0; k < s_z; k = k + x_stride_z)
  {
    for (j = 0; j < s_y; j =  j + x_stride_y)
    {
      // every sample but the last
      for (i = odd_start_offset ; i < s_x - x_stride_x; i = i + x_stride_x)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + (index - interval) * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + (index + interval) * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // last sample
      int last_index = (s_x * s_y * k) + (s_x * j) + (s_x - interval);
      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;

        // left
        memcpy (&left, wb + (last_index - interval) * bytes_for_datatype, bytes_for_datatype);

        // right
        if (positive_rank_x != -1)
          memcpy (&right, positive_buffer_x + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&right, wb + (last_index - interval) * bytes_for_datatype, bytes_for_datatype);

        // center
        memcpy (&new_val, wb + last_index * bytes_for_datatype, bytes_for_datatype);

        //if (file->idx_c->grank == 1)
        //  printf("%f %f %f\n", left, new_val, right);

        // new value
        new_val = new_val - 0.5 * (left + right);

        memcpy(wb + last_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
    }
  }

  return PIDX_success;
}



static PIDX_return_code wavelet_comm_p2p_even_x (PIDX_io file, int gi, int v, int l, int x_stride_y, int x_stride_z)
{
  int ret = 0;
  int count = 0;
  int lindex = 0;
  int req_count = 0;
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  req_count = 0;
  if (negative_rank_x != -1)
  {
    receive_x = receive_x + ((int)ceil((double)s_y/x_stride_y) * (int)ceil((double)s_z/x_stride_z) * bytes_for_datatype);
    MPI_Irecv(negative_buffer_x, ((int)ceil((double)s_y/x_stride_y) * (int)ceil((double)s_z/x_stride_z) * bytes_for_datatype), MPI_BYTE, negative_rank_x, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  if (positive_rank_x != -1)
  {
    count = 0;
    lindex = 0;
    for (k = 0; k < s_z; k = k + x_stride_z)
    {
      for (j = 0; j < s_y; j =  j + x_stride_y)
      {
        i = s_x - (int)pow(2, l);
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(positive_buffer_x + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_x = sent_x + ((int)ceil((double)s_y/x_stride_y) * (int)ceil((double)s_z/x_stride_z) * bytes_for_datatype);
    MPI_Isend(positive_buffer_x, ((int)ceil((double)s_y/x_stride_y) * (int)ceil((double)s_z/x_stride_z) * bytes_for_datatype), MPI_BYTE, positive_rank_x, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  MPI_Status status2[2];
  ret = MPI_Waitall(req_count, request, status2);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }

  return PIDX_success;
}


static PIDX_return_code wavelet_comp_even_x (PIDX_io file, int gi, int v, int l, int x_stride_x, int x_stride_y, int x_stride_z)
{
  int sample_count = 0;
  int interval = (int)pow(2,l);
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  // Even Samples
  sample_count = 0;
  for (k = 0; k < s_z; k = k + x_stride_z)
  {
    for (j = 0; j < s_y; j = j + x_stride_y)
    {
      for (i = x_stride_x; i < s_x; i = i + x_stride_x)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + (index - interval) * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + (index + interval) * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);
           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // first sample
      int first_index = (s_x * s_y * k) + (s_x * j);
      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;
        // left
        if (negative_rank_x != -1)
          memcpy (&left, negative_buffer_x + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&left, wb + (first_index + interval) * bytes_for_datatype, bytes_for_datatype);

        // right
        memcpy (&right, wb + (first_index + interval) * bytes_for_datatype, bytes_for_datatype);

        // existing value
        memcpy (&new_val, wb + first_index * bytes_for_datatype, bytes_for_datatype);

        // wavlet transformation
        new_val = new_val + 0.25 * (left + right);

        //if (file->idx_c->grank == 1)
        //  printf("%f %f %f\n", left, new_val, right);

        // new value
        memcpy(wb + first_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
    }
  }

  return PIDX_success;
}


static PIDX_return_code wavelet_comm_p2p_odd_y (PIDX_io file, int gi, int v, int y_stride_x, int y_stride_z)
{
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  //printf("My rank %d [%d + (%d/%d)] Rank I get data from %d Rank I send data to %d\n", file->idx_c->grank, file->idx_c->grank_x, file->idx_c->gnproc_x, file->idx->number_processes[0], positive_rank_x, negative_rank_x);

  int req_count = 0;
  if (positive_rank_y != -1)
  {
    receive_y = receive_y + ((int)ceil((double)s_z/y_stride_z) * (int)ceil((double)s_x/y_stride_x) * bytes_for_datatype);
    MPI_Irecv(positive_buffer_y, ((int)ceil((double)s_z/y_stride_z) * (int)ceil((double)s_x/y_stride_x) * bytes_for_datatype), MPI_BYTE, positive_rank_y, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }


#if 1
  int count = 0;
  int lindex = 0;
  if (negative_rank_y != -1)
  {
    for (k = 0; k < s_z; k = k + y_stride_z)
    {
      for (i = 0; i < s_x; i =  i + y_stride_x)
      {
        j = 0;
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(negative_buffer_y + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_y = sent_y + ((int)ceil((double)s_z/y_stride_z) * (int)ceil((double)s_x/y_stride_x) * bytes_for_datatype);
    MPI_Isend(negative_buffer_y, ((int)ceil((double)s_z/y_stride_z) * (int)ceil((double)s_x/y_stride_x) * bytes_for_datatype), MPI_BYTE, negative_rank_y, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  int ret;
  MPI_Status status[2];
  ret = MPI_Waitall(req_count, request, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif
  return PIDX_success;

}



static PIDX_return_code wavelet_comp_odd_y (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  int odd_start_offset = (int)pow(2, l);// * s_x;
  int interval = (int)pow(2,l);
  int up_index = 0;
  int down_index = 0;
  int index = 0;
  int i = 0, j = 0, k = 0;
  int sample_count = 0;

#if 1
  for (k = 0; k < s_z; k = k + stride_z)
  {
    for (i = 0 ; i < s_x; i = i + stride_x)
    {
      // every sample but the last
      for (j = odd_start_offset ; j < s_y - stride_y; j =  j + stride_y)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         up_index = (s_x * s_y * k) + (s_x * (j + interval)) + i;
         down_index = (s_x * s_y * k) + (s_x * (j - interval)) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + up_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + down_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // last sample
      int last_index = (s_x * s_y * k) + (s_x * (s_y - interval)) + i;
      int last_index_down = (s_x * s_y * k) + (s_x * (s_y - interval - interval)) + i;
      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;

        // left
        memcpy (&left, wb + last_index_down * bytes_for_datatype, bytes_for_datatype);

        // right
        if (positive_rank_y != -1)
          memcpy (&right, positive_buffer_y + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&right, wb + last_index_down * bytes_for_datatype, bytes_for_datatype);

        // center
        memcpy (&new_val, wb + last_index * bytes_for_datatype, bytes_for_datatype);

        // new value
        new_val = new_val - 0.5 * (left + right);

        memcpy(wb + last_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
    }
  }
#endif

  return PIDX_success;
}



static PIDX_return_code wavelet_comm_p2p_even_y (PIDX_io file, int gi, int v, int l, int stride_x, int stride_z)
{
  int ret = 0;
  int count = 0;
  int lindex = 0;
  int req_count = 0;
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  req_count = 0;
  if (negative_rank_y != -1)
  {
    receive_y = receive_y + ((int)ceil((double)s_z/stride_z) * (int)ceil((double)s_x/stride_x) * bytes_for_datatype);
    MPI_Irecv(negative_buffer_y, ((int)ceil((double)s_z/stride_z) * (int)ceil((double)s_x/stride_x) * bytes_for_datatype), MPI_BYTE, negative_rank_y, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  if (positive_rank_y != -1)
  {
    count = 0;
    lindex = 0;
    for (k = 0; k < s_z; k = k + stride_z)
    {
      for (i = 0; i < s_x; i =  i + stride_x)
      {
        j = s_y - (int)pow(2, l);
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(positive_buffer_y + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_y = sent_y + ((int)ceil((double)s_z/stride_z) * (int)ceil((double)s_x/stride_x) * bytes_for_datatype);
    MPI_Isend(positive_buffer_y, ((int)ceil((double)s_z/stride_z) * (int)ceil((double)s_x/stride_x) * bytes_for_datatype), MPI_BYTE, positive_rank_y, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  MPI_Status status2[2];
  ret = MPI_Waitall(req_count, request, status2);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }

  return PIDX_success;
}


static PIDX_return_code wavelet_comp_even_y (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  int sample_count = 0;
  int interval = (int)pow(2,l);
  int up_index = 0;
  int down_index = 0;
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  // Even Samples
  sample_count = 0;
  for (k = 0; k < s_z; k = k + stride_z)
  {
    for (i = 0; i < s_x; i = i + stride_x)
    {
      for (j = stride_y; j < s_y; j = j + stride_y)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         up_index = (s_x * s_y * k) + (s_x * (j + interval)) + i;
         down_index = (s_x * s_y * k) + (s_x * (j - interval)) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + up_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + down_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);
           //if (file->idx_c->grank == 0)
           //  printf("[%d %d] -> %f %f %f\n", i, j, left, right, new_val);
           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // first sample
      int first_index = (s_x * s_y * k) + (s_x * 0) + i;
      int up_first_index = (s_x * s_y * k) + (s_x * interval) + i;

      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;
        // left
        if (negative_rank_y != -1)
          memcpy (&left, negative_buffer_y + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&left, wb + up_first_index * bytes_for_datatype, bytes_for_datatype);

        // right
        memcpy (&right, wb + up_first_index * bytes_for_datatype, bytes_for_datatype);

        // existing value
        memcpy (&new_val, wb + first_index * bytes_for_datatype, bytes_for_datatype);

        //if (file->idx_c->grank == 2)
        //  printf("%f %f %f\n", left, new_val, right);

        // wavlet transformation
        new_val = new_val + 0.25 * (left + right);

        // new value
        memcpy(wb + first_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
    }
  }

  return PIDX_success;
}


static PIDX_return_code wavelet_comm_p2p_odd_z (PIDX_io file, int gi, int v, int stride_x, int stride_y)
{
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  //int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  //printf("My rank %d [%d + (%d/%d)] Rank I get data from %d Rank I send data to %d\n", file->idx_c->grank, file->idx_c->grank_x, file->idx_c->gnproc_x, file->idx->number_processes[0], positive_rank_x, negative_rank_x);

  int req_count = 0;
  if (positive_rank_z != -1)
  {
    receive_z = receive_z + ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype);
    MPI_Irecv(positive_buffer_z, ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype), MPI_BYTE, positive_rank_z, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }


#if 1
  int count = 0;
  int lindex = 0;

  if (negative_rank_z != -1)
  {
    for (j = 0; j < s_y; j = j + stride_y)
    {
      for (i = 0; i < s_x; i =  i + stride_x)
      {
        k = 0;
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(negative_buffer_z + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_z = sent_z + ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype);
    MPI_Isend(negative_buffer_z, ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype), MPI_BYTE, negative_rank_z, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  int ret;
  MPI_Status status[2];
  ret = MPI_Waitall(req_count, request, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_wavelet;
  }
#endif
  return PIDX_success;

}



static PIDX_return_code wavelet_comp_odd_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  int odd_start_offset = (int)pow(2, l);// * s_x;
  int interval = (int)pow(2,l);
  int up_index = 0;
  int down_index = 0;
  int index = 0;
  int i = 0, j = 0, k = 0;
  int sample_count = 0;

#if 1
  for (j = 0; j < s_y; j = j + stride_y)
  {
    for (i = 0 ; i < s_x; i = i + stride_x)
    {
      // every sample but the last
      for (k = odd_start_offset ; k < s_z - stride_z; k =  k + stride_z)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         up_index = (s_x * s_y * (k + interval)) + (s_x * j) + i;
         down_index = (s_x * s_y * (k - interval)) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           //if (file->idx_c->grank == 0)
           //  printf ("[%d %d] -> %d %d %d\n", i, j, index, up_index, down_index);
           float left, right, new_val;
           memcpy (&left, wb + up_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + down_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);

           new_val = new_val - 0.5 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // last sample
      int last_index = (s_x * s_y * (s_z - interval)) + (s_x * j) + i;
      int last_index_down = (s_x * s_y * (s_z - interval - interval)) + (s_x * j) + i;
      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;

        // left
        //memcpy (&left, wb + (last_index - interval) * bytes_for_datatype, bytes_for_datatype);
        memcpy (&left, wb + last_index_down * bytes_for_datatype, bytes_for_datatype);

        // right
        if (positive_rank_z != -1)
          memcpy (&right, positive_buffer_z + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&right, wb + last_index_down * bytes_for_datatype, bytes_for_datatype);

        // center
        memcpy (&new_val, wb + last_index * bytes_for_datatype, bytes_for_datatype);

        // new value
        //if (file->idx_c->grank == 0)
        //  printf ("[%d %d] -> %f %f %f\n", i, j, new_val, left, right);
        new_val = new_val - 0.5 * (left + right);

        memcpy(wb + last_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
    }
  }
#endif

  return PIDX_success;
}



// TODO: Check before if it is even possible to do wavelet transform in the direction.
static PIDX_return_code wavelet_comm_p2p_even_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y)
{
  int ret = 0;
  int count = 0;
  int lindex = 0;
  int req_count = 0;
  MPI_Request request[2];
  int i = 0, j = 0, k = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  req_count = 0;
  if (negative_rank_z != -1)
  {
    receive_z = receive_z + ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype);
    MPI_Irecv(negative_buffer_z, ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype), MPI_BYTE, negative_rank_z, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  if (positive_rank_z != -1)
  {
    count = 0;
    lindex = 0;
    for (j = 0; j < s_y; j = j + stride_y)
    {
      for (i = 0; i < s_x; i =  i + stride_x)
      {
        k = s_z - (int)pow(2, l);
        lindex = (s_x * s_y * k) + (s_x * j) + i;
        memcpy(positive_buffer_z + count * bytes_for_datatype, wb + lindex * bytes_for_datatype, bytes_for_datatype);
        count++;
      }
    }

    sent_z = sent_z + ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype);
    MPI_Isend(positive_buffer_z, ((int)ceil((double)s_x/stride_x) * (int)ceil((double)s_y/stride_y) * bytes_for_datatype), MPI_BYTE, positive_rank_z, 1234, file->idx_c->global_comm, &request[req_count]);
    req_count++;
  }

  MPI_Status status2[2];
  ret = MPI_Waitall(req_count, request, status2);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }

  return PIDX_success;
}


static PIDX_return_code wavelet_comp_even_z (PIDX_io file, int gi, int v, int l, int stride_x, int stride_y, int stride_z)
{
  int sample_count = 0;
  int interval = (int)pow(2,l);
  int up_index = 0;
  int down_index = 0;
  int index = 0;
  int i = 0, j = 0, k = 0;

  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var = var_grp->variable[v];
  int chunk_size = file->idx->chunk_size[0] * file->idx->chunk_size[1] * file->idx->chunk_size[2];
  int bytes_for_datatype = ((var_grp->variable[v]->bpv / 8) * chunk_size * var_grp->variable[v]->vps) / file->idx->compression_factor;

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[0];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[0];
  unsigned char* wb = var->rst_patch_group[0]->reg_patch->buffer;

  // Even Samples
  sample_count = 0;
  for (j = 0; j < s_y; j = j + stride_y)
  {
    for (i = 0; i < s_x; i = i + stride_x)
    {
      for (k = stride_z; k < s_z; k = k + stride_z)
      {
         index = (s_x * s_y * k) + (s_x * j) + i;
         up_index = (s_x * s_y * (k + interval)) + (s_x * j) + i;
         down_index = (s_x * s_y * (k - interval)) + (s_x * j) + i;
         if (bytes_for_datatype == 4)
         {
           float left, right, new_val;
           memcpy (&left, wb + up_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&right, wb + down_index * bytes_for_datatype, bytes_for_datatype);
           memcpy (&new_val, wb + index * bytes_for_datatype, bytes_for_datatype);
           //if (file->idx_c->grank == 0)
           //  printf("[%d %d] -> %f %f %f\n", i, j, left, right, new_val);
           new_val = new_val + 0.25 * (left + right);
           memcpy(wb + index * bytes_for_datatype, &new_val, bytes_for_datatype);
         }
      }

      // first sample
      int first_index = (s_x * s_y * 0) + (s_x * j) + i;
      int up_first_index = (s_x * s_y * interval) + (s_x * j) + i;

      if (bytes_for_datatype == 4)
      {
        float left, right, new_val;
        // left
        if (negative_rank_z != -1)
          memcpy (&left, negative_buffer_z + sample_count * bytes_for_datatype, bytes_for_datatype);
        else
          memcpy (&left, wb + up_first_index * bytes_for_datatype, bytes_for_datatype);

        // right
        memcpy (&right, wb + up_first_index * bytes_for_datatype, bytes_for_datatype);

        // existing value
        memcpy (&new_val, wb + first_index * bytes_for_datatype, bytes_for_datatype);

        // wavlet transformation
        new_val = new_val + 0.25 * (left + right);

        // new value
        memcpy(wb + first_index * bytes_for_datatype, &new_val, bytes_for_datatype);
      }
      sample_count++;
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

  int s_x = (int) var->rst_patch_group[0]->reg_patch->size[0] / file->idx->chunk_size[0];
  int s_y = (int) var->rst_patch_group[0]->reg_patch->size[1] / file->idx->chunk_size[1];
  int s_z = (int) var->rst_patch_group[0]->reg_patch->size[2] / file->idx->chunk_size[2];

  int o_x = (int) var->rst_patch_group[0]->reg_patch->offset[0] / file->idx->chunk_size[0];
  int o_y = (int) var->rst_patch_group[0]->reg_patch->offset[1] / file->idx->chunk_size[1];
  int o_z = (int) var->rst_patch_group[0]->reg_patch->offset[2] / file->idx->chunk_size[2];

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

  if (file->idx_c->grank == 0)
  {
    printf("STENCIL\n");
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
          //printf("%f\t", v);
          if (v != 0)
            printf("%d %d %d - %f\n", i, j, k, v);
        }
        //printf("\n");
      }
      //printf("\n");
    }
    printf("\n");
  }

  if (file->idx_c->grank == 0)
    free(global_buffer);

  exit:

  return PIDX_success;
}

