#include "../PIDX_io.h"


PIDX_return_code create_hz_buffers(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index)
{
  int rank = 0, nprocs = 1;
  int ret = 0;

#if PIDX_HAVE_MPI
  if (file->idx_d->parallel_mode == 1)
  {
    MPI_Comm_size(file->comm,  &nprocs);
    MPI_Comm_rank(file->comm,  &rank);
  }
#endif

  //PIDX_time time = file->idx_d->time;
  int start_index = 0, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    //time->startup_start[start_index] = PIDX_get_time();
    end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);


    /*------------------------------------Create ALL the IDs [start]---------------------------------------*/
    /* Create the chunking ID */
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the compression ID */
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, start_var_index, start_index, end_index);

    /* Create the HZ encoding ID */
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, start_var_index, start_index, end_index);
    /*------------------------------------Create ALL the IDs [end]-------------------------------------------*/



    /*------------------------------------Adding communicator [start]----------------------------------------*/
#if PIDX_HAVE_MPI
    if (file->idx_d->parallel_mode)
    {
      /* Attaching the communicator to the chunking phase */
      ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }

      /* Attaching the communicator to the compression phase */
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }

      /* Attaching the communicator to the HZ encodig phase phase */
      PIDX_hz_encode_set_global_communicator(file->hz_id, file->comm);
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }
#endif
    /*------------------------------------Adding communicator [end]------------------------------------------*/

    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_hz_encode_set_resolution(file->hz_id, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    //time->startup_end[start_index] = PIDX_get_time();


    /*----------------------------------------Chunking [start]------------------------------------------------*/
    //time->chunk_start[start_index] = PIDX_get_time();

    /* Creating the buffers required for chunking */
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

    /* Perform Chunking */
    if (file->idx_dbg->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }
    }

    //time->chunk_end[start_index] = PIDX_get_time();
    /*-----------------------------------------Chunking [end]------------------------------------------------*/


    /*----------------------------------------Compression [start]--------------------------------------------*/
    //time->compression_start[start_index] = PIDX_get_time();

    /* Perform Compression */
    if (file->idx_dbg->debug_do_compress == 1)
    {
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }
    }

    //time->compression_end[start_index] = PIDX_get_time();
    /*------------------------------------------Compression [end]--------------------------------------------*/



    /*---------------------------------------------HZ [start]------------------------------------------------*/
    //time->hz_start[start_index] = PIDX_get_time();

    /* Creating the buffers required for HZ encoding */
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    /* Perform HZ encoding */
    if (file->idx_dbg->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    /* Verify the HZ encoding */
    if(file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    /* Destroy buffers allocated during chunking phase */
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

    //time->hz_end[start_index] = PIDX_get_time();
    /*----------------------------------------------HZ [end]-------------------------------------------------*/
  }

  return PIDX_success;
}


PIDX_return_code destroy_hz_buffers(PIDX_hybrid_idx_io file, int start_var_index, int end_var_index)
{
  PIDX_return_code ret;

  PIDX_time time = file->idx_d->time;
  int start_index = 0;//, end_index = 0;
  for (start_index = start_var_index; start_index < end_var_index; start_index = start_index + (file->idx_d->var_pipe_length + 1))
  {
    time->startup_start[start_index] = PIDX_get_time();
    //end_index = ((start_index + file->idx_d->var_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (start_index + file->idx_d->var_pipe_length);

    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[start_index] = PIDX_get_time();

    /* Destroy buffers allocated during HZ encoding phase */
    ret = PIDX_hz_encode_buf_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    ret = PIDX_hz_encode_meta_data_destroy(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    ret = PIDX_chunk_meta_data_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    /* Deleting the compression ID */
    PIDX_compression_finalize(file->comp_id);

    /* Deleting the chunking ID */
    PIDX_chunk_finalize(file->chunk_id);

    time->finalize_end[start_index] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);
  }

  return PIDX_success;
}
