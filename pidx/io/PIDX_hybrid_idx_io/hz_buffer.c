#include "../PIDX_io.h"


PIDX_return_code create_hz_buffers(PIDX_hybrid_idx_io file, int svi, int evi)
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

  int si = 0, ei = 0;
  for (si = svi; si < evi; si = si + (file->idx_d->var_pipe_length + 1))
  {
    ei = ((si + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->var_pipe_length);


    // Create the chunking ID
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, svi, si, ei);

    // Create the compression ID
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, svi, si, ei);

    // Create the HZ encoding ID
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, svi, si, ei);


    if (file->idx_d->parallel_mode)
    {
      // Attaching the communicator to the chunking phase
      ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }

      // Attaching the communicator to the compression phase
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }

      // Attaching the communicator to the HZ encodig phase phase
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // metadata for chunking phase
    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    // resolution for HZ encoding
    ret = PIDX_hz_encode_set_resolution(file->hz_id, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    // metadata for hz encoding phase
    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }


    // Creating the buffers required for chunking
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

    // Perform Chunking
    if (file->idx_dbg->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }
    }


    // Perform Compression
    if (file->idx_dbg->debug_do_compress == 1)
    {
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }
    }


    // Creating the buffers required for HZ encoding
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    // Perform HZ encoding
    if (file->idx_dbg->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Verify the HZ encoding
    if(file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Destroy buffers allocated during chunking phase
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
  }

  return PIDX_success;
}


PIDX_return_code populate_hz_buffers(PIDX_hybrid_idx_io file, int svi, int evi)
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

  int si = 0, ei = 0;
  for (si = svi; si < evi; si = si + (file->idx_d->var_pipe_length + 1))
  {
    ei = ((si + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->var_pipe_length);

    // Verify the HZ encoding
    if(file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Perform HZ encoding
    if (file->idx_dbg->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_read(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Perform Compression
    if (file->idx_dbg->debug_do_compress == 1)
    {
      ret = PIDX_decompression(file->comp_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }
    }

    // Perform Chunking
    if (file->idx_dbg->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }
    }


    // Destroy buffers allocated during chunking phase
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
  }

  return PIDX_success;
}


PIDX_return_code setup_hz_buffers(PIDX_hybrid_idx_io file, int svi, int evi)
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

  int si = 0, ei = 0;
  for (si = svi; si < evi; si = si + (file->idx_d->var_pipe_length + 1))
  {
    ei = ((si + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->var_pipe_length);


    // Create the chunking ID
    file->chunk_id = PIDX_chunk_init(file->idx, file->idx_d, svi, si, ei);

    // Create the compression ID
    file->comp_id = PIDX_compression_init(file->idx, file->idx_d, svi, si, ei);

    // Create the HZ encoding ID
    file->hz_id = PIDX_hz_encode_init(file->idx, file->idx_d, svi, si, ei);


    if (file->idx_d->parallel_mode)
    {
      // Attaching the communicator to the chunking phase
      ret = PIDX_chunk_set_communicator(file->chunk_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }

      // Attaching the communicator to the compression phase
      ret = PIDX_compression_set_communicator(file->comp_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }

      // Attaching the communicator to the HZ encodig phase phase
      ret = PIDX_hz_encode_set_communicator(file->hz_id, file->comm);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // metadata for chunking phase
    ret = PIDX_chunk_meta_data_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    // resolution for HZ encoding
    ret = PIDX_hz_encode_set_resolution(file->hz_id, file->idx_d->reduced_res_from, file->idx_d->reduced_res_to);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }

    // metadata for hz encoding phase
    ret = PIDX_hz_encode_meta_data_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_rst;
    }


    // Creating the buffers required for chunking
    ret = PIDX_chunk_buf_create(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }

    /*
    // Perform Chunking
    if (file->idx_dbg->debug_do_chunk == 1)
    {
      ret = PIDX_chunk(file->chunk_id, PIDX_WRITE);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_chunk;
      }
    }


    // Perform Compression
    if (file->idx_dbg->debug_do_compress == 1)
    {
      ret = PIDX_compression(file->comp_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_compress;
      }
    }
    */


    // Creating the buffers required for HZ encoding
    ret = PIDX_hz_encode_buf_create(file->hz_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_hz;
    }

    /*
    // Perform HZ encoding
    if (file->idx_dbg->debug_do_hz == 1)
    {
      ret = PIDX_hz_encode_write(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Verify the HZ encoding
    if(file->idx_dbg->debug_hz == 1)
    {
      ret = HELPER_Hz_encode(file->hz_id);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_hz;
      }
    }

    // Destroy buffers allocated during chunking phase
    ret = PIDX_chunk_buf_destroy(file->chunk_id);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_chunk;
    }
    */
  }

  return PIDX_success;
}




PIDX_return_code destroy_hz_buffers(PIDX_hybrid_idx_io file, int svi, int evi)
{
  PIDX_return_code ret;

  PIDX_time time = file->idx_d->time;
  int si = 0;//, ei = 0;
  for (si = svi; si < evi; si = si + (file->idx_d->var_pipe_length + 1))
  {
    time->startup_start[si] = PIDX_get_time();
    //ei = ((si + file->idx_d->var_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx_d->var_pipe_length);

    /*-------------------------------------------finalize [start]---------------------------------------------*/
    time->finalize_start[si] = PIDX_get_time();

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

    time->finalize_end[si] = PIDX_get_time();
    /*-----------------------------------------finalize [end]--------------------------------------------------*/

    /* Deleting the HZ encoding ID */
    PIDX_hz_encode_finalize(file->hz_id);
  }

  return PIDX_success;
}
