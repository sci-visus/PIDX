#include "../PIDX_inc.h"

 /// 0 for patch per process and 1 for multi patch per process
static int rst_case_type = 0;


/// Initialiazation and creation of buffers for restructuring phase
/// 1. Find if this is a patch per process problem or  multi patch per process problem
/// 2. If the restructuring box size is provided then use that else find restructuring box size
/// 3.
PIDX_return_code restructure_init(PIDX_io file, int gi, int svi, int evi)
{
  int ret = 0;
  int start_index = 0, end_index = 0;
  int pipe_length = file->idx->variable_count;
  int factor = 1;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  int grp_count = var0->patch_group_count;
  int max_grp_count = 0;
  int patch_count = var0->sim_patch_count;
  int max_patch_count = 0;

  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->comm);
  if (max_patch_count > 1)
    rst_case_type = 1;

  if (rst_case_type == 0)
  {
    for (start_index = svi; start_index < evi; start_index = start_index + (pipe_length + 1))
    {
      factor = 1;
      end_index = ((start_index + pipe_length) >= (evi)) ? (evi - 1) : (start_index + pipe_length);

      /// Initialize the restructuring phase
      file->rst_id = PIDX_rst_init(file->idx, file->idx_d, svi, start_index, end_index);

      /// attach communicator to the restructuring phase
      ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

      if (file->idx->reg_box_set == 1)
      {
        /// Finding the regular power two box dimensions
        if (!(file->idx->reg_patch_size[0] == 0 && file->idx->reg_patch_size[1] == 0 && file->idx->reg_patch_size[2] == 0))
        {
          PIDX_rst_set_reg_patch_size(file->rst_id, file->idx->reg_patch_size);

          /// Creating the metadata to perform retructuring
          ret = PIDX_rst_meta_data_create(file->rst_id);
          if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
      }
      else
      {
        restructure_tag:
        PIDX_rst_auto_set_reg_patch_size(file->rst_id, factor);

        /// Creating the metadata to perform retructuring
        ret = PIDX_rst_meta_data_create(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        grp_count = var0->patch_group_count;
        max_grp_count = 0;
        MPI_Allreduce(&grp_count, &max_grp_count, 1, MPI_INT, MPI_MAX, file->comm );
        if (max_grp_count > 1)
        {
          ret = PIDX_rst_meta_data_destroy(file->rst_id);
          if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

          factor = factor * 2;

          /// Ensusring that the partitioning step works fine by
          /// making sure that number of boxes after restructuring is 1
          goto restructure_tag;
        }
      }

      /// Saving the metadata info needed for reading back the data.
      /// Especially when number of cores is different from number of cores
      /// used to create the dataset
      ret = PIDX_rst_meta_data_write(file->rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

      //// Creating the buffers required for restructurig
      ret = PIDX_rst_buf_create(file->rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

      /// Aggregating the aligned small buffers after restructuring into one single buffer
      ret = PIDX_rst_aggregate_buf_create(file->rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    }
  }

  /// Case for more than one patch per process
  else
  {
    for (start_index = svi; start_index < evi; start_index = start_index + (pipe_length + 1))
    {
      factor = 1;
      end_index = ((start_index + pipe_length) >= (evi)) ? (evi - 1) : (start_index + pipe_length);

      /// Initialize the restructuring phase
      file->multi_patch_rst_id = PIDX_multi_patch_rst_init(file->idx, file->idx_d, svi, start_index, end_index);

      /// attach communicator to the restructuring phase
      ret = PIDX_multi_patch_rst_set_communicator(file->multi_patch_rst_id, file->comm);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}


      if (file->idx->reg_box_set == 1)
      {
        /// Finding the regular power two box dimensions
        if (!(file->idx->reg_patch_size[0] == 0 && file->idx->reg_patch_size[1] == 0 && file->idx->reg_patch_size[2] == 0))
        {
            PIDX_multi_patch_rst_set_reg_patch_size(file->multi_patch_rst_id, file->idx->reg_patch_size);

            /// Creating the metadata to perform retructuring
            ret = PIDX_multi_patch_rst_meta_data_create(file->multi_patch_rst_id);
            if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
        }
      }
      else
      {
        multi_patch_restructure_tag:
        PIDX_multi_patch_rst_auto_set_reg_patch_size(file->multi_patch_rst_id, factor);

        /// Creating the metadata to perform retructuring
        ret = PIDX_multi_patch_rst_meta_data_create(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        grp_count = var0->patch_group_count;
        max_grp_count = 0;
        MPI_Allreduce(&grp_count, &max_grp_count, 1, MPI_INT, MPI_MAX, file->comm );
        if (max_grp_count > 1)
        {
          ret = PIDX_multi_patch_rst_meta_data_destroy(file->multi_patch_rst_id);
          if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

          factor = factor * 2;

          /// Ensusring that the partitioning step works fine by
          /// making sure that number of boxes after restructuring is 1
          goto multi_patch_restructure_tag;
        }
      }

      /// Saving the metadata info needed for reading back the data.
      /// Especially when number of cores is different from number of cores
      /// used to create the dataset
      ret = PIDX_multi_patch_rst_meta_data_write(file->multi_patch_rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

      //// Creating the buffers required for restructurig
      ret = PIDX_multi_patch_rst_buf_create(file->multi_patch_rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

      /// Aggregating the aligned small buffers after restructuring into one single buffer
      ret = PIDX_multi_patch_rst_aggregate_buf_create(file->multi_patch_rst_id);
      if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure(PIDX_io file, int mode)
{
  int ret = 0;

  if (rst_case_type == 0)
  {
    if (mode == PIDX_WRITE)
    {
      if (file->idx_dbg->debug_do_rst == 1)
      {
        /// Perform data restructuring
        ret = PIDX_rst_staged_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        /// Aggregating in memory restructured buffers into one large buffer
        ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_WRITE);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        /// Destroying the restructure buffers (as they are now combined into one large buffer)
        ret = PIDX_rst_buf_destroy(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
    else if (mode == PIDX_READ)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_rst_buf_aggregate(file->rst_id, PIDX_READ);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        HELPER_rst(file->rst_id);

        ret = PIDX_rst_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        ret = PIDX_rst_buf_destroy(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
  }

  /// Case for more than one patch per process
  else
  {
    if (mode == PIDX_WRITE)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_multi_patch_rst_staged_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_WRITE);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
    else if (mode == PIDX_READ)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_multi_patch_rst_buf_aggregate(file->multi_patch_rst_id, PIDX_READ);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        ret = PIDX_multi_patch_rst_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

        ret = PIDX_multi_patch_rst_buf_destroy(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure_io(PIDX_io file, int mode)
{
  int ret = 0;

  if (rst_case_type == 0)
  {
    if (mode == PIDX_WRITE)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_rst_buf_aggregated_write(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
    else if (mode == PIDX_READ)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_rst_buf_aggregated_read(file->rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
  }

  /// Case for more than one patch per process
  else
  {
    if (mode == PIDX_WRITE)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_multi_patch_rst_buf_aggregated_write(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
    else if (mode == PIDX_READ)
    {
      /* Perform data restructuring */
      if (file->idx_dbg->debug_do_rst == 1)
      {
        ret = PIDX_multi_patch_rst_buf_aggregated_read(file->multi_patch_rst_id);
        if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
      }
    }
  }

  return PIDX_success;
}



PIDX_return_code restructure_cleanup(PIDX_io file, int gi)
{
  int ret = 0;

  if (rst_case_type == 0)
  {
    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_rst_aggregate_buf_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_rst_meta_data_destroy(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    /* Deleting the restructuring ID */
    PIDX_rst_finalize(file->rst_id);
  }

  /// Case for more than one patch per process
  else
  {
    /* Destroy buffers allocated during restructuring phase */
    ret = PIDX_multi_patch_rst_aggregate_buf_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_multi_patch_rst_meta_data_destroy(file->multi_patch_rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    /* Deleting the restructuring ID */
    PIDX_multi_patch_rst_finalize(file->multi_patch_rst_id);
  }

  return PIDX_success;
}


PIDX_return_code restructure_forced_read(PIDX_io file, int svi, int evi)
{
  int ret = 0;
  int start_index = 0, end_index = 0;
  int pipe_length = file->idx->variable_count;

  for (start_index = svi; start_index < evi; start_index = start_index + (pipe_length + 1))
  {
    end_index = ((start_index + pipe_length) >= (evi)) ? (evi - 1) : (start_index + pipe_length);

    /// Initialize the restructuring phase
    file->rst_id = PIDX_rst_init(file->idx, file->idx_d, svi, start_index, end_index);

    /// attach communicator to the restructuring phase
    ret = PIDX_rst_set_communicator(file->rst_id, file->comm);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_rst_forced_raw_read(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}

    ret = PIDX_rst_finalize(file->rst_id);
    if (ret != PIDX_success) {fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__); return PIDX_err_rst;}
  }

  return PIDX_success;
}
