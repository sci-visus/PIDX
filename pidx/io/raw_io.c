#include "../PIDX_inc.h"

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);
static PIDX_return_code dump_debug_data_init(PIDX_io file);
static PIDX_return_code dump_debug_data_finalie (PIDX_io file);
static PIDX_return_code dump_process_extent(PIDX_io file);

/// Raw Write Steps
/********************************************************
*  Step 0: Setup Group related meta-data                *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform data Restructuring                   *
*  Step 3: Perform actual file IO                       *
*  Step 4: cleanup for Steps 1                          *
*********************************************************/

PIDX_return_code PIDX_raw_write(PIDX_io file, int gi, int svi, int evi)
{  
  int si = 0, ei = 0;
  PIDX_return_code ret;

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    // Step 1: Setup restructuring buffers
    ret = restructure_setup(file, gi, si, ei, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 2: Perform data restructuring
    ret = restructure(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 3: Write out restructured data
    ret = restructure_io(file, PIDX_WRITE);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }


    // Step 4: Cleanup all buffers and ids
    ret = restructure_cleanup(file);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
  }

  dump_process_extent(file);

  return PIDX_success;
}



/// Raw Read Steps
/********************************************************
*  Step 0: Setup Group related meta-data                *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform actual file IO                       *
*  Step 3: Perform data Restructuring                   *
*  Step 4: cleanup for Steps 1                          *
*********************************************************/

PIDX_return_code PIDX_raw_read(PIDX_io file, int gi, int svi, int evi)
{
  int si = 0, ei = 0;
  PIDX_return_code ret;

  // Step 0
  ret = group_meta_data_init(file, gi, svi, evi, PIDX_READ);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  int rst_case_type = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  PIDX_variable var0 = var_grp->variable[svi];
  int patch_count = var0->sim_patch_count;
  int max_patch_count = 0;

  MPI_Allreduce(&patch_count, &max_patch_count, 1, MPI_INT, MPI_MAX, file->idx_c->global_comm);
  if (max_patch_count > 1)
    rst_case_type = 1;

  file->idx->variable_pipe_length = file->idx->variable_count;
  if (file->idx_d->data_core_count == file->idx_c->gnprocs && rst_case_type == 0)
  {
    for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      // Step 1: Setup restructuring buffers
      ret = restructure_setup(file, gi, si, ei, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 2: Write out restructured data
      ret = restructure_io(file, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 3: Perform data restructuring
      ret = restructure(file, PIDX_READ);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }

      // Step 4: Cleanup all buffers and ids
      ret = restructure_cleanup(file);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
  }
  else
  {
    for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
    {
      ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
      file->idx->variable_grp[gi]->variable_tracker[si] = 1;

      ret = restructure_forced_read(file, si, ei);
      if (ret != PIDX_success)
      {
        fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_file;
      }
    }
  }

  return PIDX_success;
}


static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode)
{
  int ret;
  PIDX_time time = file->idx_d->time;

  time->set_reg_box_start = MPI_Wtime();
  ret = set_rst_box_size(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->set_reg_box_end = MPI_Wtime();

  if (mode == PIDX_WRITE)
  {
    time->header_io_start = PIDX_get_time();
    // Creates the file heirarchy and writes the header info for all binary files
    ret = init_raw_headers_layout(file, gi, svi, evi, file->idx->filename);
    if (ret != PIDX_success)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_file;
    }
    time->header_io_end = PIDX_get_time();
  }

  return PIDX_success;
}


static PIDX_return_code dump_debug_data_init(PIDX_io file)
{
  if (file->idx_dbg->dump_process_state == 1)
  {
    int ret;
    char io_file_name[1024];
    ret = mkdir(file->idx_dbg->process_state_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in aggregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, file->idx_dbg->io_dump_dir_name);
      return PIDX_err_io;
    }

    MPI_Barrier(file->idx_c->global_comm);

    sprintf(io_file_name, "%s/rank_%d", file->idx_dbg->process_state_dump_dir_name, file->idx_c->grank);
    file->idx_dbg->process_state_dump_fp = fopen(io_file_name, "a+");
    if (!file->idx_dbg->process_state_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] io_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, io_file_name);
      return PIDX_err_io;
    }
  }

  return PIDX_success;
}

static PIDX_return_code dump_debug_data_finalie (PIDX_io file)
{

  if (file->idx_dbg->dump_process_state == 1)
  {
    fclose(file->idx_dbg->process_state_dump_fp);
  }

  return PIDX_success;
}

static PIDX_return_code dump_process_extent(PIDX_io file)
{
  int i, j, k;
  dump_debug_data_init(file);
  for (i = 0; i < file->idx->variable_group_count; i++)
  {
    PIDX_variable_group var_grp = file->idx->variable_grp[i];
    for (j = 0; j < var_grp->variable_count; j++)
    {
      for (k = 0; k < var_grp->variable[j]->sim_patch_count; k++)
      {
        if (file->idx_dbg->dump_process_state == 1)
        {
          fprintf(file->idx_dbg->process_state_dump_fp, "[%d] [%d] %d %d %d %d %d %d\n", j, k, (int)var_grp->variable[j]->sim_patch[k]->offset[0], (int)var_grp->variable[j]->sim_patch[k]->offset[1], (int)var_grp->variable[j]->sim_patch[k]->offset[2], (int)var_grp->variable[j]->sim_patch[k]->size[0], (int)var_grp->variable[j]->sim_patch[k]->size[1], (int)var_grp->variable[j]->sim_patch[k]->size[2]);
        }
      }
      if (file->idx_dbg->dump_process_state == 1)
        fprintf(file->idx_dbg->process_state_dump_fp, "\n");
    }
  }
  dump_debug_data_finalie(file);

  return PIDX_success;
}
