#include "../PIDX_inc.h"

static PIDX_return_code group_meta_data_init(PIDX_io file, int gi, int svi, int evi, int mode);

// IDX Write Steps
/********************************************************
*  Step 0: Setup Group and IDX related meta-data        *
*                                                       *
*  Step 1: Setup Restructuring Phase                    *
*  Step 2: Perform data Restructuring                   *
*  Step 3: Setup HZ encoding Phase                      *
*  Step 4: Perform HZ encoding                          *
*  Step 5: Setup aggregation buffers                    *
*  Step 6: Perform data aggregation                     *
*  Step 7: Perform actual file IO                       *
*  Step 8: cleanup for Steps 1, 3, 5                    *
*                                                       *
*  Step 9: Cleanup the group and IDX related meta-data  *
*********************************************************/

PIDX_return_code PIDX_idx_insitu(PIDX_io file, int gi, int svi, int evi)
{
  pa1 = MPI_Wtime();
  PIDX_return_code ret;
  PIDX_time time = file->idx_d->time;

  // Step 0
  time->set_reg_box_start = PIDX_get_time();
  ret = set_rst_box_size(file, gi, svi);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  time->set_reg_box_end = MPI_Wtime();


  pa2 = MPI_Wtime();
  // Step 1: Setup restructuring buffers
  ret = restructure_setup(file, gi, svi, evi - 1, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }


  pa3 = MPI_Wtime();
  // Step 2: Perform data restructuring
  ret = restructure(file, PIDX_WRITE);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }

  pa4 = MPI_Wtime();
  int si = 0, ei = 0;
  file->idx->variable_pipe_length = file->idx->variable_count;
  for (si = svi; si < evi; si = si + (file->idx->variable_pipe_length + 1))
  {
    ei = ((si + file->idx->variable_pipe_length) >= (evi)) ? (evi - 1) : (si + file->idx->variable_pipe_length);
    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    file->insitu_id = PIDX_in_situ_init(file->idx, file->idx_d, file->idx_c, file->idx_dbg, gi, svi, evi);
    PIDX_in_situ_perform(file->insitu_id);
    PIDX_in_situ_finalize(file->insitu_id);
  }


  pa5 = MPI_Wtime();
  ret = restructure_cleanup(file);
  if (ret != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_file;
  }
  pa6 = MPI_Wtime();

  //printf("[%d] BOX %f SETUP %f RST %f INSITU %f CLEANUP %f\n", file->idx_c->grank, (a2 - a1), (a3 - a2), (a4 - a3), (a5 - a4), (a6 - a5));

  return PIDX_success;
}
