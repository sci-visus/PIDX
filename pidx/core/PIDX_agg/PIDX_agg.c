#include "../../PIDX_inc.h"



PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int fi, int li)
{
  PIDX_agg_id id;

  id = malloc(sizeof (*id));
  memset(id, 0, sizeof (*id));

  id->idx = idx_meta_data;
  id->idx_d = idx_d;
  id->idx_c = idx_c;

  id->gi = 0;
  id->fi = fi;
  id->li = li;
  //id->lvi = lvi;

  //fprintf(stderr, "fi : li : lvi :: %d : %d : %d\n", fi, li, lvi);

  return id;
}



PIDX_return_code PIDX_agg_finalize(PIDX_agg_id id)
{
  free(id);
  id = 0;

  return PIDX_success;
}
