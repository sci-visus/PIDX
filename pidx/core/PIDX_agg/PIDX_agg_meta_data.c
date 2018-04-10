#include "../../PIDX_inc.h"


PIDX_return_code PIDX_agg_meta_data_create(PIDX_agg_id id, Agg_buffer ab, PIDX_block_layout lbl)
{
  int i = 0, j = 0, j1 = 0, d = 0;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  ab->aggregator_interval = id->idx_c->lnprocs / ((id->li - id->fi + 1) * lbl->efc * ab->agg_f);
  assert(ab->aggregator_interval != 0);

  ab->buffer_size = 0;
  ab->sample_number = -1;
  ab->var_number = -1;
  ab->file_number = -1;

  id->agg_r = malloc(lbl->efc * sizeof (int**));
  memset(id->agg_r, 0, lbl->efc * sizeof (int**));
  for (i = 0; i < lbl->efc; i++)
  {
    id->agg_r[i] = malloc((id->li - id->fi + 1) * sizeof (int*));
    memset(id->agg_r[i], 0, (id->li - id->fi + 1) * sizeof (int*));
    for (j = id->fi; j <= id->li; j++)
    {
      j1 = j - id->fi;
      PIDX_variable var = var_grp->variable[j];
      id->agg_r[i][j1] = malloc(var->vps * sizeof (int) * ab->agg_f);
      memset(id->agg_r[i][j1], 0, var->vps * sizeof (int) * ab->agg_f);
      for (d = 0 ; d < var->vps * ab->agg_f; d++)
        id->agg_r[i][j1][d] = -1;
    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_agg_meta_data_destroy(PIDX_agg_id id, PIDX_block_layout lbl)
{
  int i = 0, j = 0;
  for (i = 0; i < lbl->efc; i++)
  {
    for (j = id->fi; j <= id->li; j++)
    {
      free(id->agg_r[i][j - id->fi]);
      id->agg_r[i][j - id->fi] = 0;
    }
    free(id->agg_r[i]);
  }
  free(id->agg_r);
  id->agg_r = 0;

  return PIDX_success;
}
