#include "../../../PIDX_inc.h"


PIDX_return_code PIDX_particle_write(PIDX_io file, int gi, int svi, int evi)
{

  int si = 0, p = 0, pt = 0, s = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];

  for (si = svi; si < evi; si++)
  {
    PIDX_variable var = var_grp->variable[si];

    int sample_count;
    int bits_per_sample;
    PIDX_get_datatype_details(var->type_name, &sample_count, &bits_per_sample);

    int bytes_for_datatype = ((var->bpv / 8) * var->vps);

    file->idx->variable_grp[gi]->variable_tracker[si] = 1;

    for (p = 0; p < var->sim_patch_count; p++)
    {
        printf("[Inside PIDX] Rank %d has %d particles [%lld %lld %lld - %lld %lld %lld]\n", file->idx_c->grank, var->sim_patch[p]->particle_count, var->sim_patch[p]->offset[0], var->sim_patch[p]->offset[1], var->sim_patch[p]->offset[2], var->sim_patch[p]->size[0], var->sim_patch[p]->size[1], var->sim_patch[p]->size[2]);
        for (pt = 0; pt < var->sim_patch[p]->particle_count; pt++)
        {
          //for (s = 0; s < sample_count; s++)
          //{
            //if ()
          //}
        }

    }
  }

  return PIDX_success;
}


PIDX_return_code PIDX_particle_read(PIDX_io file, int gi, int svi, int evi)
{

  return PIDX_success;
}
