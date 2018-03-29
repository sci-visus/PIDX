#include "../PIDX.h"

/// Function to populate file descriptor when opening an existing IDX file
/// with the requested version
PIDX_return_code PIDX_metadata_parse(FILE *fp, PIDX_file* file, int version)
{
  if(version == 7)
    return PIDX_metadata_parse_v7(fp, file);
  else if(version == 6)
    return PIDX_metadata_parse_v6(fp, file);
  
  fprintf(stderr, "Metadata version %d not supported\n", version);
  
  return PIDX_err_metadata;
}
