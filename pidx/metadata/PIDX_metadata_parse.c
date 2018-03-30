#include "../PIDX.h"

/// Function to populate file descriptor when opening an existing IDX file
/// with the requested version
PIDX_return_code PIDX_metadata_parse(FILE *fp, PIDX_file* file, char* version)
{
  if(strcmp(version, "6.1") == 0)
    return PIDX_metadata_parse_v6_1(fp, file);
  else if(strcmp(version, "6") == 0)
    return PIDX_metadata_parse_v6_0(fp, file);
  
  fprintf(stderr, "Metadata version %s not supported\n", version);
  
  return PIDX_err_metadata;
}
