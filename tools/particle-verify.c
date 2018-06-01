#if 0

sort -k1,1 -k2n wfile > swfile
sort -k1,1 -k2n rfile > srfile
diff swfile srfile

#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

int main(int argc, char **argv)
{

  char write_file[1024];
  char read_file[1024];

  if ( argc != 4 )
  {
    fprintf(stderr, "Missing arguments (provide a filename template and number of write files and number of read files)\n");
    return 0;
  }

  char wfile[1024];
  sprintf(wfile, "%s", "wfile");
  int write_file_count = atoi(argv[2]);
  FILE *wfp = fopen(wfile, "w");
  int wcount = 0;
  for (int wc = 0; wc < write_file_count; wc++)
  {
    sprintf(write_file, "%s_w_%d", argv[1], wc);
    FILE *fp = fopen(write_file, "r");

    char line [512];
    while (fgets(line, sizeof (line), fp) != NULL)
    {
      fprintf(wfp, "%s", line);
      wcount++;
    }

    fclose(fp);
  }
  fclose(wfp);


  char rfile[1024];
  sprintf(rfile, "%s", "rfile");
  int read_file_count = atoi(argv[3]);
  FILE *rfp = fopen(rfile, "w");
  int rcount = 0;
  for (int wc = 0; wc < read_file_count; wc++)
  {
    sprintf(read_file, "%s_r_%d", argv[1], wc);
    FILE *fp = fopen(read_file, "r");

    char line [512];
    while (fgets(line, sizeof (line), fp) != NULL)
    {
      fprintf(rfp, "%s", line);
      rcount++;
    }

    fclose(fp);
  }
  fclose(rfp);

  printf("Number of particles W%d R%d\n", wcount, rcount);

  assert(wcount == rcount);


  /*
  char** wstring;
  wstring = malloc(sizeof(*wstring) * rcount);
  for (int pc = 0; pc < rcount; pc++)
  {
    wstring[pc] = malloc(sizeof(*wstring[pc] * 1024));
    memset(wstring[pc], 0, sizeof(*wstring[pc] * 1024));
  }

  char** rstring;

  char wtemp[1024];
  wfp = fopen(wfile, "r");
  char line [512];
  int wcnt = 0;
  while (fgets(line, sizeof (line), wfp) != NULL)
  {
    sprintf(wstring[wcnt], "%s", line);
    printf("%s", wstring[wcnt]);
    wcnt++;
  }
  fclose(wfp);
  */

  return 0;
}
