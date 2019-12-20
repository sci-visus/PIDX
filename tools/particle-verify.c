/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */
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

  char write_file[PIDX_FILE_PATH_LENGTH];
  char read_file[PIDX_FILE_PATH_LENGTH];

  if ( argc != 4 )
  {
    fprintf(stderr, "Missing arguments (provide a filename template and number of write files and number of read files)\n");
    return 0;
  }

  char wfile[PIDX_FILE_PATH_LENGTH];
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


  char rfile[PIDX_FILE_PATH_LENGTH];
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
