/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

#ifndef PIDXTEST_H
#define PIDXTEST_H

#include <PIDX.h>

/* Kind of test to run */
enum Kind { DEFAULT = 0, SERIAL_READER, PARALLEL_READER, SERIAL_WRITER, PARALLEL_WRITER, PARALLEL_MULTI_PATCH_WRITER};

/*kindToStr*/
char* kindToStr(enum Kind k);

/*strToKind*/
enum Kind strToKind(const char *str);

/*Args*/
struct Args
{
  /* type of test to run */
  enum Kind kind;

  /* global dimensions of 3D volume */
  int extents[5];

  /* per-process dimensions of each sub-block */
  int count_local[5];

  /* Number of time-steps */
  int time_step;

  /* Number of Variables */
  int variable_count;
  
  /* output IDX file Name Template*/
  char output_file_template[512];
  char *output_file_name;
};

/*main*/
int main(int argc, char **argv);

/*parse_args*/
int parse_args(struct Args *args, int argc, char **argv);

/*print_error*/
int print_error(char *error_message, char* file, int line);

#endif