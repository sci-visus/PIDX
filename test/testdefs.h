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

#ifndef TESTDEFS_H
#define TESTDEFS_H

// int test_writer(struct Args args, int rank, int nprocs);
// int usage_writer();

int test_reader(struct Args args, int rank, int nprocs);
int usage_reader();

// int test_one_var_writer(struct Args args, int rank, int nprocs);
// int usage_one_var_writer();

// int test_multi_var_writer(struct Args args, int rank, int nprocs);
// int usage_multi_var_writer();

int serial_writer(struct Args args);
int usage_serial();

int serial_reader(struct Args args);
int usage_serial_reader();

// int test_multi_patch_writer(struct Args args, int rank, int nprocs);
// int usage_multi_patch_writer();

int test_multi_idx_writer(struct Args args, int rank, int nprocs);
int usage_multi_idx_writer();

#endif
