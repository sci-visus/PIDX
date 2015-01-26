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

#ifndef __PIDX_FILE_NAME_H
#define __PIDX_FILE_NAME_H

int generate_file_name(int blocks_per_file, char* filename_template, int file_number, char* filename, int maxlen);

int generate_file_name_template(int maxh, int bits_per_block, char* filename, int current_time_step, char* filename_template);

void adjust_file_name(char* bin_file, char* adjusted_name);

void mira_create_folder_name(char* bin_file, char* folder_name);

#endif //__PIDX_FILE_NAME_H
