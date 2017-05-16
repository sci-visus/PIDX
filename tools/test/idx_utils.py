#/*****************************************************
# **  PIDX Parallel I/O Library                      **
# **  Copyright (c) 2010-2014 University of Utah     **
# **  Scientific Computing and Imaging Institute     **
# **  72 S Central Campus Drive, Room 3750           **
# **  Salt Lake City, UT 84112                       **
# **                                                 **
# **  PIDX is licensed under the Creative Commons    **
# **  Attribution-NonCommercial-NoDerivatives 4.0    **
# **  International License. See LICENSE.md.         **
# **                                                 **
# **  For information about this project see:        **
# **  http://www.cedmav.com/pidx                     **
# **  or contact: pascucci@sci.utah.edu              **
# **  For support: PIDX-support@visus.net            **
# **                                                 **
# *****************************************************/

import os
import sys, getopt

vars_file = "VARS"
write_executable = "../build/tutorial/single_buffer_raw_write.app/Contents/MacOS/single_buffer_raw_write"
read_executable = "../build/tutorial/single_buffer_single_variable_read.app/Contents/MacOS/single_buffer_single_variable_read"

# These are the procs configuration that will be used for different core counts
# (first is the default)
procs_conf = dict()
procs_conf[4] = [(2,2,1)]
procs_conf[6] = [(3,2,1)]
procs_conf[8] = [(2,2,2), (4,2,1)]
procs_conf[9] = [(1,3,3)]
procs_conf[10] = [(5,2,1)]
procs_conf[12] = [(3,2,2),(4,3,1)]
procs_conf[14] = [(7,2,1)]
procs_conf[16] = [(2,2,4), (4,4,1)]
procs_conf[32] = [(2,4,4), (8,4,1)]
procs_conf[64] = [(4,4,4), (8,4,2)]
procs_conf[128] = [(4,4,8), (8,8,2), (16,4,2)]

# Append to the profile file the profile output
def append_profile(filename, pfile):
  file = open(filename, "r")
  prof_file = open(pfile, "a")

  for line in file:
    if 'IRPIWCCHHAI' in line:
      prof_file.write(line)

  prof_file.close()

  return 0

# Verify output file of read executable
def verify_read(filename, count):
  file = open(filename, "r")

  correct = 0
  incorrect = -1

  for line in file:
    if 'Correct Sample Count' in line:
      words = line.split()
      correct = int(words[3])
      incorrect = int(words[7])

  if correct == count and incorrect == 0:
    file.close()
    #print "Test SUCCESSED"
    return 0
  
  print "Test FAILED: count "+str(count)+" correct " + str(correct)+ " incorrect "+ str(incorrect)

  return -1

# Generate variables list
def generate_vars(n_vars, type):
  file = open(vars_file,"w") 
  
  file.write("(fields)\n")

  for i in range(n_vars):
    file.write("var_"+str(i)+" "+type+" +\n")

  file.write("(end)")
  file.close() 

