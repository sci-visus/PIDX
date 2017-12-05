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
import os.path
from test_config import *

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
def generate_vars(n_vars, *type):
  try:
    file = open(vars_file,"w") 
  
    file.write("(fields)\n")

    if len(type) > 1:
      for i in range(n_vars):
        file.write("var_"+str(i)+" "+type[0]+" +\n")
    else:
      for i in range(n_vars):
        if i % 2 == 0:
          file.write("var_"+str(i)+" "+type[0]+" +\n")
        else:
          file.write("var_"+str(i)+" "+type[1]+" +\n")

    file.write("(end)")
    file.close() 
  except IOError, inst:
    print 'ERROR: generating vars file:', inst.errno, inst.strerror

