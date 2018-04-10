##
## BSD 3-Clause License
## 
## Copyright (c) 2010-2018 ViSUS L.L.C., 
## Scientific Computing and Imaging Institute of the University of Utah
## 
## ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
## University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
##  
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## * Redistributions of source code must retain the above copyright notice, this
## list of conditions and the following disclaimer.
## 
## * Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
## 
## * Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## For additional information about this project contact: pascucci@acm.org
## For support: support@visus.net
## 
##

import os
import sys, getopt
import os.path
from test_config import *

# Append to the profile file the profile output
def append_travis(line):
  file = open("travis_tests.sh", "a")
  file.write(line+"\n")
  file.close()
  return 0

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

