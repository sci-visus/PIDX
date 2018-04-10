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
import platform
from idx_utils import *
from test_config import *

#
#   Test and profile PIDX
#   

if platform.system() == "Darwin":
  write_executable = write_executable+".app/Contents/MacOS/idx_write"
  read_executable = read_executable+".app/Contents/MacOS/idx_read"

if not(os.path.isfile(write_executable)):
  print "ERROR: write executable not found!", write_executable
  sys.exit(2)

if not(os.path.isfile(read_executable)):
  print "ERROR: read executable not found!", read_executable
  sys.exit(2)

profiling = 0
profile_file_write = "pidx_write.prof"
profile_file_read = "pidx_read.prof"

vars_file = "VARS"

def execute_test(n_cores, n_cores_read, g_box_n, l_box_n, r_box_n, n_ts, n_vars, var_type):

  g_box = "%dx%dx%d" % (g_box_n[0], g_box_n[1], g_box_n[2])
  l_box = "%dx%dx%d" % (l_box_n[0], l_box_n[1], l_box_n[2])
  r_box = "%dx%dx%d" % (r_box_n[0], r_box_n[1], r_box_n[2])

  l_box_read = l_box
  l_box_read_n = l_box_n

  #pconf = procs_conf[n_cores_read][0]

  for pconf in procs_conf[n_cores_read]:

    if n_cores != n_cores_read:  
      l_box_read_n = (g_box_n[0]/pconf[0], g_box_n[1]/pconf[1], g_box_n[2]/pconf[2])
      l_box_read = "%dx%dx%d" % (l_box_read_n[0], l_box_read_n[1], l_box_read_n[2])

    if (int(g_box_n[0]/l_box_read_n[0])*int(g_box_n[1]/l_box_read_n[1])*int(g_box_n[2]/l_box_read_n[2])) != n_cores_read \
      or (g_box_n[0]%l_box_read_n[0]) != 0 or (g_box_n[1]%l_box_read_n[1]) != 0 or (g_box_n[2]%l_box_read_n[2]) != 0:
      print "INVALID test configuration g_box ", g_box, "read l_box ", l_box_read, " g_box/l_box != ", pconf
      print "Try to change the patch size to get an integer value for g_box/l_box"
      return 1

    g_box_v = g_box.split('x')
    el_count = int(g_box_v[0])*int(g_box_v[1])*int(g_box_v[2])

    n_comp = int(var_type.split('*')[0])
    data_count = g_box_n[0]*g_box_n[1]*g_box_n[2]*n_comp

    if("8" in var_type or "16" in var_type):
      print "WARNING: testing ", var_type, "this datatype can be tested only on small domains"

    success = 0
    n_tests = 0

    generate_vars(n_vars, var_type, var_type)
    test_str = mpirun+" -np "+str(n_cores)+" "+write_executable+" -g "+g_box+" -l "+l_box+" -t "+str(n_ts)+" -v "+vars_file+" -f data"
    #test_str = mpirun+" -np "+str(n_cores)+" "+write_executable+" -g "+g_box+" -l "+l_box+" -r "+r_box+" -t "+str(n_ts)+" -v "+vars_file+" -f data"
    
    if(debug_print>0):
      print "EXECUTE write:", test_str

    if(travis_mode > 0):
      append_travis(test_str)
    else:
      os.popen(test_str+" >> _out_write.txt 2>&1")

    if profiling > 0:
      append_profile("_out_write.txt", prof_file_write)

    for t in range(1, n_ts+1):
      for vr in range(0, n_vars):
        n_tests = n_tests + 1

        test_str= mpirun+" -np "+str(n_cores_read)+" "+read_executable+" -g "+g_box+" -l "+l_box_read+" -t "+str(t-1)+" -v "+str(vr)+" -f data"

        if(debug_print>0):
          print "EXECUTE read:", test_str

        if(travis_mode > 0):
          append_travis(test_str)
          res = 0
        else:
          os.popen(test_str+" >> _out_read.txt 2>&1")
          res = verify_read("_out_read.txt", data_count)

        if profiling > 0:
          append_profile("_out_read.txt", prof_file_read)

        if res < 0:
          print "Test t="+str(t)+" v="+str(vr)+" type="+var_type + " FAILED"
        else:
          success = success + 1
          #os.popen("rm -R _out_read.txt")
          #os.popen("rm -R _out_write.txt")

    #os.popen("rm -R data*")
  if(travis_mode == 0):
    print "Success %d/%d" % (success, n_tests)

  return n_tests - success

def pow_2(n_cores, n_cores_read, var_type, n_vars, n_ts):
  print "---RUN TESTS---"
  
  #even_factor = int(n_cores ** (1. / 3))

  pconfs = procs_conf[n_cores]

  for pconf in pconfs:
    g_box = (pconf[0]*patch_size[0], pconf[1]*patch_size[1], pconf[2]*patch_size[2])
    l_box = patch_size

    r_box = l_box
    #print "r == l", r_box
    
    succ = execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  #r_box = (l_box[0]/2, l_box[1]/2, l_box[2]/2)
  #print "r < l", r_box
  #succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  #r_box = (l_box[0]*2, l_box[1]*2, l_box[2]*2)
  #print "r > l", r_box
  #succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  if succ == 0 and travis_mode == 0:
    print "TEST PASSED"
  
  return succ

def non_pow_2(n_cores, n_cores_read, var_type, n_vars, n_ts):
  print "---NON-POW 2 TESTS---"

  pconf = procs_conf[n_cores][0] #int(n_cores ** (1. / 3))

  g_box = (124, 48, 36)
  l_box = (g_box[0]/pconf[0], g_box[1]/pconf[1], g_box[2]/pconf[2])

  r_box = (64, 32, 32)
  print "r == l", r_box
  succ = execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  r_box = (64, 32, 16)
  print "r < l", r_box
  succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  r_box = (l_box[0]*2, l_box[1]*2, l_box[2]*2)
  print "r > l", r_box
  succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  if succ == 0:
    print "TEST PASSED"
  
  return succ

def print_usage():
  print 'test.py -w <wcores> -r <rcores> -p <profilefile> -m <mpirun>'

def main(argv):
  succ = 0
  global profiling
  global prof_file_write
  global prof_file_read
  global var_types
  global mpirun
  global travis_mode

  # defaults
  n_cores = 8
  n_cores_read = 8

  try:
    opts, args = getopt.getopt(argv,"h:w:r:m:p:t",["pfile="])
  except getopt.GetoptError:
    print_usage()
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print 'test_idx.py -w <wcores> -r <rcores> -p <profilefile> -m <mpirun>'
      sys.exit()
    elif opt in ("-w", "--wcores"):
      n_cores = int(arg)
    elif opt in ("-r", "--rcores"):
      n_cores_read = int(arg)
    elif opt in ("-m", "--mpirun"):
      mpirun = arg
    elif opt in ("-p", "--pfile"):
      prof_file = arg
      prof_file_write = prof_file+"_write.prof"
      prof_file_read = prof_file+"_read.prof"
      profiling = 1
    elif opt in ("-t", "--travismode"):
      file = open("travis_tests.sh", "w")
      file.write("#!/bin/sh\n\n")
      file.close()
      travis_mode = 1
      print "----RUNNING IN TRAVIS TEST MODE----"

  if profiling > 0:
    os.popen("rm "+prof_file+"*.prof")

  if not(n_cores in procs_conf) or not(n_cores_read in procs_conf):
    print "Procs configuration not available use one these: ", sorted(procs_conf.keys()), "\nOr add a new configuration to idx_utils.py"
    sys.exit(2)
  
  for var in var_types:
    succ = succ + pow_2(n_cores, n_cores_read, var, 1, 1)
  #  succ = succ + non_pow_2(n_cores, n_cores_read, var, 1, 1)

  #print "latest outputs:"
  #os.popen("cat _out_write.txt")
  #os.popen("cat _out_read.txt")

  sys.exit(succ)

if __name__ == "__main__":
  main(sys.argv[1:])

