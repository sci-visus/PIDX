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

  pconf = procs_conf[n_cores_read][0]

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
  test_str = mpirun+" -np "+str(n_cores)+" "+write_executable+" -g "+g_box+" -l "+l_box+" -r "+r_box+" -t "+str(n_ts)+" -v "+vars_file+" -f data"
  
  if(debug_print>0):
    print "EXECUTE write:", test_str
  os.popen(test_str+" 2&> _out_write.txt")

  if profiling > 0:
    append_profile("_out_write.txt", prof_file_write)

  for t in range(1, n_ts+1):
    for vr in range(0, n_vars):
      n_tests = n_tests + 1

      test_str= mpirun+" -np "+str(n_cores_read)+" "+read_executable+" -g "+g_box+" -l "+l_box_read+" -t "+str(t-1)+" -v "+str(vr)+" -f data"

      #print "Testing t="+str(t)+" v="+str(vr)+" type="+var_type
      os.popen(test_str+" 2&> _out_read.txt")
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

  print "Success %d/%d" % (success, n_tests)

  return n_tests - success

def pow_2(n_cores, n_cores_read, var_type, n_vars, n_ts):
  print "---POW 2 TESTS---"
  
  #even_factor = int(n_cores ** (1. / 3))

  pconf = procs_conf[n_cores][0]

  g_box = (pconf[0]*patch_size[0], pconf[1]*patch_size[1], pconf[2]*patch_size[2])
  l_box = patch_size

  r_box = l_box
  print "r == l", r_box
  succ = execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  r_box = (l_box[0]/2, l_box[1]/2, l_box[2]/2)
  print "r < l", r_box
  succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  r_box = (l_box[0]*2, l_box[1]*2, l_box[2]*2)
  print "r > l", r_box
  succ = succ + execute_test(n_cores, n_cores_read, g_box, l_box, r_box, n_ts, n_vars, var_type)

  if succ == 0:
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

  # defaults
  n_cores = 8
  n_cores_read = 8

  try:
    opts, args = getopt.getopt(argv,"h:w:r:m:p",["pfile="])
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
      prof_file_write = arg+"_write.prof"
      prof_file_read = arg+"_read.prof"
      profiling = 1

  if profiling > 0:
    os.popen("rm "+prof_file+"*.prof")

  if not(n_cores in procs_conf) or not(n_cores_read in procs_conf):
    print "Procs configuration not available use one these: ", sorted(procs_conf.keys()), "\nOr add a new configuration to idx_utils.py"
    sys.exit(2)
  
  for var in var_types:
    succ = succ + pow_2(n_cores, n_cores_read, var, 1, 1)
    succ = succ + non_pow_2(n_cores, n_cores_read, var, 1, 1)

  print "latest outputs:"
  os.popen("cat _out_write.txt")
  os.popen("cat _out_read.txt")

  sys.exit(succ)

if __name__ == "__main__":
  main(sys.argv[1:])

