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
write_idx_executable = "../../build/examples/idx_write"
read_idx_executable = "../../build/examples/idx_read"
write_compressed_executable = "../../build/examples/idx_write_compressed"
write_partitioned_executable = "../../build/examples/idx_write_partitioned"

mpirun="mpirun --oversubscribe"

# enable debug prints 
debug_print=1

# travis mode does not execute any test
# but it creates a bash script to run them 
travis_mode = 0

var_types = ["1*float32", "1*int32", "1*float64", 
             "2*float32", "2*int32", "2*float64", 
             "3*float32", "3*int32", "3*float64",
             "1*int64"  , "2*int64", "3*int64"]

var_types_compression = ["1*float32", "1*float64", 
                         "2*float32", "2*float64", 
                         "3*float32", "3*float64"]


vars_file = "VARS"

patch_size = (15, 33, 71)

# These are the procs configuration that will be used for different core counts
# (first is the default)
procs_conf = dict()
procs_conf[2] = [(2,1,1)]
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
