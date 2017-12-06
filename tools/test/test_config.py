# User settings
# Set write_executable and read_executable according to your build directory
write_executable = "../../build/examples/idx_write"
read_executable = "../../build/examples/idx_read"
mpirun="mpirun"

# enable debug prints 
debug_print=1

# travis mode does not execute any test
# but it creates a bash script to run them 
travis_mode = 0

#var_types = ["1*float32", "1*int32", "1*float64", 
#             "2*float32", "2*int32", "2*float64", 
#             "3*float32", "3*int32", "3*float64"]

var_types = ["1*float32"]

vars_file = "VARS"

patch_size = (12, 12, 12)

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
