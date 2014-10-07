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
#Simulation Style IO
# 4 Variables : 
#var1 : 1 sample
#var2 : 1 sample
#var3 : 3 sample
#var4 : 11 sample

#Will Create 8 Files (128 Aggregators)
mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 8FILES_Benchmark1 -t 10 -b 512

#Will Create 16 Files (256 Aggregators)
mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 16FILES_Benchmark1 -t 10 -b 256

#Will Create 32 Files (512 Aggregators)
mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 32FILES_Benchmark1 -t 10 -b 128

#Will Create 64 Files (1024 Aggregators)
mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 64FILES_Benchmark1 -t 10 -b 64


#Regular IO 
# 1 Variabe
#var1 : 1 sample

#Will Create 128 Files (128 Aggregators)
mpirun -n 4096 ./k-benchmark-2 -g 512x512x512 -l 32x32x32 -f 128FILES_Benchmark2 -t 10 -b 32

#Will Create 64 Files (64 Aggregators)
mpirun -n 4096 ./k-benchmark-2 -g 512x512x512 -l 32x32x32 -f 64FILES_Benchmark2 -t 10 -b 64

#Will Create 32 Files (32 Aggregators)
mpirun -n 4096 ./k-benchmark-2 -g 512x512x512 -l 32x32x32 -f 32FILES_Benchmark2 -t 10 -b 128

#Will Create 16 Files (16 Aggregators)
mpirun -n 4096 ./k-benchmark-2 -g 512x512x512 -l 32x32x32 -f 16FILES_Benchmark2 -t 10 -b 256

#Will Create 8 Files (8 Aggregators)
mpirun -n 4096 ./k-benchmark-2 -g 512x512x512 -l 32x32x32 -f 8FILES_Benchmark2 -t 10 -b 512
