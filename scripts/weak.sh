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

#Less Aggregator
mpirun -n 1024 ./k-benchmark-1 -g 512x256x256 -l 32x32x32 -f 1K_A -t 15 -b 512

mpirun -n 2048 ./k-benchmark-1 -g 512x512x256 -l 32x32x32 -f 2K_A -t 15 -b 512

mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 4K_A -t 15 -b 512

mpirun -n 8192 ./k-benchmark-1 -g 1024x512x512 -l 32x32x32 -f 8K_A -t 15 -b 512

mpirun -n 16384 ./k-benchmark-1 -g 1024x1024x512 -l 32x32x32 -f 16K_A -t 15 -b 512

mpirun -n 32768 ./k-benchmark-1 -g 1024x1024x1024 -l 32x32x32 -f 32K_A -t 15 -b 512

mpirun -n 65536 ./k-benchmark-1 -g 2048x1024x1024 -l 32x32x32 -f 64K_A -t 15 -b 512




#More Aggregator
mpirun -n 1024 ./k-benchmark-1 -g 512x256x256 -l 32x32x32 -f 1K_B -t 15 -b 64

mpirun -n 2048 ./k-benchmark-1 -g 512x512x256 -l 32x32x32 -f 2K_B -t 15 -b 64

mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 4K_B -t 15 -b 64

mpirun -n 8192 ./k-benchmark-1 -g 1024x512x512 -l 32x32x32 -f 8K_B -t 15 -b 64

mpirun -n 16384 ./k-benchmark-1 -g 1024x1024x512 -l 32x32x32 -f 16K_B -t 15 -b 64

mpirun -n 32768 ./k-benchmark-1 -g 1024x1024x1024 -l 32x32x32 -f 32K_B -t 15 -b 64

mpirun -n 65536 ./k-benchmark-1 -g 2048x1024x1024 -l 32x32x32 -f 64K_B -t 15 -b 64

