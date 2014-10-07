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

mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 4K_1 -t 15 -b 256

mpirun -n 4096 ./k-benchmark-1 -g 512x512x512 -l 32x32x32 -f 4K_2 -t 15 -b 512

mpirun -n 8192 ./k-benchmark-1 -g 1024x512x512 -l 32x32x32 -f 8K_1 -t 15 -b 256

mpirun -n 8192 ./k-benchmark-1 -g 1024x512x512 -l 32x32x32 -f 8K_2 -t 15 -b 512

