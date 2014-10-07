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

mpirun -n 8 ./PIDX-analysis-1 -g 64x64x64 -l 32x32x32 -f A1 -t 5 -b 256

mpirun -n 16 ./PIDX-analysis-1 -g 128x64x64 -l 32x32x32 -f B1 -t 5 -b 256

mpirun -n 32 ./PIDX-analysis-1 -g 128x128x64 -l 32x32x32 -f C1 -t 5 -b 256

mpirun -n 64 ./PIDX-analysis-1 -g 128x128x128 -l 32x32x32 -f D1 -t 5 -b 256


mpirun -n 8 ./PIDX-analysis-2 -g 64x64x64 -l 32x32x32 -f A2 -t 5 -b 256

mpirun -n 16 ./PIDX-analysis-2 -g 128x64x64 -l 32x32x32 -f B2 -t 5 -b 256

mpirun -n 32 ./PIDX-analysis-2 -g 128x128x64 -l 32x32x32 -f C2 -t 5 -b 256

mpirun -n 64 ./PIDX-analysis-2 -g 128x128x128 -l 32x32x32 -f D2 -t 5 -b 256
