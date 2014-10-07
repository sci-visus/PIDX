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

#Creating IDX file with 3 time steps with global volume 128x128x128 with 8 processes each with bounds 64x64x64


#One file Per Process (Therefore 8 files in all)
mpirun -n 8 ./test-PIDX-writer -g 128x128x128 -l 64x64x64 -f One_File_Per_Process -t 3 -r 0 -a 1


#One File In all (Therefore 1 file)
mpirun -n 8 ./test-PIDX-writer -g 128x128x128 -l 64x64x64 -f One_File_In_All -t 3 -r 0 -a 1
