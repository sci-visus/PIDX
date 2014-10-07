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


###########################################
# PIDX_ADD_EXECUTABLE
###########################################

MACRO(PIDX_ADD_EXECUTABLE targetname files)
  #MESSAGE("Adding executable " ${targetname} " from sources: " ${files})                        
  ADD_EXECUTABLE(${targetname} "MACOSX_BUNDLE" ${files})
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE C)
  INSTALL(TARGETS ${targetname} 
                  RUNTIME DESTINATION bin
                  BUNDLE DESTINATION  bin)
ENDMACRO()


###########################################
# PIDX_ADD_LIBRARY
###########################################

MACRO(PIDX_ADD_LIBRARY targetname files)
  #MESSAGE("Adding library " ${targetname} " from sources: " ${files})                        
  ADD_LIBRARY(${targetname} ${files})
  SET_TARGET_PROPERTIES(${targetname} PROPERTIES LINKER_LANGUAGE C)
  INSTALL(TARGETS ${targetname} 
                  ARCHIVE DESTINATION lib
                  LIBRARY DESTINATION lib)
ENDMACRO()
