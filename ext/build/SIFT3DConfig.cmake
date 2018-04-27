################################################################################
# Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
################################################################################
# SIFT3D package configuration file. Defines the following variables:
#   SIFT3D_INCLUDE_DIRS - Include directories
#   SIFT3D_LIBRARIES - Library targets
################################################################################

# Import the targets, if not already imported
if (NOT TARGET reg)
        include (${CMAKE_CURRENT_LIST_DIR}/SIFT3D-targets.cmake)
endif() 

# Get the libraries
set (SIFT3D_LIBRARIES reg sift3D imutil)
set (SIFT3D_INCLUDE_DIRS include/sift3d)
