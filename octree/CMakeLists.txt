################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. Licensed under the MIT License. See LICENSE file in the
# project root for full license information.
################################################################################

project(OCTREE)

set(OCTREE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# include common test header files
include_directories(${OCTREE_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})


add_subdirectory(tests)
#####add_subdirectory(examples)

## Any other builds needed for this project would go here
