################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# This program and the accompanying materials are made
# available under the terms of the Eclipse Public License 2.0
# which is available at https://www.eclipse.org/legal/epl-2.0/
#
# See LICENSE.md file in the project root for full license information.
#
# SPDX-License-Identifier: EPL-2.0
#
# Contributors:
#    David D. Marshall - initial code and implementation
################################################################################

project(OCTREE)

set(OCTREE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# include common test header files
include_directories(${OCTREE_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})


add_subdirectory(tests)
#####add_subdirectory(examples)

## Any other builds needed for this project would go here
