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

project(OCTREE_TESTS)

# set what libraries need to be used
set(_OCTREE_TESTS_LIBRARIES_LIST ${GTEST_BOTH_LIBRARIES})
set(_OCTREE_TESTS_DEPENDENCIES_LIST ${GTEST_DEPENDENCY_NAME} ${EIGEN3_DEPENDENCY_NAME})

# add the octree test
set (OCTREE_TEST_EXE octree_test)
set (OCTREE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/octree_test.cpp)
set (OCTREE_TEST_HEADER ${OCTREE_INCLUDE_DIR}/member.h
                        ${OCTREE_INCLUDE_DIR}/node.h
                        ${OCTREE_INCLUDE_DIR}/octree.h)

add_executable(${OCTREE_TEST_EXE} ${OCTREE_TEST_SOURCE} ${OCTREE_TEST_HEADER})
target_link_libraries(${OCTREE_TEST_EXE} ${_OCTREE_TESTS_LIBRARIES_LIST})
add_dependencies(${OCTREE_TEST_EXE} ${_OCTREE_TESTS_DEPENDENCIES_LIST})
add_test(NAME ${OCTREE_TEST_EXE} COMMAND $<TARGET_FILE:${OCTREE_TEST_EXE}> --text)

