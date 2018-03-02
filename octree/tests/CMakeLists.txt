################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. Licensed under the MIT License. See LICENSE file in the
# project root for full license information.
################################################################################

project(OCTREE_TESTS)

# set what libraries need to be used
set(_OCTREE_TESTS_LIBRARIES_LIST gtest_main)

# add the octree test
set (OCTREE_TEST_EXE octree_test)
set (OCTREE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/octree_test.cpp)
set (OCTREE_TEST_HEADER ${OCTREE_INCLUDE_DIR}/member.h
                        ${OCTREE_INCLUDE_DIR}/node.h
                        ${OCTREE_INCLUDE_DIR}/octree.h)

add_executable(${OCTREE_TEST_EXE} ${OCTREE_TEST_SOURCE} ${OCTREE_TEST_HEADER})
target_link_libraries(${OCTREE_TEST_EXE} ${_OCTREE_TESTS_LIBRARIES_LIST})
add_test(NAME ${OCTREE_TEST_EXE} COMMAND $<TARGET_FILE:${OCTREE_TEST_EXE}> --text)

