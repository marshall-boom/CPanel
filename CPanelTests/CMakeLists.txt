################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. Licensed under the MIT License. See LICENSE file in the
# project root for full license information.
################################################################################

project(CPANEL_TESTS)

# set what libraries need to be used
set(_CPANEL_TESTS_LIBRARIES_LIST gtest_main)

# add the geometry test
set (GEOMETRY_TEST_EXE geometry_test)
set (GEOMETRY_TEST_SOURCE ${PROJECT_SOURCE_DIR}/geometry_test.cpp)
set (GEOMETRY_TEST_HEADER ${CPANEL_INCLUDE_DIR}/geometry.h
                          ${CPANEL_INCLUDE_DIR}/surface.h
                          ${CPANEL_INCLUDE_DIR}/liftingSurf.h
                          ${CPANEL_INCLUDE_DIR}/wake.h
                          ${CPANEL_INCLUDE_DIR}/bodyPanel.h
                          ${CPANEL_INCLUDE_DIR}/wakePanel.h)

add_executable(${GEOMETRY_TEST_EXE} ${GEOMETRY_TEST_SOURCE} ${GEOMETRY_TEST_HEADER})
target_link_libraries(${GEOMETRY_TEST_EXE} ${_CPANEL_TESTS_LIBRARIES_LIST})
add_dependencies(${GEOMETRY_TEST_EXE} ${_CPANEL_TESTS_LIBRARIES_LIST})
add_test(NAME ${GEOMETRY_TEST_EXE} COMMAND $<TARGET_FILE:${GEOMETRY_TEST_EXE}>)

# add the influence test
set (INFLUENCE_TEST_EXE influence_test)
set (INFLUENCE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/influence_test.cpp)
set (INFLUENCE_TEST_HEADER ${OCTREE_INCLUDE_DIR}/member.h)

add_executable(${INFLUENCE_TEST_EXE} ${INFLUENCE_TEST_SOURCE} ${INFLUENCE_TEST_HEADER})
target_link_libraries(${INFLUENCE_TEST_EXE} ${_CPANEL_TESTS_LIBRARIES_LIST})
add_dependencies(${INFLUENCE_TEST_EXE} ${_CPANEL_TESTS_LIBRARIES_LIST})
add_test(NAME ${INFLUENCE_TEST_EXE} COMMAND $<TARGET_FILE:${INFLUENCE_TEST_EXE}>)

