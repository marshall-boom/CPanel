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

project(CPANEL_TESTS)

# set what libraries need to be used
set(_CPANEL_TESTS_LIBRARIES_LIST ${GTEST_BOTH_LIBRARIES})
set(_CPANEL_TESTS_DEPENDENCIES_LIST ${GTEST_DEPENDENCY_NAME} ${EIGEN3_DEPENDENCY_NAME})

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
add_dependencies(${GEOMETRY_TEST_EXE} ${_CPANEL_TESTS_DEPENDENCIES_LIST})
add_test(NAME ${GEOMETRY_TEST_EXE} COMMAND $<TARGET_FILE:${GEOMETRY_TEST_EXE}>)

# add the influence test
set (INFLUENCE_TEST_EXE influence_test)
set (INFLUENCE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/influence_test.cpp)
set (INFLUENCE_TEST_HEADER ${OCTREE_INCLUDE_DIR}/member.h)

add_executable(${INFLUENCE_TEST_EXE} ${INFLUENCE_TEST_SOURCE} ${INFLUENCE_TEST_HEADER})
target_link_libraries(${INFLUENCE_TEST_EXE} ${_CPANEL_TESTS_LIBRARIES_LIST})
add_dependencies(${INFLUENCE_TEST_EXE} ${_CPANEL_TESTS_DEPENDENCIES_LIST})
add_test(NAME ${INFLUENCE_TEST_EXE} COMMAND $<TARGET_FILE:${INFLUENCE_TEST_EXE}>)

