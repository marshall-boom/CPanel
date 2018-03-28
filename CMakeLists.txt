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

cmake_minimum_required(VERSION 3.0)

project(CPANEL
        VERSION 0.2.1
        LANGUAGES CXX)

# prevent in source builds
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif ()

# add the project specific cmake files to module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/cmake-util
                      ${PROJECT_SOURCE_DIR}/cmake)

################################################################################
# Setup CPanel with version numbers
################################################################################
math(EXPR CPANEL_NUMBER_VERSION "10000*${CPANEL_VERSION_MAJOR}
                                +100*${CPANEL_VERSION_MINOR}
                                +${CPANEL_VERSION_PATCH}")

configure_file(${CPANEL_SOURCE_DIR}/cmake/cpanel_versions.hpp.in
               ${CPANEL_BINARY_DIR}/include/cpanel/versions.h @ONLY)

set(CPANEL_INCLUDE_DIR ${CPANEL_SOURCE_DIR}/CPanel)

include_directories(${CPANEL_INCLUDE_DIR})
include_directories(${CPANEL_BINARY_DIR}/include)

################################################################################
# Compiler Shenanigans to Ensure IEEE
#list(APPEND CMAKE_C_FLAGS  "-std=c99 -mfpmath=sse -msse2")
################################################################################

################################################################################
# Find Dependencies
################################################################################
find_package(Eigen3_Local REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(GTest_Local 1.8.0 REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include(CTest)

################################################################################
# Add subdirectories
################################################################################

add_subdirectory(octree)
set(OCTREE_INCLUDE_DIR ${CPANEL_SOURCE_DIR}/octree/include)
include_directories(${OCTREE_INCLUDE_DIR})

add_subdirectory(CPanelTests)
add_subdirectory(CPanel)

################################################################################
# Configure doxygen
################################################################################
find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(${CPANEL_SOURCE_DIR}/cmake/Doxyfile.in
                 ${CPANEL_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    COMMAND ${DOXYGEN_EXECUTABLE} ${CPANEL_BINARY_DIR}/Doxyfile
                    SOURCES ${CPANEL_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CPANEL_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen")
endif()

