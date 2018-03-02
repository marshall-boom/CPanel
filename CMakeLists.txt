################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. Licensed under the MIT License. See LICENSE file in the
# project root for full license information.
################################################################################

cmake_minimum_required(VERSION 3.0)

project(CPANEL
        VERSION 1.0.0
        LANGUAGES CXX)

# prevent in source builds
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif ()

# add the project specific cmake files to module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

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
####find_package(<package-name>)

################################################################################
# Find and configure Eigen3
################################################################################

# figure out what the user wants
if (DEFINED  EIGEN3_VERSION_REQ)
  CHECK_VERSION_REQ(_eli_eigen3_use ${ELI_EIGEN3_VERSION_REQ})
else()
  set(_eli_eigen3_use 0)
endif()

# get the information for CPPTest
if (_eli_eigen3_use EQUAL 0)
  find_package(Eigen3 REQUIRED)
  add_custom_target(Eigen3) # using system version means need to add empty target for dependencies
elseif (_eli_eigen3_use EQUAL 1)
  find_package(Eigen3_Local REQUIRED)
elseif (_eli_eigen3_use EQUAL 2)
  find_package(Eigen3_Local ${ELI_EIGEN3_VERSION_REQ} REQUIRED)
else()
  message(FATAL_ERROR " Unknown Eigen flag ->${ELI_EIGEN3_VERSION_REQ}<- requested.")
endif()
message(STATUS "Eigen3 Variables:\n\tEIGEN3_FOUND=${EIGEN3_FOUND}\n\tEIGEN3_INCLUDE=${EIGEN3_INCLUDE_DIR}\n\tEIGEN3_VERSION=${EIGEN3_VERSION}\n\tEIGEN3_VERSION_OK=${EIGEN3_VERSION_OK}")

# add the Eigen3 include directories
include_directories(${EIGEN3_INCLUDE_DIR})

################################################################################
# Setup Google Test
################################################################################
find_package(GTest_Local 1.8.0 REQUIRED)
include(CTest)


#########include(${CPANEL_SOURCE_DIR}/cmake/googletest-setup.cmake)

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

