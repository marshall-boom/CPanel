################################################################################
# Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. Licensed under the MIT License. See LICENSE file in the
# project root for full license information.
################################################################################

project(CPANEL_APP)

set(CPANEL_LIBRARIES_LIST boost_system boost_filesystem)

# add the CPanel application
set (CPANEL_EXE cpanel)
set (CPANEL_SOURCE ${PROJECT_SOURCE_DIR}/cpanel.cpp
                   ${PROJECT_SOURCE_DIR}/bodyPanel.cpp
                   ${PROJECT_SOURCE_DIR}/cpFile.cpp
                   ${PROJECT_SOURCE_DIR}/VTUfile.cpp
                   ${PROJECT_SOURCE_DIR}/bodyStreamline.cpp
                   ${PROJECT_SOURCE_DIR}/cpNode.cpp
                   ${PROJECT_SOURCE_DIR}/octreeFile.cpp
                   ${PROJECT_SOURCE_DIR}/wake.cpp
                   ${PROJECT_SOURCE_DIR}/chtlsnd.cpp
                   ${PROJECT_SOURCE_DIR}/edge.cpp
                   ${PROJECT_SOURCE_DIR}/panel.cpp
                   ${PROJECT_SOURCE_DIR}/wakeLine.cpp
                   ${PROJECT_SOURCE_DIR}/convexHull.cpp
                   ${PROJECT_SOURCE_DIR}/geometry.cpp
                   ${PROJECT_SOURCE_DIR}/panelOctree.cpp
                   ${PROJECT_SOURCE_DIR}/wakePanel.cpp
                   ${PROJECT_SOURCE_DIR}/CPanelMgr.cpp
                   ${PROJECT_SOURCE_DIR}/inputParams.cpp
                   ${PROJECT_SOURCE_DIR}/streamline.cpp
                   ${PROJECT_SOURCE_DIR}/cpCase.cpp
                   ${PROJECT_SOURCE_DIR}/liftingSurf.cpp
                   ${PROJECT_SOURCE_DIR}/surface.cpp)
set (OCTREE_TEST_HEADER ${PROJECT_SOURCE_DIR}/bodyPanel.h
                        ${PROJECT_SOURCE_DIR}/cpCase.h
                        ${PROJECT_SOURCE_DIR}/inputParams.h
                        ${PROJECT_SOURCE_DIR}/streamline.h
                        ${PROJECT_SOURCE_DIR}/wakePanel.h
                        ${PROJECT_SOURCE_DIR}/bodyStreamline.h
                        ${PROJECT_SOURCE_DIR}/cpFile.h
                        ${PROJECT_SOURCE_DIR}/liftingSurf.h
                        ${PROJECT_SOURCE_DIR}/surface.h
                        ${PROJECT_SOURCE_DIR}/chtlsnd.h
                        ${PROJECT_SOURCE_DIR}/cpNode.h
                        ${PROJECT_SOURCE_DIR}/octreeFile.h
                        ${PROJECT_SOURCE_DIR}/VTUfile.h
                        ${PROJECT_SOURCE_DIR}/convexHull.h
                        ${PROJECT_SOURCE_DIR}/edge.h
                        ${PROJECT_SOURCE_DIR}/panel.h
                        ${PROJECT_SOURCE_DIR}/wake.h
                        ${PROJECT_SOURCE_DIR}/CPanelMgr.h
                        ${PROJECT_SOURCE_DIR}/geometry.h
                        ${PROJECT_SOURCE_DIR}/panelOctree.h
                        ${PROJECT_SOURCE_DIR}/wakeLine.h
                        ${OCTREE_INCLUDE_DIR}/member.h
                        ${OCTREE_INCLUDE_DIR}/node.h
                        ${OCTREE_INCLUDE_DIR}/octree.h)

add_executable(${CPANEL_EXE} ${CPANEL_SOURCE} ${CPANEL_HEADER})
target_link_libraries(${CPANEL_EXE} ${CPANEL_LIBRARIES_LIST})
#add_dependencies(${CPANEL_EXE} ${CPANEL_LIBRARIES_LIST})
