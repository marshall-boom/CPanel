/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Chris Satterwhite - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel__panelOctree__
#define __CPanel__panelOctree__

#include <iostream>
#include "octree.h"
#include "panel.h"

class panel;

class panelOctree : public octree<panel>
{
    
public:
    panelOctree() : octree() {}
    
    Eigen::Vector3d findRefPoint(const panel &obj);
    
//    panel* getClosestPanel(const Eigen::Vector3d &pnt);
};

#endif /* defined(__CPanel__panelOctree__) */
