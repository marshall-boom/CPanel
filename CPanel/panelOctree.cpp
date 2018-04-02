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

#include "panelOctree.h"

Eigen::Vector3d panelOctree::findRefPoint(const panel &obj)
{
    return obj.getCenter();
}

//panel* panelOctree::getClosestPanel(const Eigen::Vector3d &pnt)
//{
//    // Check panel is within octree
//    if (isInsideOctree(pnt))
//    {
//        std::vector<panel*> pans = findNodeContainingPnt(pnt)->getMembers();
//        
//        panel* closestPan;
//        double dist;
//        double minDist = 1000000;
//        for (int i=0; i<pans.size(); i++)
//        {
//            dist = (pans[i]->getCenter()-pnt).norm();
//            if (dist < minDist)
//            {
//                closestPan = pans[i];
//                minDist = dist;
//            }
//        }
//        if (closestPan->isOnPanel(
//    }
//    
//}
