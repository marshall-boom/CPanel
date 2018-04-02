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

#ifndef __CPanel___Unstructured_Panel_Code__bodyStreamline__
#define __CPanel___Unstructured_Panel_Code__bodyStreamline__

#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"
#include "cpNode.h"

class geometry;

class bodyStreamline
{
    std::vector<Eigen::Vector3d> pnts;
    std::vector<Eigen::Vector3d> velocities;
    Eigen::Vector3d Vinf;
    double Vmag;
    double PG;
    int marchDir;
    
//    Eigen::Vector3d trailingEdgePnt(bodyPanel* p);
    edge* edgeIntersection(bodyPanel* pan,const Eigen::Vector3d &pnt, double pntPot, Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge, double maxAngle, edge* lastEdge, bool &stagFlag);
    
    bool stagnationPnt(const Eigen::Vector3d vel, const Eigen::Vector3d &velOld, double maxAngle);
    
    double safeInvCos(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
public:
    bodyStreamline(Eigen::Vector3d startPnt, bodyPanel* startPan, const Eigen::Vector3d &Vinf, double PG, geometry *geom, int pntsPerPanel, bool marchFwd);
    
    std::vector<Eigen::Vector3d> getPnts() {return pnts;}
    std::vector<Eigen::Vector3d> getVelocities() {return velocities;}
    
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__bodyStreamline__) */
