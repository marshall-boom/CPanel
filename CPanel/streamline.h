/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
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

#ifndef __CPanel___Unstructured_Panel_Code__streamline__
#define __CPanel___Unstructured_Panel_Code__streamline__

#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include "bodyPanel.h"
#include "wakePanel.h"
#include "geometry.h"

class streamline
{
    std::vector<Eigen::Vector3d> pnts;
    std::vector<Eigen::Vector3d> velocities;
    Eigen::Vector3d Vinf;
    double PG;
    geometry* geom;
    
    Eigen::Matrix<double,1,6> coeff5,coeff4;
    
    
    Eigen::Vector3d rkf(const Eigen::Vector3d &x0,double dt,double &error);
    
public:
    streamline(const Eigen::Vector3d &startPnt, double xMax, double tol, const Eigen::Vector3d &Vinf, double PG, geometry* geom);
    
    std::vector<Eigen::Vector3d> getPnts() {return pnts;}
    std::vector<Eigen::Vector3d> getVelocities() {return velocities;}
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__streamline__) */
