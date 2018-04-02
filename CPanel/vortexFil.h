/*******************************************************************************
 * Copyright (c) 2016 Connor Sousa
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
 *    Connor Sousa - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel___Unstructured_Panel_Code__vortexFil__
#define __CPanel___Unstructured_Panel_Code__vortexFil__

#include <iostream>
#include "wakePanel.h"
#include "geometry.h"
#include <Eigen/Dense>

class vortexFil{
    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    double strength;
    wakePanel* parentPan;
    
public:
    vortexFil(Eigen::Vector3d p1, Eigen::Vector3d p2, double strength, wakePanel* parentPan);
    Eigen::Vector3d velInfl(Eigen::Vector3d POI);
    void setStrength(double fStrength){strength = fStrength;};
    double getStrength(){return strength;};
    Eigen::Vector3d getP1(){return p1;};
    Eigen::Vector3d getP2(){return p2;};
    
    Eigen::Vector3d moveNode(Eigen::Vector3d pos, double dt, std::vector<double> bodyKin);
    void moveFilament(std::vector<double> bodyKin, double dt);
};






#endif /* defined(__CPanel___Unstructured_Panel_Code__vortexFil__) */
