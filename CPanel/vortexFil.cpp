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

#include "vortexFil.h"

vortexFil::vortexFil(Eigen::Vector3d pp1, Eigen::Vector3d pp2, double sstrength, wakePanel* pparentPan)
  : p1(pp1), p2(pp2), strength(sstrength), parentPan(pparentPan) {}

Eigen::Vector3d vortexFil::velInfl(Eigen::Vector3d POI){
    
    // VSAero doublet velocity influence formulation
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(3);
    
    Eigen::Vector3d a,b,s;
    
    a = POI-p1;
    b = POI-p2;
    s = p2-p1;
    
    vel = parentPan->vortexV(a,b,s);
    
    return vel*strength/(4*M_PI);
}

//Eigen::Vector3d vortexFil::moveNode(Eigen::Vector3d pos, double dt, std::vector<double> bodyKin){
//    
//    
//    Eigen::Vector3d localMovement;
//    
//    // U = U3 + (-q*z + r*y)
//    localMovement.x() = bodyKin[0] - bodyKin[4]*pos.z() + bodyKin[5]*pos.y();
//    
//    // V = V3 + (-r*x + p*z)
//    localMovement.y() = bodyKin[1] - bodyKin[5]*pos.x() + bodyKin[3]*pos.z();
//    
//    // W = W3 + (-p*y + q*x)
//    localMovement.z() = bodyKin[2] - bodyKin[3]*pos.y() + bodyKin[4]*pos.x();
//    
//    return  pos + localMovement*dt;
//    
//}
//
//void vortexFil::moveFilament(std::vector<double> bodyKin, double dt){
//    
//    Eigen::Vector3d newP1 = moveNode(p1, dt, bodyKin);
//    Eigen::Vector3d newP2 = moveNode(p2, dt, bodyKin);
//
//    p1 = newP1;
//    p2 = newP2;
//}
