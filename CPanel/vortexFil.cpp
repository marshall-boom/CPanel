//
//  vortexFil.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 9/13/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#include "vortexFil.h"

vortexFil::vortexFil(Eigen::Vector3d p1, Eigen::Vector3d p2, double strength, wakePanel* parentPan) : p1(p1), p2(p2), strength(strength), parentPan(parentPan) {};

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
