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
