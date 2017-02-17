//
//  vortexFil.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 9/13/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#include "vortexFilament.h"

vortexFilament::vortexFilament(Eigen::Vector3d p1, Eigen::Vector3d p2, double strength, wakePanel* parentPan) : p1(p1), p2(p2), strength(strength), parentPan(parentPan) {};

Eigen::Vector3d vortexFilament::velInfl(Eigen::Vector3d POI){
    
    // VSAero doublet velocity influence formulation
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(3);
    
    Eigen::Vector3d a,b,s;
    
    a = POI-p1;
    b = POI-p2;
    s = p2-p1;
    
//    vel = parentPan->vortexV(a,b,s); /// This was an attempt to reuse the code in 'panel.cpp' but the linking error was impossibly frustrating. Uncomment in you're up for a challenge (or know more than the basics).

    double core = 0.05;
    vel = (a.cross(b)*(a.norm()+b.norm()))/(a.norm()*b.norm()*((a.norm()*b.norm())+a.dot(b))+(pow(core,2)));
    
    return vel*strength/(4*M_PI);
}
