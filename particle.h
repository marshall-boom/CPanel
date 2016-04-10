//
//  particle.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 3/5/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__particle__
#define __CPanel___Unstructured_Panel_Code__particle__

#include <iostream>
#include <Eigen/Dense>

class particle {
    
    int nParticles;
    Eigen::Vector3d pos;
    Eigen::Vector3d strength;
    double radius;
    double h;
    
public:
    particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius);
    
    // Model after cpNode
    Eigen::Vector3d getPos() {return pos;};
    void setPos(Eigen::Vector3d newPos) {pos = newPos;};
    
    Eigen::Vector3d getStrength() {return strength;};
    Eigen::Vector3d partVelInfl(const Eigen::Vector3d &POI);
    void partStretching(particle* part);
};





#endif /* defined(__CPanel___Unstructured_Panel_Code__particle__) */
