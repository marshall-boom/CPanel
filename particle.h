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
    
    double Uinf = 0.1;
    double dt = 0.1; //GET FROM INPUT
    double h = Uinf*dt; // Average distance between particles is taken to be the starting distance between released particles. Should modify this to find the smaller/larger of either particles shed in the stream direction or the dist between particles on the panel.
    
    //      dt = K*chord/Uinf
    
    
public:
    particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius);
    
    // Model after cpNode
    Eigen::Vector3d getPos() {return pos;};
    void setPos(Eigen::Vector3d nPos) {pos = nPos;};
    
    Eigen::Vector3d getStrength() {return strength;};
    
    Eigen::Vector3d partVelInfl(const Eigen::Vector3d &POI);
    double partPotInfl(const Eigen::Vector3d &POI);
    void stepParticle(Eigen::Vector3d &velOnPart);
    
};





#endif /* defined(__CPanel___Unstructured_Panel_Code__particle__) */
