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
//#include <math.h>

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
    void setStrength(Eigen::Vector3d newStrength) {strength = newStrength;};
    Eigen::Vector3d partVelInfl(const Eigen::Vector3d &POI);
    Eigen::Vector3d partVelInflGaussian(const Eigen::Vector3d &POI);
    Eigen::Vector3d partStretching(particle* part);
    Eigen::Matrix3d partStretchingGaussian(particle* part);
};





#endif /* defined(__CPanel___Unstructured_Panel_Code__particle__) */
