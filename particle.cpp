//
//  particle.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 3/5/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#include "particle.h"
#include <Eigen/Geometry>
#include "math.h"


particle::particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius) : pos(pos), strength(strength), radius(radius) {};

Eigen::Vector3d particle::partVelInfl(const Eigen::Vector3d &POI){
    
    
    double coreOverlap = 1.3; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
    // Research more into what's going on here
    double sigma = coreOverlap*radius;
    
    Eigen::Vector3d dist = POI-pos;
    
    // ** Regularized Influence ** //
    return -1/(4*M_PI)*(pow(dist.norm(),2) + 2.5*pow(sigma,2))/(pow(pow(dist.norm(),2) + pow(sigma,2),2.5)) * dist.cross(this->strength);
    
    // ** Singular Influence ** //
//    return -1/(4*M_PI)*(1/pow(dist.norm(),3))*dist.cross(strength);
    
};

void particle::partStretching(particle* part){
    
//    std::cout << "Pre stretching: " << this->strength.x() << ", " << this->strength.y() << ", " << this->strength.z() << std::endl;
    double sigma = 1.3*radius; // coreOverlap of 1.3;
    
    double volP   =  4*M_PI/3*pow(this->radius,3);
    double volThis = 4*M_PI/3*pow(this->radius,3);

    Eigen::Vector3d dist = part->pos - this->pos;
    double nu = 1.983e-5;

    Eigen::Vector3d strengthChange = 1/(4*M_PI)*(((dist.norm()*dist.norm() + 2.5*sigma*sigma)/(pow(dist.norm()*dist.norm() + sigma*sigma,2.5))*part->strength).cross(this->strength) +
                                                 
                                                 3*((dist.norm()*dist.norm() + 3.5*sigma*sigma)/(pow(dist.norm()*dist.norm() + sigma*sigma,3.5)))*((part->strength).dot(dist.cross(this->strength)))*dist)+
    
                                                 105*nu*(pow(sigma,4)/(pow(dist.norm()*dist.norm()+sigma*sigma,4.5))*(volP*this->strength - volThis*part->getStrength()));
    
//    std::cout << "Stretching influence: " << strengthChange.x() << ", " << strengthChange.y() << ", " << strengthChange.z() << std::endl;

    this->strength += strengthChange;
}



