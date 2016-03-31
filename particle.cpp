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
    
    
    Eigen::Vector3d velInfl;
    double coreOverlap = 1.3; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
    h = radius; // Research more into what's going on here
    double sigma = coreOverlap*h;
    
    Eigen::Vector3d dist = pos-POI;
    
    // ** Regularized Influence ** //
    //velInfl= -1/(4*M_PI)*(pow(dist.norm(),2) + 2.5*pow(sigma,2))/(pow(pow(dist.norm(),2) + pow(sigma,2),2.5)) * dist.cross(strength);
    
    // ** Singular Influence ** //
    velInfl = -1/(4*M_PI)*(1/pow(dist.norm(),3))*dist.cross(strength);
    
    return velInfl;
};

double particle::partPotInfl(const Eigen::Vector3d &POI){
    
    double pot = (1/(4*M_PI))*(POI-pos).norm()/1000;
    
    return pot;
}


void particle::stepParticle(Eigen::Vector3d &velOnPart){
    
    pos += velOnPart*dt;
    
}


