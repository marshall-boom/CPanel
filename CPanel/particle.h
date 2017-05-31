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
#include "bodyPanel.h"
#include "wakePanel.h"
//#include <math.h>

class wakePanel;

class particle {
    
    Eigen::Vector3d previousVelInfl, previousStrengthUpdate; //used for Adams-Bashforth stepping routine
    double coreOverlap = 1.5; // More info: calebretta pp. 47
    
public:
    
    // Make public for faster access
    Eigen::Vector3d pos, strength;
    double radius;
    Eigen::Vector3d velOn; // Used for FMM.
    int shedTime;
    wakePanel* parentPanel;
    
    particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius, Eigen::Vector3d previousVelInfl, Eigen::Vector3d previousStrengthUpdate, int shedTime);
    
    void setPrevVelInfl(Eigen::Vector3d vel) {previousVelInfl = vel;};
    Eigen::Vector3d getPrevVelInfl() {return previousVelInfl;};
    
    void setprevStrengthUpdate(Eigen::Vector3d strengthUpdate) {previousStrengthUpdate = strengthUpdate;};
    Eigen::Vector3d getprevStrengthUpdate() {return previousStrengthUpdate;};
    
    void setPos(Eigen::Vector3d newPos) {pos = newPos;};
    void setStrength(Eigen::Vector3d newStrength) {strength = newStrength;};
    void setParentWake(wakePanel* pan){parentPanel = pan;};
    
    
    Eigen::Vector3d velInflAlgSmooth(const Eigen::Vector3d &POI);
    
    Eigen::Vector3d velInfl(particle* part);
    Eigen::Vector3d velInfl(const Eigen::Vector3d &POI);
    
    Eigen::Vector3d vortexStretching(particle* part);
    Eigen::Vector3d viscousDiffusion(particle* part);
    

//    Eigen::Vector3d partStrengthUpdate(particle* part);

    
//    Eigen::Vector3d partVelInfl(const Eigen::Vector3d &POI);
//    Eigen::Vector3d partStrengthUpdate(particle* part);
};





#endif /* defined(__CPanel___Unstructured_Panel_Code__particle__) */
