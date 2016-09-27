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
    double coreOverlap = 1.5; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
public:
    double radius;
    Eigen::Vector3d velOn; // Used for FMM. 

    wakePanel* parentPanel;
    
    particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius, Eigen::Vector3d previousVelInfl, Eigen::Vector3d previousStrengthUpdate);
    
    void setPrevVelInfl(Eigen::Vector3d vel) {previousVelInfl = vel;};
    Eigen::Vector3d getPrevVelInfl() {return previousVelInfl;};
    
    void setprevStrengthUpdate(Eigen::Vector3d strengthUpdate) {previousStrengthUpdate = strengthUpdate;};
    Eigen::Vector3d getprevStrengthUpdate() {return previousStrengthUpdate;};
    
    // Make public for faster access
    Eigen::Vector3d pos, strength;
    void setPos(Eigen::Vector3d newPos) {pos = newPos;};
    void setStrength(Eigen::Vector3d newStrength) {strength = newStrength;};
    
    Eigen::Vector3d partVelInfl(const Eigen::Vector3d &POI);
    Eigen::Vector3d partVelInflGaussian(particle* part);
    Eigen::Vector3d partVelInflGaussian(const Eigen::Vector3d &POI);
    Eigen::Vector3d partStrengthUpdate(particle* part);
    Eigen::Vector3d partStrengthUpdate(bodyPanel* bPan);
    
    Eigen::Vector3d vortexStretchingGaussian(particle* part);
    Eigen::Vector3d viscousDiffusionGaussian(particle* part);
    
    Eigen::Vector3d ptSourceStretching(bodyPanel* bPan);
    Eigen::Vector3d ptDoubletStretching(bodyPanel* bPan);
    
    void setParentWake(wakePanel* pan){parentPanel = pan;};

};





#endif /* defined(__CPanel___Unstructured_Panel_Code__particle__) */
