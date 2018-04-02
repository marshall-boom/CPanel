/*******************************************************************************
 * Copyright (c) 2016 Connor Sousa
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Connor Sousa - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

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
