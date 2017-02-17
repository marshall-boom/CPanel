//
//  vortexFilament.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 9/13/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__vortexFilament__
#define __CPanel___Unstructured_Panel_Code__vortexFilament__

#include <iostream>
#include "wakePanel.h"
#include <Eigen/Dense>

class wakePanel;

class vortexFilament{
    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    double strength;
    wakePanel* parentPan;
    
public:
    vortexFilament(Eigen::Vector3d p1, Eigen::Vector3d p2, double strength, wakePanel* parentPan);
    Eigen::Vector3d velInfl(Eigen::Vector3d POI);
    void setStrength(double fStrength){strength = fStrength;};
    double getStrength(){return strength;};
    Eigen::Vector3d getP1(){return p1;};
    Eigen::Vector3d getP2(){return p2;};
};






#endif /* defined(__CPanel___Unstructured_Panel_Code__vortexFilament__) */
