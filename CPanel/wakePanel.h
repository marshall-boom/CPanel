//
//  wakePanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wakePanel__
#define __CPanel__wakePanel__

#include <iostream>
#include "panelOctree.h"
#include "wakeLine.h"
#include "particle.h"
#include "vortexFilament.h"
#include "panel.h"

class wake;
class bodyPanel;
class edge;
class panel;
class vortexFilament;
class particle;

class wakePanel : public panel
{
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    wake* parentWake;
    
    // VP mods
    wakePanel* bufferParent = nullptr;
    double prevStrength = 0;
    vortexFilament* vortFil;
    
public:
    wakePanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake, int surfID);
    
    wakePanel(const wakePanel &copy);
    double getPrevStrength() {return prevStrength;};
    void setPrevStrength(double previousStrength) {prevStrength = previousStrength;};
    
    void setTEpanel();
    void setUpper(bodyPanel* up);
    void setLower(bodyPanel* lp);
    void setParentPanels(bodyPanel* upper, bodyPanel* lower);
    void setParentWake(wake* w) {parentWake = w;}
//    wakePanel* makeVortexSheet();
    void interpPanels(std::vector<bodyPanel*> &interpP, double &interpC);
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    void setMu();
    void setMu(double value){doubletStrength = value;}
    void setStrength();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
    edge* getTE();
    bool isTEpanel() {return TEpanel;}
    
    void setBufferParent( wakePanel* wPan ){ bufferParent = wPan;}
    wakePanel* getBufferParent(){return bufferParent;};
    bool isSecondRow = false;

    std::vector<edge*> edgesInOrder();
    std::vector<cpNode*> pointsInOrder();
    
    vortexFilament* getVortFil(){return vortFil;};
    void setVortFil(vortexFilament* filament){vortFil = filament;}

    double getPartRadius(Eigen::Vector3d Vinf, double &dt);
    Eigen::Vector3d partStretching(particle* part);


};

#endif /* defined(__CPanel__wakePanel__) */
