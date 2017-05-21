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
#include <numeric>
#include "panel.h"
#include "panelOctree.h"
#include "wakeLine.h"

class wake;
class bodyPanel;
class particle;
class edge;
class vortexFil;

class wakePanel : public panel
{
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    wake* parentWake = nullptr;
    double prevStrength = 0;
    vortexFil* vortFil;
    wakePanel* bufferParent;
    
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
    void setMu(double strength); //2BW
    void setStrength();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
    edge* getTE();
    bool isTEpanel() {return TEpanel;}
    
    Eigen::Vector3d panToPartStrengthT1();
    Eigen::Vector3d panToPartStrength();
    
    //    std::vector<Eigen::Vector3d> vortexRingVectors();
    double getPartRadius(Eigen::Vector3d Vinf, double &dt);
    bool isSecondRow = false;
    
    
    std::vector<int> sort_indexes(std::vector<double> &v);
    Eigen::Vector3d partStretching(particle* part);
    
    std::vector<cpNode*> pointsInOrder();
    std::vector<edge*> edgesInOrder();
    
    vortexFil* getVortFil(){return vortFil;};
    void setVortFil(vortexFil* filament){vortFil = filament;};
    
    Eigen::Vector3d partSeedPt(Eigen::Vector3d &Vinf, double &dt);
    void setBufferParent(wakePanel* pan){bufferParent = pan;};
    wakePanel* getBufferParent(){return bufferParent;};
    
    Eigen::Vector3d edgeStrength( edge* curEdge, int edgeNum);
    
};

#endif /* defined(__CPanel__wakePanel__) */
