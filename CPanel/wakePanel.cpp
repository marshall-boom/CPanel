//
//  wakePanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakePanel.h"
#include "wake.h"
#include "bodyPanel.h"
#include "edge.h"


wakePanel::wakePanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake,int surfID)
  : panel(nodes,pEdges,bezNorm,surfID), upperPan(nullptr), lowerPan(nullptr),
	TEpanel(false), parentWake(parentWake), vortFil(nullptr), bufferParent(nullptr)
{
    for (size_t i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addWakePan(this);
    }
}

wakePanel::wakePanel(const wakePanel &copy)
  : panel(copy), upperPan(copy.upperPan), lowerPan(copy.lowerPan), TEpanel(copy.TEpanel),
	parentWake(copy.parentWake), vortFil(copy.vortFil), bufferParent(copy.bufferParent) {}

void wakePanel::setTEpanel()
{
    TEpanel = true;
    parentWake->addTEPanel(this);
}
void wakePanel::setUpper(bodyPanel* up) {upperPan = up;}
void wakePanel::setLower(bodyPanel* lp) {lowerPan = lp;}

void wakePanel::interpPanels(std::vector<bodyPanel*> &interpPans, double &interpCoeff)
{
    wakeLine* wl1 = nullptr;
    wakeLine* wl2 = nullptr;
    std::vector<wakeLine*> wakeLines = parentWake->getWakeLines();
    if (center(1) < wakeLines[1]->getY())
    {
        wl1 = wakeLines[0];
        wl2 = wakeLines[1];
    }
    else if (center(1) >= wakeLines.end()[-1]->getY())
    {
        wl1 = wakeLines.end()[-2]; //Second to last wakeline
        wl2 = wakeLines.end()[-1]; //Last wakeline
    }
    else
    {
        for (size_t i=1; i<wakeLines.size()-1; i++)
        {
            if ((wakeLines[i]->getY() <= center(1) && wakeLines[i+1]->getY() > center(1)))
            {
                wl1 = wakeLines[i];
                wl2 = wakeLines[i+1];
            }
        }
    }
    assert(wl1 != nullptr || wl2 != nullptr);
    interpCoeff = (center(1)-wl1->getY())/(wl2->getY()-wl1->getY());
    
    interpPans[0] = wl1->getUpper();
    interpPans[1] = wl1->getLower();
    interpPans[2] = wl2->getUpper();
    interpPans[3] = wl2->getLower();
}


double wakePanel::panelPhi(const Eigen::Vector3d &POI)
{
    return -doubletStrength*dubPhiInf(POI);
}

Eigen::Vector3d wakePanel::panelV(const Eigen::Vector3d &POI)
{
    return doubletStrength*dubVInf(POI);
}

void wakePanel::setMu()
{
    std::vector<bodyPanel*> interpPans(4);
    double interpCoeff;
    interpPanels(interpPans, interpCoeff);
    doubletStrength = (1-interpCoeff)*interpPans[0]->getMu() + (interpCoeff-1)*interpPans[1]->getMu() + interpCoeff*interpPans[2]->getMu() - interpCoeff*interpPans[3]->getMu();
}

void wakePanel::setMu(double strength){
    doubletStrength = strength;
}

void wakePanel::setStrength()
{
    doubletStrength = upperPan->getMu()-lowerPan->getMu();
}

void wakePanel::setParentPanels(bodyPanel* upper, bodyPanel* lower)
{
    // Set flags used in gathering surrounding panels on same side of discontinuity in velocity calculation.
    upper->setUpper();
    lower->setLower();
    if (!TEpanel)
    {
        setTEpanel();
    }
    
    if (upperPan == nullptr)
    {
        // Parents have not yet been set
        upperPan = upper;
        lowerPan = lower;
    }
    else
    {
        // Parents already set, wake panel is at wing body joint. Choose panels further upstream
        
        if (upper->getCenter()(0) < center(0))
        {
            upperPan = upper;
            lowerPan = lower;
        }
    }
    
    
    
    wakeLine* wLine = new wakeLine(upperPan,lowerPan,normal);
    parentWake->addWakeLine(wLine);
    
}

edge* wakePanel::getTE()
{
    for (edges_index_type i=0; i<pEdges.size(); i++)
    {
        if (pEdges[i]->isTE())
        {
            return pEdges[i];
        }
    }
    return nullptr;
}

//wakePanel* wakePanel::makeVortexSheet()
//{
//    Eigen::Vector3d p1,p2,p3,p4,deltaVec;
//    Eigen::VectorXi sheetVerts = Eigen::VectorXi::Zero(4);
//    double length = 100;
//
//    p1(1) = -1000000;
//
//    // Find points on trailing edge;
//    Eigen::Vector3i neighbVerts = upperPan->getVerts();
//    bool breakFlag = false;
//    for (int i=0; i<verts.size(); i++)
//    {
//        for (int j=0; j<neighbVerts.size(); j++)
//        {
//            if (nodes->row(verts(i)) == nodes->row(neighbVerts(j)))
//            {
//                Eigen::Vector3d pnt = nodes->row(verts(i));
//                if (pnt(1) > p1(1))
//                {
//                    sheetVerts(1) = sheetVerts(0);
//                    p2 = p1;
//
//                    sheetVerts(0) = verts(i);
//                    p1 = pnt;
//
//                }
//                else
//                {
//                    sheetVerts(1) = verts(i);
//                    p2 = pnt;
//                    breakFlag = true;
//                    break;
//                }
//            }
//        }
//        if (breakFlag)
//        {
//            break;
//        }
//    }
//    sheetVerts(2) = (int)nodes->rows();
//    sheetVerts(3) = sheetVerts(2) + 1;
//
//    deltaVec = length*normal.cross((p2-p1).normalized());
//    p3 = p2+deltaVec;
//    p4 = p1+deltaVec;
//
//    nodes->conservativeResize(nodes->rows()+2, 3);
//    nodes->row(sheetVerts(2)) = p3;
//    nodes->row(sheetVerts(3)) = p4;
//
//    wakePanel* sheet = new wakePanel(sheetVerts,nodes,normal,ID,parentWake);
//    sheet->setTEpanel();
//    sheet->setUpper(upperPan);
//    sheet->setLower(lowerPan);
//    return sheet;
//
//
//}


double wakePanel::getPartRadius(Eigen::Vector3d Vinf, double &dt){
    // Particle radius will be the average distance between the spanwise particle seed points. Quackenbush, eq. 9
    
    Eigen::Vector3d currPnt = this->getCenter() + dt*Vinf;
    
    wakePanel* neighbor1 = this->getEdges()[1]->getOtherWakePan(this);
    wakePanel* neighbor2 = this->getEdges()[3]->getOtherWakePan(this);
    
    std::vector<double> dist;
    if(neighbor1) // If there is a neighbor, then use it
    {
        Eigen::Vector3d neighbor1Pnt = neighbor1->getCenter() + dt*Vinf;
        dist.push_back(std::abs((currPnt-neighbor1Pnt).norm()));
    }
    if(neighbor2)
    {
        Eigen::Vector3d neighbor2Pnt = neighbor2->getCenter() + dt*Vinf;
        dist.push_back(std::abs((currPnt-neighbor2Pnt).norm()));
    }
    
    return std::accumulate(dist.begin(), dist.end(), 0.0)/dist.size();
};


std::vector<int> wakePanel::sort_indexes(std::vector<double> &v) {
    
    // if function errors out, its because I replaced all 'std::size_t' with 'int' except for vector<int>
    
    // initialize original index locations
    std::vector<int> idx(v.size());
    for (size_t i = 0; i !=idx.size(); ++i) idx[i] = i;
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),[v](int i1, int i2) {return v[i1] < v[i2];});
    
    return idx;
}

Eigen::Vector3d wakePanel::partStretching(particle* part){
    Eigen::Vector3d POI = part->pos;
    Eigen::Vector3d partStretching;
    double dist2panel = (POI-center).norm();
    bool isFarField = false;
    Eigen::MatrixXd velGradMat = Eigen::Matrix3d::Zero();
    
    if(dist2panel/longSide > 5){
        isFarField = true;
    }
    
    if(isFarField)
    {
        velGradMat += velocityGradientPointDoublet(POI);
    }
    else
    {
        velGradMat += velocityGradientDoublet(POI);
    }
    
    partStretching = velGradMat*part->strength;
    
    return partStretching;
}


std::vector<cpNode*> wakePanel::pointsInOrder(){
    // will use edges in order to find the correct points
    
    //           0
    //       1------0       --> y
    //       |      |      |
    //      1|      |3      V
    //       |      |
    //       2------3      x
    //           2
    
    
    std::vector<edge*> pEdges = this->edgesInOrder();
    std::vector<cpNode*> ptsIO(4);
    
    // Look at points on the trailing edge first
    cpNode* n1 = pEdges[0]->getN1();
    cpNode* n2 = pEdges[0]->getN2();
    
    // If n1 shares a node with edge3, then it is node0
    if(pEdges[3]->containsNode(n1))
    {
        ptsIO[0] = n1;
        ptsIO[1] = n2;
    }
    else{
        ptsIO[0] = n2;
        ptsIO[1] = n1;
    }
    
    // Now look at other edge
    n1 = pEdges[2]->getN1();
    n2 = pEdges[2]->getN2();
    
    // If n1 shares a node with edge1 then it is node 2
    if(pEdges[1]->containsNode(n1))
    {
        ptsIO[2] = n1;
        ptsIO[3] = n2;
    }
    else
    {
        ptsIO[2] = n2;
        ptsIO[3] = n1;
    }
    
    
    return ptsIO;
}


std::vector<edge*> wakePanel::edgesInOrder(){
    // Can put in the constructor at a later time
    
    // Function will find the correct panel edges that look like this for a flat wake, but will still find the correct edges when the panel is flipped
    
    //       ----0---       --> y
    //       |      |      |
    //       1      3      V
    //       |      |      x
    //       ----2---
    //
    
    std::vector<edge *> edgesIO(4);
    std::vector<edge*> pEdges = this->getEdges();
    
    // The trailing edge is built first and the parallel downstream edge is always built third
    edgesIO[0] = pEdges[0];
    edgesIO[2] = pEdges[2];
    
    
    // Find the vectors in the drawing
    Eigen::Vector3d a = pEdges[0]->getMidPoint() - this->getCenter();
    Eigen::Vector3d b = this->getNormal();
    
    Eigen::Vector3d c = a.cross(b);
    Eigen::Vector3d d = center+c;
    
    // See if edge1 midpoint or edge 3 midpoint is closer to the point d;
    double distToE1 = (pEdges[1]->getMidPoint() - d).norm();
    double distToE3 = (pEdges[3]->getMidPoint() - d).norm();
    
    if(distToE3 < distToE1)
    {
        edgesIO[3] = pEdges[3];
        edgesIO[1] = pEdges[1];
    }
    else
    {
        edgesIO[3] = pEdges[1];
        edgesIO[1] = pEdges[3];
    }
    
    
    return edgesIO;
}


Eigen::Vector3d wakePanel::partSeedPt(Eigen::Vector3d &Vinf, double &dt){
    // Using the trailing edge and first wake panel to build off so it can be modified if second wake panel is implemented
    //  1---0
    //  |   |  spanwise orientation is not accounted for and might be opposite
    //  2---3  of the diagram but average still works
    
    Eigen::Vector3d n0,n1,n2,n3;
    edge* edge0 = this->getEdges()[0]; // Node proj. works from trailing edge only
    edge* edge2 = this->getEdges()[2]; // Downstream edge
    
    n0 = edge2->getN1()->getPnt();
    n1 = edge2->getN2()->getPnt();
    
    n2 = edge0->getN1()->secProjNode(dt, Vinf.norm());
    n3 = edge0->getN2()->secProjNode(dt, Vinf.norm());
    
    return (n0+n1+n2+n3)/4;
}

Eigen::Vector3d wakePanel::edgeStrength( edge* curEdge, int edgeNum){
    
    
    Eigen::Vector3d strength;
    std::vector<cpNode*> ptsIO = this->pointsInOrder();
    
    if(edgeNum == 0){
        std::cout << "Don't try to collapse the upstream edge" << std::endl;
        std::exit(0);
    }
    if(edgeNum == 2)
    {
        Eigen::Vector3d Rj = this->pointsInOrder()[2]->getPnt();
        Eigen::Vector3d Ri = this->pointsInOrder()[3]->getPnt();
        strength = (this->getMu()-this->getPrevStrength())*(Ri-Rj);
    }
    else if(edgeNum == 1)
    {
        wakePanel* otherPan = curEdge->getOtherWakePan(this);
        Eigen::Vector3d Rj = ptsIO[1]->getPnt();
        Eigen::Vector3d Ri = ptsIO[2]->getPnt();
        
        if(otherPan) // Panel has neighbor
        {
            strength = (this->getMu()-otherPan->getMu())*(Ri-Rj);
        }else{
            strength = this->getMu()*(Ri-Rj);
        }
    }
    else // Is edge 3
    {
        wakePanel* otherPan = curEdge->getOtherWakePan(this);
        Eigen::Vector3d Rj = ptsIO[3]->getPnt();
        Eigen::Vector3d Ri = ptsIO[0]->getPnt();
        
        if(otherPan)
        {
            strength = (this->getMu()-otherPan->getMu())*(Ri-Rj);
        }else{
            strength = this->getMu()*(Ri-Rj);
        }
    }
    return strength;
}
