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


wakePanel::wakePanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake,int surfID) : panel(nodes,pEdges,bezNorm,surfID), TEpanel(false), parentWake(parentWake), upperPan(nullptr), lowerPan(nullptr)
{
    for (int i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addWakePan(this);
    }
}

wakePanel::wakePanel(const wakePanel &copy) : panel(copy), upperPan(copy.upperPan), lowerPan(copy.lowerPan), TEpanel(copy.TEpanel), parentWake(copy.parentWake) {}

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
        for (int i=1; i<wakeLines.size()-1; i++)
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

void wakePanel::manuallySetMu(double strength) //VPP 
{
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
    for (int i=0; i<pEdges.size(); i++)
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

std::vector<cpNode*> wakePanel::pointsInOrder(){
    wakePanel* pan = this;
    
    std::vector<edge*> panEdges = pan->getEdges();
    std::vector<cpNode*> ptsIO;
    ptsIO.push_back(panEdges[0]->getN1());
    ptsIO.push_back(panEdges[0]->getN2());
    
    int pt = panEdges[0]->getN2()->getIndex();
    int edge1n1 = panEdges[0]->getN1()->getIndex();
    int edge1n2 = panEdges[0]->getN2()->getIndex();
    edge* nextEdge;
    
    for(int i=0; i<panEdges.size(); i++){
        int n1 = panEdges[i]->getN1()->getIndex();
        int n2 = panEdges[i]->getN2()->getIndex();
        
        bool samepanel =false;
        if((n1==edge1n1 && n2==edge1n2) || (n1==edge1n2 && n2==edge1n1)){
            samepanel = true;
        }
        if(n1==pt && samepanel==false){
            nextEdge = panEdges[i];
            pt = nextEdge->getN2()->getIndex();
            ptsIO.push_back(panEdges[i]->getN2());
            break;
        }
        if(n2==pt && samepanel == false){
            nextEdge = panEdges[i];
            pt = nextEdge->getN1()->getIndex();
            ptsIO.push_back(panEdges[i]->getN1());
            break;
        }
    }
    
    edge1n1 = nextEdge->getN1()->getIndex();
    edge1n2 = nextEdge->getN2()->getIndex();
    
    for(int i=0; i<panEdges.size(); i++){
        int n1 = panEdges[i]->getN1()->getIndex();
        int n2 = panEdges[i]->getN2()->getIndex();
        
        bool samepanel = false;
        if((n1==edge1n1 && n2==edge1n2) || (n1==edge1n2 && n2==edge1n1)){
            samepanel = true;
        }
        if(n1==pt && samepanel==false){
            ptsIO.push_back(panEdges[i]->getN2());
        }
        if(n2==pt && samepanel == false){
            ptsIO.push_back(panEdges[i]->getN1());
        }
    }
    
    return ptsIO;
}

//std::vector<Eigen::Vector3d> wakePanel::VPshedPts(){ //BW2
//    wakePanel* pan = this;
//    std::vector<Eigen::Vector3d> shedPts;
//    Eigen::Vector3d a,b,c,d,cen;
//    
//    a = pan->pointsInOrder()[0]->getPnt();
//    b = pan->pointsInOrder()[1]->getPnt();
//    c = pan->pointsInOrder()[2]->getPnt();
//    d = pan->pointsInOrder()[3]->getPnt();
//    cen = pan->getCenter();
//    
//    shedPts.push_back((a+((a+b)/2)+cen+((d+a)/2))/4);
//    shedPts.push_back((b+((b+c)/2)+cen+((a+b)/2))/4);
//    shedPts.push_back((c+((c+d)/2)+cen+((b+c)/2))/4);
//    shedPts.push_back((d+((d+a)/2)+cen+((c+d)/2))/4);
//    
//    return shedPts;
//}

Eigen::Vector3d wakePanel::partSeedPt(Eigen::Vector3d &Vinf, double &dt){ //vpp
    
//    Eigen::Vector3d projAngle = (this->getUpper()->getNormal() + this->getLower()->getNormal())/2;
    edge* edge2 = this->getEdges()[this->sortedEdgeInd()[1]];
    edge* edge4 = this->getEdges()[this->sortedEdgeInd()[3]]; // node proj. works from trailing edge only
    Eigen::Vector3d TEn1, TEn2, projn1, projn2;
    TEn1 = edge2->getN1()->getPnt();
    TEn2 = edge2->getN2()->getPnt();

    projn1 = edge4->getN1()->secProjNode(dt, Vinf.norm());
    projn2 = edge4->getN2()->secProjNode(dt, Vinf.norm());
    
    // average points for all of these. that should be it.
    
    
    Eigen::Vector3d seedPt = (TEn1+ TEn2+ projn1+ projn2)/4;
    return seedPt;
}

Eigen::Vector3d wakePanel::panToPartStrengthT1(){
    //                      --> y
    //       |      |      |
    //      3|      |1     V
    //       |______|      x
    //          2
    // For the first time step, there is no downstream vorticity so edge 2 needs to be accounted for. Edge 4 is never included because that edge is forced to have the same strength as the body edge in order to enforce the kutta condition. Collapse is illustrated on pg. 25
    
    Eigen::Vector3d edge1strength, edge2strength, edge3strength;
    std::vector<int> sortedEdgeInd = this->sortedEdgeInd();
    std::vector<Eigen::Vector3d> ringVecs = this->vortexRingVectors();
    
    
    // Edge 1
    wakePanel* neighbPan = this->getEdges()[sortedEdgeInd[0]]->getOtherWakePan(this);
    
    double neighbPanMu;
    if(neighbPan){
        neighbPanMu = neighbPan->getMu();
        edge1strength = (this->getMu() - neighbPanMu)/2*ringVecs[0]; // If neighbor panel, have to share strength with it.
    }else{ // If no wake pan, use whole strength
        edge1strength = (this->getMu())*ringVecs[0];
    }

    // Edge 2
    edge2strength = (this->getMu())*ringVecs[1];
    
    // Edge 3
    neighbPan = this->getEdges()[sortedEdgeInd[2]]->getOtherWakePan(this);
    
    if(neighbPan){
        neighbPanMu = neighbPan->getMu();
        edge3strength = (neighbPanMu - this->getMu())/2*ringVecs[2];
    }else{ // If no wake pan,
        edge3strength = (this->getMu())*ringVecs[2];
    }

    return -edge1strength + edge2strength +edge3strength;
}

Eigen::Vector3d wakePanel::panToPartStrength(){

    Eigen::Vector3d edge1strength, edge3strength;
    std::vector<int> sortedEdgeInd = this->sortedEdgeInd();
    std::vector<Eigen::Vector3d> ringVecs = this->vortexRingVectors();
    
    
    // Edge 1
    wakePanel* neighbPan = this->getEdges()[sortedEdgeInd[0]]->getOtherWakePan(this);
    
    double neighbPanMu;
    if(neighbPan){
        neighbPanMu = neighbPan->getMu();
        edge1strength = (this->getMu() - neighbPanMu)/2*ringVecs[0]; // If neighbor panel, have to share strength with it.
    }else{ // If no wake pan, use whole strength
        edge1strength = (this->getMu())*ringVecs[0];
    }
    
    // Edge 3
    neighbPan = this->getEdges()[sortedEdgeInd[2]]->getOtherWakePan(this);
    
    if(neighbPan){
        neighbPanMu = neighbPan->getMu();
        edge3strength = (neighbPanMu - this->getMu())/2*ringVecs[2];
    }else{ // If no wake pan,
        edge3strength = (this->getMu())*ringVecs[2];
    }
    
    return -edge1strength + edge3strength;
}


//std::vector<double> wakePanel::shedPanletArea(){
//    // Calculate the area of the panlets. Using an iregular polygon to calculate the area for more accuracy over trapezoidal approximation. More info here: keisan.casio.com/exec/system/1322718508
//    wakePanel* pan = this;
//    int numPanlets = 4;
//    std::vector<double> panletAreas;
//    
//    Eigen::Vector3d pt0,pt1,pt2,pt3;
//    pt0 = pan->getNodes()[0]->getPnt();
//    pt1 = pan->getNodes()[1]->getPnt();
//    pt2 = pan->getNodes()[2]->getPnt();
//    pt3 = pan->getNodes()[3]->getPnt();
//    
//    Eigen::Vector3d mid01,mid12,mid23,mid30; // Midpoint of panel edges
//    mid01 = (pt1+pt0)/2;
//    mid12 = (pt2+pt1)/2;
//    mid23 = (pt3+pt2)/2;
//    mid30 = (pt0+pt3)/2;
//    
//    std::vector<std::vector<Eigen::Vector3d>> panlets;
//    Eigen::Vector3d dummy;
//    panlets.resize(numPanlets,std::vector<Eigen::Vector3d>(pan->getNodes().size(),dummy));
//    panlets[0][0] = pt0; panlets[0][1] = mid01; panlets[0][2] = center; panlets[0][3] = mid30;
//    panlets[1][0] = pt1; panlets[1][1] = mid12; panlets[1][2] = center; panlets[1][3] = mid01;
//    panlets[2][0] = pt2; panlets[2][1] = mid23; panlets[2][2] = center; panlets[2][3] = mid12;
//    panlets[3][0] = pt3; panlets[3][1] = mid30; panlets[3][2] = center; panlets[3][3] = mid23;
//    
//    
//    for(int i=0; i<numPanlets; i++){
//        double aLen,bLen,cLen,dLen; // Vector Magnitudes
//        double th1,th2,theta,s,polyArea;
//        Eigen::Vector3d p0,p1,p2,p3;
//        Eigen::Vector3d avec,bvec,cvec,dvec;
//        
//        p0 = panlets[i][0]; p1 = panlets[i][1]; p2 = panlets[i][2]; p3 = panlets[i][3];
//        
//        avec = p1-p0;
//        bvec = p2-p1;
//        cvec = p3-p2;
//        dvec = p0-p3;
//        
//        aLen = avec.norm();
//        bLen = bvec.norm();
//        cLen = cvec.norm();
//        dLen = dvec.norm();
//        s = (aLen+bLen+cLen+dLen)/2;
//        
//        th1 = acos((avec.dot(bvec))/(aLen*bLen));
//        th2 = acos((cvec.dot(dvec))/(cLen*dLen));
//        theta = th1+th2;
//        
//        polyArea = pow(((s-aLen)*(s-bLen)*(s-cLen)*(s-dLen)-aLen*bLen*cLen*dLen*pow(cos(theta/2),2)),.5);
//        
//        panletAreas.push_back(polyArea);
//    }
//    
//    return panletAreas;
//}

std::vector<int> wakePanel::sortedEdgeInd(){
   
    // Function puts edges in order so that the proper edge can be identified when collapsing the panel.
    //        __3___        --> y
    //       |      |      |
    //      2|      |0     V
    //       |______|      x
    //          1
    
    
    
    std::vector<int> sortedEdgeInd(this->getEdges().size());
    
    std::vector<double> mpX(this->getEdges().size());
    std::vector<double> mpY(this->getEdges().size());
    
    for(int i=0; i<mpX.size();i++){
        mpX[i] = this->getEdges()[i]->getMidPoint()[0];
        mpY[i] = this->getEdges()[i]->getMidPoint()[1];
    }
    
    std::vector<int> mpXsorted = sort_indexes(mpX);
    std::vector<int> mpYsorted = sort_indexes(mpY);

    sortedEdgeInd[3] = mpXsorted[0]; // Comparing the edge midpoints to find max and min x,y for use in
    sortedEdgeInd[1] = mpXsorted[3];
    sortedEdgeInd[2] = mpYsorted[0];
    sortedEdgeInd[0] = mpYsorted[3];
   
    return sortedEdgeInd;
}

std::vector<Eigen::Vector3d> wakePanel::vortexRingVectors(){
    // Finds the unit direction vectors of the panels
    // CHANGING THIS to point in positive local x and y directions
    //        __4__>        --> y
    //       |      |      |
    //      3|      |1     V
    //       \/___>\/      x
    //          2
    
    wakePanel* pan = this;
    std::vector<Eigen::Vector3d> ringVecs;
    std::vector<int> edgeIndexes = pan->sortedEdgeInd();
    
    // Edge 1:
    Eigen::Vector3d n1 = pan->getEdges()[edgeIndexes[0]]->getN1()->getPnt();
    Eigen::Vector3d n2 = pan->getEdges()[edgeIndexes[0]]->getN2()->getPnt();
    if(n2.x() > n1.x()){
        ringVecs.push_back((n2-n1)/(n2-n1).norm()); // Return unit vector
    }else{
        ringVecs.push_back((n1-n2)/(n1-n2).norm());
    }
    // Edge 2:
    n1 = pan->getEdges()[edgeIndexes[1]]->getN1()->getPnt();
    n2 = pan->getEdges()[edgeIndexes[1]]->getN2()->getPnt();
    
    if(n2.y() > n1.y()){
        ringVecs.push_back((n2-n1)/(n2-n1).norm());
    }else{
        ringVecs.push_back((n1-n2)/(n1-n2).norm());
    }
    // Edge 3:
    n1 = pan->getEdges()[edgeIndexes[2]]->getN1()->getPnt();
    n2 = pan->getEdges()[edgeIndexes[2]]->getN2()->getPnt();
    if(n2.x() > n1.x()){ //just changed this.
        ringVecs.push_back((n2-n1)/(n2-n1).norm());
    }else{
        ringVecs.push_back((n1-n2)/(n1-n2).norm());
    }
    // Edge 4:
    n1 = pan->getEdges()[edgeIndexes[3]]->getN1()->getPnt();
    n2 = pan->getEdges()[edgeIndexes[3]]->getN2()->getPnt();
    
    if(n2.y() > n1.y()){
        ringVecs.push_back((n2-n1)/(n2-n1).norm());
    }else{
        ringVecs.push_back((n1-n2)/(n1-n2).norm());
    }
    
    return ringVecs;
}

//Eigen::Vector3d wakePanel::findPartStrength(){
//    Eigen::Vector3d strength;
//    
//    std::vector<Eigen::Vector3d> ringVecs = this->vortexRingVectors();
//    std::vector<int> sortedInd = this->sortedEdgeInd();
//    
//    
//    std::vector<double> edgeStrengths;
//    for(int i=0; i<ringVecs.size(); i++){
//        wakePanel* neighbPan = this->getEdges()[sortedEdgeInd()[i]]->getOtherWakePan(this);
//        double neighbPanMu;
//        
//        if(neighbPan){
//            neighbPanMu = neighbPan->getMu();
//        }else{
//            neighbPanMu = 0;
//        }
//        edgeStrengths.push_back(this->getMu() - neighbPanMu);
////        std::cout << "edge "+std::to_string(i)+" strength: " << edgeStrengths[i] << std::endl;
//    }
//
////    std::cout << "\n\n\n";
////    for(int i=0;i<4;i++){
////        std::cout << "plot3(" << this->getNodes()[i]->getPnt().x() << "," << this->getNodes()[i]->getPnt().y() << "," << this->getNodes()[i]->getPnt().z() << ",'g*');" << std::endl;
////    }
//    
//    // Vector Direction
//    Eigen::Vector3d vec13avg = (ringVecs[0]-ringVecs[2])/2; // Filament 3 pointing -x dir.
//    Eigen::Vector3d vec24avg = (ringVecs[3]-ringVecs[1])/2; // Filament 2 pointing -y dir.
//    
////    for(int i=0;i<this->getNodes().size();i++){
////        std::cout << "pt"+std::to_string(i)+" = [" << this->getNodes()[i]->getPnt().x() << "," << this->getNodes()[i]->getPnt().y() << ","<< this->getNodes()[i]->getPnt().z() << "];" << std::endl;
////    }
//    
////    std::cout << "streamWise = [" << vec13avg.x() << "," << vec13avg.y() << ","<< vec13avg.z() << "];" << std::endl;
////    std::cout << "spanWise = [" << vec24avg.x() << "," << vec24avg.y() << ","<< vec24avg.z() << "];" << std::endl;
////    std::cout << "center = [" << this->getCenter().x() << "," << this->getCenter().y() << ","<< this->getCenter().z() << "];" << std::endl;
//
//    strength = (edgeStrengths[2]-edgeStrengths[0])*vec24avg + (edgeStrengths[3]-edgeStrengths[1])*vec13avg;
//    
////    std::cout << "test using the plot panel function in matlab" << std::endl;
//    return strength;
//}

double wakePanel::getPartRadius(Eigen::Vector3d &Vinf, double &dt){
    //radius will be the average distance between the spanwise particle seed points a la Quackenbush et al. eq. (9)

    Eigen::Vector3d currPnt = this->partSeedPt(Vinf, dt);
    
    wakePanel* neighbor1 = this->getEdges()[sortedEdgeInd()[0]]->getOtherWakePan(this);
    wakePanel* neighbor2 = this->getEdges()[sortedEdgeInd()[2]]->getOtherWakePan(this);
    
    Eigen::Vector3d neighbor1Pnt = neighbor1->partSeedPt(Vinf, dt);
    Eigen::Vector3d neighbor2Pnt = neighbor2->partSeedPt(Vinf, dt);
    
    return 0.5*(std::abs((currPnt-neighbor1Pnt).norm()) + std::abs((currPnt-neighbor2Pnt).norm()));
    
};


std::vector<int> wakePanel::sort_indexes(std::vector<double> &v) {
    
    // initialize original index locations
    std::vector<int> idx(v.size());
    for (std::size_t i = 0; i !=idx.size(); ++i) idx[i] = i;
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),[v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});
    
    return idx;
}


