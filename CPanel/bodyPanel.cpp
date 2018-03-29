//
//  bodyPanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "bodyPanel.h"
#include "edge.h"
#include "cpNode.h"
#include "surface.h"

bodyPanel::bodyPanel(std::vector<cpNode*> nnodes, std::vector<edge*> ppEdges,
		            Eigen::Vector3d bezNorm,surface* pparentSurf, int surfID)
  : panel(nnodes,ppEdges,bezNorm,surfID), parentSurf(pparentSurf), TSorder(3), upper(false), lower(false),
	TEpanel(false), TE(nullptr), tipFlag(false), streamFlag(false), index(-1)
{
    for (size_t i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addBodyPan(this);
    }
    for (size_t i=0; i<nodes.size(); i++)
    {
        nodes[i]->addBodyPanel(this);
    }
}

void bodyPanel::addNeighbor(bodyPanel* p)
{
    neighbors.push_back(p);
}

void bodyPanel::setUpper() {upper = true;}
void bodyPanel::setLower() {lower = true;}

void bodyPanel::setTEpanel(edge* trailingEdge)
{
    TEpanel = true;
    parentSurf->setLSflag();
    TE = trailingEdge;
}

void bodyPanel::setIndex(int i) {index = i;}

void bodyPanel::setTipFlag()
{
    if (!parentSurf->isLiftingSurf())
    {
        tipFlag = false;
        return;
    }
    
    int count = 0;
    double angle;
    Eigen::Vector3d nNormal;
    for (bodyPanels_index_type i=0; i<neighbors.size(); i++)
    {
        nNormal = neighbors[i]->getNormal();
        double dot = normal.dot(nNormal)/(normal.norm()*nNormal.norm());
        if (dot > 1)
        {
            dot = 1;
        }
        else if (dot < -1)
        {
            dot = -1;
        }
        
        angle = acos(dot);
        
        if (angle < pow(10,-8) && asin(nNormal(0)) > -M_PI/12)
        {
            // asin(nNormal(0)) is the angle normal vector makes with xy plane.  Catches bug where some leading edge panels were being flagged as on the tip patch.
            count++;
        }
    }
    if (count >= 2)
    {
        tipFlag = true;
    }
    
}

void bodyPanel::setSigma(Eigen::Vector3d Vinf, double Vnorm)
{
    sourceStrength = (-Vinf.dot(normal)+Vnorm); // CS: Vnorm is eminating from the panel
}

void bodyPanel::setMu(double dubStrength)
{
    doubletStrength = dubStrength;
}

void bodyPanel::setStreamFlag()
{
    streamFlag = true;
}

double bodyPanel::panelPhi(const Eigen::Vector3d &POI)
{
    double phiSrc,phiDub;
    phiSrc = 0;
    phiDub = 0;
    
    panelPhiInf(POI,phiSrc,phiDub);
    return -sourceStrength*phiSrc-doubletStrength*phiDub;
}

Eigen::Vector3d bodyPanel::panelV(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d vSrc = Eigen::Vector3d::Zero();
    Eigen::Vector3d vDub = Eigen::Vector3d::Zero();
    
    panelVInf(POI,vSrc,vDub);
    
    return sourceStrength*vSrc+doubletStrength*vDub;
}

void bodyPanel::panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub)
{
    Eigen::Vector3d pjk = POI-center;
    bool itselfFlag = false;
    if (pjk.norm() < pow(10,-10))
    {
        itselfFlag = true;
    }
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
    {
        phiSrc = pntSrcPhi(pjk.norm());
        phiDub = pntDubPhi(PN,pjk.norm());
    }
    else
    {
        double Al,phiV;
        Eigen::Vector3d a,b,s;
        for (nodes_index_type i=0; i<nodes.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            Al = local.row(2).dot(s.cross(a));
            if (!itselfFlag)
            {
// Note: last parameter was not used
//                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1));
                phiDub += phiV;
            }
            phiSrc += srcSidePhi(PN,Al,phiV,a,b,s);
            
        }
        phiSrc /= (4*M_PI);
        if (!itselfFlag)
        {
            phiDub /= (4*M_PI);
        }
        else
        {
            phiDub = -0.5;
        }
    }
}

void bodyPanel::panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub)
{
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
    {
        vSrc = -pntSrcV(pjk); // Connor added this negative. Figure 10.17 (Katz) shows a positive velocity influence from pnt source. CPanel matches this in tests, however, sign convention for CPanel/VSaero is opposite of Katz
        vDub = pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s,l,m,n;
        l = local.row(0);
        m = local.row(1);
        n = local.row(2);
        pjk = POI-center;
        double Al;
        for (nodes_index_type i=0; i<nodes.size(); i++)
        {
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = a-b;
            Al = n.dot(s.cross(a));
            
            vDub += vortexV(a,b,s);
            vSrc += srcSideV(PN,Al,a,b,s,l,m,n);
        }
        vDub /= (4*M_PI);
        vSrc /= (4*M_PI);
    }
}

Eigen::Vector3d bodyPanel::pntVInf(const Eigen::Vector3d &POI){
    // Function should only be used for CPanel test function.
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    
    return sourceStrength*pntSrcV(pjk) + doubletStrength*pntDubV(local.row(2),pjk); //pntDubVinf
}


double bodyPanel::srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
    	GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
    return (Al*GL-PN*phiV);
}

Eigen::Vector3d bodyPanel::srcSideV(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
        GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
// NOTE: last parameter is not used
//    double CJK = vortexPhi(PN,Al,a,b,s,l,m,n);
    double CJK = vortexPhi(PN,Al,a,b,s,l,m);
    return (GL*(s.dot(m)*l-s.dot(l)*m)+CJK*n);
}

inline double bodyPanel::pntSrcPhi(const double &PJK)
{
    return area/(4*M_PI*PJK);
}

inline Eigen::Vector3d bodyPanel::pntSrcV(const Eigen::Vector3d &pjk)
{
    return area*pjk/(4*M_PI*pow(pjk.norm(),3));
}

void bodyPanel::setCluster()
{
    int dim;
    int buffer = 10;
    if (tipFlag)
    {
        dim = 2;
    }
    else
    {
        dim = 3;
    }
    
    double nObs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))-1+buffer; // Binomial Coefficient
    double nPanels = std::ceil(nObs/2);
    bodyPanels_index_type oldSize = cluster.size();
    cluster.push_back(this);
    bool upFlag = upper;
    bool lowFlag = lower;
    while (cluster.size() < nPanels)
    {
        std::vector<bodyPanel*> toAdd;
        for (bodyPanels_index_type i=oldSize; i<cluster.size(); i++)
        {
            std::vector<bodyPanel*> temp = cluster[i]->getNeighbors();
            for (bodyPanels_index_type j=0; j<temp.size(); j++)
            {
                if (clusterTest(temp[j], 5*M_PI/6,upFlag,lowFlag))
                {
                    // Do not include panels on other side of discontinuity (i.e. wake), panels already in cluster, or this panel
                    if (temp[j]->isUpper())
                    {
                        upFlag = true;
                    }
                    else if (temp[j]->isLower())
                    {
                        lowFlag = true;
                    }
                    if (std::find(toAdd.begin(), toAdd.end(), temp[j]) == toAdd.end())
                    {
                        toAdd.push_back(temp[j]);
                    }
                }
            }
        }
        
        if (toAdd.size() == 0 && cluster.size() < nPanels)
        {
            // No valid panels could be found so decrease order of taylor series
            while (cluster.size() < nPanels+1 && TSorder >= 1)
            {
                if (buffer > 0)
                {
                    buffer -= 1;
                }
                else
                {
                    TSorder -= 1;
                }
                
                // Do Something to fix error when cluster can't be found.
                
                nObs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder)) + buffer; // Binomial Coefficient
                nPanels = ceil(nObs/2);
                
            }
            
            if (cluster.size() > 0) {
                cluster.erase(cluster.begin()+nPanels+1,cluster.end());
            }
            break;
        }
        else
        {
            oldSize = cluster.size();
            for (bodyPanels_index_type i=0; i<toAdd.size(); i++)
            {
                cluster.push_back(toAdd[i]);
                if (cluster.size() == nPanels+1)
                {
                    break;
                }
            }
        }
    }
}


Eigen::Vector3d bodyPanel::pntVelocity(const Eigen::Vector3d &pnt,double pntPotential, double PG, const Eigen::Vector3d &Vinf)
{
    Eigen::Vector3d vel;
    
    if (cluster.size() == 0)
    {
        setCluster();
    }
    
    if (tipFlag)
    {
        vel = velocity2D(pnt,pntPotential);
    }
    else
    {
        vel = velocity3D(pnt,pntPotential);
        if (vel != vel)
        {
            // NaN can be returned if supporting data was all on a flat patch (essentially 2D)
            vel = velocity2D(pnt,pntPotential);
        }
    }
    
    if (vel != vel)
    {
        // Last resort approximation if CHTLS not able to compute velocity. From Kinney (CBAERO).
        vel = (bezNormal.cross(Vinf)).cross(bezNormal);
    }
    
    vel(0) /= PG; // Prandlt Glauert Correction
    return vel;
}


Eigen::Vector3d bodyPanel::velocity2D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 2;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Eigen::Vector3d pntLocal = global2local(pnt,true);
    Eigen::Vector2d X0 = pntLocal.head(2);
    Eigen::Vector2d V0 = Eigen::Vector2d::Zero();
    Xf.resize(clust.size(),3);
    Xb = Eigen::MatrixXd::Zero(0,dim);
    Vb = Eigen::MatrixXd::Zero(0,dim);
    df.resize(clust.size());
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(i) = global2local(clust[i]->getCenter(),true);
        df(i) = clust[i]->getPotential()-pntPotential;
    }
    Eigen::MatrixXd xLocal = Xf.block(0,0,clust.size(),2);
    chtlsnd tipV(X0,xLocal,TSorder,Xb,Vb,V0);
    Eigen::Vector3d vLocal;
    vLocal(0) = tipV.getF().row(0)*df;
    vLocal(1) = tipV.getF().row(1)*df;
    vLocal(2) = 0;
    vel = local2global(vLocal,false);
    
    return vel;
}

Eigen::Vector3d bodyPanel::velocity3D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 3;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Xf.resize(clust.size(),dim);
    Xb.resize(clust.size(),dim);
    Vb.resize(clust.size(),dim);
    df.resize(clust.size());
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(i) = clust[i]->getCenter();
        Vb.row(i) = clust[i]->getBezNormal();
        df(i) = clust[i]->getPotential()-pntPotential;
    }
    Xb = Xf;
    chtlsnd vWeights(pnt,Xf,TSorder,Xb,Vb,Eigen::Vector3d::Zero());
    
    vel(0) = vWeights.getF().row(0)*df;
    vel(1) = vWeights.getF().row(1)*df;
    vel(2) = vWeights.getF().row(2)*df;
    
    return vel;
}

void bodyPanel::computeCp(double Vinf)
{
    Cp = (1-pow(velocity.norm()/Vinf,2));
}

void bodyPanel::computeCp(double Vinf,double dt){
    double dPhi_dt = ( prevPotential - potential ) / dt;
    
    Cp = 1 - pow( velocity.norm()/Vinf , 2) - 2/(Vinf*Vinf) * dPhi_dt; // Katz 13.168
}

void bodyPanel::computeVelocity(double PG, const Eigen::Vector3d &Vinf)
{
    velocity = pntVelocity(center,potential,PG,Vinf);
}



Eigen::Vector3d bodyPanel::computeMoments(const Eigen::Vector3d &cg)
{
    Eigen::Vector3d r = center-cg;
    Eigen::Vector3d F = -Cp*bezNormal*area;
    return r.cross(F);
}

bool bodyPanel::clusterTest(bodyPanel* other,double angle,bool upFlag,bool lowFlag)
{
    if ((upFlag && other->isLower()) || (lowFlag && other->isUpper()))
    {
        return false;
    }
    
    if (tipFlag != other->isTipPan())
    {
        return false;
    }
    
    double dot = other->getNormal().dot(normal);
    // Floating point error can cause panels with the same normal vector to result in a dot product greater than one, causing acos to return nan.
    if (dot > 1)
    {
        dot = 1;
    }
    else if (dot < -1)
    {
        dot = -1;
    }
    return (acos(dot) < angle && std::find(cluster.begin(),cluster.end(),other)==cluster.end() && other != this);
}

void bodyPanel::setLSflag()
{
    parentSurf->setLSflag();
}

bool bodyPanel::nearTrailingEdge()
{
    // Returns true if panel is on trailing edge or a neighbor is on trailing edge
    if (parentSurf->isLiftingSurf())
    {
        if (upper || lower)
        {
            return true;
        }
        else
        {
            std::vector<bodyPanel*> neighbs = getNeighbors();
            for (bodyPanels_index_type i=0; i<neighbs.size(); i++)
            {
                if (neighbs[i]->isUpper() || neighbs[i]->isLower())
                {
                    return true;
                }
            }
        }
    }
    return false;
}

double bodyPanel::dist2Pan(bodyPanel* other)
{
    return (other->getCenter()-center).norm();
}

std::vector<bodyPanel*> bodyPanel::getRelatedPanels()
{
    std::vector<bodyPanel*> pans;
    pans.push_back(this);
    std::vector<bodyPanel*> nodePans;
    
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        nodePans = nodes[i]->getBodyPans();
        for (size_t j=0; j<nodePans.size(); j++)
        {
            if (std::find(pans.begin(),pans.end(),nodePans[j]) == pans.end())
            {
                pans.push_back(nodePans[j]);
            }
            
        }
    }
    return pans;
}

Eigen::Vector3d bodyPanel::partStretching(particle* part){
    
    Eigen::Vector3d partStretching;
    double dist2panel = (part->pos-center).norm();
    Eigen::MatrixXd velGradMat = Eigen::Matrix3d::Zero();
    
    // Far field convention usually either 3 or 5, see Katz or Chris' thesis
    bool isFarField = false;
    if(dist2panel/longSide > 5){
        isFarField = true;
    }
    
    if(isFarField)
    {
        velGradMat += velocityGradientPointSource(part->pos);
        velGradMat += velocityGradientPointDoublet(part->pos);
    }
    else if (this->getNodes().size() == 3)
    {
        velGradMat += velocityGradientTriSource(part->pos);
        velGradMat += velocityGradientDoublet(part->pos);
    }
    else
    {
        velGradMat += velocityGradientQuadSource(part->pos);
        velGradMat += velocityGradientDoublet(part->pos);
    }
    
    partStretching = velGradMat*part->strength;
    
    return partStretching;
}



Eigen::Matrix3d bodyPanel::velocityGradientPointSource(Eigen::Vector3d POI){
    // Equations in documentation
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    double x, y, z, x0, y0, z0;
    
    x = this->getCenter().x(); y = this->getCenter().y(); z = this->getCenter().z();
    
    x0 = POI.x(); y0 = POI.y(); z0 = POI.z();
    
    double sdd = pow(pow(x-x0,2) + pow(y-y0,2) + z*z,2.5); // source deriv denom
    
    double derConst = sourceStrength*area/(4*M_PI);
    velGradMat(0,0) = (-2*(x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z)/sdd;
    velGradMat(1,0) = -3*(x-x0)*(y-y0)/sdd;
    velGradMat(2,0) = -3*(x-x0)*z/sdd;
    
    velGradMat(0,1) = -3*(x-x0)*(y-y0)/sdd;
    velGradMat(1,1) = ((x-x0)*(x-x0)- (y-y0)*(y-y0) + z*z)/sdd;
    velGradMat(2,1) = -3*z*(y-y0)/sdd;
    
    velGradMat(0,2) = -3*(x-x0)*(z-z0)/sdd;
    velGradMat(1,2) = -3*(y-y0)*(z-z0)/sdd;
    velGradMat(2,2) = (-3*z*(z-z0) + ((x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z))/sdd;
    
    velGradMat *= derConst;
    
    return velGradMat;
}


Eigen::Matrix3d bodyPanel::velocityGradientQuadSource(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
    Eigen::Vector3d n4global = this->getNodes()[3]->getPnt();
    
    Eigen::Vector3d n1 = global2local(n1global, true);
    Eigen::Vector3d n2 = global2local(n2global, true);
    Eigen::Vector3d n3 = global2local(n3global, true);
    Eigen::Vector3d n4 = global2local(n4global, true);
    Eigen::Vector3d POIloc = global2local(POI, true);
    
    double x1, y1, x2, y2, x3, y3, x4, y4;
    x1 = n1.x(); x2 = n2.x(); x3 = n3.x(), x4 = n4.x();
    y1 = n1.y(); y2 = n2.y(); y3 = n3.y(), y4 = n4.y();
    
    double x, y, z;
    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
    
    
    double derConst = sourceStrength/(4*M_PI);
    
    // Terms
    double d12,d23,d34,d41;
    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
    d34 = pow((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3),0.5);
    d41 = pow((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4),0.5);
    
    double r1,r2,r3,r4;
    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    r4 = pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dx, dr2dx, dr3dx, dr4dx;
    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dx = (x-x4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dy, dr2dy, dr3dy, dr4dy;
    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dy = (y-y4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dz, dr2dz, dr3dz, dr4dz;
    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dz = (z)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double m12,m23,m34,m41;
    m12 = (y2-y1)/(x2-x1);
    m23 = (y3-y2)/(x3-x2);
    m34 = (y4-y3)/(x4-x3);
    m41 = (y1-y4)/(x1-x4);
    
    double e1,e2,e3,e4;
    e1 = (x-x1)*(x-x1) + z*z;
    e2 = (x-x2)*(x-x2) + z*z;
    e3 = (x-x3)*(x-x3) + z*z;
    e4 = (x-x4)*(x-x4) + z*z;
    
    double h1,h2,h3,h4;
    h1 = (x-x1)*(y-y1);
    h2 = (x-x2)*(y-y2);
    h3 = (x-x3)*(y-y3);
    h4 = (x-x4)*(y-y4);
    
    double de1dx,de2dx,de3dx,de4dx;
    de1dx = 2*(x-x1);
    de2dx = 2*(x-x2);
    de3dx = 2*(x-x3);
    de4dx = 2*(x-x4);
    
    double dh1dx,dh2dx,dh3dx,dh4dx;
    dh1dx = (y-y1);
    dh2dx = (y-y2);
    dh3dx = (y-y3);
    dh4dx = (y-y4);
    
    double dh1dy,dh2dy,dh3dy,dh4dy;
    dh1dy = (x-x1);
    dh2dy = (x-x2);
    dh3dy = (x-x3);
    dh4dy = (x-x4);
    
    double de1dz,de2dz,de3dz,de4dz;
    de1dz = 2*z;
    de2dz = 2*z;
    de3dz = 2*z;
    de4dz = 2*z;
    
    
    velGradMat(0,0)=(y2-y1)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dx + dr4dx)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dx + dr1dx)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(1,0)=(y2-y1)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dy + dr4dy)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dy + dr1dy)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(2,0)=(y2-y1)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dz + dr4dz)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dz + dr1dz)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(0,1)=(x1-x2)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dx + dr4dx)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dx + dr1dx)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(1,1)=(x1-x2)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dy + dr4dy)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dy + dr1dy)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(2,1)=(x1-x2)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dz + dr4dz)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dz + dr1dz)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(0,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*(m12*de1dx - dh1dx) - (m12*e1-h1)*(z*dr1dx))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*(m12*de2dx - dh2dx) - (m12*e2-h2)*(z*dr2dx))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*(m23*de2dx - dh2dx) - (m23*e2-h2)*(z*dr2dx))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*(m23*de3dx - dh3dx) - (m23*e3-h3)*(z*dr3dx))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((z*r3*(m34*de3dx - dh3dx) - (m34*e3-h3)*(z*dr3dx))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((z*r4*(m34*de4dx - dh4dx) - (m34*e4-h4)*(z*dr4dx))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((z*r4*(m41*de4dx - dh4dx) - (m41*e4-h4)*(z*dr4dx))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((z*r1*(m41*de1dx - dh1dx) - (m41*e1-h1)*(z*dr1dx))/pow(z*r1,2));
    
    velGradMat(1,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m12*e1-h1)*(z*dr1dy))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m12*e2-h2)*(z*dr2dy))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m23*e2-h2)*(z*dr2dy))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m23*e3-h3)*(z*dr3dy))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m34*e3-h3)*(z*dr3dy))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((-z*r4*dh4dy - (m34*e4-h4)*(z*dr4dy))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((-z*r4*dh4dy - (m41*e4-h4)*(z*dr4dy))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m41*e1-h1)*(z*dr1dy))/pow(z*r1,2));
    
    velGradMat(2,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*m12*de1dz - (m12*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*m12*de2dz - (m12*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*m23*de2dz - (m23*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*m23*de3dz - (m23*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((z*r3*m34*de3dz - (m34*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((z*r4*m34*de4dz - (m34*e4-h4)*(r4+z*dr4dz))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((z*r4*m41*de4dz - (m41*e4-h4)*(r4+z*dr4dz))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((z*r1*m41*de1dz - (m41*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2));
    
    velGradMat *= derConst;
    
    return velGradMat;
}

Eigen::Matrix3d bodyPanel::velocityGradientTriSource(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
    
    Eigen::Vector3d n1 = global2local(n1global, true);
    Eigen::Vector3d n2 = global2local(n2global, true);
    Eigen::Vector3d n3 = global2local(n3global, true);
    Eigen::Vector3d POIloc = global2local(POI, true);
    
    double x1, y1, x2, y2, x3, y3;
    x1 = n1.x(); x2 = n2.x(); x3 = n3.x();
    y1 = n1.y(); y2 = n2.y(); y3 = n3.y();
    
    double x, y, z;
    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
    
    double derConst = sourceStrength/(4*M_PI);
    
    // Terms
    double d12,d23,d31;
    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
    d31 = pow((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3),0.5);
    
    double r1,r2,r3;
    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dx, dr2dx, dr3dx;
    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dy, dr2dy, dr3dy;
    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dz, dr2dz, dr3dz;
    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double m12,m23,m31;
    m12 = (y2-y1)/(x2-x1);
    m23 = (y3-y2)/(x3-x2);
    m31 = (y1-y3)/(x1-x3);
    
    double e1,e2,e3;
    e1 = (x-x1)*(x-x1) + z*z;
    e2 = (x-x2)*(x-x2) + z*z;
    e3 = (x-x3)*(x-x3) + z*z;
    
    double h1,h2,h3;
    h1 = (x-x1)*(y-y1);
    h2 = (x-x2)*(y-y2);
    h3 = (x-x3)*(y-y3);
    
    double de1dx,de2dx,de3dx;
    de1dx = 2*(x-x1);
    de2dx = 2*(x-x2);
    de3dx = 2*(x-x3);
    
    double dh1dx,dh2dx,dh3dx;
    dh1dx = (y-y1);
    dh2dx = (y-y2);
    dh3dx = (y-y3);
    
    double dh1dy,dh2dy,dh3dy;
    dh1dy = (x-x1);
    dh2dy = (x-x2);
    dh3dy = (x-x3);
    
    double de1dz,de2dz,de3dz;
    de1dz = 2*z;
    de2dz = 2*z;
    de3dz = 2*z;
    
    
    velGradMat(0,0)=(y2-y1)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dx + dr1dx)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(1,0)=(y2-y1)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dy + dr1dy)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(2,0)=(y2-y1)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dz + dr1dz)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(0,1)=(x1-x2)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dx + dr1dx)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(1,1)=(x1-x2)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dy + dr1dy)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(2,1)=(x1-x2)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dz + dr1dz)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(0,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*(m12*de1dx - dh1dx) - (m12*e1-h1)*(z*dr1dx))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*(m12*de2dx - dh2dx) - (m12*e2-h2)*(z*dr2dx))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*(m23*de2dx - dh2dx) - (m23*e2-h2)*(z*dr2dx))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*(m23*de3dx - dh3dx) - (m23*e3-h3)*(z*dr3dx))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((z*r3*(m31*de3dx - dh3dx) - (m31*e3-h3)*(z*dr3dx))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((z*r1*(m31*de1dx - dh1dx) - (m31*e1-h1)*(z*dr1dx))/pow(z*r1,2));
    
    velGradMat(1,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m12*e1-h1)*(z*dr1dy))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m12*e2-h2)*(z*dr2dy))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m23*e2-h2)*(z*dr2dy))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m23*e3-h3)*(z*dr3dy))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m31*e3-h3)*(z*dr3dy))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m31*e1-h1)*(z*dr1dy))/pow(z*r1,2));
    
    velGradMat(2,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*m12*de1dz - (m12*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*m12*de2dz - (m12*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*m23*de2dz - (m23*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*m23*de3dz - (m23*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((z*r3*m31*de3dz - (m31*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((z*r1*m31*de1dz - (m31*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2));
    
    velGradMat *= derConst;
    
    return velGradMat;
}




