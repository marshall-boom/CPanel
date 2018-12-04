/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
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
 *    Chris Satterwhite - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include "bodyPanel.h"
#include "edge.h"
#include "cpNode.h"
#include "surface.h"

#define _USE_MATH_DEFINES
#include <math.h>

bodyPanel::bodyPanel(std::vector<cpNode*> nnodes, std::vector<edge*> ppEdges,
		            Eigen::Vector3d bezNorm,surface* pparentSurf, size_t surfID)
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
    //sourceStrength = (-Vinf.dot(normal)+Vnorm); // CS: Vnorm is eminating from the panel
	sourceStrength = (Vinf.dot(normal) + Vnorm);
	//sourceStrength = (-Vinf.dot(normal) + Vnorm) / supAreaCorrect;
	//sourceStrength = (-Vinf.dot(normal) + Vnorm) * supAreaCorrect;
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
		double Al;
		double phiV = 0;
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

void bodyPanel::srcPanelPhiInf(const Eigen::Vector3d &POI, double &phi)
{
	Eigen::Vector3d pjk = POI - center;
	bool itselfFlag = false;
	if (pjk.norm() < pow(10, -10))
	{
		itselfFlag = true;
	}
	Eigen::Matrix3d local = getLocalSys();
	double PN = pjk.dot(local.row(2));
	if (pjk.norm() / longSide > 5)
	{
		phi = pntSrcPhi(pjk.norm());
	}
	else
	{
		double Al;
		double phiV = 0;
		Eigen::Vector3d a, b, s;
		for (nodes_index_type i = 0; i<nodes.size(); i++)
		{
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			if (i != nodes.size() - 1)
			{
				p1 = nodes[i]->getPnt();
				p2 = nodes[i + 1]->getPnt();
			}
			else
			{
				p1 = nodes[i]->getPnt();
				p2 = nodes[0]->getPnt();
			}
			a = POI - p1;
			b = POI - p2;
			s = p2 - p1;
			Al = local.row(2).dot(s.cross(a));
			if (!itselfFlag)
			{
				// Note: last parameter was not used
				//                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
				phiV = vortexPhi(PN, Al, a, b, s, local.row(0), local.row(1));
			}
			phi += srcSidePhi(PN, Al, phiV, a, b, s);
		}
		phi /= (4 * M_PI);
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
    size_t buffer = 10;
    if (tipFlag)
    {
        dim = 2;
    }
    else
    {
        dim = 3;
    }
    
    size_t nObs = static_cast<size_t>(chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))-1)+buffer; // Binomial Coefficient
    size_t nPanels = (nObs+1)/2; // ceil(nObs/2) for integers
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
                
                nObs = static_cast<size_t>(chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))) + buffer; // Binomial Coefficient
                nPanels = (nObs+1)/2; // ceil(nObs/2) for integers
                
            }
            
            if (cluster.size() > 0) {
                cluster.erase(cluster.begin()+static_cast<bodyPanels_type::difference_type>(nPanels+1),cluster.end());
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
    Xf.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),3);
    Xb = Eigen::MatrixXd::Zero(0,dim);
    Vb = Eigen::MatrixXd::Zero(0,dim);
    df.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()));
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(static_cast<Eigen::MatrixXd::Index>(i)) = global2local(clust[i]->getCenter(),true);
        df(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getPotential()-pntPotential;
    }
    Eigen::MatrixXd xLocal = Xf.block(0,0,static_cast<Eigen::MatrixXd::Index>(clust.size()),2);
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
    
    Xf.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    Xb.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    Vb.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    df.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()));
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getCenter();
        Vb.row(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getBezNormal();
        df(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getPotential()-pntPotential;
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
	//std::cout << velocity << "\n" << std::endl;
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


void bodyPanel::linComputeVelocity(double PG, Eigen::Vector3d &Vinf)
{
	double mu_x, mu_y;
	Eigen::Vector3d vertDubStrengths, linDubConsts, panVel;
	Eigen::Matrix3d vertsMat = linVertsMatrix(true);
	vertDubStrengths = linGetDubStrengths();

	linDubConsts = vertsMat.inverse() * vertDubStrengths;
	//mu_0 = linDubConsts[0];
	mu_x = linDubConsts[1];
	mu_y = linDubConsts[2];

	panVel[0] = mu_x;
	panVel[1] = mu_y;
	panVel[2] = 0;

	velocity = local2global(panVel, false);
	velocity(0) /= PG;
}


Eigen::Vector3d bodyPanel::supComputeVelocity(Eigen::Vector3d Vinf, const double mach, bool velCorrection)
{
	double mu_x, mu_y;
	Eigen::Vector3d vertDubStrengths, linDubConsts, panPertVel, pertVelOrig, pertVel;
	Eigen::Matrix3d vertsMat = supVertsMatrix(supLocalNodes);
	vertDubStrengths = linGetDubStrengths();

	linDubConsts = vertsMat.inverse() * vertDubStrengths;
	//mu_0 = linDubConsts[0];
	mu_x = -linDubConsts[1];
	mu_y = -linDubConsts[2];

	panPertVel[0] = mu_x;
	panPertVel[1] = mu_y;
	panPertVel[2] = 0;

	//panPertVel += 2*abs(sourceStrength) * -normal;

	pertVelOrig = (supTransMat.transpose()).inverse() * panPertVel;
	//pertVelOrig.x() -= abs(sourceStrength)/2;
	//pertVelOrig -= 2*(sourceStrength * normal);

	if (velCorrection)
	{
		pertVel = supVelCorrection(pertVelOrig, mach);
	}
	else
	{
		pertVel = pertVelOrig;
	}

	velocity = Vinf + pertVel;

	std::cout << velocity.x() << "\t" << velocity.y() << std::endl;

	return pertVel;
}


Eigen::Vector3d bodyPanel::supVelCorrection(Eigen::Vector3d pertVel, const double mach)
{
	Eigen::Vector3d pertVelCorrected, pertMassFlux, totMassFlux;
	double Bmach2, normLocDensity, lessCheck, moreCheck;
	Bmach2 = pow(mach, 2) - 1; // Bmach2 = 1 w/ mach = sqrt(2)

	pertMassFlux = pertVel;
	pertMassFlux.x() *= Bmach2;
	totMassFlux = velocity + pertMassFlux;

	if (pertVel.x() < 0)
	{
		//velocity = totMassFlux / (1 - pow(mach, 2)*(pertVel.x()/Vinf.norm()));
		pertVelCorrected = pertMassFlux / (1 - pow(mach, 2)*pertVel.x());
	}

	return pertVelCorrected;
}


void bodyPanel::supComputeCp(Eigen::Vector3d Vinf, const double mach, Eigen::Vector3d pertVel)
{
	double beta2 = 1 - pow(mach,2);

	Cp = (-2 * pertVel.x() / Vinf.norm()) - ((pow(pertVel.y(), 2) + pow(pertVel.z(), 2)) / pow(Vinf.norm(), 2));
	//Cp = -2 * pertVel.x() / Vinf.norm();
}


//Eigen::Vector3d bodyPanel::linComputeVelocity2(double PG, Eigen::Vector3d &Vinf, Eigen::Vector3d &POI)
//{
//	//double mu_x, mu_y;
//	Eigen::Vector3d vertDubStrengths, linDubConsts, panVel, dubVec, POIloc;
//	Eigen::Matrix3d vertsMat, Jints;
//	double mu;
//
//	POIloc = global2local(POI, true);
//
//	vertsMat = linVertsMatrix(true);
//	vertDubStrengths = linGetDubStrengths();
//
//	linDubConsts = vertsMat.inverse() * vertDubStrengths;
//	mu = linDubConsts[0] + linDubConsts[1] * POIloc.x() + linDubConsts[2] * POIloc.y();
//	dubVec[0] = mu;
//	dubVec[1] = linDubConsts[1];
//	dubVec[2] = linDubConsts[2];
//
//	//Jints = linDubVInf(POI) / (4.0*M_PI);
//	Jints = linDubVInf(POI);
//
//	panVel = Jints * dubVec;
//
//	/*mu_0 = linDubConsts[0];
//	mu_x = linDubConsts[1];
//	mu_y = linDubConsts[2];*/
//
//	//velocity += local2global(panVel, false);
//	//velocity += panVel;
//	//velocity(0) /= PG;
//
//	return panVel;
//}


double bodyPanel::linGetTEdubStrength()
{
	edge* trailEdge;
	double mu_0, mu_x, mu_y, mu;
	Eigen::Vector3d vertDubStrengths, linDubConsts, midPntTE;

	Eigen::Matrix3d vertsMat = linVertsMatrix(true);
	vertDubStrengths = linGetOrigDubStrengths();

	linDubConsts = vertsMat.inverse() * vertDubStrengths;
	mu_0 = linDubConsts[0];
	mu_x = linDubConsts[1];
	mu_y = linDubConsts[2];

	trailEdge = getTrailingEdge();
	midPntTE = global2local(trailEdge->getMidPoint(), true);

	mu = mu_0 + mu_x * midPntTE.x() + mu_y * midPntTE.y();

	return mu_0;
}


// Don't actually need edgeFlags. Just do DOIflag, then check R1 and R2 when doing inf coeff calcs

bool bodyPanel::supDOIcheck(Eigen::Vector3d &P, const double Mach, Eigen::Vector3d &windDir)
{
	// P is POI, Q is panel center

	//Eigen::Vector3d c0(1, 0, 0); // eventually to be created from freestream direction vector

	bool DOIflag = false;
	Eigen::Vector3d Q, Pc, Qc, triDiams;
	double B, RBcntr, x0, y0, dist, triRad;
	B = sqrt(pow(Mach,2) - 1);
	Q = center;

	// Hyperbolic distance from POI to panel center
	RBcntr = pow(P.x() - Q.x(), 2) - pow(B, 2)*pow(P.y() - Q.y(), 2) - pow(B, 2)*pow(P.z() - Q.z(), 2);

	// CSYS shift for panel-to-Mach cone distance check
	Pc = P - P;
	Qc = Q - P;
	x0 = (Qc - Pc).dot(windDir);
	y0 = sqrt(pow((Qc - Pc).norm(), 2) - pow(x0, 2));
	if (y0 >= B * x0)
	{
		dist = sqrt(pow(x0 + B * y0, 2) / (1 + pow(B, 2)));
	}
	else
	{
		dist = sqrt(pow((Qc - Pc).norm(), 2));
	}

	for (size_t i = 0; i < nodes.size(); i++)
	{

		if (i != nodes.size() - 1)
		{
			triDiams[i] = (nodes[i]->getPnt() - nodes[i + 1]->getPnt()).norm();
		}
		else
		{
			triDiams[i] = (nodes[0]->getPnt() - nodes[i]->getPnt()).norm();
		}
	}
	triRad = triDiams.maxCoeff() / 2;

	if ((Pc - Qc).dot(windDir) >= 0 || dist < triRad)
	{
		if (dist > triRad) // check if P is further from Mach cone than panel radius
		{
			if (RBcntr >= 0)
			{
				return DOIflag = true;
			}
		}
		else
		{
			double RBvert;
			Eigen::Vector3d vert;
			std::vector<bool> pntFlags;
			for (size_t i = 0; i < nodes.size(); i++)
			{
				vert = nodes[i]->getPnt();
				RBvert = pow(P.x() - vert.x(), 2) - pow(B, 2)*pow(P.y() - vert.y(), 2) - pow(B, 2)*pow(P.z() - vert.z(), 2);
				if (RBvert >= 0 && (P - vert).dot(windDir) >= 0)
				{
					return DOIflag = true;
				}
			}

			// Use Mach wedge calc to check later after coord. trans.
			return DOIflag = true;

			/////////////////////////////////////////////////////////////
			//supOutputGeom(POI, true);
			/////////////////////////////////////////////////////////////

			//// Find intersection point of Mach cone and panel
			//Eigen::Vector3d interPnt;
			//interPnt = supConePanelInter(P, Mach, windDir);

			//// Check if inter. pnt. is within perimiter of panel or downstream of panel
			//Eigen::Vector3d p0, p1, p2, AB, AC, AP, myNorm;
			//double alpha, beta, gamma;
			//p0 = nodes[0]->getPnt();
			//p1 = nodes[1]->getPnt();
			//p2 = nodes[2]->getPnt();
			//AB = p1 - p0;
			//AC = p2 - p0;
			//AP = interPnt - p0;

			//// Need to use normal computed here due to sign issues
			//myNorm = AB.cross(AC);

			//gamma = myNorm.dot(AB.cross(AP)) / myNorm.dot(normal);
			//beta = myNorm.dot(AP.cross(AC)) / myNorm.dot(normal);
			//alpha = 1.0 - gamma - beta;

			//if ((alpha >= 0) && (beta >= 0) && (gamma >= 0))
			//{
			//	return DOIflag = true;
			//}
			//else
			//{
			//	if ((interPnt - center).dot(windDir) >= 0)
			//	{
			//		return DOIflag = true;
			//	}
			//}
		}
	}

	return DOIflag;
}


void bodyPanel::supTransformPanel(const double Bmach, double alpha, double beta, const double M)
{
	// alpha and beta are input in degrees the transformed to radians

	alpha *= 180 / M_PI;
	beta *= 180 / M_PI;
	supSetG2LSmatrix(Bmach, alpha, beta, M);
	//supSetG2LSmatrixPilot(Bmach, alpha, beta, M);

	/*Eigen::Vector3d windDir;
	windDir << cos(alpha)*cos(beta), -sin(beta), sin(alpha)*cos(beta);
	supTransMat = supGetLocalSys(windDir);*/

	for (size_t i = 0; i < nodes.size(); i++)
	{
		supLocalNodes.push_back(supTransMat * (nodes[i]->getPnt() - center));
	}
}


 void bodyPanel::supSetG2LSmatrix(const double Bmach, const double a, const double b, const double M)
{
	// alpha and beta are input in radians
	double B = Bmach;
	Eigen::Vector3d cRef, nRef, nWind, scaleV, nWindScaled, vRef, uRef, col0, col1, col2;
	Eigen::Matrix3d ref2wind, machScaleB, machScaleC, machScaleWind;
	double s, r;
	//Eigen::Matrix3d transMat;

	if (a != 0 || b != 0)
	{
		// Freestream direction in reference CSYS
		cRef << cos(a)*cos(b), -sin(b), sin(a)*cos(b);

		// Reference to wind transformation
		ref2wind << cos(a)*cos(b), -sin(b), sin(a)*cos(b),
			cos(a)*sin(b), cos(b), sin(a)*sin(b),
			-sin(a), 0, cos(a);
	}
	else
	{
		cRef << 1, 0, 0;
		ref2wind.setIdentity();
	}

	// Matrices to remove B (i.e. Mach) from integration
	s = -1.0;
	machScaleB << s * pow(B, 2), 0, 0,
		0, 1, 0,
		0, 0, 1;
	machScaleC << 1, 0, 0,
		0, s*pow(B, 2), 0,
		0, 0, s*pow(B, 2);

	// Combined wind and Mach scaling matrix
	machScaleWind = ref2wind.transpose() * machScaleC * ref2wind;

	// Panel normal in reference CSYS
	nRef = normal;

	// Stuff
	vRef = nRef.cross(cRef) / (nRef.cross(cRef)).norm();
	uRef = vRef.cross(nRef);
	r = 1.0;
	if (signbit(nRef.dot(machScaleB * nRef)))
	{
		r = -1.0;
	}

	col0 = (1.0 / sqrt(abs(nRef.dot(machScaleB * nRef)))) * machScaleWind * uRef;
	col1 = (r * s / B) * machScaleWind * vRef;
	col2 = B * nRef / sqrt(abs(nRef.dot(machScaleB * nRef)));

	/*col0 = machScaleWind * uRef;
	col1 = (r * s / B) * machScaleWind * vRef;
	col2 = B * nRef;*/

	supTransMat.col(0) = col0;
	supTransMat.col(1) = col1;
	supTransMat.col(2) = col2;

	supTransMat.transposeInPlace();

	//////////////////////////// Correction term just greater than 1 for 5 deg test case -> small effect on results
	/*Eigen::Vector3d POI;
	supOutputGeom(POI, false);*/
	double nxCo = (ref2wind.row(0)).dot(nRef);
	//supAreaCorrect = B / sqrt(1 - pow(M, 2)*pow(nxCo, 2));
	supAreaCorrect = B * sqrt(1 - pow(M, 2)*pow(nxCo, 2));
}


 void bodyPanel::supSetG2LSmatrixPilot(const double Bmach, const double a, const double b, const double M)
 {
	 // alpha and beta are input in radians
	 double B = Bmach;
	 Eigen::Vector3d cRef, nRef, nC, coNorm;
	 Eigen::Matrix3d ref2wind, Bmat, A1, A2;
	 double s, gam, denom;
	 //Eigen::Matrix3d transMat;

	 if (a != 0 || b != 0)
	 {
		 // Freestream direction in reference CSYS
		 cRef << cos(a)*cos(b), -sin(b), sin(a)*cos(b);

		 // Reference to wind transformation
		 ref2wind << cos(a)*cos(b), -sin(b), sin(a)*cos(b),
			 cos(a)*sin(b), cos(b), sin(a)*sin(b),
			 -sin(a), 0, cos(a);
	 }
	 else
	 {
		 cRef << 1, 0, 0;
		 ref2wind.setIdentity();
	 }

	 //// Matrices to remove B (i.e. Mach) from integration
	 s = -1.0;
	 Bmat << s * pow(B, 2), 0, 0,
		 0, 1, 0,
		 0, 0, 1;

	 // Panel normal in reference CSYS
	 nRef = normal;

	 nC << nRef.dot(ref2wind.row(0)), nRef.dot(ref2wind.row(1)), nRef.dot(ref2wind.row(2));

	 // Panel conormal
	 coNorm = Bmat * nC;
	 gam = sqrt(pow(nC.y(), 2) + pow(nC.z(), 2));
	 denom = sqrt(abs(nC.dot(coNorm)));

	 A1 << 1, 0, 0,
		 0, B*nC.z() / gam, -B * nC.y() / gam,
		 0, B*nC.y() / gam, B * nC.z() / gam;

	 A2 << gam / denom, 0, -s * B*nC.x() / denom,
		 0, 1, 0,
		 B*nC.x() / denom, 0, gam / denom;
	 //sin(nRef.dot(coNorm))

	 supTransMat = A2 * A1 * ref2wind;

	 /*Eigen::Vector3d POI;
	 supOutputGeom(POI, false);
	 std::cout << "\n" << ref2wind << std::endl;
	 std::cout << "\n" << A1 << std::endl;
	 std::cout << "\n" << A2 << std::endl;
	 std::cout << "\n" << supTransMat << std::endl;*/

	 //supTransMat.transposeInPlace();
 }


Eigen::Matrix3d bodyPanel::supGetLocalSys(Eigen::Vector3d &windDir)
{
	// Same as original local CSYS, except x is in direction of freestream
	Eigen::Matrix3d transLocalFS;
	//Eigen::Vector3d FSdir(1, 0, 0);

	transLocalFS.row(1) = windDir.cross(normal) / (windDir.cross(normal)).norm();
	transLocalFS.row(0) = normal.cross(transLocalFS.row(1));
	transLocalFS.row(2) = normal;

	return transLocalFS;
}


Eigen::Vector3d bodyPanel::supConePanelInter(const Eigen::Vector3d &POI, const double Mach, Eigen::Vector3d &windDir)
{
	// Initialize
	Eigen::Matrix3d transLocalFS;
	Eigen::Vector3d POIlocal, n, w, u, locInterPnt, interPnt;
	double h, uz, s;
	bool zFlag = false;

	// Local CSYS w/ x-aligned freestream
	transLocalFS = supGetLocalSys(windDir);

	// Convert POI to local CSYS
	POIlocal = transLocalFS * (POI - center);

	h = (POI - center).dot(transLocalFS.row(2));
	if (signbit(normal.z()) == signbit(h))
	{
		uz = -1.0 / Mach;
	}
	else if (signbit(normal.z()) != signbit(h))
	{
		uz = 1.0 / Mach;
		if (normal.x() == 0 && normal.z() == 0)
		{
			zFlag = true;
		}
	}
	else
	{
		std::cout << "uh oh" << std::endl;
	}
	u << -sqrt(1 - pow(1.0 / Mach, 2)), 0, uz;
	if (!zFlag)
	{
		u.z() = transLocalFS.row(2) * u;
	}
	u = u / u.norm();
	n << 0, 0, 1;
	w = POIlocal;
	s = -n.dot(w) / n.dot(u);

	locInterPnt = POIlocal + s * u;

	//transLocalFS.transposeInPlace();
	interPnt = transLocalFS.transpose() * locInterPnt + center;

	return interPnt;
}


void bodyPanel::supPhiInf(const Eigen::Vector3d &POI, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow, double &srcPhi, bool DOIflag, const double Mach, Eigen::Vector3d &windDir)
{
	if (DOIflag)
	{
		///////////////////////////////////////////////////////////////////////////
		//supOutputGeom(POI, true);

		// Initialize
		double x, y, z, x1, y1, x2, y2;
		double m, lam, xm, xmc, ym1, ym2, ym1c, ym2c, s1, s2, sm1, sm2, r1, r2, R1, R2;
		Eigen::Vector2d integralCoeffs; // [Q1, w0]
		double srcCoeff = 0;
		Eigen::Vector3d localPOI, dubCoeffs, dubCoeffMat, dubVertsPhi;
		Eigen::Matrix3d dubVertsMat;
		bool mFlag, zFlag;
		zFlag = false;
		double Bmach = sqrt(pow(Mach, 2) - 1);

		integralCoeffs.setZero();
		dubCoeffs.setZero();

		// if changed, don't forget about big eps
		double epsGen = 1.0e-10;

		localPOI = supTransMat * (POI - center);
		x = localPOI.x();
		y = localPOI.y();
		z = localPOI.z();

		/////////////////////////////////////////// Need to verify -- this is true, panel trans makes this true
		Eigen::Vector3d localWindDir(1, 0, 0);

		// Iterate through panel edges
		for (nodes_index_type i = 0; i < nodes.size(); i++)
		{
			integralCoeffs.setZero();
			mFlag = false;
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			double eps1;
			double eps2;
			if (i != nodes.size() - 1)
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[i + 1];

				/*eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[i + 1]->linGetCPoffset();*/
			}
			else
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[0];

				/*eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[0]->linGetCPoffset();*/
			}

			eps1 = 1.0e-8;
			eps2 = 1.0e-8;

			x1 = p1.x();
			y1 = p1.y();
			x2 = p2.x();
			y2 = p2.y();

			// Compute inf. coeff. geometric values
			m = (y2 - y1) / (x2 - x1);
			s1 = y - y1;
			s2 = y - y2;
			r1 = sqrt(pow(s1, 2) + pow(z, 2));
			r2 = sqrt(pow(s2, 2) + pow(z, 2));

			////////////////////////////////////////////////////////////////////////////////////////////////////////
			//supOutputGeom(POI, true);

			if (abs(m) < epsGen) // Edge parallel to freestream
			{
				mFlag = true;
				m = 0;
				xmc = -s1;
				ym1c = -(x - x1);
				ym2c = -(x - x2);
				sm1 = x - x1;
				sm2 = x - x2;

				xm = xmc;
			}
			else if (abs(m) > 1.0 / epsGen) // Edge perinducular to freestream
			{
				lam = 0;
				xm = (x - x1) - (y - y1)*lam;
				sm1 = xm + s1 * lam;
				sm2 = xm + s2 * lam;
				ym1 = s1 - sm1 * lam;
				ym2 = s2 - sm2 * lam;
			}
			else
			{
				lam = 1 / m;
				xm = (x - x1) - (y - y1) / m;
				sm1 = xm + s1 / m;
				sm2 = xm + s2 / m;
				ym1 = s1 - sm1 / m;
				ym2 = s2 - sm2 / m;
				xmc = m * xm;
				ym1c = m * ym1;
				ym2c = m * ym2;
			}

			R1 = 0;
			R2 = 0;
			if ((localPOI - p1).dot(localWindDir) >= 0)
			{
				R1 = sqrt(pow(sm1, 2) - pow(r1, 2));
				if (isnan(R1))
				{
					R1 = 0;
				}
			}
			if ((localPOI - p2).dot(localWindDir) >= 0)
			{
				R2 = sqrt(pow(sm2, 2) - pow(r2, 2));
				if (isnan(R2))
				{
					R2 = 0;
				}
			}

			// Check that edge intersects Mach cone and upstream of POI
			if (R1 > abs(eps1) || R2 > abs(eps2))
			{
				//////////////////////////////////////////////////////////
				//supOutputGeom(POI, true);
				//////////////////////////////////////////////////////////

				if (abs(z) > 1e-10)
				{
					if (abs(m) < 1) // subsonic edge
					{
						integralCoeffs = supEdgeInfSub(R1, R2, ym1c, ym2c, xmc, m, z, eps1, eps2, mFlag);
					}
					else if (abs(m) > 1) // supersonic edge
					{
						integralCoeffs = supEdgeInfSup(R1, R2, ym1, ym2, xm, lam, z, eps1, eps2);
					}
					else if (abs(1 - abs(m)) < epsGen) // sonic edge
					{
						integralCoeffs = supEdgeInfSon(ym1, ym2, xmc, ym1c, ym2c, R1, R2, lam, z);
					}
				}
				else // z = 0
				{
					zFlag = true;
					if (abs(1 - abs(m)) < epsGen) // sonic edge
					{

					}
					else if (abs(m) < 1) // subsonic edge
					{

					}
					else if (abs(m) > 1) // supersonic edge
					{

					}
				}

				if (isnan(integralCoeffs[0]) || isnan(integralCoeffs[1]))
				{
					double stuff = 0;
				}
			}
			else if (abs(m) > 1)
			{
				double wedgeCheck = pow(xm, 2) + (pow(z, 2) / pow(m, 2)) - pow(z, 2);
				if ((xm > 0) && (signbit(ym1) != signbit(ym2)) && (abs(ym1) > epsGen && abs(ym2) > epsGen) && (wedgeCheck >= 0))
				{
					////////////////////////////////////////////////////////////////////////////////////////////////////////
					//supOutputGeom(POI, true);

					double Q11, w01, Q12, w02, Q1, w0;

					double s1, t1;
					s1 = 1.0;
					t1 = 1.0;
					if (signbit(z * ym1))
					{
						s1 = -1.0;
					}
					if (signbit(ym1))
					{
						t1 = -1.0;
					}

					Q11 = s1 * M_PI / 2.0;
					w01 = t1 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					double s2, t2;
					s2 = 1.0;
					t2 = 1.0;
					if (signbit(z * ym2))
					{
						s2 = -1.0;
					}
					if (signbit(ym2))
					{
						t2 = -1.0;
					}
					Q12 = s2 * M_PI / 2.0;
					w02 = t2 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					Q1 = Q12 - Q11;
					w0 = w02 - w01;
					integralCoeffs[0] = Q1;
					integralCoeffs[1] = w0;

					if (isnan(integralCoeffs[0]) || isnan(integralCoeffs[1]))
					{
						double stuff = 0;
					}
				}
			}

			if (isnan(integralCoeffs[0]) || isnan(integralCoeffs[1]))
			{
				double stuff = 0;
			}

			// Sum doublet coefficients 
			// integralCoeffs = [Q1 w0]
			// dubCoeffs = [Q1 w0 w0/m]
			// srcCoeff = [xm*w0]
			dubCoeffs[0] += integralCoeffs[0];
			dubCoeffs[1] += integralCoeffs[1];
			if (!mFlag)
			{
				dubCoeffs[2] += integralCoeffs[1] / m;
			}

			// Sum source coefficient
			srcCoeff += xm * integralCoeffs[1];
			//srcCoeff += (xm * integralCoeffs[1]) / supAreaCorrect;
		}

		// Compute final source coefficient
		double Q1sum = dubCoeffs[0];
		srcPhi = (srcCoeff - z * Q1sum) / (2.0 * M_PI);
		//srcPhi = ((srcCoeff - z * Q1sum) / supAreaCorrect) / (2.0 * M_PI);
		//srcPhi = ((srcCoeff - z * Q1sum) / (2.0 * M_PI)) / supAreaCorrect;
		//srcPhi = ((srcCoeff - z * Q1sum) / (2.0 * M_PI)) * supAreaCorrect;

		// Compute final doublet coefficients
		if (zFlag)
		{
			z = 0;
		}
		dubCoeffs[0] = -dubCoeffs[0];
		dubCoeffs[1] = -z * dubCoeffs[1];
		dubCoeffs[2] = -z * dubCoeffs[2];
		dubCoeffs = dubCoeffs / (2.0 * M_PI);

		// Build matrix of doublet coefficients
		dubCoeffMat[0] = dubCoeffs[0];
		dubCoeffMat[1] = x * dubCoeffs[0] - dubCoeffs[1];
		dubCoeffMat[2] = y * dubCoeffs[0] - dubCoeffs[2];

		// Convert linear doublet equation to vertex based linear doublet equations
		dubVertsMat = supVertsMatrix(supLocalNodes);

		// Compute influence of each vertex on the field point
		dubVertsPhi = dubCoeffMat.transpose() * dubVertsMat.inverse();

		// Fill A matrix
		// Need to make this more efficient
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset()) // using orig. node should be fine
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
			else
			{
				Arow[nodes[i]->getIndex()] += dubVertsPhi(i);
			}
		}
	}
	else
	//if (!DOIflag)
	{
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset()) // using orig. node should be fine
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////// Funcion for testing against MATLAB
void bodyPanel::supOutputGeom(const Eigen::Vector3d &POI, bool outPOI)
{
	std::ofstream fid;
	std::string	myFile = "D:\\Desktop\\Thesis\\Code\\MATLAB\\Supersonic\\stuff.csv";
	fid.open(myFile);

	fid << nodes[0]->getPnt().x() << "," << nodes[0]->getPnt().y() << "," << nodes[0]->getPnt().z() << "\n";
	fid << nodes[1]->getPnt().x() << "," << nodes[1]->getPnt().y() << "," << nodes[1]->getPnt().z() << "\n";
	fid << nodes[2]->getPnt().x() << "," << nodes[2]->getPnt().y() << "," << nodes[2]->getPnt().z() << "\n";
	if (outPOI)
	{
		fid << POI.x() << "," << POI.y() << "," << POI.z();
		fid << "\n";
	}
	fid.close();
}


Eigen::Vector2d bodyPanel::supEdgeInfSon(const double ym1, const double ym2, const double xmc, const double ym1c, const double ym2c, const double R1, const double R2, const double lam, const double z)
{
	Eigen::Vector2d integralCoeffs;



	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSub(const double R1, const double R2, const double ym1c, const double ym2c, const double xmc, const double m, const double z, const double eps1, const double eps2, bool mFlag)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1cNeg, ym2cNeg;
	ym1cNeg = -ym1c;
	ym2cNeg = -ym2c;

	if (R1 > eps1 && R2 > eps2)
	{
		Q1 = atan2(z*xmc * (ym1cNeg*R2 - ym2cNeg * R1), pow(z, 2) * ym1cNeg * ym2cNeg + pow(xmc, 2)*R1*R2);
		if (mFlag)
		{
			//w0 = (1.0 / sqrt(1 - pow(m, 2))) * log((ym2cNeg + R2 * sqrt(1 - pow(m, 2))) / (ym1cNeg + R1 * sqrt(1 - pow(m, 2))));
			w0 = log((ym2cNeg + R2) / (ym1cNeg + R1));
			if (isnan(w0))
			{
				double stuff = 0;
			}
		}
		else
		{
			w0 = (m / sqrt(1 - pow(m, 2))) * log((ym2cNeg + R2 * sqrt(1 - pow(m, 2))) / (ym1cNeg + R1 * sqrt(1 - pow(m, 2))));
		}
	}
	else
	{
		double Q11, Q12, w01, w02, s;
		s = 1.0;
		if (signbit(z))
		{
			s = -1.0;
		}

		if (mFlag)
		{
			/*Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
			w01 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));

			Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
			w02 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));*/

			if (R1 > eps1)
			{
				Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
				w01 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q11 = 0;
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
				w02 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q12 = 0;
				w02 = 0;
			}
		}
		else
		{
			/*Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
			Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
			if (abs(ym1c) > eps1)
			{
				w01 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				w01 = 0;
			}
			if (abs(ym2c) > eps1)
			{
				w02 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				w02 = 0;
			}*/
			
			/*Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
			w01 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));

			Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
			w02 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));*/

			if (R1 > eps1)
			{
				Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
				w01 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q11 = 0;
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
				w02 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q12 = 0;
				w02 = 0;
			}
		}

		Q1 = Q12 - Q11;
		w0 = w02 - w01;

		if (isnan(Q1) || isnan(w0))
		{
			double stuff = 0;
		}

	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSup(const double R1, const double R2, const double ym1, const double ym2, const double xm, const double lam, const double z, const double eps1, const double eps2)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1Neg, ym2Neg;
	ym1Neg = -ym1;
	ym2Neg = -ym2;

	if (R1 > eps1 && R2 > eps2)
	{
		Q1 = atan2(z * xm * (ym1Neg*R2 - ym2Neg * R1), pow(z,2)*ym2Neg*ym1Neg + pow(xm,2)*R1*R2);
		w0 = (1 / sqrt(1 - pow(lam, 2))) * atan2(sqrt(1 - pow(lam, 2)) * (ym1Neg*R2 - ym2Neg * R1), ym1Neg*ym2Neg + (1 - pow(lam, 2))*R1*R2);
	}
	else
	{
		double Q11, Q12, w01, w02;
		if (R1 > eps1)
		{
			Q11 = atan2(z*ym1, xm*R1);
			w01 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym1, R1*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double s1, t1;
			s1 = 1.0;
			t1 = 1.0;
			if (signbit(z * ym1))
			{
				s1 = -1.0;
			}
			if (signbit(ym1))
			{
				t1 = -1.0;
			}

			Q11 = s1 * M_PI / 2;
			w01 = t1 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}
		if (R2 > eps2)
		{
			Q12 = atan2(z*ym2, xm*R2);
			w02 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym2, R2*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double s2, t2;
			s2 = 1.0;
			t2 = 1.0;
			if (signbit(z * ym2))
			{
				s2 = -1.0;
			}
			if (signbit(ym2))
			{
				t2 = -1.0;
			}
			Q12 = s2 * M_PI / 2;
			w02 = t2 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}
		Q1 = Q12 - Q11;
		w0 = w02 - w01;
	}

	if (isnan(Q1) || isnan(w0))
	{
		double stuff;
	}

	if (isnan(integralCoeffs[0]) || isnan(integralCoeffs[1]))
	{
		double stuff = 0;
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


//Eigen::Matrix3d bodyPanel::supG2LSmatrix(const double Bmach)
//{
//	double a, b;	// represent AoA and sideslip for now
//	a = 0;
//	b = 0;
//
//	double B = Bmach;
//	Eigen::Vector3d cRef, nRef, nWind, scaleV, nWindScaled, vRef, uRef;
//	Eigen::Matrix3d ref2wind, machScale, machScaleWind;
//	double s, r;
//	Eigen::Matrix3d transMat;
//
//	if (a != 0 || b != 0)
//	{
//		// Convert AoA and sideslip to radians
//		a = a * M_PI / 180.0;
//		a = b * M_PI / 180.0;
//
//		// Freestream direction in reference CSYS
//		cRef << cos(a)*cos(b), -sin(b), sin(a)*cos(b);
//
//		// Reference to wind transformation
//		ref2wind << cos(a)*cos(b), -sin(b), sin(a)*cos(b),
//			cos(a)*sin(b), cos(b), sin(a)*sin(b),
//			-sin(a), 0, cos(a);
//	}
//	else
//	{
//		cRef << 1, 0, 0;
//		ref2wind.setIdentity();
//	}
//
//	// Matrix to remove B (i.e. Mach) from integration
//	s = -1.0;
//	machScale << 1, 0, 0,
//		0, s*pow(B, 2), 0,
//		0, 0, s*pow(B, 2);
//
//	machScaleWind = ref2wind.transpose() * machScale * ref2wind;
//	//std::cout << machScaleWind << std::endl;
//
//	// Panel normal in reference CSYS
//	nRef = normal;
//
//	// Panel normal in wind CSYS
//	nWind = ref2wind * nRef;
//
//	// Panel normal in wind CSYS and scaled
//	nWindScaled = nWind;
//	nWindScaled[0] = nWind[0] * -pow(B, 2);
//	/*std::cout << nWindScaled << std::endl;*/
//
//	vRef = nRef.cross(cRef) / (nRef.cross(cRef)).norm();
//	uRef = vRef.cross(nRef);
//	/*std::cout << vRef << std::endl;
//	std::cout << uRef << std::endl;*/
//	r = 1.0;
//	if (signbit(nRef.dot(nRef)))
//	{
//		r = -1.0;
//	}
//
//	transMat << (1 / sqrt(abs(nRef.dot(nRef)))) * machScaleWind * uRef,
//		r*s / B * machScaleWind * vRef,
//		B*nRef / sqrt(abs(nRef.dot(nRef)));
//
//	transMat.transposeInPlace();
//	//std::cout << "\n" << transMat << std::endl;
//
//	return transMat;
//}



//std::vector<Eigen::Vector3d> bodyPanel::supTransformPanel(Eigen::Vector3d &localPOI, double Bmach, double alpha, double beta)
//{
//	std::vector<Eigen::Vector3d> localNodes;
//	alpha *= 180 / M_PI;
//	beta *= 180 / M_PI;
//	supSetG2LSmatrix(Bmach, alpha, beta);
//	//Eigen::Matrix3d transMat = supG2LSmatrix(Bmach);
//
//	/*for (size_t i = 0; i < nodes.size(); i++)
//	{
//		localNodes.push_back(transMat * (nodes[i]->getPnt() - center));
//	}
//	localPOI = transMat * (localPOI - center);*/
//
//	for (size_t i = 0; i < nodes.size(); i++)
//	{
//		localNodes.push_back(supTransMat * (nodes[i]->getPnt() - center));
//	}
//	localPOI = supTransMat * (localPOI - center);
//
//	return localNodes;
//}


//
//bool bodyPanel::supEdgeCheck(edge* myEdge, Eigen::Vector3d &P, const double B)
//{
//	bool edgeFlag = false;
//	double x, y, z, xi1, eta1, zeta1, xi2, eta2, zeta2, m, n, l, x1, y1, z1, x2, y2, z2;
//	/*Eigen::Vector3d pnt1, pnt2;
//	Eigen::Vector3d c0(1, 0, 0);*/
//
//	// Just finished this function. Need to write last condition and test in MATLAB, then test here  
//
//	x = P.x();
//	y = P.y();
//	z = P.z();
//	xi1 = myEdge->getN1()->getPnt().x();
//	eta1 = myEdge->getN1()->getPnt().y();
//	zeta1 = myEdge->getN1()->getPnt().z();
//	xi2 = myEdge->getN2()->getPnt().x();
//	eta2 = myEdge->getN2()->getPnt().y();
//	zeta2 = myEdge->getN2()->getPnt().z();
//
//	m = (eta2 - eta1) / (xi2 - xi1);
//	n = (zeta2 - zeta1) / (eta2 - eta1);
//	l = (zeta2 - zeta1) / (xi2 - xi1);
//
//	if (xi1 == xi2)
//	{
//		x1 = xi1;
//		x2 = xi1;
//		if (eta1 == eta2)
//		{
//			y1 = eta1;
//			y2 = eta1;
//			z1 = z + sqrt(pow((xi1 - x), 2) / pow(B, 2) - pow((eta1 - y), 2));
//			z2 = z - sqrt(pow((xi1 - x), 2) / pow(B, 2) - pow((eta1 - y), 2));
//		}
//		else if (zeta1 == zeta2)
//		{
//			z1 = zeta1;
//			z2 = zeta1;
//			y1 = y + sqrt(pow((xi1 - x), 2) / pow(B, 2) - pow((zeta1 - z), 2));
//			y2 = y - sqrt(pow((xi1 - x), 2) / pow(B, 2) - pow((zeta1 - z), 2));
//		}
//		else
//		{
//			z1 = (B*zeta1 + n * sqrt(-pow(B, 2) * pow(eta1, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(n, 2) * y - 2 * pow(B, 2) * eta1*n*z + 2 * pow(B, 2) * eta1*n*zeta1 - pow(B, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * n*y*z - 2 * pow(B, 2) * n*y*zeta1 - pow(B, 2) * pow(z, 2) + 2 * pow(B, 2) * z*zeta1 - pow(B, 2) * pow(zeta1, 2) + pow(n, 2) * pow(x, 2) - 2 * pow(n, 2) * x*xi1 + pow(n, 2) * pow(xi1, 2) + pow(x, 2) - 2 * x*xi1 + pow(xi1, 2)) - B * eta1*n + B * n*y + B * pow(n, 2) * z) / (B*pow(n, 2) + B);
//			y1 = (z1 - zeta1) / n + eta1;
//
//			z2 = (B*zeta1 - n * sqrt(-pow(B, 2) * pow(eta1, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(n, 2) * y - 2 * pow(B, 2) * eta1*n*z + 2 * pow(B, 2) * eta1*n*zeta1 - pow(B, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * n*y*z - 2 * pow(B, 2) * n*y*zeta1 - pow(B, 2) * pow(z, 2) + 2 * pow(B, 2) * z*zeta1 - pow(B, 2) * pow(zeta1, 2) + pow(n, 2) * pow(x, 2) - 2 * pow(n, 2) * x*xi1 + pow(n, 2) * pow(xi1, 2) + pow(x, 2) - 2 * x*xi1 + pow(xi1, 2)) - B * eta1*n + B * n*y + B * pow(n, 2) * z) / (B*pow(n, 2) + B);
//			y2 = (z2 - zeta1) / n + eta1;
//		}
//
//		if ((!isnan(y1) && !isnan(z1)) || (!isnan(y2) && !isnan(z2)))
//		{
//			edgeFlag = true;
//		}
//	}
//	else if (eta1 == eta2)
//	{
//		y1 = eta1;
//		y2 = eta1;
//		if (zeta1 == zeta2)
//		{
//			z1 = zeta1;
//			z2 = zeta2;
//			x1 = x + sqrt(pow(B, 2) * pow((eta1 - y), 2) + pow((zeta1 - z), 2));
//			x2 = x - sqrt(pow(B, 2) * pow((eta1 - y), 2) + pow((zeta1 - z), 2));
//		}
//		else
//		{
//			z1 = (l*xi1 - zeta1 + l * sqrt(-pow(B, 2) * pow(eta1, 2) * pow(l, 2) + pow(B, 2) * pow(eta1, 2) + 2 * pow(B, 2) * eta1*pow(l, 2) * y - 2 * pow(B, 2) * eta1*y - pow(B, 2) * pow(l, 2) * pow(y, 2) + pow(B, 2) * pow(y, 2) + pow(l, 2) * pow(xi1, 2) + 2 * l*xi1*z - 2 * l*xi1*zeta1 + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) + pow(l, 2) * z) / (pow(l, 2) - 1);
//			x1 = (z1 - zeta1) / l + xi1;
//
//			z2 = -(zeta1 - l * xi1 + l * sqrt(-pow(B, 2) * pow(eta1, 2) * pow(l, 2) + pow(B, 2) * pow(eta1, 2) + 2 * pow(B, 2) * eta1*pow(l, 2) * y - 2 * pow(B, 2) * eta1*y - pow(B, 2) * pow(l, 2) * pow(y, 2) + pow(B, 2) * pow(y, 2) + pow(l, 2) * pow(xi1, 2) + 2 * l*xi1*z - 2 * l*xi1*zeta1 + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) - pow(l, 2) * z) / (pow(l, 2) - 1);
//			x2 = (z2 - zeta1) / l + xi1;
//		}
//
//		if ((!isnan(x1) && !isnan(z1)) || (!isnan(x2) && !isnan(z2)))
//		{
//			edgeFlag = true;
//		}
//	}
//	else if (zeta1 == zeta2)
//	{
//		z1 = zeta1;
//		z2 = zeta1;
//
//		x1 = (B*sqrt(-pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) - x + pow(B, 2) * pow(m, 2) * xi1 - pow(B, 2) * eta1*m + pow(B, 2) * m*y) / (pow(B, 2) * pow(m, 2) - 1);
//		y1 = m * (x1 - xi1) + eta1;
//
//		x2 = -(x + B * sqrt(-pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) - pow(B, 2) * pow(m, 2) * xi1 + pow(B, 2) * eta1*m - pow(B, 2) * m*y) / (pow(B, 2) * pow(m, 2) - 1);
//		y2 = m * (x2 - xi1) + eta1;
//
//		if ((!isnan(x1) && !isnan(y1)) || (!isnan(x2) && !isnan(y2)))
//		{
//			edgeFlag = true;
//		}
//	}
//	else
//	{
//		// Intersection 1
//		x1 = (y - eta1 + m * xi1 + (y - eta1 - m * x + m * xi1 + B * m*sqrt(-pow(B, 2) * pow(eta1, 2) * pow(m, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) * y - 2 * pow(B, 2) * eta1*pow(m, 2) * n*z + 2 * pow(B, 2) * eta1*pow(m, 2) * n*zeta1 - pow(B, 2) * pow(m, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * pow(m, 2) * n*y*z - 2 * pow(B, 2) * pow(m, 2) * n*y*zeta1 - pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(n, 2) * pow(x, 2) - 2 * pow(m, 2) * pow(n, 2) * x*xi1 + pow(m, 2) * pow(n, 2) * pow(xi1, 2) + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*n*x*z + 2 * m*n*x*zeta1 + 2 * m*n*xi1*z - 2 * m*n*xi1*zeta1 - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) + pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) - pow(B, 2) * pow(m, 2) * pow(n, 2) * y + pow(B, 2) * pow(m, 2) * n*z - pow(B, 2) * pow(m, 2) * n*zeta1) / (pow(B, 2) * pow(m, 2) * pow(n, 2) + pow(B, 2) * pow(m, 2) - 1)) / m;
//		z1 = zeta1 - eta1 * n + n * y + (n*(y - eta1 - m * x + m * xi1 + B * m*sqrt(-pow(B, 2) * pow(eta1, 2) * pow(m, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) * y - 2 * pow(B, 2) * eta1*pow(m, 2) * n*z + 2 * pow(B, 2) * eta1*pow(m, 2) * n*zeta1 - pow(B, 2) * pow(m, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * pow(m, 2) * n*y*z - 2 * pow(B, 2) * pow(m, 2) * n*y*zeta1 - pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(n, 2) * pow(x, 2) - 2 * pow(m, 2) * pow(n, 2) * x*xi1 + pow(m, 2) * pow(n, 2) * pow(xi1, 2) + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*n*x*z + 2 * m*n*x*zeta1 + 2 * m*n*xi1*z - 2 * m*n*xi1*zeta1 - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) + pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) - pow(B, 2) * pow(m, 2) * pow(n, 2) * y + pow(B, 2) * pow(m, 2) * n*z - pow(B, 2) * pow(m, 2) * n*zeta1)) / (pow(B, 2) * pow(m, 2) * pow(n, 2) + pow(B, 2) * pow(m, 2) - 1);
//		y1 = m * (x1 - xi1) + eta1;
//
//		// Intersection 2
//		x2 = -(eta1 - y - m * xi1 + (eta1 - y + m * x - m * xi1 + B * m*sqrt(-pow(B, 2) * pow(eta1, 2) * pow(m, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) * y - 2 * pow(B, 2) * eta1*pow(m, 2) * n*z + 2 * pow(B, 2) * eta1*pow(m, 2) * n*zeta1 - pow(B, 2) * pow(m, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * pow(m, 2) * n*y*z - 2 * pow(B, 2) * pow(m, 2) * n*y*zeta1 - pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(n, 2) * pow(x, 2) - 2 * pow(m, 2) * pow(n, 2) * x*xi1 + pow(m, 2) * pow(n, 2) * pow(xi1, 2) + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*n*x*z + 2 * m*n*x*zeta1 + 2 * m*n*xi1*z - 2 * m*n*xi1*zeta1 - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) - pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) + pow(B, 2) * pow(m, 2) * pow(n, 2) * y - pow(B, 2) * pow(m, 2) * n*z + pow(B, 2) * pow(m, 2) * n*zeta1) / (pow(B, 2) * pow(m, 2) * pow(n, 2) + pow(B, 2) * pow(m, 2) - 1)) / m;
//		z2 = zeta1 - eta1 * n + n * y + (n*(y - eta1 - m * x + m * xi1 + B * m*sqrt(-pow(B, 2) * pow(eta1, 2) * pow(m, 2) * pow(n, 2) + 2 * pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) * y - 2 * pow(B, 2) * eta1*pow(m, 2) * n*z + 2 * pow(B, 2) * eta1*pow(m, 2) * n*zeta1 - pow(B, 2) * pow(m, 2) * pow(n, 2) * pow(y, 2) + 2 * pow(B, 2) * pow(m, 2) * n*y*z - 2 * pow(B, 2) * pow(m, 2) * n*y*zeta1 - pow(B, 2) * pow(m, 2) * pow(z, 2) + 2 * pow(B, 2) * pow(m, 2) * z*zeta1 - pow(B, 2) * pow(m, 2) * pow(zeta1, 2) + pow(eta1, 2) + 2 * eta1*m*x - 2 * eta1*m*xi1 - 2 * eta1*y + pow(m, 2) * pow(n, 2) * pow(x, 2) - 2 * pow(m, 2) * pow(n, 2) * x*xi1 + pow(m, 2) * pow(n, 2) * pow(xi1, 2) + pow(m, 2) * pow(x, 2) - 2 * pow(m, 2) * x*xi1 + pow(m, 2) * pow(xi1, 2) - 2 * m*n*x*z + 2 * m*n*x*zeta1 + 2 * m*n*xi1*z - 2 * m*n*xi1*zeta1 - 2 * m*x*y + 2 * m*xi1*y + pow(y, 2) + pow(z, 2) - 2 * z*zeta1 + pow(zeta1, 2)) + pow(B, 2) * eta1*pow(m, 2) * pow(n, 2) - pow(B, 2) * pow(m, 2) * pow(n, 2) * y + pow(B, 2) * pow(m, 2) * n*z - pow(B, 2) * pow(m, 2) * n*zeta1)) / (pow(B, 2) * pow(m, 2) * pow(n, 2) + pow(B, 2) * pow(m, 2) - 1);
//		y2 = m * (x2 - xi1) + eta1;
//
//		if ((!isnan(x1) && !isnan(y1) && !isnan(z1)) || (!isnan(x2) && !isnan(y2) && !isnan(z2)))
//		{
//			edgeFlag = true;
//		}
//	}
//
//	return edgeFlag;
//
//	/*pnt1[0] = x1;
//	pnt1[1] = y1;
//	pnt1[2] = z1;
//	pnt2[0] = x2;
//	pnt2[1] = y2;
//	pnt2[2] = z2;*/
//
//	// Pretty sure I just need to check isnan()
//
//	//// fix this stuff. need to account for z != 0 in minPnt calc
//	// dot c0 to P, not pntMin
//
//
//	//Eigen::Vector3d pntMin(x - sqrt(pow(B, 2)*pow(z, 2)), y, z);
//	//if (abs(eta2 - y1) < abs(eta2 - eta1) && (pntMin - pnt1).dot(c0) >= 0)
//	//{
//	//	edgeFlag = true;
//	//}
//	//else if (abs(eta2 - y2) < abs(eta2 - eta1) && (pntMin - pnt2).dot(c0) >= 0)
//	//{
//	//	edgeFlag = true;
//	//}
//}


//// Don't actually need edgeFlags. Just do DOIflag, then check R1 and R2 when doing inf coeff calcs
//
//bool bodyPanel::supDOIcheck(Eigen::Vector3d &P,const double M)
//{
//	// P is POI, Q is panel center
//
//	Eigen::Vector3d c0(1,0,0); // eventually to be created from freestream direction vector
//
//	bool DOIflag = false;
//	Eigen::Vector3d Q, Pc, Qc, triDiams;
//	double B, RBcntr, x0, y0, dist, triRad;
//	B = sqrt(pow(M, 2) - 1);
//	Q = center;
//
//	// Hyperbolic distance from POI to panel center
//	RBcntr = pow(P.x() - Q.x(), 2) - pow(B, 2)*pow(P.y() - Q.y(), 2) - pow(B, 2)*pow(P.z() - Q.z(), 2);
//
//	// CSYS shift for panel-to-Mach cone distance check
//	Pc = P - P;
//	Qc = Q - P;
//	x0 = (Qc - Pc).dot(c0);
//	y0 = sqrt(pow((Qc - Pc).norm(), 2) - pow(x0, 2));
//	if (y0 >= B * x0)
//	{
//		dist = sqrt(pow(x0+B*y0,2) / (1+pow(B,2)));
//	}
//	else
//	{
//		dist = sqrt(pow((Qc-Pc).norm(),2));
//	}
//	
//	for (size_t i = 0; i < nodes.size(); i++)
//	{
//		
//			if (i != nodes.size() - 1)
//			{
//				triDiams[i] = (nodes[i]->getPnt() - nodes[i + 1]->getPnt()).norm();
//			}
//			else
//			{
//				triDiams[i] = (nodes[0]->getPnt() - nodes[i]->getPnt()).norm();
//			}
//	}
//	triRad = triDiams.maxCoeff()/2;
//
//	if ((Pc - Qc).dot(c0) >= 0 || dist < triRad)
//	{
//		if (dist > triRad)
//		{
//			if (RBcntr > 0)
//			{
//				DOIflag = true;
//				edgeFlags.push_back(true);
//				edgeFlags.push_back(true);
//				edgeFlags.push_back(true);
//			}
//		}
//		else
//		{
//			double RBvert;
//			Eigen::Vector3d vert;
//			std::vector<bool> pntFlags;
//			size_t count = 0;
//			for (size_t i = 0; i < nodes.size(); i++)
//			{
//				vert = nodes[i]->getPnt();
//				RBvert = pow(P.x() - vert.x(), 2) - pow(B, 2)*pow(P.y() - vert.y(), 2) - pow(B, 2)*pow(P.z() - vert.z(), 2);
//				if (RBvert > 0 && (P - vert).dot(c0))
//				{
//					count += 1;
//					pntFlags.push_back(true);
//				}
//				else
//				{
//					pntFlags.push_back(false);
//				}
//			}
//
//			if (count > 1)
//			{
//				DOIflag = true;
//				edgeFlags.push_back(true);
//				edgeFlags.push_back(true);
//				edgeFlags.push_back(true);
//			}
//			else if (count == 1)
//			{
//				DOIflag = true;
//				if (pntFlags[0])
//				{
//					edgeFlags.push_back(true);
//					edgeFlags.push_back(supEdgeCheck(getEdges()[1], P, B));
//					edgeFlags.push_back(true);
//				}
//				else if (pntFlags[1])
//				{
//					edgeFlags.push_back(true);					
//					edgeFlags.push_back(true);
//					edgeFlags.push_back(supEdgeCheck(getEdges()[2], P, B));
//				}
//				else if (pntFlags[2])
//				{
//					edgeFlags.push_back(supEdgeCheck(getEdges()[0], P, B));
//					edgeFlags.push_back(true);
//					edgeFlags.push_back(true);
//				}
//			}
//			else
//			{
//				edgeFlags.push_back(supEdgeCheck(getEdges()[0], P, B));
//				edgeFlags.push_back(supEdgeCheck(getEdges()[1], P, B));
//				edgeFlags.push_back(supEdgeCheck(getEdges()[2], P, B));
//
//				std::vector<bool> edgeFlagsFalse;
//				edgeFlagsFalse.push_back(false);
//				edgeFlagsFalse.push_back(false);
//				edgeFlagsFalse.push_back(false);
//				if (edgeFlags == edgeFlagsFalse)
//				{
//					DOIflag = false;
//				}
//			}
//		}
//	}
//	else
//	{
//		DOIflag = false;
//		edgeFlags.push_back(false);
//		edgeFlags.push_back(false);
//		edgeFlags.push_back(false);
//	}
//	
//	return DOIflag;
//}