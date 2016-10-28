
//  panel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panel.h"
#include "cpNode.h"
#include "edge.h"
#include "surface.h"

panel::panel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, int surfID) : ID(surfID), nodes(nodes), pEdges(pEdges), bezNormal(bezNorm)
{
    setGeom();
}

//panel::panel(const panel &copy) : ID(copy.ID), nodes(copy.nodes), pEdges(copy.pEdges)
//{
//    setGeom();
//}

void panel::setGeom()
{
    longSide = 0;
    
    for (int i=0; i<pEdges.size(); i++)
    {
        double l = pEdges[i]->length();
        if (l > longSide)
        {
            longSide = l;
        }
    }
    
    if (pEdges.size() == 3)
    {
        Eigen::Vector3d p0,p1,p2;
        Eigen::Vector3d a,b;
        p0 = nodes[0]->getPnt();
        p1 = nodes[1]->getPnt();
        p2 = nodes[2]->getPnt();
        a = p1-p0;
        b = p2-p0;
        center = (p0+p1+p2)/3;
        
        double theta = acos(a.dot(b)/(a.norm()*b.norm()));
        area = 0.5*a.norm()*b.norm()*sin(theta);
        normal = a.cross(b);
        normal.normalize();
        
        // If normals weren't included in input, set bezNormal to calculated normal
        if (bezNormal.isZero())
        {
            bezNormal = normal;
        }
        else
        {
            bezNormal.normalize();
        }
    }
    else if (pEdges.size() == 4)
    {
        Eigen::Vector3d p0,p1,p2,p3,m1,m2;
        Eigen::Vector3d a,b,c,d,p,q;
        p0 = nodes[0]->getPnt();
        p1 = nodes[1]->getPnt();
        p2 = nodes[2]->getPnt();
        p3 = nodes[3]->getPnt();
        a = p1-p0;
        b = p2-p1;
        c = p3-p2;
        d = p0-p3;
        p = b+c;
        q = a+b;
        m1 = p1+0.5*p;
        m2 = p0+0.5*q;
        center = 0.5*(m1+m2);
        area = 0.5*p.cross(q).norm();
        normal = a.cross(b);
        normal.normalize();
        
        // If normals weren't included in input, set bezNormal to calculated normal
        if (bezNormal.isZero())
        {
            bezNormal = normal;
        }
        else
        {
            bezNormal.normalize();
        }
    }
}

void panel::setPotential(Eigen::Vector3d Vinf)
{
    prevPotential = potential;
    
    potential = Vinf.dot(center) - doubletStrength; // Katz 11.74, 13.157
}

bool panel::inPanelProjection(const Eigen::Vector3d &POI, Eigen::Vector3d &projectedPnt)
{
    // Returns true if point is contained in extrusion of panel infinitely in normal direction
    Eigen::MatrixXd points(nodes.size()+1,3);
    std::vector<Eigen::Vector3d> nodesLocal;
    for (int i=0; i<nodes.size(); i++)
    {
        nodesLocal.push_back(global2local(nodes[i]->getPnt(), true));
        points.row(i) = nodesLocal[i];
    }
    points.row(nodes.size()) = global2local(POI,true);
    
    convexHull hull(points,true);
    
    if (hull.compareNodes(nodesLocal))
    {
        Eigen::Vector3d vec = POI-center;
        Eigen::Vector3d projVec = vec-(vec.dot(normal))*normal;
        projectedPnt = center + projVec;
//        projectedPnt = points.row(nodes.size());
//        projectedPnt(2) = 0; // Get point in panel plane.
//        projectedPnt = local2global(projectedPnt, true);
        return true;
    }
    
    projectedPnt = POI;
    return false;
}

bool panel::onPanelCheck(const Eigen::Vector3d &POI){
    // check to see if a point of interest lies on a panel. Will be using my completely arbitrary definition of 5% of the panel's longest side as the cutoff for 'in the panel'
    
    Eigen::Vector3d dummy = Eigen::Vector3d::Zero(); //isPanelProjection returns the projected point which is not needed for this application.
    Eigen::Vector3d pointInLocal = global2local(POI, true);
    
    // see if is in projection AND if point in local coords is within the distance of panel
    if(inPanelProjection(POI, dummy) && (pointInLocal.z() < 0.05*this->longSide))
    {
        return true;
    }
    return false;
    
}


Eigen::Matrix3d panel::getLocalSys()
{
    // Local Coordinate System
    // X : Points from center of panel to first vertex
    // Y : Normal crossed with X to obtain right hand coordinate system
    // Z : Normal to the panel
    Eigen::Matrix3d local = Eigen::Matrix3d::Zero();
    local.row(0) = getUnitVector(center,nodes[0]->getPnt());
    local.row(1) = normal.cross(local.row(0));
    local.row(2) = normal;
    
    return local;
}


Eigen::Vector3d panel::getUnitVector(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2)
{
    //Returns unit vector pointing from p1 to p2;
    Eigen::Vector3d unit;
    for (int i=0; i<3; i++)
    {
        unit(i) = p2(i)-p1(i);
    }
    unit.normalize();
    return unit;
}

Eigen::Vector3d panel::global2local(const Eigen::Vector3d &globalVec, bool translate)
{
    Eigen::Vector3d toTrans = globalVec;
    if (translate)
    {
        toTrans = globalVec-center;
    }
    return getLocalSys()*toTrans;
}

Eigen::Vector3d panel::local2global(const Eigen::Vector3d &localVec, bool translate)
{
    Eigen::Matrix3d transMat = getLocalSys();
    transMat.transposeInPlace();
    if (translate)
    {
        return transMat*localVec+center;
    }
    else
    {
        return transMat*localVec;
    }
    
}

double panel::dubPhiInf(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    
    if (pjk.norm() < 0.0000001)
    {
        return -0.5;
    }
    
    if (pjk.norm()/longSide > 5)
    {
        return pntDubPhi(PN,pjk.norm());
    }
    else
    {
        double phi = 0;
        double Al;
        Eigen::Vector3d a,b,s;
        for (int i=0; i<nodes.size(); i++)
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
            phi += vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
        }
        return phi/(4*M_PI);
    }
}

Eigen::Vector3d panel::dubVInf(const Eigen::Vector3d &POI)
{
    // VSAero doublet velocity influence formulation
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    //    if (pjk.norm() < 0.0000001)
    //    {
    //        vel << 0,0,0;
    //        return vel;
    //    }
    if (pjk.norm()/longSide > 5)
    {
        return pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s;
        int i1,i2;
        for (int i=0; i<nodes.size(); i++)
        {
            if (i!=nodes.size()-1)
            {
                i1 = i;
                i2 = i+1;
            }
            else
            {
                i1 = i;
                i2 = 0;
            }
            p1 = nodes[i1]->getPnt();
            p2 = nodes[i2]->getPnt();
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            
            vel += vortexV(a,b,s);
            
        }
        return vel/(4*M_PI);
    }
}


Eigen::Vector3d panel::vortexV(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s)
{
    
    return (a.cross(b)*(a.norm()+b.norm()))/(a.norm()*b.norm()*((a.norm()*b.norm())+a.dot(b))+pow(core*s.norm(),2)); // Connor: s is side length and was not included before. Excluding side length makes for an arbitrary core size for different geometry
}

double panel::vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    double eps = pow(10, -15);
    double PA,PB,num,denom;
    
    PA = a.dot(l.cross(a.cross(s)));
    PB = PA-Al*s.dot(m);
    num = s.dot(m)*PN*(b.norm()*PA-a.norm()*PB);
    denom = PA*PB+pow(PN,2)*a.norm()*b.norm()*pow(s.dot(m),2);
    if (denom == 0 && std::abs(PN) < eps)
    {
        // Point is on edge.
        if (PN >= 0)
        {
            return 0.5*M_PI;
        }
        else
        {
            return -0.5*M_PI;
        }
    }
    return atan2(num,denom);
}

double panel::pntDubPhi(const double &PN, const double &PJK)
{
    return PN*area/(4*M_PI*pow(PJK,3));
}

Eigen::Vector3d panel::pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk)
{
    return area*(3*pjk.dot(n)*pjk-pow(pjk.norm(),2)*n)/(4*M_PI*pow(pjk.norm(),5));
}

Eigen::VectorXi panel::getVerts()
{
    Eigen::VectorXi verts(nodes.size());
    for (int i=0; i<nodes.size(); i++)
    {
        verts(i) = nodes[i]->getIndex();
    }
    return verts;
}

std::vector<Eigen::Vector3d> panel::pntsAroundPnt(int nPnts,const Eigen::Vector3d &POI,double r)
{
    std::vector<Eigen::Vector3d> pnts;
    double theta;
    Eigen::Vector3d pnt;
    for (int i=0; i<nPnts; i++)
    {
        theta = (double)i/nPnts*(2*M_PI);
        pnt(0) = r*cos(theta);
        pnt(1) = r*sin(theta);
        pnt(2) = 0;
        pnt = local2global(pnt,true)+(POI-center);
        pnts.push_back(pnt);
    }
    
    return pnts;
}

Eigen::Vector3d panel::pntNearEdge(edge* e)
{
    Eigen::Vector3d pnt = center+0.95*(e->getMidPoint()-center);
    return pnt;
}


Eigen::Matrix3d panel::velocityGradientPointDoublet(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz

    double x, y, z, x0, y0, z0;
    x = this->getCenter().x(); y = this->getCenter().y(); z = this->getCenter().z();
    x0 = POI.x(); y0 = POI.y(); z0 = POI.z();
    
    double ddd = pow(pow(x-x0,2) + pow(y-y0,2) + z*z,3.5); // doublet deriv denom
    double uvConst = 3*doubletStrength*area/(4*M_PI);
    double wConst = -doubletStrength*area/(4*M_PI);
    
    velGradMat(0,0) = uvConst*z*(-4*pow(x-x0,2)+pow(y-y0,2)+z*z)/ddd;
    velGradMat(1,0) = uvConst*(-5)*z*(y-y0)*(x-x0)/ddd;
    velGradMat(2,0) = uvConst*(x-x0)*(pow(x-x0,2)+pow(y-y0,2)-4*z*z)/ddd;
    
    velGradMat(0,1) = uvConst*(-5)*(x-x0)*(y-y0)*z/ddd;
    velGradMat(1,1) = uvConst*z*(pow(x-x0,2) - 4*pow(y-y0,2) + z*z)/ddd;
    velGradMat(2,1) = uvConst*(y-y0)*(pow(x-x0,2) + pow(y-y0,2) - 4*z*z)/ddd;
    
    velGradMat(0,2) = wConst*(-3)*(x-x0)*(pow(x-x0,2) + pow(y-y0,2) - 4*z*z)/ddd;
    velGradMat(1,2) = wConst*(-3)*(y-y0)*(pow(x-x0,2) + pow(y-y0,2) - 4*z*z)/ddd;
    velGradMat(2,2) = wConst*(-3)*z*(3*pow(x-x0,2) + 3*pow(y-y0,2) - 2*z*z)/ddd;

    return velGradMat;
}

Eigen::Matrix3d panel::velocityGradientDoublet(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat = Eigen::Matrix3d::Zero();
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    Eigen::Vector3d p1,p2,a,b,s;
    int i1,i2;
    for (int i=0; i<nodes.size(); i++)
    {
        if (i!=nodes.size()-1)
        {
            i1 = i;
            i2 = i+1;
        }
        else
        {
            i1 = i;
            i2 = 0;
        }
        p1 = nodes[i1]->getPnt();
        p2 = nodes[i2]->getPnt();
        a = POI-p1;
        b = POI-p2;
        s = p2-p1;
        
        bool isOnPanel = onPanelCheck(POI);
        bool isNearFilament = nearFilamentCheck(p1,p2,POI);
        
        //If the POI is within the core distance of the filament (isNearFilament), the stretching influence is not accounted for
        if(isNearFilament == false && isOnPanel == false)
        {
            velGradMat += gradDoub(a,b,s);
        }
        else if(isNearFilament == false && isOnPanel == true)
        {
            // Only add the w gradint values, as u and v are zero on the panel. Katz
            velGradMat.col(2) += gradDoub(a,b,s).col(2);
            std::cout << "particle is on Panel" << std::endl; // Also, convet the Katz values to VSaero by dividing(or mult) by 4pi. Do for both source and doublet.
        }
    }
    
    return velGradMat*doubletStrength/(4*M_PI);
}

Eigen::Matrix3d panel::gradDoub(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s){
    Eigen::Matrix3d velGradMat = Eigen::Matrix3d::Zero();
    
    Eigen::Vector3d top = (a.norm() + b.norm())*(a.cross(b));
    double bot = a.norm()*b.norm()*(a.norm()*b.norm()+(a.dot(b)))+pow(core*s.norm(),2);
    
    Eigen::Vector3d dadx = {1,0,0};
    Eigen::Vector3d dbdx = {1,0,0};
    Eigen::Vector3d dady = {0,1,0};
    Eigen::Vector3d dbdy = {0,1,0};
    Eigen::Vector3d dadz = {0,0,1};
    Eigen::Vector3d dbdz = {0,0,1};
    
    double dandx = a.x()/a.norm();
    double dbndx = b.x()/b.norm();
    double dandy = a.y()/a.norm();
    double dbndy = b.y()/b.norm();
    double dandz = a.z()/a.norm();
    double dbndz = b.z()/b.norm();
    
    Eigen::Vector3d dtopdx = (a.norm()+b.norm())*(dadx.cross(b)+a.cross(dbdx)) + a.cross(b)*(dandx + dbndx);
    Eigen::Vector3d dtopdy = (a.norm()+b.norm())*(dady.cross(b)+a.cross(dbdy)) + a.cross(b)*(dandy + dbndy);
    Eigen::Vector3d dtopdz = (a.norm()+b.norm())*(dadz.cross(b)+a.cross(dbdz)) + a.cross(b)*(dandz + dbndz);
    
    double dbotdx = dandx*b.norm()*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*dbndx*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*b.norm()*(dandx*b.norm()+a.norm()*dbndx+(dadx).dot(b)+a.dot(dbdx));
    
    double dbotdy = dandy*b.norm()*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*dbndy*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*b.norm()*(dandy*b.norm()+a.norm()*dbndy+(dady).dot(b)+a.dot(dbdy));
    
    double dbotdz = dandz*b.norm()*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*dbndz*(a.norm()*b.norm()+a.dot(b))+
                    a.norm()*b.norm()*(dandz*b.norm()+a.norm()*dbndz+(dadz).dot(b)+a.dot(dbdz));
    
    velGradMat.row(0) = (bot*dtopdx - top*dbotdx)/(bot*bot);
    velGradMat.row(1) = (bot*dtopdy - top*dbotdy)/(bot*bot);
    velGradMat.row(2) = (bot*dtopdz - top*dbotdz)/(bot*bot);

    return velGradMat;
}

bool panel::nearFilamentCheck(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &POI){
    // Function checks if a point of interest is near a filament. First, the POI (point of interest) is compared with the filament endpoints to see if it is within the core radius value. Next, the cross product is taken between the two vectors going from each endpoint to the POI. If this cross product is within the core value times the edge length, then it is near.
    
    
    bool isNear = false;
    double edgeLength = (p2-p1).norm();
    double p1d = (p1-POI).norm();
    double p2d = (p2-POI).norm();
    
    // Check to see if POI is within the core radius of end points
    if(p1d < core*edgeLength || p2d < core*edgeLength)
    {
        isNear = true;
        return isNear;
    }
    
    // Check to see if point is also outside of the edges
    bool outsideEdges = false;
    if(p1d > edgeLength || p2d > edgeLength){
        outsideEdges = true;
    }
    
    // Check to see if point is within core distance height
    bool isWithinHeight = false;
    double dist = ((p1-POI).cross(p2-POI)).norm()/edgeLength;
    
    if(dist < core*edgeLength)
    {
        isWithinHeight = true;
    }
    
    if(isWithinHeight && !outsideEdges){
        isNear = true;
    }

    return isNear;
}


























































//Eigen::Matrix3d panel::velocityGradientTriDoublet(Eigen::Vector3d POI){
//    Eigen::Matrix3d velGradMat = Eigen::Matrix3d::Zero();
//    // dudx  dvdx  dwdx
//    // dudy  dvdy  dwdy
//    // dudz  dvdz  dwdz
//    
//    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
//    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
//    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
//    
//    Eigen::Vector3d n1 = global2local(n1global, true);
//    Eigen::Vector3d n2 = global2local(n2global, true);
//    Eigen::Vector3d n3 = global2local(n3global, true);
//    Eigen::Vector3d POIloc = global2local(POI, true);
//    
//    double x1, y1, x2, y2, x3, y3;
//    x1 = n1.x(); x2 = n2.x(); x3 = n3.x();
//    y1 = n1.y(); y2 = n2.y(); y3 = n3.y();
//    
//    double x, y, z;
//    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
//    
//    double derConst = doubletStrength/(4*M_PI);
//    
//    double d12,d23,d31;
//    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
//    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
//    d31 = pow((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3),0.5);
//    
//    double r1,r2,r3;
//    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    
//    double dr1dx, dr2dx, dr3dx;
//    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    
//    double dr1dy, dr2dy, dr3dy;
//    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    
//    double dr1dz, dr2dz, dr3dz;
//    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    
//    double m12,m23,m31;
//    m12 = (y2-y1)/(x2-x1);
//    m23 = (y3-y2)/(x3-x2);
//    m31 = (y1-y3)/(x1-x3);
//    
//    double e1,e2,e3;
//    e1 = (x-x1)*(x-x1) + z*z;
//    e2 = (x-x2)*(x-x2) + z*z;
//    e3 = (x-x3)*(x-x3) + z*z;
//    
//    double h1,h2,h3;
//    h1 = (x-x1)*(y-y1);
//    h2 = (x-x2)*(y-y2);
//    h3 = (x-x3)*(y-y3);
//    
//    double de1dx,de2dx,de3dx;
//    de1dx = 2*(x-x1);
//    de2dx = 2*(x-x2);
//    de3dx = 2*(x-x3);
//    
//    double dh1dx,dh2dx,dh3dx;
//    dh1dx = (y-y1);
//    dh2dx = (y-y2);
//    dh3dx = (y-y3);
//    
//    double dh1dy,dh2dy,dh3dy;
//    dh1dy = (x-x1);
//    dh2dy = (x-x2);
//    dh3dy = (x-x3);
//    
//    double de1dz,de2dz,de3dz;
//    de1dz = 2*z;
//    de2dz = 2*z;
//    de3dz = 2*z;
//    
//    double beta1,beta2,beta3;
//    beta1 = (r1*r2 - ((x-x1)*(x-x2) + (y-y1)*(y-y2) + z*z));
//    beta2 = (r2*r3 - ((x-x2)*(x-x3) + (y-y2)*(y-y3) + z*z));
//    beta3 = (r3*r1 - ((x-x3)*(x-x1) + (y-y3)*(y-y1) + z*z));
//    
//    double bot1,bot2,bot3;
//    bot1 = (r1*r2*(r1*r2 - ((x-x1)*(x-x2) + (y-y1)*(y-y2) + z*z)));
//    bot2 = (r2*r3*(r2*r3 - ((x-x2)*(x-x3) + (y-y2)*(y-y3) + z*z)));
//    bot3 = (r3*r1*(r3*r1 - ((x-x3)*(x-x1) + (y-y3)*(y-y1) + z*z)));
//    
//    double top1, top2, top3;
//    
//    // w: Will do w first because in the case that the particle winds up on the panel, then the u and v terms will be zero and the function call can end.
//    top1 = (((x-x2)*(y-y1) - (x-x1)*(y-y2))*(r1+r2));
//    top2 = (((x-x3)*(y-y2) - (x-x2)*(y-y3))*(r2+r3));
//    top3 = (((x-x1)*(y-y3) - (x-x3)*(y-y1))*(r3+r1));
//    
//    velGradMat(0,2)=
//    (bot1*(((y-y1)-(y-y2))*(r1+r2)+((x-x2)*(y-y1)-(x-x1)*(y-y2))*(dr1dx+dr2dx))-top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2+r1*dr2dx-(x-x2)-(x-x1))))/((bot1*bot1)) +
//    (bot2*(((y-y2)-(y-y3))*(r2+r3)+((x-x3)*(y-y2)-(x-x2)*(y-y3))*(dr2dx+dr3dx))-top2*(dr2dx*r3*beta2+r2*dr3dx*beta2+r2*r3*(dr2dx*r3+r2*dr3dx-(x-x3)-(x-x2))))/((bot2*bot2)) +
//    (bot3*(((y-y3)-(y-y1))*(r3+r1)+((x-x1)*(y-y3)-(x-x3)*(y-y1))*(dr3dx+dr1dx))-top3*(dr3dx*r1*beta3+r3*dr1dx*beta3+r3*r1*(dr3dx*r1+r3*dr1dx-(x-x1)-(x-x3))))/((bot3*bot3));
//    
//    velGradMat(1,2)=
//    (bot1*(((x-x2)-(x-x1))*(r1+r2)+((x-x2)*(y-y1)-(x-x1)*(y-y2))*(dr1dy+dr2dy))-top1*(dr1dy*r2*beta1+r1*dr2dy*beta1+r1*r2*(dr1dy*r2+r1*dr2dy-(y-y1)-(y-y2))))/((bot1*bot1)) +
//    (bot2*(((x-x3)-(x-x2))*(r2+r3)+((x-x3)*(y-y2)-(x-x2)*(y-y3))*(dr2dy+dr3dy))-top2*(dr2dy*r3*beta2+r2*dr3dy*beta2+r2*r3*(dr2dy*r3+r2*dr3dy-(y-y2)-(y-y3))))/((bot2*bot2)) +
//    (bot3*(((x-x1)-(x-x3))*(r3+r1)+((x-x1)*(y-y3)-(x-x3)*(y-y1))*(dr3dy+dr1dy))-top3*(dr3dy*r1*beta3+r3*dr1dy*beta3+r3*r1*(dr3dy*r1+r3*dr1dy-(y-y3)-(y-y1))))/((bot3*bot3));
//    
//    velGradMat(2,2)=(bot1*(((x-x2)*(y-y1)-(x-x1)*(y-y2))*(dr1dz+dr2dz)) - top1*(dr1dz*r2*beta1+r1*dr2dz*beta1 + r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/((bot1*bot1)) +
//    (bot2*(((x-x3)*(y-y2)-(x-x2)*(y-y3))*(dr2dz+dr3dz))-top2*(dr2dz*r3*beta2 + r2*dr3dz*beta2 + r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/((bot2*bot2)) +
//    (bot3*(((x-x1)*(y-y3)-(x-x3)*(y-y1))*(dr3dz+dr1dz))-top3*(dr3dz*r1*beta3 + r3*dr1dz*beta3 + r3*r1*(dr3dz*r1+r3*dr1dz-2*z)))/((bot3*bot3));
//    
//    
//    // u
//    top1 = (z*(y1-y2)*(r1+r2));
//    top2 = (z*(y2-y3)*(r2+r3));
//    top3 = (z*(y3-y1)*(r3+r1));
//
//    velGradMat(0,0)=(bot1*(z*(y1-y2)*(dr1dx+dr2dx))-top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2+r1*dr2dx-(x-x2)-(x-x1))))/(bot1*bot1) +
//                    (bot2*(z*(y2-y3)*(dr2dx+dr3dx))-top2*(dr2dx*r3*beta2+r2*dr3dx*beta2+r2*r3*(dr2dx*r3+r2*dr3dx-(x-x3)-(x-x2))))/(bot2*bot2) +
//                    (bot3*(z*(y3-y1)*(dr3dx+dr1dx))-top3*(dr3dx*r1*beta3+r3*dr1dx*beta3+r3*r1*(dr3dx*r1+r3*dr1dx-(x-x1)-(x-x3))))/(bot3*bot3);
//
//    velGradMat(1,0)=(bot1*(z*(y1-y2)*(dr1dy+dr2dy))-top1*(dr1dy*r2*beta1+r1*dr2dy*beta1+r1*r2*(dr1dy*r2+r1*dr2dy-(y-y2)-(y-y1))))/(bot1*bot1) +
//                    (bot2*(z*(y2-y3)*(dr2dy+dr3dy))-top2*(dr2dy*r3*beta2+r2*dr3dy*beta2+r2*r3*(dr2dy*r3+r2*dr3dy-(y-y3)-(y-y2))))/(bot2*bot2) +
//                    (bot3*(z*(y3-y1)*(dr3dy+dr1dy))-top3*(dr3dy*r1*beta3+r3*dr1dy*beta3+r3*r1*(dr3dy*r1+r3*dr1dy-(y-y1)-(y-y3))))/(bot3*bot3);
//
//    velGradMat(2,0)=(bot1*((y1-y2)*(z*(dr1dz+dr2dz)+(r1+r2)))-top1*(dr1dz*r2*beta1+r1*dr2dz*beta1+r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/(bot1*bot1) +
//                    (bot2*((y2-y3)*(z*(dr2dz+dr3dz)+(r2+r3)))-top2*(dr2dz*r3*beta2+r2*dr3dz*beta2+r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/(bot2*bot2) +
//                    (bot3*((y3-y1)*(z*(dr3dz+dr1dz)+(r3+r1)))-top3*(dr3dz*r1*beta3+r3*dr1dz*beta3+r3*r1*(dr3dz*r1+r3*dr1dz-2*z)))/(bot3*bot3);
//    
//    // v
//    top1 = z*(x2-x1)*(r1+r2);
//    top2 = z*(x3-x2)*(r2+r3);
//    top3 = z*(x1-x3)*(r3+r1);
//    
//    velGradMat(0,1)=(bot1*(z*(x2-x1)*(dr1dx+dr2dx))-top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2+r1*dr2dx-(x-x2)-(x-x1))))/(bot1*bot1) +
//                    (bot2*(z*(x3-x2)*(dr2dx+dr3dx))-top2*(dr2dx*r3*beta2+r2*dr3dx*beta2+r2*r3*(dr2dx*r3+r2*dr3dx-(x-x3)-(x-x2))))/(bot2*bot2) +
//                    (bot3*(z*(x1-x3)*(dr3dx+dr1dx))-top3*(dr3dx*r1*beta3+r3*dr1dx*beta3+r3*r1*(dr3dx*r1+r3*dr1dx-(x-x1)-(x-x3))))/(bot3*bot3);
//
//    velGradMat(1,1)=(bot1*(z*(x2-x1)*(dr1dy+dr2dy))-top1*(dr1dy*r2*beta1+r1*dr2dy*beta1+r1*r2*(dr1dy*r2+r1*dr2dy-(y-y2)-(y-y1))))/(bot1*bot1) +
//                    (bot2*(z*(x3-x2)*(dr2dy+dr3dy))-top2*(dr2dy*r3*beta2+r2*dr3dy*beta2+r2*r3*(dr2dy*r3+r2*dr3dy-(y-y3)-(y-y2))))/(bot2*bot2) +
//                    (bot3*(z*(x1-x3)*(dr3dy+dr1dy))-top3*(dr3dy*r1*beta3+r3*dr1dy*beta3+r3*r1*(dr3dy*r1+r3*dr1dy-(y-y1)-(y-y3))))/(bot3*bot3);
//    
//    velGradMat(2,1)=(bot1*((x2-x1)*(z*(dr1dz+dr2dz)+(r1+r2)))-top1*(dr1dz*r2*beta1+r1*dr2dz*beta1+r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/(bot1*bot1) +
//                    (bot2*((x3-x2)*(z*(dr2dz+dr3dz)+(r2+r3)))-top2*(dr2dz*r3*beta2+r2*dr3dz*beta2+r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/(bot2*bot2) +
//                    (bot3*((x1-x3)*(z*(dr3dz+dr1dz)+(r3+r1)))-top3*(dr3dz*r1*beta3+r3*dr1dz*beta3+r3*r1*(dr3dz*r1+r3*dr1dz-2*z)))/(bot3*bot3);
//    
//
//    velGradMat *= derConst;
//    
//    std::cout << velGradMat << "\n\n" << std::endl;
//
//    return velGradMat;
//    
//}
//
//Eigen::Matrix3d panel::velocityGradientQuadDoublet(Eigen::Vector3d POI){
//    Eigen::Matrix3d velGradMat;
//    // dudx  dvdx  dwdx
//    // dudy  dvdy  dwdy
//    // dudz  dvdz  dwdz
//    
//    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
//    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
//    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
//    Eigen::Vector3d n4global = this->getNodes()[3]->getPnt();
//    
//    Eigen::Vector3d n1 = global2local(n1global, true);
//    Eigen::Vector3d n2 = global2local(n2global, true);
//    Eigen::Vector3d n3 = global2local(n3global, true);
//    Eigen::Vector3d n4 = global2local(n4global, true);
//    Eigen::Vector3d POIloc = global2local(POI, true);
//    
//    bool onPanel = inPanelProjection(<#const Eigen::Vector3d &POI#>, <#Eigen::Vector3d &projectedPnt#>)
//    double x1, y1, x2, y2, x3, y3, x4, y4;
//    x1 = n1.x(); x2 = n2.x(); x3 = n3.x(), x4 = n4.x();
//    y1 = n1.y(); y2 = n2.y(); y3 = n3.y(), y4 = n4.y();
//    
//    double x, y, z;
//    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
//    
//    double derConst = doubletStrength/(4*M_PI);
//    
//    double d12,d23,d34,d41;
//    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
//    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
//    d34 = pow((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3),0.5);
//    d41 = pow((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4),0.5);
//    
//    double r1,r2,r3,r4;
//    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    r4 = pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
//    
//    double dr1dx, dr2dx, dr3dx, dr4dx;
//    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    dr4dx = (x-x4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
//    
//    double dr1dy, dr2dy, dr3dy, dr4dy;
//    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    dr4dy = (y-y4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
//    
//    double dr1dz, dr2dz, dr3dz, dr4dz;
//    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
//    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
//    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
//    dr4dz = (z)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
//    
//    double m12,m23,m34,m41;
//    m12 = (y2-y1)/(x2-x1);
//    m23 = (y3-y2)/(x3-x2);
//    m34 = (y4-y3)/(x4-x3);
//    m41 = (y1-y4)/(x1-x4);
//    
//    double e1,e2,e3,e4;
//    e1 = (x-x1)*(x-x1) + z*z;
//    e2 = (x-x2)*(x-x2) + z*z;
//    e3 = (x-x3)*(x-x3) + z*z;
//    e4 = (x-x4)*(x-x4) + z*z;
//    
//    double h1,h2,h3,h4;
//    h1 = (x-x1)*(y-y1);
//    h2 = (x-x2)*(y-y2);
//    h3 = (x-x3)*(y-y3);
//    h4 = (x-x4)*(y-y4);
//    
//    double de1dx,de2dx,de3dx,de4dx;
//    de1dx = 2*(x-x1);
//    de2dx = 2*(x-x2);
//    de3dx = 2*(x-x3);
//    de4dx = 2*(x-x4);
//    
//    double dh1dx,dh2dx,dh3dx,dh4dx;
//    dh1dx = (y-y1);
//    dh2dx = (y-y2);
//    dh3dx = (y-y3);
//    dh4dx = (y-y4);
//    
//    double dh1dy,dh2dy,dh3dy,dh4dy;
//    dh1dy = (x-x1);
//    dh2dy = (x-x2);
//    dh3dy = (x-x3);
//    dh4dy = (x-x4);
//    
//    double de1dz,de2dz,de3dz,de4dz;
//    de1dz = 2*z;
//    de2dz = 2*z;
//    de3dz = 2*z;
//    de4dz = 2*z;
//    
//    double beta1,beta2,beta3,beta4;
//    beta1 = (r1*r2 - ((x-x1)*(x-x2) + (y-y1)*(y-y2) + z*z));
//    beta2 = (r2*r3 - ((x-x2)*(x-x3) + (y-y2)*(y-y3) + z*z));
//    beta3 = (r3*r4 - ((x-x3)*(x-x4) + (y-y3)*(y-y4) + z*z));
//    beta4 = (r4*r1 - ((x-x4)*(x-x1) + (y-y4)*(y-y1) + z*z));
//    
//    double bot1,bot2,bot3,bot4;
//    bot1 = (r1*r2*(r1*r2 - ((x-x1)*(x-x2) + (y-y1)*(y-y2) + z*z)));
//    bot2 = (r2*r3*(r2*r3 - ((x-x2)*(x-x3) + (y-y2)*(y-y3) + z*z)));
//    bot3 = (r3*r4*(r3*r4 - ((x-x3)*(x-x4) + (y-y3)*(y-y4) + z*z)));
//    bot4 = (r4*r1*(r4*r1 - ((x-x4)*(x-x1) + (y-y4)*(y-y1) + z*z)));
//    
//    double top1, top2, top3, top4;
//    
//    // u
//    top1 = (z*(y1-y2)*(r1+r2));
//    top2 = (z*(y2-y3)*(r2+r3));
//    top3 = (z*(y3-y4)*(r3+r4));
//    top4 = (z*(y4-y1)*(r4+r1));
//    
//    velGradMat(0,0)=(bot1*(z*(y1-y2)*(dr1dx+dr2dx))-top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2+r1*dr2dx-(x-x2)-(x-x1))))/(bot1*bot1) +
//    (bot2*(z*(y2-y3)*(dr2dx+dr3dx))-top2*(dr2dx*r3*beta2+r2*dr3dx*beta2+r2*r3*(dr2dx*r3+r2*dr3dx-(x-x3)-(x-x2))))/(bot2*bot2) +
//    (bot3*(z*(y3-y4)*(dr3dx+dr4dx))-top3*(dr3dx*r4*beta3+r3*dr4dx*beta3+r3*r4*(dr3dx*r4+r3*dr4dx-(x-x4)-(x-x3))))/(bot3*bot3) +
//    (bot4*(z*(y4-y1)*(dr4dx+dr1dx))-top4*(dr4dx*r1*beta4+r4*dr1dx*beta4+r4*r1*(dr4dx*r1+r4*dr1dx-(x-x1)-(x-x4))))/(bot4*bot4);
//    
//    velGradMat(1,0)=(bot1*(z*(y1-y2)*(dr1dy+dr2dy))-top1*(dr1dy*r2*beta1+r1*dr2dy*beta1+r1*r2*(dr1dy*r2+r1*dr2dy-(y-y2)-(y-y1))))/(bot1*bot1) +
//    (bot2*(z*(y2-y3)*(dr2dy+dr3dy))-top2*(dr2dy*r3*beta2+r2*dr3dy*beta2+r2*r3*(dr2dy*r3+r2*dr3dy-(y-y3)-(y-y2))))/(bot2*bot2) +
//    (bot3*(z*(y3-y4)*(dr3dy+dr4dy))-top3*(dr3dy*r4*beta3+r3*dr4dy*beta3+r3*r4*(dr3dy*r4+r3*dr4dy-(y-y4)-(y-y3))))/(bot3*bot3) +
//    (bot4*(z*(y4-y1)*(dr4dy+dr1dy))-top4*(dr4dy*r1*beta4+r4*dr1dy*beta4+r4*r1*(dr4dy*r1+r4*dr1dy-(y-y1)-(y-y4))))/(bot4*bot4);
//    
//    velGradMat(2,0)=(bot1*((y1-y2)*(z*(dr1dz+dr2dz)+(r1+r2)))-top1*(dr1dz*r2*beta1+r1*dr2dz*beta1+r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/(bot1*bot1) +
//    (bot2*((y2-y3)*(z*(dr2dz+dr3dz)+(r2+r3)))-top2*(dr2dz*r3*beta2+r2*dr3dz*beta2+r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/(bot2*bot2) +
//    (bot3*((y3-y4)*(z*(dr3dz+dr4dz)+(r3+r4)))-top3*(dr3dz*r4*beta3+r3*dr4dz*beta3+r3*r4*(dr3dz*r4+r3*dr4dz-2*z)))/(bot3*bot3) +
//    (bot4*((y4-y1)*(z*(dr4dz+dr1dz)+(r4+r1)))-top4*(dr4dz*r1*beta4+r4*dr1dz*beta4+r4*r1*(dr4dz*r1+r4*dr1dz-2*z)))/(bot4*bot4);
//    
//    // v
//    top1 = z*(x2-x1)*(r1+r2);
//    top2 = z*(x3-x2)*(r2+r3);
//    top3 = z*(x4-x3)*(r3+r4);
//    top4 = z*(x1-x4)*(r4+r1);
//    
//    velGradMat(0,1)=(bot1*(z*(x2-x1)*(dr1dx+dr2dx))-top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2+r1*dr2dx-(x-x2)-(x-x1))))/(bot1*bot1) +
//    (bot2*(z*(x3-x2)*(dr2dx+dr3dx))-top2*(dr2dx*r3*beta2+r2*dr3dx*beta2+r2*r3*(dr2dx*r3+r2*dr3dx-(x-x3)-(x-x2))))/(bot2*bot2) +
//    (bot3*(z*(x4-x3)*(dr3dx+dr4dx))-top3*(dr3dx*r4*beta3+r3*dr4dx*beta3+r3*r4*(dr3dx*r4+r3*dr4dx-(x-x4)-(x-x3))))/(bot3*bot3) +
//    (bot4*(z*(x1-x4)*(dr4dx+dr1dx))-top4*(dr4dx*r1*beta4+r4*dr1dx*beta4+r4*r1*(dr4dx*r1+r4*dr1dx-(x-x1)-(x-x4))))/(bot4*bot4);
//    
//    velGradMat(1,1)=(bot1*(z*(x2-x1)*(dr1dy+dr2dy))-top1*(dr1dy*r2*beta1+r1*dr2dy*beta1+r1*r2*(dr1dy*r2+r1*dr2dy-(y-y2)-(y-y1))))/(bot1*bot1) +
//    (bot2*(z*(x3-x2)*(dr2dy+dr3dy))-top2*(dr2dy*r3*beta2+r2*dr3dy*beta2+r2*r3*(dr2dy*r3+r2*dr3dy-(y-y3)-(y-y2))))/(bot2*bot2) +
//    (bot3*(z*(x4-x3)*(dr3dy+dr4dy))-top3*(dr3dy*r4*beta3+r3*dr4dy*beta3+r3*r4*(dr3dy*r4+r3*dr4dy-(y-y4)-(y-y3))))/(bot3*bot3) +
//    (bot4*(z*(x1-x4)*(dr4dy+dr1dy))-top4*(dr4dy*r1*beta4+r4*dr1dy*beta4+r4*r1*(dr4dy*r1+r4*dr1dy-(y-y1)-(y-y4))))/(bot4*bot4);
//    
//    velGradMat(2,1)=(bot1*((x2-x1)*(z*(dr1dz+dr2dz)+(r1+r2)))-top1*(dr1dz*r2*beta1+r1*dr2dz*beta1+r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/(bot1*bot1) +
//    (bot2*((x3-x2)*(z*(dr2dz+dr3dz)+(r2+r3)))-top2*(dr2dz*r3*beta2+r2*dr3dz*beta2+r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/(bot2*bot2) +
//    (bot3*((x4-x3)*(z*(dr3dz+dr4dz)+(r3+r4)))-top3*(dr3dz*r4*beta3+r3*dr4dz*beta3+r3*r4*(dr3dz*r4+r3*dr4dz-2*z)))/(bot3*bot3) +
//    (bot4*((x1-x4)*(z*(dr4dz+dr1dz)+(r4+r1)))-top4*(dr4dz*r1*beta4+r4*dr1dz*beta4+r4*r1*(dr4dz*r1+r4*dr1dz-2*z)))/(bot4*bot4);
//    
//    // w
//    top1 = (((x-x2)*(y-y1) - (x-x1)*(y-y2))*(r1+r2));
//    top2 = (((x-x3)*(y-y2) - (x-x2)*(y-y3))*(r2+r3));
//    top3 = (((x-x4)*(y-y3) - (x-x3)*(y-y4))*(r3+r4));
//    top4 = (((x-x1)*(y-y4) - (x-x4)*(y-y1))*(r4+r1));
//    
//    velGradMat(0,2)=(bot1*(((y-y1)-(y-y2))*(r1+r2) + ((x-x2)*(y-y1) - (x-x1)*(y-y2))*(dr1dx + dr2dx)) - top1*(dr1dx*r2*beta1+r1*dr2dx*beta1+r1*r2*(dr1dx*r2 + r1*dr2dx - (x-x2) - (x-x1))))/((bot1*bot1)) +
//    (bot2*(((y-y2)-(y-y3))*(r2+r3) + ((x-x3)*(y-y2) - (x-x2)*(y-y3))*(dr2dx + dr3dx)) - top2*(dr2dx*r3*beta2+r2*dr3dx*beta2 + r2*r3*(dr2dx*r3 + r2*dr3dx - (x-x3) - (x-x2))))/((bot2*bot2)) +
//    (bot3*(((y-y3)-(y-y4))*(r3+r4) + ((x-x4)*(y-y3) - (x-x3)*(y-y4))*(dr3dx + dr4dx)) - top3*(dr3dx*r4*beta3+r3*dr4dx*beta3 + r3*r4*(dr3dx*r4 + r3*dr4dx - (x-x4) - (x-x3))))/((bot3*bot3)) +
//    (bot4*(((y-y4)-(y-y1))*(r4+r1) + ((x-x1)*(y-y4) - (x-x4)*(y-y1))*(dr4dx + dr1dx)) - top4*(dr4dx*r1*beta4+r4*dr1dx*beta4 + r4*r1*(dr4dx*r1 + r4*dr1dx - (x-x1) - (x-x4))))/((bot4*bot4));
//    
//    velGradMat(1,2)=(bot1*(((x-x2)-(x-x1))*(r1+r2)+((x-x2)*(y-y1)-(x-x1)*(y-y2))*(dr1dy + dr2dy))-top1*(dr1dy*r2*beta1 + r1*dr2dy*beta1 + r1*r2*(dr1dy*r2 + r1*dr2dy - (y-y1) - (y-y2))))/((bot1*bot1)) +
//    (bot2*(((x-x3)-(x-x2))*(r2+r3)+((x-x3)*(y-y2)-(x-x2)*(y-y3))*(dr2dy + dr3dy))-top2*(dr2dy*r3*beta2 + r2*dr3dy*beta2 + r2*r3*(dr2dy*r3 + r2*dr3dy - (y-y2) - (y-y3))))/((bot2*bot2)) +
//    (bot3*(((x-x4)-(x-x3))*(r3+r4)+((x-x4)*(y-y3)-(x-x3)*(y-y4))*(dr3dy + dr4dy))-top3*(dr3dy*r4*beta3 + r3*dr4dy*beta3 + r3*r4*(dr3dy*r4 + r3*dr4dy - (y-y3) - (y-y4))))/((bot3*bot3)) +
//    (bot4*(((x-x1)-(x-x4))*(r4+r1)+((x-x1)*(y-y4)-(x-x4)*(y-y1))*(dr4dy + dr1dy))-top4*(dr4dy*r1*beta4 + r4*dr1dy*beta4 + r4*r1*(dr4dy*r1 + r4*dr1dy - (y-y4) - (y-y1))))/((bot4*bot4));
//    
//    velGradMat(2,2)=(bot1*(((x-x2)*(y-y1)-(x-x1)*(y-y2))*(dr1dz+dr2dz)) - top1*(dr1dz*r2*beta1 + r1*dr2dz*beta1 + r1*r2*(dr1dz*r2+r1*dr2dz-2*z)))/((bot1*bot1)) +
//    (bot2*(((x-x3)*(y-y2)-(x-x2)*(y-y3))*(dr2dz+dr3dz))-top2*(dr2dz*r3*beta2 + r2*dr3dz*beta2 + r2*r3*(dr2dz*r3+r2*dr3dz-2*z)))/((bot2*bot2)) +
//    (bot3*(((x-x4)*(y-y3)-(x-x3)*(y-y4))*(dr3dz+dr4dz))-top3*(dr3dz*r4*beta3 + r3*dr4dz*beta3 + r3*r4*(dr3dz*r4+r3*dr4dz-2*z)))/((bot3*bot3)) +
//    (bot4*(((x-x1)*(y-y4)-(x-x4)*(y-y1))*(dr4dz+dr1dz))-top4*(dr4dz*r1*beta4 + r4*dr1dz*beta4 + r4*r1*(dr4dz*r1+r4*dr1dz-2*z)))/((bot4*bot4));
//    
//    velGradMat *= derConst;
//    
//    return velGradMat;
//}
//
