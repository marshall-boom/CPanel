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

#include "panel.h"
#include "cpNode.h"
#include "edge.h"
#include "surface.h"

#define _USE_MATH_DEFINES
#include <math.h>

panel::panel(std::vector<cpNode*> nnodes, std::vector<edge*> ppEdges, Eigen::Vector3d bbezNorm, size_t surfID)
  : nodes(nnodes), pEdges(ppEdges), bezNormal(bbezNorm), ID(surfID)
{
    setGeom();
}


void panel::setGeom()
{
    longSide = 0;
    
    for (edges_index_type i=0; i<pEdges.size(); i++)
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

        /*normal = a.cross(b);*/
		normal = b.cross(a);
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
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        nodesLocal.push_back(global2local(nodes[i]->getPnt(), true));
        points.row(static_cast<Eigen::Vector3d::Index>(i)) = nodesLocal[i];
    }
    points.row(static_cast<Eigen::Vector3d::Index>(nodes.size())) = global2local(POI,true);
    
    convexHull hull(points,true);
    
    if (hull.compareNodes(nodesLocal))
    {
        Eigen::Vector3d vec = POI-center;
        Eigen::Vector3d projVec = vec-(vec.dot(normal))*normal;
        projectedPnt = center + projVec;

        return true;
    }
    
    projectedPnt = POI;
    return false;
}

bool panel::onPanelCheck(const Eigen::Vector3d &POI){
    // Check to see if a point of interest lies on a panel. Will be using my completely arbitrary definition of 5% of the panel's longest side as the cutoff for 'in the panel'
    
    Eigen::Vector3d dummy = Eigen::Vector3d::Zero(); // Projected point which is not needed for this application.
    Eigen::Vector3d pointInLocal = global2local(POI, true);
    
    // See if is in projection AND if point in local coords is within the distance of panel
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
// NOTE: last paremeter is not used
//            phi += vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
            phi += vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1));
        }
        return phi/(4*M_PI);
    }
}

// Linear doublet influence
void panel::linDubPhiInf(const Eigen::Vector3d &POI, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow)
{
	// Define and get stuff
	Eigen::Vector3d pjk = POI - center;
	Eigen::Matrix3d local = getLocalSys();
	Eigen::Vector3d POIloc, vertsPhi;
	Eigen::Matrix<size_t, Eigen::Dynamic, 1> verts = getVerts();					//////// MIGHT NEED TO CHANGE TO getBodyVerts()
	double PN = pjk.dot(local.row(2)); // normal dist from panel to field point
	POIloc = global2local(POI, true);

	///////// FAR FIELD NOT YET WORKING, EXCLUDED FROM CALCULATION ///////////

	// Check far field condition
	//if (pjk.norm() / longSide > 5)
	//{
		//std::cout << "\tFar Field" << std::endl;
		//vertsPhi = linPntDubPhi(PN, pjk.norm(), POIloc);
		//std::cout << "Far: " << vertsPhiFar << "\n" << std::endl;
	//}
	//else
	//{
		// Define stuff
		double Al;
		Eigen::Vector3d a, b, s;
		Eigen::Vector3d Iints, Imat;
		Eigen::VectorXd Hints;
		Hints.resize(Eigen::Index(3));
		Hints.setZero();
		Eigen::Matrix3d vertsMat;
		double l1, l2, c1, c2, s1, s2, g, nuXi, nuEta, myAl;

		// Iterate through panel edges and get H integrals
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

			nuXi = (p2.y() - p1.y()) / s.norm();
			nuEta = (p2.x() - p1.x()) / s.norm();
			myAl = (POI.x() - p1.x())*nuXi - (POI.y() - p1.y())*nuEta;
			g = sqrt(pow(myAl, 2) + pow(PN, 2));
			l1 = (p1.x() - POI.x())*nuEta + (p1.y() - POI.y())*nuXi;
			l2 = (p2.x() - POI.x())*nuEta + (p2.y() - POI.y())*nuXi;
			double m = (p2.y() - p1.y()) / (p2.x() - p1.x());
			double mt = -nuXi / nuEta;

			s1 = sqrt(pow(l1,2) + pow(g,2));
			s2 = sqrt(pow(l2, 2) + pow(g, 2));
			c1 = pow(g, 2) + abs(PN)*s1;
			c2 = pow(g, 2) + abs(PN)*s2;

			linPhiHintegrals(Hints,g,myAl,Al,l1,l2,c1,c2,nuEta,nuXi,PN,a,b,s, local.row(0), local.row(1));
		}

		// Get final H integrals and compute I integrals

		Hints[0] = Hints[0] / PN;
		Hints[1] = -Hints[1];
		Hints[2] = -Hints[2];				// Hints = [H113 H213 H123]
		Iints = PN * Hints / (4.0 * M_PI);	// Iints = [I11 I21 I12]
		
		// Build I matrix
		Imat[0] = Iints[0];
		Imat[1] = POIloc.x()*Iints[0] + Iints[1];
		Imat[2] = POIloc.y()*Iints[0] + Iints[2];
		
		// Convert linear doublet equation to vertex based linear doublet equations
		vertsMat = linVertsMatrix(true);

		// Compute influence of each vertex on the field point
		vertsPhi = Imat.transpose() * vertsMat.inverse();

	//}
	
	// Fill A matrix
	// Need to make this more efficient
	for (size_t i = 0; i < nodes.size(); i++)
	{
		// Check if influencing node is the same as the influenced node
		if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset())
		{
			Arow[nodes[i]->getIndex()] = -0.5;
		}
		else
		{
			Arow[nodes[i]->getIndex()] += vertsPhi(i);
		}
	}
}

// Linear doublet influence
Eigen::Matrix3d panel::linDubVInf(const Eigen::Vector3d &POI)
{
	// Define and get stuff
	Eigen::Vector3d pjk = POI - center;
	Eigen::Matrix3d local = getLocalSys();
	Eigen::Vector3d POIloc, vertsPhi;
	Eigen::Matrix<size_t, Eigen::Dynamic, 1> verts = getVerts();					//////// MIGHT NEED TO CHANGE TO getBodyVerts()
	double PN = pjk.dot(local.row(2)); // normal dist from panel to field point
	POIloc = global2local(POI, true);

	Eigen::Matrix3d Jints;
	Jints.setZero();

	if (pjk.norm() < 0.0000001)
	{
		return Jints;
	}
	else
	{
		///////// FAR FIELD NOT YET WORKING, EXCLUDED FROM CALCULATION ///////////

		// Check far field condition
		//if (pjk.norm() / longSide > 5)
		//{
		//std::cout << "\tFar Field" << std::endl;
		//vertsPhi = linPntDubPhi(PN, pjk.norm(), POIloc);
		//std::cout << "Far: " << vertsPhiFar << "\n" << std::endl;
		//}
		//else
		//{
		// Define stuff
		double Al;
		Eigen::Vector3d a, b, s;
		//Eigen::Vector3d Iints, Imat;
		Eigen::VectorXd Hints;
		Hints.resize(Eigen::Index(9));
		Hints.setZero();
		Eigen::Matrix3d vertsMat, Jints;
		Jints.setZero();
		double l1, l2, c1, c2, g, nuXi, nuEta;

		double F123 = 0;

		double s1, s2, Etest;

		// Iterate through panel edges and get H integrals
		for (nodes_index_type i = 0; i<nodes.size(); i++)
		{
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			Eigen::Vector3d Eints;
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

			nuXi = s.dot(local.row(1)) / s.norm();
			nuEta = s.dot(local.row(0)) / s.norm();
			l1 = -a.dot(local.row(1))*nuXi - a.dot(local.row(0))*nuEta;
			l2 = -b.dot(local.row(1))*nuXi - b.dot(local.row(0))*nuEta;
			g = sqrt(pow(Al, 2) + pow(PN, 2));

			c1 = a.norm();
			c2 = b.norm();

			/*s1 = sqrt(pow(l1, 2) + pow(g, 2));
			s2 = sqrt(pow(l2, 2) + pow(g, 2));
			c1 = pow(g, 2) + abs(PN)*s1;
			c2 = pow(g, 2) + abs(PN)*s2;*/

			Eints[0] = 1.0 / c2 - 1.0 / c1;		// E111
			Etest = -a.dot(local.row(0))/c2 + b.dot(local.row(0))/c1;
			Eints[1] = (p2[0] - POI[0]) / c2 - (p1[0] - POI[0]) / c1;
			Eints[2] = (p2[1] - POI[1]) / c2 - (p1[1] - POI[1]) / c1;

			//linPhiHintegrals(Hints, g, Al, l1, l2, c1, c2, nuEta, nuXi, PN, a, b, s, local.row(0), local.row(1));
			linVelHintegrals(Hints, F123, Eints, g, Al, nuEta, nuXi);
		}

		// Get final H integrals
		Hints[0] = Hints[0] / PN;	// H113
		Hints[1] = -Hints[1];		// H213
		Hints[2] = -Hints[2];		// H123

		Hints[3] = (1.0 / 3.0 / pow(PN, 2)) * (Hints[0] + Hints[3]);	// H115
		Hints[4] = (1.0 / 3.0) * -Hints[4];								// H215
		Hints[5] = (1.0 / 3.0) * -Hints[5];								// H125

		Hints[6] = (1.0 / 3.0) * -Hints[6];						// H225
		Hints[7] = (1.0 / 3.0) * (-Hints[0] - F123);			// H135
		Hints[8] = -Hints[7] - pow(PN, 2)*Hints[3] + Hints[0];	// H315

		// Get J integrals
		Jints = linVelJintegrals(Hints, PN);

		//Iints = PN * Hints / (4.0 * M_PI);	// Iints = [I11 I21 I12]

		//									// Build I matrix
		//Imat[0] = Iints[0];
		//Imat[1] = POIloc.x()*Iints[0] + Iints[1];
		//Imat[2] = POIloc.y()*Iints[0] + Iints[2];

		//// Convert linear doublet equation to vertex based linear doublet equations
		//vertsMat = linVertsMatrix();

		//// Compute influence of each vertex on the field point
		//vertsPhi = Imat.transpose() * vertsMat.inverse();

		return Jints;

		//}
	}
}

void panel::setPanelVel(Eigen::Vector3d Vel)
{
	velocity = Vel;
}

void panel::linPhiHintegrals(Eigen::VectorXd &Hints, const double g, const double myAl, const double Al, const double l1, const double l2, const double c1, const double c2, const double nuEta, const double nuXi, const double &PN, const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l, const Eigen::Vector3d &m)
{
	// Define stuff
	double F111, num, denom;	// num and denom not used anymore, along with c1 and c2
	double Htest = 0;

	double Htan;

	// Compute F111 integral
	if (l1 >= 0 && l2 >= 0)
	{
		F111 = log((sqrt(pow(l2, 2) + pow(g, 2)) + l2) / (sqrt(pow(l1, 2) + pow(g, 2)) + l1));
	}
	else if (l1 < 0 && l2 < 0)
	{
		F111 = log((sqrt(pow(l1, 2) + pow(g, 2)) - l1) / (sqrt(pow(l2, 2) + pow(g, 2)) - l2));
	}
	else if (l1 < 0 && l2 >= 0)
	{
		F111 = log(((sqrt(pow(l1, 2) + pow(g, 2)) - l1) * (sqrt(pow(l2, 2) + pow(g, 2)) + l2)) / pow(g,2));
	}

	/*num = myAl * (l2*c1 - l1 * c2);
	denom = c1 * c2 + pow(myAl, 2)*l1*l2;
	Htan = atan2(num, denom);*/

	//Hints[0] += atan2(num, denom);

	Hints[0] += vortexPhi(PN, Al, a, b, s, l, m);	// H113 integral
	Hints[1] += F111 * nuXi;	// H213 integral
	Hints[2] += F111 * nuEta;	// H123 integral
}


void panel::linVelHintegrals(Eigen::VectorXd &Hints, double F123, Eigen::Vector3d &Eints, const double g, const double Al, const double nuEta, const double nuXi)
{
	// Define stuff
	double E111, E211, E121, F113, F123_temp;

	E111 = Eints[0];
	E211 = Eints[1];
	E121 = Eints[2];

	// Compute F integrals   check sign of nuEta if not working
	F113 = (1.0 / pow(g, 2)) * (-nuEta*E211 + nuXi*E121);
	F123_temp = nuEta * Al*F113 - nuXi * E111;
	F123 += nuEta * F123_temp;

	// Compute H integrals
	Hints[3] += Al * F113;			// H115
	Hints[4] += nuXi * F113;		// H215
	Hints[5] += nuEta * F113;		// H125
	Hints[6] += nuXi * F123_temp;	// H225
}


Eigen::Matrix3d panel::linVelJintegrals(Eigen::VectorXd &Hints, const double &PN)
{
	// Build matrix of J integrals
	//		[ J11x J21x J12x;
	//		  J11y J21y J12y;
	//		  J11z J21z J12z ]

	Eigen::Matrix3d Jints;

	Jints(0, 0) = 3.0 * PN * Hints[4];							// J11x
	Jints(0, 1) = -(1.0 / PN) * Hints[0] + 3.0 * PN*Hints[8];	// J21x
	Jints(0, 2) = 3.0 * PN * Hints[6];							// J12x

	Jints(1, 0) = 3.0 * PN * Hints[5];		// J11y
	Jints(1, 1) = 3.0 * PN * Hints[6];		// J21y
	Jints(1, 2) = 3.0 * PN * Hints[7];		// J12y

	Jints(2, 0) = Hints[0] - 3.0 * pow(PN, 2) * Hints[3];	// J11z
	Jints(2, 1) = -Hints[1] - 3.0 * pow(PN, 2) * Hints[4];	// J21z
	Jints(2, 2) = -Hints[2] - 3.0 * pow(PN, 2) * Hints[5];	// J12z

	Jints = (1.0 / 4.0 / M_PI) * Jints;

	return Jints;
}


Eigen::Matrix3d panel::linVertsMatrix(bool translate)
{
	Eigen::Vector3d vertLoc;
	Eigen::Matrix3d vertsMat;

	for (nodes_index_type i = 0; i < nodes.size(); i++)
	{
		vertLoc = global2local(nodes[i]->getPnt(), translate);
		vertsMat(i, 0) = 1.0;
		vertsMat(i, 1) = vertLoc.x();
		vertsMat(i, 2) = vertLoc.y();
	}
	return vertsMat;
}


Eigen::Matrix3d panel::supVertsMatrix(std::vector<Eigen::Vector3d> &supPnts)
{
	//Eigen::Vector3d vertLoc;
	Eigen::Matrix3d vertsMat;

	for (size_t i = 0; i < supPnts.size(); i++)
	{
		//vertLoc = global2local(nodes[i]->getPnt(), translate);
		vertsMat(i, 0) = 1.0;
		vertsMat(i, 1) = supPnts[i].x();
		vertsMat(i, 2) = supPnts[i].y();
	}
	return vertsMat;
}


void panel::linComputeVelocity(double PG,Eigen::Vector3d &Vinf)
{
	double mu_x, mu_y;
	Eigen::Vector3d vertDubStrengths, linDubConsts, panVel;
	Eigen::Matrix3d vertsMat = linVertsMatrix(false);
	vertDubStrengths = linGetDubStrengths();
	
	linDubConsts = vertsMat.inverse() * vertDubStrengths;
	//mu_0 = linDubConsts[0];
	mu_x = linDubConsts[1];
	mu_y = linDubConsts[2];

	panVel[0] = -mu_x;
	panVel[1] = -mu_y;
	panVel[2] = 0.0;
	
	velocity = local2global(panVel, false);
	velocity(0) /= PG;
}


//void panel::linGetConstDubStrength()
//{
//	double mu_0; //mu_x, mu_y;
//	Eigen::Vector3d vertDubStrengths, linDubConsts;
//	Eigen::Matrix3d vertsMat = linVertsMatrix(true);
//	vertDubStrengths = linGetDubStrengths();
//
//	linDubConsts = vertsMat.inverse() * vertDubStrengths;
//	mu_0 = linDubConsts(0);
//	/*mu_x = linDubConsts(1);
//	mu_y = linDubConsts(2);*/
//
//	doubletStrength = mu_0;
//	//doubletStrength = linDubConsts[0];
//}

//
//Eigen::Vector3d panel::supGetDubDiffs()
//{
//	Eigen::Vector3d vertDubDiffs;	// [mu1 mu2 mu3]
//
//	for (nodes_index_type i = 0; i < nodes.size(); i++)
//	{
//		//vertDubStrengths[i] = nodes[i]->linGetMu();
//		vertDubDiffs[i] = nodes[i]->linGetMu();
//	}
//
//	return vertDubDiffs;
//}


Eigen::Vector3d panel::linGetDubStrengths()
{
	Eigen::Vector3d vertDubStrengths;	// [mu1 mu2 mu3]

	for (nodes_index_type i = 0; i < nodes.size(); i++)
	{
		//vertDubStrengths[i] = nodes[i]->linGetMu();
		vertDubStrengths[i] = nodes[i]->linGetPotential();
	}

	return vertDubStrengths;
}


Eigen::Vector3d panel::linGetOrigDubStrengths()
{
	Eigen::Vector3d vertDubStrengths;	// [mu1 mu2 mu3]

	for (nodes_index_type i = 0; i < nodes.size(); i++)
	{
		vertDubStrengths[i] = nodes[i]->linGetMu();
		//vertDubStrengths[i] = nodes[i]->linGetPotential();
	}

	return vertDubStrengths;
}


Eigen::Vector3d panel::dubVInf(const Eigen::Vector3d &POI)
{
    // VSAero doublet velocity influence formulation
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();

    if (pjk.norm()/longSide > 5)
    {
        return pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s;
        nodes_index_type i1,i2;
        for (nodes_index_type i=0; i<nodes.size(); i++)
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

double panel::vortexPhi(const double &PN, const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l, const Eigen::Vector3d &m)
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

//Eigen::Vector3d panel::linPntDubPhi(const double &PN, const double &PJK, const Eigen::Vector3d &POIloc)
//{
//	double unitInf;
//	Eigen::Vector3d farImat, vertsPhi;
//	Eigen::Matrix3d vertsMat;
//		
//	unitInf = pntDubPhi(PN, PJK);
//	//farImat[0] = unitInf;
//	//farImat[1] = unitInf * POIloc.x();
//	//farImat[2] = unitInf * POIloc.y();
//	vertsMat = linVertsMatrix();
//	farImat[0] = unitInf;
//	farImat[1] = unitInf;
//	farImat[2] = unitInf;
//
//	vertsPhi = (farImat.transpose() * vertsMat.inverse());
//	//vertsPhi = farImat;
//
//	//double unitInfTest, stuff;
//	//Eigen::Vector3d vertsPhiTest;
//	////stuff = vertsPhi.mean();
//	//unitInfTest = pntDubPhi(PN, PJK)/3.0;
//	////vertsPhiTest = unitInfTest * vertsMat;
//	////stuff = vertsPhiTest.mean();
//	//vertsPhiTest[0] = unitInfTest;
//	//vertsPhiTest[1] = unitInfTest;
//	//vertsPhiTest[2] = unitInfTest;
//	
//	return vertsPhi;
//}

Eigen::Vector3d panel::pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk)
{
    return area*(3*pjk.dot(n)*pjk-pow(pjk.norm(),2)*n)/(4*M_PI*pow(pjk.norm(),5));
}

Eigen::Matrix<size_t, Eigen::Dynamic, 1> panel::getVerts()
{
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> verts(nodes.size());
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        verts(static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Index>(i)) = nodes[i]->getIndex();
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

    double x, y, z, x0, y0;
    x = this->getCenter().x(); y = this->getCenter().y(); z = this->getCenter().z();
    x0 = POI.x(); y0 = POI.y();
    
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
    nodes_index_type i1,i2;
    for (nodes_index_type i=0; i<nodes.size(); i++)
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
            std::cout << "Warning: particle is on panel..." << std::endl; // Not necessarily a bad thing if this panel is a wake panel. If it is a body panel, then a divergence free velocity field doesn't hold. Try increasing mesh resolution...
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
