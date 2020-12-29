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

#ifndef __CPanel__panel__
#define __CPanel__panel__

#include <iostream>
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include "math.h"
#include "convexHull.h"
#include "panelOctree.h"
#include "chtlsnd.h"
//#include "cpNode.h"
//#include "edge.h"

class panelOctree;
class edge;
class cpNode;
class surface;

class panel
{    
protected:
	using nodes_type = std::vector<cpNode *>;
	using nodes_index_type = nodes_type::size_type;
	using edges_type = std::vector<edge *>;
	using edges_index_type = edges_type::size_type;

    nodes_type nodes;
    edges_type pEdges;
    Eigen::Vector3d center;
    Eigen::Vector3d normal;
    Eigen::Vector3d bezNormal; //Used in derivative calculation
    double area;
    double longSide;

	//lin
	nodes_type bodyNodes;
	Eigen::Vector3d linVelocity;
	std::vector<Eigen::Vector3d> linLocalNodes;

    double doubletStrength = 0;
    double potential = 0;
    double prevPotential = 0;
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    double Cp;
    size_t ID;
    double core = 0.05;
    
    double vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m);
    Eigen::Vector3d getUnitVector(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);
    Eigen::Matrix3d getLocalSys();
    
    double pntDubPhi(const double &PN, const double &PJK);
    
    Eigen::Vector3d pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk);
    
    Eigen::Matrix3d velocityGradientPointDoublet(Eigen::Vector3d POI);
    Eigen::Matrix3d velocityGradientDoublet(Eigen::Vector3d POI);
    Eigen::Matrix3d gradDoub(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s);

public:
    panel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, size_t surfID);
    
    virtual ~panel() {}
    
    void setGeom();
        
    void setPotential(Eigen::Vector3d Vinf);

    Eigen::Vector3d vortexV(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s); // Public so vortFil can access it.

    bool inPanelProjection(const Eigen::Vector3d &POI, Eigen::Vector3d &projectedPnt);
    bool onPanelCheck(const Eigen::Vector3d &POI);

    Eigen::Vector3d global2local(const Eigen::Vector3d &globalVec,bool translate);
    Eigen::Vector3d local2global(const Eigen::Vector3d &localVec,bool translate);
    
    double dubPhiInf(const Eigen::Vector3d &POI);
    Eigen::Vector3d dubVInf(const Eigen::Vector3d &POI);
    Eigen::Vector3d pntDubVInf(const Eigen::Vector3d &POI);

    virtual double panelPhi(const Eigen::Vector3d &POI) = 0;
    virtual Eigen::Vector3d panelV(const Eigen::Vector3d &POI) = 0;
    
    std::vector<Eigen::Vector3d> pntsAroundPnt(int nPnts,const Eigen::Vector3d &POI, double r);
    Eigen::Vector3d pntNearEdge(edge* e);
    
    size_t getID() {return ID;}
    Eigen::Vector3d getCenter() const {return center;}
    Eigen::Vector3d getNormal() const {return normal;}
    Eigen::Vector3d getBezNormal() const {return bezNormal;}
    std::vector<cpNode*> getNodes() {return nodes;}
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> getVerts();
    std::vector<edge*> getEdges() {return pEdges;}
    double getLongSide() {return longSide;}
    double getArea() {return area;}
    double getMu() {return doubletStrength;}
    double getPotential() {return potential;}
    
    bool nearFilamentCheck(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &POI);

	void linDubPhiInf(const Eigen::Vector3d &POI, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow);
	void linPhiHintegrals(Eigen::VectorXd &Hints, const double g, const double myAl, const double Al, const double l1, const double l2, const double c1, const double c2, const double nuEta, const double nuXi, const double &PN, const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l, const Eigen::Vector3d &m);
	//Eigen::Vector3d linPntDubPhi(const double &PN, const double &PJK, const Eigen::Vector3d &POIloc);
	Eigen::Matrix3d linVertsMatrix(bool translate);

	Eigen::Vector3d linGetDubStrengths();
	void linGetConstDubStrength();
	void linComputeVelocity(double PG,Eigen::Vector3d &Vinf);

	Eigen::Vector3d getVel() { return velocity; }

	Eigen::Matrix3d linDubVInf(const Eigen::Vector3d &POI);

	void linVelHintegrals(Eigen::VectorXd &Hints, double F123, Eigen::Vector3d &Eints, const double g, const double Al, const double nuEta, const double nuXi);
	Eigen::Matrix3d linVelJintegrals(Eigen::VectorXd &Hints, const double &PN);
	void setPanelVel(Eigen::Vector3d Vel);
	Eigen::Vector3d linGetOrigDubStrengths();

	Eigen::Matrix3d supVertsMatrix(std::vector<Eigen::Vector3d> &supNodes);
	Eigen::Vector3d supGetDubDiffs();


	void supOutputGeom(const Eigen::Vector3d &POI, bool outPOI);

};

#endif /* defined(__CPanel__panel__) */
