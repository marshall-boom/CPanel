//
//  cpNode.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 2/2/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "cpNode.h"
#include "edge.h"
#include "bodyPanel.h"


cpNode::cpNode(Eigen::Vector3d pnt,int index) : pnt(pnt), index(index), TEnode(false) {}

//cpNode::cpNode(const cpNode& copy) : pnt(copy.pnt), index(copy.index),TEnode(copy.TEnode)
//{
//    
//}

Eigen::Vector3d cpNode::operator-=(const cpNode &rhs)
{
    return pnt-rhs.getPnt();
}

Eigen::Vector3d operator-(cpNode lhs, const cpNode &rhs)
{
    return lhs -= rhs;
}

Eigen::Vector3d cpNode::operator+=(const cpNode &rhs)
{
    return pnt+rhs.getPnt();
}

Eigen::Vector3d operator+(cpNode lhs, const cpNode &rhs)
{
    return lhs += rhs;
}

void cpNode::addEdge(edge* e)  {edges.push_back(e);}

void cpNode::addBodyPanel(bodyPanel* p) {bodyPans.push_back(p);}

edge* cpNode::getTE(edge* exception)
{
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != exception)
        {
            if (edges[i]->getN2() == this)
            {
                edges[i]->flipDir();
            }
            return edges[i];
        }
    }
    return nullptr;
}

void cpNode::setTE() {TEnode = true;}

void cpNode::setIndex(int i) {index = i;}

edge* cpNode::getOtherTrailEdge(edge* current)
{
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != current)
        {
            return edges[i];
        }
    }
    return nullptr;
}


std::vector<edge*> cpNode::getTrailingEdges()
{
    std::vector<edge*> trailingEdges;
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE())
        {
            trailingEdges.push_back(edges[i]);
        }
    }
    return trailingEdges;
}


double cpNode::nodeWakeProjAngle(cpNode* tePoint){
    std::vector<edge*> tedges = tePoint->getTrailingEdges();
    
    double bisect = 0;
    for (int j = 0; j < tedges.size(); j++) {
        std::vector<bodyPanel*> tempPan = tedges[j]->getBodyPans();
        
        Eigen::Vector3d nVec1 = tempPan[0]->getNormal();
        Eigen::Vector3d nVec2 = tempPan[1]->getNormal();
        
        // unit vec to angle
        double thetaTop = acos(nVec1.x()/pow(nVec1.x()*nVec1.x()+nVec1.y()*nVec1.y()+nVec1.z()*nVec1.z(),.5));
        double thetaBot = acos(nVec2.x()/pow(nVec2.x()*nVec2.x()+nVec2.y()*nVec2.y()+nVec2.z()*nVec2.z(),.5));
        bisect += (thetaTop+thetaBot)/2-1.570796;
        
    }
    bisect = bisect/tedges.size();
    
    return bisect;
}


Eigen::Vector3d cpNode::firstProjNode(cpNode* TEnode, double dt, double c_w, double inputV){
    
    double bisectAngle = TEnode->nodeWakeProjAngle(TEnode);
    Eigen::Vector3d proj1 = TEnode->getPnt(); // Temporarily set it to the node to build off of it
    proj1.x() += c_w*dt*inputV*cos(bisectAngle);
    proj1.z() += c_w*dt*inputV*sin(bisectAngle);
    
    // Project node out perpendicular to panel
    std::vector<edge*> tedges = TEnode->getTrailingEdges();
    double edgeAngle=0; //Angle perpedicular to edge to find projected node y coord.
    for (int j = 0; j < tedges.size(); j++) {
        Eigen::Vector3d node1 = tedges[j]->getN1()->getPnt();
        Eigen::Vector3d node2 = tedges[j]->getN2()->getPnt();
        
        edgeAngle += atan((node2.x()-node1.x())/(node2.y()-node1.y())); // Aircraft coordinates are different from traditional y/x tangent
    }
    edgeAngle = edgeAngle/tedges.size();
    proj1.y() -= c_w*dt*inputV*sin(edgeAngle);
    
    return proj1;
}



Eigen::Vector3d cpNode::secProjNode(cpNode* TEnode, double dt, double c_w, double inputV){
    
    double Uinf = 10;
    
    double bisectAngle = TEnode->nodeWakeProjAngle(TEnode);
    Eigen::Vector3d proj2 = TEnode->getPnt();
    proj2.x() += (c_w+1)*dt*inputV*cos(bisectAngle);
    proj2.z() += (c_w+1)*dt*inputV*sin(bisectAngle);
    
    // Project node out perpendicular to panel
    std::vector<edge*> tedges = TEnode->getTrailingEdges();
    
    double edgeAngle=0;
    for (int j = 0; j < tedges.size(); j++) {
        Eigen::Vector3d node1 = tedges[j]->getN1()->getPnt();
        Eigen::Vector3d node2 = tedges[j]->getN2()->getPnt();
        
        edgeAngle += atan((node2.x()-node1.x())/(node2.y()-node1.y()));
    }
    edgeAngle = edgeAngle/tedges.size();
    proj2.y() -= (c_w+1)*dt*inputV*sin(edgeAngle);
    
    return proj2;
}

