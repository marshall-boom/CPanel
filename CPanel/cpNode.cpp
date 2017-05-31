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


std::vector<edge*> cpNode::getTrailingEdges(){
    
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

Eigen::Vector3d cpNode::nodeWakeProjAngle(){
    
    std::vector<edge*> tedges = this->getTrailingEdges();
    Eigen::Vector3d avgNorm = Eigen::Vector3d::Zero();

    for (int j = 0; j < tedges.size(); j++) {
        std::vector<bodyPanel*> bPans = tedges[j]->getBodyPans();
        
        // In the case that a side patch is included with the top and bottom panels,
        // the displacement in the z direction will work, but not the y. CPanel seems
        // to exclude these from the node values though...
        for (int j=0; j<bPans.size(); j++) {
            avgNorm += bPans[j]->getNormal();
        }
        
        avgNorm = avgNorm/avgNorm.size(); // Sum to average
    }
    
    return avgNorm.normalized();
}


Eigen::Vector3d cpNode::firstProjNode(double dt, double inputV){
    
    Eigen::Vector3d bisectVec = this->nodeWakeProjAngle();
    
    return this->getPnt() += c_w*dt*inputV*bisectVec;

}


Eigen::Vector3d cpNode::secProjNode(double dt, double inputV){
    
    Eigen::Vector3d bisectVec = this->nodeWakeProjAngle();

    return this->getPnt() += (c_w+1)*dt*inputV*bisectVec;

}



