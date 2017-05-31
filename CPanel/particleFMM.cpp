//
//  particleFMM.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 5/10/16.
//  Copyright (c) 2016 Connor Sousa. All rights reserved.
//

#include "particleFMM.h"
#include "particle.h"
#include "cpCase.h"

void particleFMM::build(particleOctree* tree)
{

    partTree = tree;
    levels = tree->numTreeLevels();

    computeMultExp();
}

void particleFMM::computeMultExp()
{
    for(int i=levels+1; i>1; i--)
    { // Ignore root node. 1 based indexing for tree
        std::vector<node<particle> *> nodes = partTree->getLevelNodes(i);
        
        for (int j=0; j<nodes.size(); j++)
        {
            if (nodes[j]->isLeafNode())
            {
                // Simply collect members and find expansion
                std::vector<particle*> parts = nodes[j]->getMembers();
                Eigen::Vector3d pos = findExpPos(nodes[j]);
                Eigen::Vector3d strength = findExpStrength(nodes[j]);
                double radius = nodes[j]->getHalfDimension().norm()*2*1.5;
                
                particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, 0); // dummy timestep.
                nodes[j]->setMultExp(p);
            }
            else
            {
                // Calculate it from children
                std::vector<particle*> childExpans = nodes[j]->getChildExpans();
                
                Eigen::Vector3d pos = findExpPos(childExpans);
                Eigen::Vector3d strength = findExpStrength(childExpans);
                double radius = nodes[j]->getHalfDimension().norm()*2*1.5;//findExpRadius(childExpans);
                
                particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, 0);
                nodes[j]->setMultExp(p);
            }
        }
    }
}


Eigen::Vector3d particleFMM::barnesHutVel(Eigen::Vector3d POI)
{
    // Ignore the root node. Root node wouldn't be beneficial in implementation since octree is always near POI. Instead, start at level 2
    
    Eigen::Vector3d velInfl = Eigen::Vector3d::Zero();
    std::vector<node<particle>*> rootChildren = partTree->getRootNode()->getChildren();
    
    for (int i=0; i<rootChildren.size(); i++)
    {
        velInfl+=rootChildren[i]->calcVel(POI);
    }

    return velInfl;
}

Eigen::Vector3d particleFMM::barnesHutVel(particle* part)
{
    // Ignore the root node. Root node wouldn't be beneficial in implementation since octree is always near POI. Instead, start at level 2
    
    Eigen::Vector3d velInfl = Eigen::Vector3d::Zero();
    std::vector<node<particle>*> rootChildren = partTree->getRootNode()->getChildren();
    
    for (int i=0; i<rootChildren.size(); i++)
    {
        velInfl += rootChildren[i]->calcVel(part);
    }
    
    return velInfl;
}

Eigen::Vector3d particleFMM::barnesHutStretch(particle* part)
{
    Eigen::Vector3d stretchInfl = Eigen::Vector3d::Zero();
    std::vector<node<particle>*> rootChildren = partTree->getRootNode()->getChildren();
    
    for (int i=0; i<rootChildren.size(); i++)
    {
        stretchInfl += rootChildren[i]->calcStretch(part);
    }
    
    return stretchInfl;
}

Eigen::Vector3d particleFMM::barnesHutDiff(particle* part)
{
    Eigen::Vector3d diffInfl = Eigen::Vector3d::Zero();
    std::vector<node<particle>*> rootChildren = partTree->getRootNode()->getChildren();
    
    for (int i=0; i<rootChildren.size(); i++)
    {
        diffInfl += rootChildren[i]->calcDiff(part);
    }
    
    return diffInfl;
}




Eigen::Vector3d particleFMM::findExpPos(node<particle>* thisNode)
{
    // Weighted position of the strength magnitude
    
    std::vector<particle*> nParts = thisNode->getMembers();
    
    Eigen::Vector3d weightedPos = Eigen::Vector3d::Zero();
    double totStrength = 0;
    
    for (int i=0; i<nParts.size(); i++) {
        weightedPos += (nParts[i]->pos) * (nParts[i]->strength).norm();
        totStrength += (nParts[i]->strength).norm();
    }
    
    return weightedPos/totStrength;
}

Eigen::Vector3d particleFMM::findExpPos(std::vector<particle*> parts)
{
    // Weighted position of the strength magnitude
    
    Eigen::Vector3d weightedPos = Eigen::Vector3d::Zero();
    double totStrength = 0;
    
    for (int i=0; i<parts.size(); i++) {
        weightedPos += (parts[i]->pos) * (parts[i]->strength).norm();
        totStrength += (parts[i]->strength).norm();
    }
    
    return weightedPos/totStrength;
}


Eigen::Vector3d particleFMM::findExpStrength(node<particle>* thisNode){
    std::vector<particle*> nodeParticles = thisNode->getMembers();
    
    Eigen::Vector3d nodeStren = Eigen::Vector3d::Zero();
    
    for (int i=0; i<nodeParticles.size(); i++) {
        nodeStren+=nodeParticles[i]->strength;
    }
    
    return nodeStren;
}

Eigen::Vector3d particleFMM::findExpStrength(std::vector<particle*> parts){
    
    Eigen::Vector3d nodeStren = Eigen::Vector3d::Zero();
    
    for (int i=0; i<parts.size(); i++) {
        nodeStren+=parts[i]->strength;
    }
    
    return nodeStren;
}

double particleFMM::findExpRadius(node<particle>* thisNode){
    
    std::vector<particle*> nodeParticles = thisNode->getMembers();
    
    double nodeRad = 0;
    
    for (int i=0; i<nodeParticles.size(); i++) {
        nodeRad+=nodeParticles[i]->radius;
    }
    
    return nodeRad;
}

double particleFMM::findExpRadius(std::vector<particle*> parts){
    // Overloaded function does same as above, but intakes a vector of particles
    
    double nodeRad = 0;
    
    for (int i=0; i<parts.size(); i++) {
        nodeRad+=parts[i]->radius;
    }
    
    return nodeRad;
}


//write first here, then make a general case in a new folder above using a general type.








