//
//  node.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__node__
#define __CPanel__node__

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <assert.h>
#include "member.h"


template<typename type>
class node
{
    node<type>* parent;
    node<type>* children[8];
    Eigen::Vector3d origin;
    Eigen::Vector3d halfDimension;
    std::vector<member<type>> members;
    short level;
    short maxMembers;
    double maxTheta;

    void createChild(int childNumber)
    {
        // Child    0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
            Eigen::Vector3d tempOrigin = origin;
            Eigen::Vector3d tempHalfDimension;
            tempOrigin[0] += halfDimension[0] * (childNumber&4 ? 0.5 : -0.5);
            tempOrigin[1] += halfDimension[1] * (childNumber&2 ? 0.5 : -0.5);
            tempOrigin[2] += halfDimension[2] * (childNumber&1 ? 0.5 : -0.5);
            for (int i=0; i<3; i++)
            {
                tempHalfDimension[i] = 0.5*halfDimension[i];
            }
            children[childNumber] = new node<type>(this,tempOrigin,tempHalfDimension,level,maxMembers,maxTheta);
    }

    void pushMember(const member<type> &member)
    {
        // Checks if child is already created. If not, creates the right child before adding member.
        int child = getChildContainingMember(member);
        if (children[child] != NULL)
        {
            children[child]->addMember(member);
        }
        else
        {
            createChild(child);
            children[child]->addMember(member);
        }
    }

    void addChild(node<type>* child, int childNumber)
    {
        children[childNumber] = child;
    }

    void addLevel()
    {
        //Recursively adds one to the level of each node if parent is added
        level++;
        if (!isLeafNode())
        {
            for (int i=0; i<8; i++)
            {
                if(children[i]){
                    children[i]->addLevel();
                }
            }
        }
    }

    std::vector<type*> membersToObjects()
    {
        std::vector<type*> objects;
        for (int i=0; i<members.size(); i++)
        {
            type* obj = members[i].getObj();
            objects.push_back(obj);
        }
        return objects;
    }

public:
    node(node<type>* parent_ptr,Eigen::Vector3d origin,Eigen::Vector3d halfDimension, short parent_level,  short maxMembers, double maxTheta)
      : parent(parent_ptr),origin(origin),halfDimension(halfDimension),maxMembers(maxMembers),maxTheta(maxTheta),multExp(nullptr)
    {
        for (int i=0; i<8; i++)
        {
            children[i] = NULL;
        }

        level = parent_level+1;
    }

    node(const node<type>& copy)
      : parent(copy.parent), origin(copy.origin), halfDimension(copy.halfDimension),
		members(copy.members), level(copy.level), maxMembers(copy.maxMembers),
		maxTheta(copy.maxTheta), multExp(copy.multExp)
    {
        for (int i=0; i<8; i++)
        {
            children[i] = new node<type>(*copy.children[i]);
        }
    }

    ~node()
    {
        if (!isLeafNode())
        {
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL)
                {
                    delete children[i];
                }
            }
        }
    }

    node<type> operator=(node<type> other)
    {
        parent = other.parent;
        origin = other.origin;
        halfDimension = other.halfDimension;
        members = other.members;
        level = other.level;
        maxMembers = other.maxMembers;
        for (int i=0; i<8; i++)
        {
            children[i] = new node<type>(*other.children[i]);
        }
        return this;
    }


    void setMaxMembers(const int &max)
    {
        maxMembers = max;
    }

    void createParent(int corner)
    {
        // Corner   0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
        Eigen::Vector3d tempOrigin(origin);
        Eigen::Vector3d tempHalfDimension;
        tempOrigin[0] = origin[0]+halfDimension[0] * (corner&4 ? 1 : -1);
        tempOrigin[1] = origin[1]+halfDimension[1] * (corner&2 ? 1 : -1);
        tempOrigin[2] = origin[2]+halfDimension[2] * (corner&1 ? 1 : -1);
        for (int i=0; i<3; i++)
        {
            tempHalfDimension[i] = 2*halfDimension[i];
        }
        parent = new node<type>(NULL,tempOrigin,tempHalfDimension,level-1,maxMembers,maxTheta);

        int child = 0;
        for (int i=0; i<3; i++)
        {
            child = corner ^ 1<<i; //Flips the bits set for the corner with bitwise XOR
        }
        parent->addChild(this,child); //Adds current node to parent
        addLevel();
    }

    void addMember(const member<type> &member)
    {
        if (isLeafNode())  //If no children exist, add member to current node.  If they do, add member to proper child.
        {
            members.push_back(member);
            if (members.size() > maxMembers)
            {
                for (int i=0; i<members.size(); i++)
                {
                    pushMember(members[i]);
                }
                members.clear();
            }
        }
        else
        {
            pushMember(member);
        }
    }

    bool isLeafNode()
    {
        bool flag = true;
        for (int i=0; i<8; i++)
        {
            if (children[i] != NULL)
            {
                flag = false;
                break;
            }
        }
        return flag;
    }

    bool isSameLevel(node<type>* otherNode)
    {
        bool isSame = this->getLevel()==otherNode->getLevel()?true:false;
        return isSame;
    }

    Eigen::Vector3d calcVel(Eigen::Vector3d POI){

        Eigen::Vector3d velInfl = Eigen::Vector3d::Zero();
        if(this->isLeafNode())
        {
            std::vector<type*> nParts = this->getMembers();
            for(int i=0; i<nParts.size(); i++)
            {
                velInfl+=nParts[i]->velInflAlgSmooth(POI); // uses faster velocity formulation
            }
        }
        else
        {
            if(this->isFarField(POI))
            {
                velInfl = this->multExp->velInflAlgSmooth(POI);
            }
            else
            {
                std::vector<node<type>*> children = this->getChildren();
                for(int i=0; i<children.size() ; i++)
                {
                    velInfl += children[i]->calcVel(POI);
                }
            }
        }

        return velInfl;
    }

    Eigen::Vector3d calcVel(type* part){

        Eigen::Vector3d velInfl = Eigen::Vector3d::Zero();
        if(this->isLeafNode())
        {
            std::vector<type*> nParts = this->getMembers();
            for(int i=0; i<nParts.size(); i++)
            {
                if(part != nParts[i])
                {
                    velInfl += nParts[i]->velInfl(part);
                }
            }
        }
        else
        {
            if(this->isFarField(part->pos))
            {
                velInfl = this->multExp->velInfl(part);
            }
            else
            {
                std::vector<node<type>*> children = this->getChildren();
                for(int i=0; i<children.size() ; i++)
                {
                    velInfl += children[i]->calcVel(part);
                }
            }
        }

        return velInfl;
    }

    Eigen::Vector3d calcStretch(type* part){

        Eigen::Vector3d stretchInfl = Eigen::Vector3d::Zero();
        if(this->isLeafNode())
        {
            std::vector<type*> nParts = this->getMembers();
            for(int i=0; i<nParts.size(); i++)
            {
                if(nParts[i] != part) // Kroneger delta func.
                {
                    stretchInfl += nParts[i]->vortexStretching(part);
                }
            }
        }
        else
        {
            if(this->isFarField(part->pos))
            {
                stretchInfl = this->multExp->vortexStretching(part);
            }
            else
            {
                std::vector<node<type>*> children = this->getChildren();
                for(int i=0; i<children.size() ; i++)
                {
                    stretchInfl += children[i]->calcStretch(part);
                }
            }
        }

        return stretchInfl;
    }


    Eigen::Vector3d calcDiff(type* part){

        Eigen::Vector3d diffInfl = Eigen::Vector3d::Zero();
        if(this->isLeafNode())
        {
            std::vector<type*> nParts = this->getMembers();
            for(size_t i=0; i<nParts.size(); i++)
            {
                if(nParts[i] != part) // Kroneger delta func.
                {
                    diffInfl += nParts[i]->viscousDiffusion(part);
                }
            }
        }
        else
        {
            if(this->isFarField(part->pos))
            {
                diffInfl = this->multExp->viscousDiffusion(part);
            }
            else
            {
                std::vector<node<type>*> children = this->getChildren();
                for(int i=0; i<children.size() ; i++)
                {
                    diffInfl += children[i]->calcDiff(part);
                }
            }
        }

        return diffInfl;
    }



//    bool isFarField(node<type>* otherNode)
//    {
//        if(otherNode->getLevel() == 2){
//            return false;
//        }
//        double centerToCenterDist = (this->getOrigin()-otherNode->getOrigin()).norm();
//        bool sameLevel = isSameLevel(otherNode);
//
//        if(std::abs(centerToCenterDist) > 2*(this->getHalfDimension()).norm() && sameLevel){
//            return true;
//        }else{
//            return false;
//        }
//    };

    bool isFarField(node<type>* otherNode){
        // Follows documentation in FTM chapter of Connor's thesis

        // Width of Node
        double s = 2*this->halfDimension.x(); // The 'x' can be used because the implementation of the octree keeps each node cubic

        // Particle distance
        double d = (this->multExp->pos - otherNode->multExp->pos).norm();

        // theta = s/d;
        if(s/d < maxTheta)
            // maxTheta was put in the particle class to allow for different values for panels if BH is implemented in the future
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    bool isFarField(Eigen::Vector3d POI){

        // Width of Node
        double s = 2*this->halfDimension.x(); // The 'x' can be used because the implementation of the octree keeps each node cubic

        // Particle distance
        double d = (this->multExp->pos - POI).norm();

        // theta = s/d;
        if(s/d < maxTheta)
        {
            return true;
        }
        else
        {
            return false;
        }
    };


    int getChildContainingMember(const member<type> &member)
    {
        return getChildContainingPnt(member.getRefPoint());
    }

    int getChildContainingPnt(const Eigen::Vector3d &pnt)
    { // Finds child quadrant
        int child = 0;
        if (pnt(0) > origin[0])
        {
            child |= 4;
        }
        if (pnt(1) > origin[1])
        {
            child |= 2;
        }
        if (pnt(2) > origin[2])
        {
            child |= 1;
        }

        return child;
    }

    std::vector<type*> getMembers()
    {
        if (isLeafNode())
        {
            return membersToObjects();
        }
        else
        {
            std::vector<type*> recursiveMembers;
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL)
                {
                    std::vector<type*> temp = children[i]->getMembers();
                    for (int j=0; j<temp.size(); j++)
                    {
                        recursiveMembers.push_back(temp[j]);
                    }
                }
            }
            return recursiveMembers;
        }

    }

    std::vector<node<type>*> getSubNodes()
    {
        if (isLeafNode())
        {
            std::vector<node<type>*> temp;
            temp.push_back(this);
            return temp;
        }
        else
        {
            std::vector<node<type>*> recursiveNodes;
            recursiveNodes.push_back(this);
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL)
                {
                    std::vector<node<type>*> temp = children[i]->getSubNodes();
                    for (int j=0; j<temp.size(); j++)
                    {
                        recursiveNodes.push_back(temp[j]);
                    }
                }
            }
            return recursiveNodes;
        }
    }

    std::vector<type*> getMembers(node<type>* exception)
    {
        // Used if searching tree starting at bottom to avoid searching the same node twice.
        if (isLeafNode() && this != exception)
        {
            return membersToObjects();
        }
        else if (isLeafNode() && this == exception)
        {
            std::vector<type*> empty;
            return empty;
        }
        else
        {
            std::vector<type*> recursiveMembers;
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL && children[i] != exception)
                {
                    std::vector<type*> temp = children[i]->getMembers(exception);
                    for (int j=0; j<temp.size(); j++)
                    {
                        recursiveMembers.push_back(temp[j]);
                    }
                }
            }
            return recursiveMembers;
        }
    }

    bool hasChildren()
    {
        for(int i=0; i<8; i++)
        {
            node<type>* subNode = this->getChild(i);
            if(subNode)
            {
                return true;
            }
        }
        return false;
    }

    std::vector<node<type>*> getChildren()
    {
        std::vector<node<type>*> children;
        for(int i=0; i<8; i++)
        {
            node<type>* subNode = this->getChild(i);
            if(subNode)
            {
                children.push_back(subNode);
            }
        }
        return children;
    }

    std::vector<type*> getChildExpans()
    {
        std::vector<node<type>*> children = this->getChildren();
        std::vector<type*> childExps;

        for(int i=0; i<children.size(); i++)
        {
            childExps.push_back(children[i]->multExp);
        }

        return childExps;
    }

    type* multExp; // Public because will access a lot
    void setMultExp(type* exp) {multExp = exp;};

    Eigen::Vector3d getOrigin() {return origin;}
    Eigen::Vector3d getHalfDimension() {return halfDimension;}
    node<type>* getChild(int childNumber) {return children[childNumber];}
    short getLevel() {return level;}
    node<type>* getParent() {return parent;}
};

#endif /* defined(__CPanel__node__) */




