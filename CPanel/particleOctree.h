//
//  particleOctree.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 5/9/16.
//  Copyright (c) 2016 Connor Sousa. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__particleOctree__
#define __CPanel___Unstructured_Panel_Code__particleOctree__

#include <iostream>
#include "octree.h"
#include "particle.h"

class particle;

class particleOctree : public octree<particle>
{
    
public:
    particleOctree() : octree() {}
    
    Eigen::Vector3d findRefPoint(const particle &obj);
    Eigen::Vector3d sumPartStrengths(node<particle*>* node);
};



#endif /* defined(__CPanel___Unstructured_Panel_Code__particleOctree__) */
